"""
    apeNodeNumbers!(net::HybridNetwork)

Change numbers of internal nodes of `net` to satisfy conditions assumed by the `ape` R 
package: leaves are 1-n, the root is n+1, and internal nodes are higher consecutive 
integers.

# Examples

```julia-repl
julia> net = readTopology("(A,(B,(C,D)));")
julia> directEdges!(net)
julia> preorder!(net)
julia> apeNodeNumbers!(net)
julia> net.node
7-element Array{PhyloNetworks.Node,1}:
 PhyloNetworks.Node:
 number:4
 name:A
 leaf node
 attached to 1 edges, numbered: 1

 PhyloNetworks.Node:
 number:3
 name:B
 leaf node
 attached to 1 edges, numbered: 2

 PhyloNetworks.Node:
 number:2
 name:C
 leaf node
 attached to 1 edges, numbered: 3

 PhyloNetworks.Node:
 number:1
 name:D
 leaf node
 attached to 1 edges, numbered: 4

 PhyloNetworks.Node:
 number:7
 attached to 3 edges, numbered: 3 4 5
               
 PhyloNetworks.Node:
 number:6
 attached to 3 edges, numbered: 2 5 6
               
 PhyloNetworks.Node:
 number:5
 attached to 2 edges, numbered: 1 6
```
"""

function apeNodeNumbers!(net::HybridNetwork)
    lnum = 1 # first number for leaves
    inum = length(net.leaf) + 1 # first number of internal nodes
    #ensure internal node numbers do not overlap with leaf numbers
    #root will be labeled ntips +1, consistent with ape "phylo" class
    for n in net.nodes_changed # topological (pre)order
        if n.leaf
            n.number = lnum
            lnum += 1
        else
            n.number = inum # root will be ntips + 1, because pre-order
            inum += 1
        end
    end
end

"""
    generateMajorEdge(net::HybridNetwork)

Generate matrix of major edges from `net` where edge[i,1]  is the number of the
parent node of edge i and edge[i,2] is the number of the child node of edge i.
Assume `nodes_changed` was updated, to list edges in pre-order.

# Examples

```julia-repl
julia> net = readTopology("(A,(B,(C,D)));")
julia> directEdges!(net)
julia> preorder!(net)
julia> apeNodeNumbers!(net)
julia> generateMajorEdge(net)
6×2 Array{Int64,2}:
 5  4
 5  6
 6  3
 6  7
 7  2
 7  1
```
"""

function generateMajorEdge(net::HybridNetwork)  
    edge = Matrix{Int}(length(net.edge)-length(net.hybrid), 2) # major edges
    i = 1 #row index for edge matrix
    for n in net.nodes_changed # topological pre-order
        !n.leaf || continue # skip leaves: associate node with children edges
        for e in n.edge
            if e.node[e.isChild1 ? 2 : 1] == n && e.isMajor
                #exclude parent edge and minor hybrid edges
                edge[i,1] = n.number # parent
                edge[i,2] = e.node[e.isChild1 ? 1 : 2].number # child
                i += 1 #increase index value
            end
        end
    end
    return edge
end

"""
    generateMajorLength(net::HybridNetwork)

Generate vector of edge lengths of major `net` edges organized in the same order 
as the `edge` matrix created via `generateMajorEdge`. Replace values of `1.0` with `#NULL` values
recognized by the `ape` library for the `R` programming language.


# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
julia> directEdges!(net)
julia> preorder!(net)
julia> apeNodeNumbers!(net)
julia> generateMajorLength(net)
8-element NullableArrays.NullableArray{Float64,1}:
#NULL
#NULL
#NULL
#NULL
#NULL
#NULL
#NULL
#NULL
```
""" #"

function generateMajorLength(net::HybridNetwork) 
    edgeLength = Array{Float64}(length(net.edge)-length(net.hybrid)) 
    i=1
    for n in net.nodes_changed #traverse tree in topological order excluding leaves
        if !n.leaf
            for e in n.edge #iterate over each edge attatched to a node
                #exclude node that is parent of current edge and minor hybrid edge
                if e.node[e.isChild1 ? 2 : 1] == n && e.isMajor 
                    edgeLength[i] = e.length
                    i=i+1
                end
            end
        end
    end
    #-1.0 interpreted as missing
    edgeLength = NullableArray(edgeLength, map(x -> x==-1.0, edgeLength))
    return edgeLength
end

"""
    generateMinorReticulation(net::HybridNetwork)

Generate a matrix of minor hybrid edges from `net` where edge[i,1] represents 
the parent node of edge i and edge[i,2] represents the child node of edge i.

Assume nodes have been organized for preorder traversal and node numbering satisfies 
the conditions assumed by the `ape` package in R.

# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
julia> directEdges!(net)
julia> preorder!(net)
julia> apeNodeNumbers!(net)
julia> generateMinorReticulation(net)
1×2 Array{Int64,2}:
 7  9
```
""" #"

function generateMinorReticulation(net::HybridNetwork)    
    reticulation = Matrix{Int}(length(net.hybrid), 2) #initialize reticulation matrix 
    #fill reticulation matrix
    j = 1 #row index for reticulation matrix
    for e in net.edge #iterate over each edge attatched to a hybrid node
        if !e.isMajor #find minor hybrid edge
            #push parent node to first column
            reticulation[j,1] = e.node[e.isChild1 ? 2 : 1].number
            #push child node to second column
            reticulation[j,2] = e.node[e.isChild1 ? 1 : 2].number 
            j += 1 #increase index value
        end
    end
    return reticulation
end

"""
    generateMinorReticulationLength(net::HybridNetwork)

Generate vector of minor edge lengths organized in the same order as the 
`reticulation` matrix created via `generateMinorReticulation`. 
Replace values of `1.0` with `#NULL` values recognized by the `ape` library for 
the `R` programming language.

# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
julia> directEdges!(net)
julia> preorder!(net)
julia> apeNodeNumbers!(net)
julia> generateApeReticulationLength(net)
1-element NullableArrays.NullableArray{Float64,1}:
#NULL
```
""" #"

function generateMinorReticulationLength(net::HybridNetwork) 
    reticulationLength = Vector{Float64}(0) #initialize 
    for e in net.edge #iterate over each edge attatched to a hybrid node
        if !e.isMajor #find minor hybrid edge
            push!(reticulationLength, e.length)
        end
    end
    reticulationLength = NullableArray(reticulationLength, map(x -> x==-1.0, reticulationLength))
    return reticulationLength
end

doc"""
    apeRExport(net::HybridNetwork)

Export `net` to R as a `evonet` or `phylo` object (depending on degree of hybridization) 
recognized by the `ape` library in R.

# Arguments

- useEdgeLength: if true, export edge lengths from `net`.
- mainTree: if true, minor hybrid nodes omitted.


# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
julia> phy = RExport(net)
```

This creates an RObject called `phy` that is essentially an `phylo` (or `evonet`) 
class S3 object. The resulting RObject can be evaluated using the tools
available in the `ape` library in R. For example, we can visualize the network using 
the `plot` function:

```julia-repl
julia> using RCall
R> library(ape)
R> plot.evonet($phy)
```
""" #" 

function apeRExport(net::HybridNetwork; mainTree::Bool=false, useEdgeLength::Bool=true) 
    if mainTree == true && net.numHybrids > 0
        net = majorTree(net)
    end
    directEdges!(net)
    preorder!(net) # create field nodes_changed
    apeNodeNumbers!(net)
    ntips = length(net.leaf) 
    totalnodes = length(net.node) 
    Nnode = totalnodes - ntips 
    o = sortperm([n.number for n in net.leaf])
    tipLabel = [net.leaf[i].name for i in o]
    edge = generateMajorEdge(net)
    R"""
    phy = list(Nnode = $Nnode, edge = $edge, tip.label = $tipLabel)
    """
    if useEdgeLength == true
        edgeLength = generateMajorLength(net)
        R"""
        phy[['edge.length']] = $edgeLength
        """
    end
    if net.numHybrids > 0
        reticulation = generateMinorReticulation(net)
        R"""
        phy[['reticulation']] = $reticulation
        class(phy) <- c("evonet", "phylo")
        """
        if useEdgeLength == true # extract minor edge lengths
            reticulationLength = generateMinorReticulationLength(net)
            R"""
            phy[['reticulation.length']] = $reticulationLength
            """
        end
    elseif net.numHybrids == 0
        R"""
        class(phy) = "phylo"
        """
    end
    phy = reval("phy")
    return(phy)
end

doc"""
    function sexp(net::HybridNetwork)

Sexp method to export HybridNework objects to the R language
as either `phylo` or `evonet` object (depending on degree of hybridization)
recognized by the R package `ape`.

Inspired by https://github.com/richardreeve/Phylo.jl/blob/master/src/rcall.jl

# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
julia> sobj = sexp(net)
```

This creates an RObject called `phy` that is essentially an `phylo` (or `evonet`)
class S3 object. The resulting RObject can be evaluated using the tools
available in the `ape` library in R. For example, we can visualize the network using
the `plot` function:

```julia-repl
julia> using RCall
R> library(ape)
R> plot.evonet($sobj)
```
""" #"

function sexp(net::HybridNetwork)
    preorder!(net) #organize nodes for preorder traversal
    directEdges!(net)
    apeNodeNumbers!(net)
    ntips = length(net.leaf)
    totalnodes = length(net.node)
    Nnode = totalnodes - ntips
    o = sortperm([n.number for n in net.leaf])
    tipLabel = [net.leaf[i].name for i in o]
    edge = generateMajorEdge(net) #generate edge matrix
    #create an object with fields for $Nnode, $tipLabel, and $edge
    phy = Dict{Symbol, Any}()
    phy[:Nnode] = Nnode
    phy[Symbol("tip.label")] = tipLabel
    phy[:edge] = edge
    edgeLength = generateMajorLength(net) #major edges only
    nBL = sum([isnull(e) for e in edgeLength])
    if nBL>0
        phy[Symbol("edge.length")] = edgeLength
    end
    if net.numHybrids > 0
        reticulation = generateMinorReticulation(net) #minor edges only
        reticulationLength = generateMinorReticulationLength(net)
        phy[:reticulation] = reticulation
        nBL = sum([isnull(e) for e in reticulationLength])
        if nBL >0
            phy[Symbol("reticulation.length")] = reticulationLength
        end
        # fixit: export gamma if available
    end
    sobj = protect(sexp(phy))
    if net.numHybrids == 0
        setclass!(sobj, sexp("phylo"))
    else
        setclass!(sobj, sexp(["evonet", "phylo"]))
    end
    unprotect(1)
    return(sobj) #export RObject
end
