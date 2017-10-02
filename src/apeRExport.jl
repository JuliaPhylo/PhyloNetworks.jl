"""

    apeNodeNumbers(net::HybridNetwork)

Change numbers of internal nodes of `net` to satisfy conditions assumed by the `ape` R 
package: leaves are 1-n, the root is n+1, and internal nodes are higher consecutive 
integers.

# Examples

```julia-repl
julia> net = readTopology("(A,(B,(C,D)));")
julia> directEdges!(net)
julia> preorder!(net)
julia> apeNodeNumbers(net)
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

function apeNodeNumbers(net::HybridNetwork)  
    lnum = 1
    inum = length(net.leaf) + 1 
    #ensure internal node numbers do not overlap with leaf numbers
    #root will be labeled ntips +1, consistent with ape "phylo" class
    for n in net.nodes_changed #traverse tree in topological order excluding leaves
        if n.leaf
            n.number = lnum
            lnum += 1
        else
            n.number = inum #change node number
            inum += 1
        end
    end
end

"""

    generateMajorEdge(net::HybridNetwork)

Generate matrix of major edges from `net` where edge[i,1] 
represents the parent node of edge i and edge[i,2] represents the child node of edge i. 

Assume nodes have been organized for preorder traversal 
and node numbering satisfies the conditions assumed by the `ape` package in R.

# Examples

```julia-repl
julia> net = readTopology("(A,(B,(C,D)));")
julia> directEdges!(net)
julia> preorder!(net)
julia> apeNodeNumbers(net)
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
    edge = Matrix{Int}(length(net.edge)-length(net.hybrid), 2) #initialize edge matrix 
    #fill edge matrix
    i = 1 #row index for edge matrix
    for n in net.nodes_changed #traverse tree in topological order excluding leaves
        if !n.leaf
            for e in n.edge #iterate over each edge attatched to a node
                if e.node[e.isChild1 ? 1 : 2] != n && e.gamma > 0.5 
                #exclude node that is parent of current edge and minor hybrid edge
                    #push parent node to first column
                    edge[i,1] = e.node[e.isChild1 ? 2 : 1].number 
                    #push child node to second column
                    edge[i,2] = e.node[e.isChild1 ? 1 : 2].number 
                    i += 1 #increase index value
                end
            end
        end
    end
    return edge
end

"""
    generateMajorLength(net::HybridNetwork)

Generate vector of edge lengths of major `net` edges organized in the same order 
as the `edge` matrix created via `generateApeEdge`.

# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
julia> directEdges!(net)
julia> preorder!(net)
julia> apeNodeNumbers(net)
julia> generateMajorLength(net)
8-element Array{Float64,1}:
 -1.0
 -1.0
 -1.0
 -1.0
 -1.0
 -1.0
 -1.0
 -1.0
```

""" #"

function generateMajorLength(net::HybridNetwork) 
    #initialize edge length vector 
    edgeLength = Array{Float64}(length(net.edge)-length(net.hybrid)) 
    #fill edge length vector
    i=1
    for n in net.nodes_changed #traverse tree in topological order excluding leaves
        if !n.leaf
            for e in n.edge #iterate over each edge attatched to a node
                #exclude node that is parent of current edge and minor hybrid edge
                if e.node[e.isChild1 ? 1 : 2] != n && e.gamma > 0.5 
                    edgeLength[i] = e.length
                    i=i+1
                end
            end
        end
    end
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
julia> apeNodeNumbers(net)
julia> generateMinorReticulation(net)
1×2 Array{Int64,2}:
 7  9
 
```

""" #"

function generateMinorReticulation(net::HybridNetwork)    
    reticulation = Matrix{Int}(length(net.hybrid), 2) #initialize reticulation matrix 
    #fill reticulation matrix
    j = 1 #row index for reticulation matrix
    for n in net.hybrid
        for e in n.edge #iterate over each edge attatched to a hybrid node
            if e.gamma < 0.5 #find minor hybrid edge
                #push parent node to first column
                reticulation[j,1] = e.node[e.isChild1 ? 2 : 1].number
                #push child node to second column
                reticulation[j,2] = e.node[e.isChild1 ? 1 : 2].number 
                j += 1 #increase index value
            end
        end
    end
    return reticulation
end

"""

    generateMinorReticulationLength(net::HybridNetwork)

Generate vector of minor edge lengths organized in the same order as the `edge` matrix 
created via `generateApeReticulation`.

# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
julia> directEdges!(net)
julia> preorder!(net)
julia> apeNodeNumbers(net)
julia> generateApeReticulationLength(net)
1-element Array{Float64,1}:
 -1.0
```

""" #"

function generateMinorReticulationLength(net::HybridNetwork) 
    edgeLength = Array{Float64}(length(net.hybrid)) #initialize edge length vector   
    #fill edge length vector
    i=1
    for n in net.nodes_changed #traverse tree in topological order excluding leaves
        if !n.leaf
            for e in n.edge #iterate over each edge attatched to a node
                #exclude node that is parent of current edge and major hybrid edge
                if e.node[e.isChild1 ? 1 : 2] != n && e.gamma < 0.5 
                    edgeLength[i] = e.length
                    i=i+1
                end
            end
        end
    end
    return edgeLength
end

"""
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
R> plot($phy)
```

""" #"

function apeRExport(net::HybridNetwork; mainTree::Bool=false, useEdgeLength::Bool=false) 
    directEdges!(net)
    preorder!(net) #organize nodes for preorder traversal 
    if mainTree == true
        net = majorTree(net)
        preorder!(net) #organize nodes for preorder traversal   
    end 
    apeNodeNumbers(net) #number `julia` network like a `phylo` variable 
    ntips = length(net.leaf) #number of leaves for Nnode 
    totalnodes = length(net.node) #total number of nodes for edge matrix and Nnode 
    Nnode = totalnodes - ntips #integer for number of internal nodes
    tipLabel = [node.name for node in net.leaf] #tiplabel vector
    edge = generateMajorEdge(net) #generate edge matrix 
    if useEdgeLength == true
        edgeLength = generateMajorLength(net) #generate vector of major edge lengths 
    elseif useEdgeLength == false
        #fill edge length vector with null lengths
        edgeLength = fill(-1.0, length(net.edge)-length(net.hybrid)) 
    end 
    if net.numHybrids > 0
        reticulation = generateMinorReticulation(net) #generate reticulation matrix 
        if useEdgeLength == true
            #generate vector of minor edge lengths
            reticulationLength = generateMinorReticulationLength(net) 
        elseif useEdgeLength == false
            #fill minor edge length vector with null lengths
            reticulationLength = fill(-1.0, length(net.hybrid)) 
        end 
        #push Nnode integer, edge and reticulation matricies and lengths, tiplabel vector 
        #to RObject $phy and assign "phylo" and "evonet" class attribute and 
        #export as an RObject
        R"""
        edgeLength = $edgeLength
        edgeLength[edgeLength==-1.0]=NA
        reticulationLength = $reticulationLength
        reticulationLength[reticulationLength==-1.0]=NA
        phy = list(Nnode = $Nnode, edge = $edge, tip.label = $tipLabel, edge.length = edgeLength, reticulation = $reticulation, reticulation.length = reticulationLength)
        class(phy) <- c("evonet", "phylo")
        """
        phy = reval("phy") 
        return(phy) 
    elseif net.numHybrids == 0
        #push Nnode integer, major edge matrix and lengths, tiplabel vector to 
        #RObject $phy and assign "phylo" class attribute and export as an RObject
        R"""
        edgeLength = $edgeLength
        edgeLength[edgeLength==-1.0]=NA
        phy = list(Nnode = $Nnode, edge = $edge, tip.label = $tipLabel, edge.length = edgeLength)
        class(phy) = "phylo"
        """
        phy = reval("phy") 
        return(phy)
    end
end

"""

    function sexp(net::HybridNetwork; numHybrid::Int64=0)

Sexp method to export HybridNework objects generated in the 
Julia package `PhyloNetworks.jl` to the R language as either `phylo` or
`evonet` object (depending on degree of hybridization) recognized by the R package `ape`.

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
R> plot($sobj)

""" #"

function sexp(net::HybridNetwork) 
    preorder!(net) #organize nodes for preorder traversal
    directEdges!(net) 
    apeNodeNumbers(net) #number `julia` network like a `phylo` variable
    ntips = length(net.leaf) #number of leaves for Nnode
    totalnodes = length(net.node) #total number of nodes for edge matrix and Nnode
    Nnode = totalnodes - ntips #integer for number of internal nodes
    tipLabel = [node.name for node in net.leaf] #tip label vector 
    edge = generateMajorEdge(net) #generate edge matrix
    #create an object with fields for $Nnode, $tipLabel, and $edge
    phy = Dict{Symbol, Any}()
    phy[:Nnode] = Nnode
    phy[Symbol("tip.label")] = tipLabel 
    phy[:edge] = edge 
    if net.numHybrids == 0
        edgeLength = generateMajorLength(net) #generate vector of major edge lengths
        #have `R` interpret edge lengths of `-1.0` as `NA` 
        edgeLength = [reinterpret(Float64, 0x7ff00000000007a2) for edge in edgeLength if edge == -1.0]
        phy[Symbol("edge.length")] = edgeLength #add $edgeLength field to $phy 
        #assign "phylo" class to $phy
        sobj = protect(sexp(phy))
        setclass!(sobj, sexp("phylo"))
        unprotect(1)
    elseif net.numHybrids > 0
        edgeLength = generateMajorLength(net)
        edgeLength = [reinterpret(Float64, 0x7ff00000000007a2) for edge in edgeLength if edge == -1.0]
        phy[Symbol("edge.length")] = edgeLength
        reticulation = generateMinorReticulation(net) #generate reticulation matrix
        #generate vector of minor edge lengths
        reticulationLength = generateMinorReticulationLength(net)
        #have R interpret edge lengths of `-1.0` as `NA`
        reticulationLength = [reinterpret(Float64, 0x7ff00000000007a2) for reticulation in reticulationLength if reticulation == -1.0]
        phy[:reticulation] = reticulation #add $reticulation field to $phy
        #add $reticulationLength field to $phy
        phy[Symbol("reticulation.length")] = reticulationLength 
        #assign "phylo" and "evonet" classes to $phy
        sobj = protect(sexp(phy))
        setclass!(sobj, sexp(["phylo", "evonet"]))
        unprotect(1)
    end
    return(sobj) #export RObject
end
