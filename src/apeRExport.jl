"""
    apeNodeNumbers!(net::HybridNetwork)

Change internal node numbers of `net` to satisfy the conditions assumed by the `ape` R 
package: leaves are 1 to n, the root is n+1, and internal nodes are higher consecutive 
integers. Assume `nodes_changed` was updated, to list nodes in pre-order.

# Examples

```julia-repl
julia> net = readTopology("(A,(B,(C,D)));");
julia> directEdges!(net); preorder!(net)
julia> PhyloNetworks.apeNodeNumbers!(net)
julia> printNodes(net)
Node    In Cycle        isHybrid        hasHybEdge      Node label      isLeaf  Edges numbers
4       -1              false           false           A               true    1
3       -1              false           false           B               true    2
2       -1              false           false           C               true    3
1       -1              false           false           D               true    4
7       -1              false           false                           false   3       4       5
6       -1              false           false                           false   2       5       6
5       -1              false           false                           false   1       6
```
"""

function apeNodeNumbers!(net::HybridNetwork)
    lnum = 1 # first number for leaves
    inum = length(net.leaf) + 1 # first number of internal nodes
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

Generate matrix of major edges from `net` where edge[i,1] is the number of the
parent node of edge i and edge[i,2] is the number of the child node of edge i.
Assume `nodes_changed` was updated, to list nodes in pre-order.

# Examples

```julia-repl
julia> net = readTopology("(A,(B,(C,D)));");
julia> directEdges!(net); preorder!(net)
julia> PhyloNetworks.apeNodeNumbers!(net)
julia> PhyloNetworks.generateMajorEdge(net)
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
                i += 1
            end
        end
    end
    return edge
end

"""
    generateMajorLength(net::HybridNetwork)

Generate vector of edge lengths of major `net` edges organized in the same order
as the `edge` matrix created via `generateMajorEdge`. Replace values of `-1.0` with
null values recognized by the `ape` library for the `R` programming language.
Output is a NullableArray.
Assume `nodes_changed` was updated, to list nodes in pre-order.

# Examples

```julia-repl
julia> net = readTopology("(((A:3.1,(B:0.2)#H1:0.3::0.9),(C,#H1:0.3::0.1):1.1),D:0.7);");
julia> directEdges!(net); preorder!(net)
julia> PhyloNetworks.generateMajorLength(net)
8-element NullableArrays.NullableArray{Float64,1}:
 #NULL
 0.7
 #NULL
 1.1
 #NULL
 3.1
 0.3
 0.2
```
""" #"

function generateMajorLength(net::HybridNetwork) 
    edgeLength = Array{Float64}(length(net.edge)-length(net.hybrid)) 
    i=1
    for n in net.nodes_changed # topological pre-order
        if !n.leaf
            for e in n.edge # for parent major edge below
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
the number of the parent node of edge i and edge[i,2] represents the number
of the child node of edge i. (node numbers may be negative, unless they were
modified by `apeNodeNumbers!`).

# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);");
julia> PhyloNetworks.generateMinorReticulation(net)
1×2 Array{Int64,2}:
 -6  3
```
""" #"

function generateMinorReticulation(net::HybridNetwork)    
    reticulation = Matrix{Int}(length(net.hybrid), 2) # initialize
    j = 1 # row index, row = reticulate edge
    for e in net.edge
        if !e.isMajor # minor (hybrid) edges only
            reticulation[j,1] = e.node[e.isChild1 ? 2 : 1].number # parent
            reticulation[j,2] = e.node[e.isChild1 ? 1 : 2].number # child
            j += 1
        end
    end
    return reticulation
end

"""
    generateMinorReticulationLength(net::HybridNetwork)

Generate vector of minor edge lengths organized in the same order as the 
`reticulation` matrix created via `generateMinorReticulation`. 
Replace values of `-1.0` with null values recognized by the `ape` library for 
the `R` programming language. Output is a NullableArray.

# Examples

```julia-repl
julia> net = readTopology("(((A:3.1,(B:0.2)#H1:0.3::0.9),(C,#H1:0.3::0.1):1.1),D:0.7);");
julia> PhyloNetworks.generateMinorReticulationLength(net)
1-element NullableArrays.NullableArray{Float64,1}:
 0.3
```
""" #"

function generateMinorReticulationLength(net::HybridNetwork) 
    reticulationLength = Vector{Float64}(0) # initialize
    for e in net.edge
        if !e.isMajor #find minor hybrid edge
            push!(reticulationLength, e.length)
        end
    end
    reticulationLength = NullableArray(reticulationLength, map(x -> x==-1.0, reticulationLength))
    return reticulationLength
end

"""
    generateMinorReticulationGamma(net::HybridNetwork)

Generate vector of minor edge gammas (inheritance probabilities) organized in the
same order as the `reticulation` matrix created via `generateMinorReticulation`.
Replace values of `-1.0` with null to be recognized as NA in `R` programming.
Output is a NullableArray.

# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);");
julia> PhyloNetworks.generateMinorReticulationGamma(net)
1-element Array{Float64,1}:
 0.1
 ```
 """ #"

function generateMinorReticulationGamma(net::HybridNetwork)    
    reticulationGamma = Vector{Float64}(0) #initialize  
    for e in net.edge
        if !e.isMajor # minor hybrid edges only
            push!(reticulationGamma, e.gamma)
        end
    end
    reticulationGamma = NullableArray(reticulationGamma, map(x -> x==-1.0, reticulationGamma))
    return reticulationGamma
end

doc"""
    apeRExport(net::HybridNetwork; mainTree=false, useEdgeLength=true)

Create an RObject of class `phylo` (and `evonet` depending on the number
of hybridizations) recognized by the `ape` library in R (S3 object). This
RObject can be evaluated using the tools available in the `ape` library in R.
For example, we can visualize the network using `ape`'s `plot` function.

# Arguments

- useEdgeLength: if true, export edge lengths from `net`.
- mainTree: if true, minor hybrid edges are omitted, but minor hybrid nodes

# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);");
julia> phy = apeRExport(net)
RCall.RObject{RCall.VecSxp}
$Nnode
[1] 5

$edge
     [,1] [,2]
[1,]    5    6
[2,]    5    1
[3,]    6    8
[4,]    6    7
[5,]    7    2
[6,]    8    4
[7,]    8    9
[8,]    9    3

$tip.label
[1] "D" "C" "B" "A"

$reticulation
     [,1] [,2]
[1,]    7    9

$reticulation.gamma
[1] 0.1

attr(,"class")
[1] "evonet" "phylo"

julia> using RCall

julia> R"library(ape)"

julia> phy
RCall.RObject{RCall.VecSxp}

    Evolutionary network with 1 reticulation

               --- Base tree ---
Phylogenetic tree with 4 tips and 5 internal nodes.

Tip labels:
[1] "D" "C" "B" "A"

Rooted; no branch lengths.

R> phy

Evolutionary network with 1 reticulation

               --- Base tree ---
Phylogenetic tree with 4 tips and 5 internal nodes.

Tip labels:
[1] "D" "C" "B" "A"

Rooted; no branch lengths.

R> str(phy)
List of 5
$ Nnode             : int 5
$ edge              : int [1:8, 1:2] 5 5 6 6 7 8 8 9 6 1 ...
$ tip.label         : chr [1:4] "D" "C" "B" "A"
$ reticulation      : int [1, 1:2] 7 9
$ reticulation.gamma: num 0.1
- attr(*, "class")= chr [1:2] "evonet" "phylo"
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
        if sum([!isnull(e) for e in edgeLength])>0
            R"""
            phy[['edge.length']] = $edgeLength
            """
        end
    end
    if net.numHybrids > 0
        reticulation = generateMinorReticulation(net)
        reticulationGamma = generateMinorReticulationGamma(net)
        R"""
        phy[['reticulation']] = $reticulation
        class(phy) <- c("evonet", "phylo")
        """
        if sum([!isnull(e) for e in reticulationGamma])>0
            R"phy[['reticulation.gamma']] = $reticulationGamma"
        end
        if useEdgeLength # extract minor edge lengths
            reticulationLength = generateMinorReticulationLength(net)
            if sum([!isnull(e) for e in reticulationLength])>0
                R"""
                phy[['reticulation.length']] = $reticulationLength
                """
            end
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

Inspired from [Phylo.jl](https://github.com/richardreeve/Phylo.jl/blob/master/src/rcall.jl)

# Examples

```julia-repl
julia> net = readTopology("(((A:.2,(B:.1)#H1:.1::0.9):.1,(C:.11,#H1:.01::0.1):.19):.1,D:.4);");
julia> using RCall
R> library(ape)
R> $net

Evolutionary network with 1 reticulation

               --- Base tree ---
Phylogenetic tree with 4 tips and 5 internal nodes.

Tip labels:
[1] "D" "C" "B" "A"

Rooted; includes branch lengths.
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
    edge = generateMajorEdge(net)
    phy = Dict{Symbol, Any}() # dictionary exported with regular sexp at the end
    phy[:Nnode] = Nnode
    phy[Symbol("tip.label")] = tipLabel
    phy[:edge] = edge
    edgeLength = generateMajorLength(net)
    nBL = sum([!isnull(e) for e in edgeLength]) # number of non-missing branch lengths
    if nBL>0
        phy[Symbol("edge.length")] = edgeLength
    end
    if net.numHybrids > 0
        reticulation = generateMinorReticulation(net)
        reticulationGamma = generateMinorReticulationGamma(net)
        reticulationLength = generateMinorReticulationLength(net)
        phy[:reticulation] = reticulation
        if sum([!isnull(e) for e in reticulationGamma]) > 0
            phy[Symbol("reticulation.gamma")] = reticulationGamma
        end
        if sum([!isnull(e) for e in reticulationLength]) > 0
            phy[Symbol("reticulation.length")] = reticulationLength
        end
    end
    sobj = protect(sexp(phy)) # RObject
    if net.numHybrids == 0
        setclass!(sobj, sexp("phylo"))
    else
        setclass!(sobj, sexp(["evonet", "phylo"]))
    end
    unprotect(1)
    return(sobj)
end
