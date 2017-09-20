using RCall
using RCall: protect, unprotect, RClass
import RCall.sexp

"""

    apeNodeNumbers(tree::HybridNetwork)

Change numbers of internal nodes of `tree` to satisfy conditions assumed by the `ape` R package: leaves are 1-n, the root is n+1, 
and internal nodes are higher consecutive integers.

# Examples

```julia-repl
julia> tree = readTopology("(A,(B,(C,D)));")
julia> directEdges!(tree)
julia> preorder!(tree)
julia> apeNodeNumbers(tree)
julia> tree.node
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

function apeNodeNumbers(tree::HybridNetwork)  
    lnum = 1
    inum = length(tree.leaf) + 1 #ensures internal node numbers do not overlap with leaf numbers
    #root will be labeled ntips +1, consistent with ape "phylo" class
    for n in tree.nodes_changed #traverse tree in topological order excluding leaves
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

    generateEdge(tree::HybridNetwork)

Generate matrix from `tree` where edge[i,1] represents the parent node of edge i and edge[i,2] represents the child node of edge i.
Assumes nodes have been organized for preorder traversal and node numbering satisfies the conditions assumed by the `ape` package in R.

#Examples

```julia-repl
julia> tree = readTopology("(A,(B,(C,D)));")
julia> directEdges!(tree)
julia> preorder!(tree)
julia> apeNodeNumbers(tree)
6×2 Array{Int64,2}:
 5  4
 5  6
 6  3
 6  7
 7  2
 7  1
```

"""

function generateEdge(tree::HybridNetwork) 
    
    edge = Matrix{Int}(length(tree.edge), 2) #initialize edge matrix
    
    #fill edge matrix
    i = 1 #row index for edge matrix
    for n in tree.nodes_changed #traverse tree in topological order excluding leaves
        if !n.leaf
            for e in n.edge #iterate over each edge attatched to a node
                if e.node[e.isChild1 ? 1 : 2] != n #exclude node that is parent of current edge
                    edge[i,1] = e.node[e.isChild1 ? 2 : 1].number #push parent node to first column
                    edge[i,2] = e.node[e.isChild1 ? 1 : 2].number #push child node to second column
                    i += 1 #increase index value
                end
            end
        end
    end
    return edge
end

"""
    generateLength(tree::HybridNetwork)

generate vector of edge lengths of `tree` edges organized in the same order as the `edge` matrix created via `generate edge`

"""

function generateLength(tree::HybridNetwork)
    edgeLength = Array{Float64}(length(tree.edge)) #initialize edge length vector
    
    #fill edge length vector
    i=1
    for n in tree.nodes_changed #traverse tree in topological order excluding leaves
        if !n.leaf
            for e in n.edge #iterate over each edge attatched to a node
                if e.node[e.isChild1 ? 1 : 2] != n #exclude node that is parent of current edge
                    edgeLength[i] = e.length
                    i=i+1
                end
            end
        end
    end
    return edgeLength
end

function generateApeReticulation(net::HybridNetwork) 
    
    reticulation = Matrix{Int}(length(net.hybrid), 2) #initialize reticulation matrix
    
    #fill reticulation matrix
    j = 1 #row index for reticulation matrix
    for n in net.hybrid
        for e in n.edge #iterate over each edge attatched to a hybrid node
            if e.gamma < 0.5 #find minor hybrid edge
                reticulation[j,1] = e.node[e.isChild1 ? 2 : 1].number #push parent node to first column
                reticulation[j,2] = e.node[e.isChild1 ? 1 : 2].number #push child node to second column
                j += 1 #increase index value
            end
        end
    end
    return reticulation
end

"""
    RExport(net::HybridNetwork)

Export `net` to R as a `evonet` object recognized by the `ape` library in R.

#Arguments

-numhybrids: Must be a non-negative integer. If > 0, will will return `evonet` object. If = 0, will return `phylo` object.
-useEdgeLength: if true, export edge lengths from net.
-mainTree: if true, minor hybrid nodes omitted.


#Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
julia> RExport(net)
```

This creates an RObject called `phy` that is essentially an `phylo` class S3 object. The resulting RObject can be evaluated using the tools
available in the `ape` library in R. For example, we can visualize the network using the `plot` function:

```julia-repl
julia> using RCall
R> library(ape)
R> plot(phy)
```

"""

function RExport(net::HybridNetwork; numHybrids::Int64=0, mainTree::Bool=false, useEdgeLength::Bool=false)
    
    directEdges!(net)
    preorder!(net) #organize nodes for preorder traversal
    if numHybrids > 0

        if mainTree == true
            net = majorTree(net)
        end
    
        apeNodeNumbers(net)
        
        ntips = length(net.leaf) #number of leaves for Nnode
        
        totalnodes = length(net.node) #total number of nodes for edge matrix and Nnode
        
        Nnode = totalnodes - ntips #integer for number of internal nodes

        tipLabel = [node.name for node in net.leaf] #tiplabel vector

        edge = generateEdge(net)
        
        reticulation = generateApeReticulation(net)
    
        #push Nnode integer, edge and reticulation matricies, tiplabel vector to RObject $phy and assign "phylo" and "evonet" class attribute and export as an RObject
        R"""
        phy = list(Nnode = $Nnode, edge = $edge, tip.label = $tipLabel, reticulation = $reticulation)
        class(phy) <- c("evonet", "phylo")
        """
        phy = reval("phy")
    
        return(phy)
        
    elseif numHybrids == 0
    
        preorder!(net) #organize nodes for preorder traversal
    
        apeNodeNumbers(net)
        ntips = length(net.leaf) #number of leaves for Nnode
        totalnodes = length(net.node) #total number of nodes for edge matrix and Nnode
        Nnode = totalnodes - ntips #integer for number of internal nodes

        tipLabel = [node.name for node in net.leaf] #tiplabel vector

        edge = generateEdge(net)
    
        if useEdgeLength == true
            edgeLength = generateLength(net)
        elseif useEdgeLength == false
            edgeLength = fill(-1.0, length(net.edge)) 
        end

        #push Nnode integer, edge matrix, tiplabel vector to RObject $phy and assign "phylo" class attribute and export as an RObject
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

Modifies existing sexp function to work with HybridNework objects generated in the Julia package PhyloNetworks.jl

#Arguments

-numhybrid: Must be a non-negative integer. If > 0, will will return `evonet` object. If = 0, will return `phylo` object.

"""


function sexp(net::HybridNetwork; numHybrid::Int64=0)
    
    preorder!(net) #organize nodes for preorder traversal
    directEdges!(net)
    
    apeNodeNumbers(net)
    ntips = length(net.leaf) #number of leaves for Nnode
    totalnodes = length(net.node) #total number of nodes for edge matrix and Nnode
    Nnode = totalnodes - ntips #integer for number of internal nodes

    tipLabel = [node.name for node in net.leaf] #tiplabel vector
     
    edge = generateEdge(net)
    
    phy = Dict{Symbol, Any}()
    phy[:Nnode] = Nnode
    phy[Symbol("tip.label")] = tipLabel
    
    phy[:edge] = edge
    
    if numHybrid == 0
        edgeLength = generateLength(net)
        phy[Symbol("edge.length")] = edgeLength
        sobj = protect(sexp(phy))
        setclass!(sobj, sexp("phylo"))
        unprotect(1)
    elseif numHybrid > 0
        reticulation = generateApeReticulation(net)
        phy[:reticulation] = reticulation
        sobj = protect(sexp(phy))
        setclass!(sobj, sexp("evonet"))
        unprotect(1)
    end
    return(sobj)
end

