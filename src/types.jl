# circularity: a node has a vector of edges, and an edge has a vector of nodes

"""
    ANode

Abstract node. An object of type [`EdgeT`](@ref) has a `node` attribute,
which is an vector of 2 objects of some subtype of `ANode`.
The concrete type [`Node`](@ref) is a subtype of `ANode`,
and has an `edge` attribute, which is vector of [`Edge`](@ref PhyloNetworks.EdgeT) objects
(where Edge is an alias for EdgeT{Node}).
"""
abstract type ANode end

"""
    EdgeT{node type}
    Edge = EdgeT{Node}
    Edge(number, length=1.0)

Data structure for an edge and its various attributes. Most notably:
- `number` (integer): serves as unique identifier;
  remains unchanged when the network is modified,
  with a nearest neighbor interchange for example
- `node`: vector of [`Node`](@ref)s, normally just 2 of them
- `isChild1` (boolean): `true` if `node[1]` is the child node of the edge,
  false if `node[1]` is the parent node of the edge
- `length`: branch length
- `hybrid` (boolean): whether the edge is a tree edge or a hybrid edge
  (in which case `isChild1` is important, even if the network is semi-directed)
- `gamma`: proportion of genetic material inherited by the child node via the edge;
  1.0 for a tree edge
- `isMajor` (boolean): whether the edge is the major path to the child node;
  `true` for tree edges, since a tree edge is the only path to its child node;
  normally true if `gamma>0.5`.

and other fields, used very internally
"""
mutable struct EdgeT{T<:ANode}
    number::Int
    length::Float64 #default 1.0
    hybrid::Bool
    y::Float64 # exp(-t), cannot set in constructor for congruence
    z::Float64 # 1-y , cannot set in constructor for congruence
    gamma::Float64 # set to 1.0 for tree edges, hybrid?gamma:1.0
    node::Vector{T}
    isChild1::Bool # used for hybrid edges to set the direction (default true)
    isMajor::Bool  # major edge treated as tree edge for network traversal
                   # true if gamma>.5, or if it is the original tree edge
    inCycle::Int # = Hybrid node number if this edge is part of a cycle created by such hybrid node
    # -1 if not part of cycle. used to add new hybrid edge. updated after edge is part of a network
    containRoot::Bool # true if this edge can contain a root given the direction of hybrid edges
                      # used to add new hybrid edge. updated after edge is part of a network
    istIdentifiable::Bool # true if the parameter t (length) for this edge is identifiable as part of a network
                          # updated after part of a network
    fromBadDiamondI::Bool # true if the edge came from deleting a bad diamond I hybridization
end
# outer constructors: ensure congruence among (length, y, z) and (gamma, hybrid, isMajor), and size(node)=2
function EdgeT{T}(number::Int, length::Float64=1.0) where {T<:ANode}
    y = exp(-length)
    EdgeT{T}(number,length,false,y,1.0-y,1.,T[],true,true,-1,true,true,false)
end
function EdgeT{T}(number::Int, length::Float64,hybrid::Bool,gamma::Float64,isMajor::Bool=(!hybrid || gamma>0.5)) where {T<:ANode}
    y = exp(-length)
    EdgeT{T}(number,length,hybrid,y,1.0-y, hybrid ? gamma : 1.,T[],true,isMajor,-1,!hybrid,true,false)
end
function EdgeT{T}(number::Int, length::Float64,hybrid::Bool,gamma::Float64,node::Vector{T}) where {T<:ANode}
    size(node,1) == 2 || error("vector of nodes must have exactly 2 values")
    y = exp(-length)
    Edge{T}(number,length,hybrid,y,1.0-y,
        hybrid ? gamma : 1., node,true, !hybrid || gamma>0.5,
        -1,!hybrid,true,false)
end
function EdgeT{T}(number::Int, length::Float64,hybrid::Bool,gamma::Float64,node::Vector{T},isChild1::Bool, inCycle::Int, containRoot::Bool, istIdentifiable::Bool) where {T<:ANode}
    size(node,1) == 2 || error("vector of nodes must have exactly 2 values")
    y = exp(-length)
    Edge{T}(number,length,hybrid,y,1.0-y,
        hybrid ? gamma : 1., node,isChild1, !hybrid || gamma>0.5,
        inCycle,containRoot,istIdentifiable,false)
end

# warning: gammaz, inCycle, isBadTriangle/Diamond updated until the node is part of a network
"""
    Node(number, leaf)
    Node(number, leaf, hybrid)

Data structure for a node and its various attributes. Most notably:

- `number` (integer): serves as unique identifier;
  remains unchanged when the network is modified,
  with a nearest neighbor interchange for example
- `leaf` (boolean): whether the node is a leaf (with data typically) or an
  internal node (no data typically)
- `name` (string): taxon name for leaves; internal node may or may not have a name
- `edge`: vector of [`Edge`](@ref PhyloNetworks.EdgeT)s that the node is attached to;
  1 if the node is a leaf, 2 if the node is the root, 3 otherwise, and
  potentially more if the node has a polytomy
- `hybrid` (boolean): whether the node is a hybrid node (with 2 or more parents)
  or a tree node (with a single parent)

Other more internal attributes include:

- `isBadDiamondI` and `isBadDiamondII` (booleans): whether the node is a
  hybrid node where the reticulation forms a cycle of 4 nodes (diamond),
  and where both parents of the hybrid nodes are connected to a leaf.
  In a bad diamond of type I, the hybrid node itself is also connected
  to a leaf but the common neighbor of the 2 hybrid's parents is not connected
  to a leaf.
  In a bad diamond of type II, the hybrid node has an internal node as child,
  and the common neighbor of the 2 hybrid's parents is connected to a leaf.
- `isBadTriangle`, `isVeryBadTriangle` and `isExtBadTriangle` (booleans):
  true if the reticulation forms a cycle of 3 nodes (triangle) and
  depending on the number of leaves attached these 3 nodes. The triangle means
  that the 2 parents of the hybrid node are directly related:
  one is the child of the other. `isBadTriangle` is true if the triangle is
  "good", as per Solís-Lemus & Ané (2016), that is, if all 3 nodes in the cycle
  are not connected to any leaves (the reticulation is detectable from quartet
  concordance factors, even though all branch lengths are not identifiable).
  `isVeryBadTriangle` is true if 2 (or all) of the 3 nodes are connected to a
  leaf, in which case the reticulation is undetectable from unrooted gene tree
  topologies (thus it's best to exclude these reticulations from a search).
  `isBadTriangle` is true if exactly 1 of the 3 nodes is connected to a leaf.

For details see Solís-Lemus & Ané (2016, doi:10.1371/journal.pgen.1005896)
"""
mutable struct Node <: ANode
    number::Int
    leaf::Bool
    hybrid::Bool
    gammaz::Float64  # notes file for explanation. gammaz if tree node, gamma2z if hybrid node.
                     # updated after node is part of network with updateGammaz!
    edge::Vector{EdgeT{Node}}
    hasHybEdge::Bool #is there a hybrid edge in edge? only needed when hybrid=false (tree node)
    isBadDiamondI::Bool # for hybrid node, is it bad diamond case I, update in updateGammaz!
    isBadDiamondII::Bool # for hybrid node, is it bad diamond case II, update in updateGammaz!
    isExtBadTriangle::Bool # for hybrid node, is it extremely bad triangle, udpate in updateGammaz!
    isVeryBadTriangle::Bool # for hybrid node, is it very bad triangle, udpate in updateGammaz!
    isBadTriangle::Bool # for hybrid node, is it very bad triangle, udpate in updateGammaz!
    inCycle::Int # = hybrid node if this node is part of a cycle created by such hybrid node, -1 if not part of cycle
    prev::Union{Nothing,Node} # previous node in cycle, used in updateInCycle. set to "nothing" to begin with
    k::Int # num nodes in cycle, only stored in hybrid node, updated after node becomes part of network
           # default -1
    typeHyb::Int8 # type of hybridization (1,2,3,4, or 5), needed for quartet network only. default -1
    name::AbstractString
end

const Edge = EdgeT{Node}

Node() = Node(-1,false,false,-1.,Edge[],false,false,false,false,false,false,-1,nothing,-1,-1,"")
Node(number::Int, leaf::Bool, hybrid::Bool=false) = Node(number,leaf,hybrid,-1.,Edge[],hybrid,false,false,false,false,false,-1.,nothing,-1,-1,"")
# set hasHybEdge depending on edge:
Node(number::Int, leaf::Bool, hybrid::Bool, edge::Vector{Edge}) = Node(number,leaf,hybrid,-1.,edge,any(e->e.hybrid,edge),false,false,false,false,false,-1.,nothing,-1,-1,"")
Node(number::Int, leaf::Bool, hybrid::Bool,gammaz::Float64, edge::Vector{Edge}) = Node(number,leaf,hybrid,gammaz,edge,any(e->e.hybrid, edge),false,false,false,false,false,-1.,nothing,-1,-1,"")

# partition type
mutable struct Partition
    cycle::Vector{Int} #hybrid node number for cycle (or cycles)
    edges::Vector{Edge} #edges in partition
end

abstract type Network end

# warning: no attempt to make sure the direction of edges matches with the root
# warning: no check if it is network or tree, node array can have no hybrids
# warning: nodes and edges need to be defined and linked before adding to a network
"""
    HybridNetwork

Subtype of abstract `Network` type.
Explicit network or tree with the following attributes:

- numTaxa (taxa are tips, i.e. nodes attached to a single edge)
- numNodes (total number of nodes: tips and internal nodes)
- numEdges
- numHybrids (number of hybrid nodes)
- edge (array of Edges)
- node (array of Nodes)
- root (index of root in vector 'node'. May be artificial, for printing and traversal purposes only.)
- hybrid (array of Nodes: those are are hybrid nodes)
- leaf (array of Nodes: those that are leaves)
- loglik (score after fitting network to data, i.e. negative log pseudolik for SNaQ)
- isRooted (true or false)
"""
mutable struct HybridNetwork <: Network
    numTaxa::Int  # cannot set in constructor for congruence
    numNodes::Int
    numEdges::Int
    node::Array{Node,1}
    edge::Array{Edge,1}
    root::Int # node[root] is the root node, default 1
    names::Array{String,1} # translate table for taxon names --but also includes hybrid names...
    hybrid::Array{Node,1} # array of hybrid nodes in network
    numHybrids::Int # number of hybrid nodes
    cladewiseorder_nodeIndex::Vector{Int} # index in 'node' for "cladewise" preorder in main tree
    visited::Array{Bool,1} # reusable array of booleans
    edges_changed::Array{Edge,1} # reusable array of edges
    nodes_changed::Array{Node,1} # reusable array of nodes. used for preorder traversal
    leaf::Array{Node,1} # array of leaves
    ht::Vector{Float64} # vector of parameters to optimize
    numht::Vector{Int} # vector of number of the hybrid nodes and edges in ht e.g. [3,6,8,...], 2 hybrid nodes 3,6, and edge 8 is the 1st identifiable
    numBad::Int # number of bad diamond I hybrid nodes, set as 0
    hasVeryBadTriangle::Bool # true if the network has extremely/very bad triangles that should be ignored
    index::Vector{Int} #index in net.edge, net.node of elements in net.ht to make updating easy
    loglik::Float64 # value of the min -loglik after optBL
    blacklist::Vector{Int} # reusable array of integers, used in afterOptBL
    partition::Vector{Partition} # to choose edges from a partition only to avoid intersecting cycles
    cleaned::Bool # attribute to know if the network has been cleaned after readm default false
    isRooted::Bool # to know if network is rooted, e.g. after directEdges! (which updates isChild1 of each edge)
    # inner constructor
    function HybridNetwork(node::Array{Node,1},edge::Array{Edge,1})
        hybrid=Node[];
        leaf=Node[];
        for n in node
            if n.hybrid push!(hybrid,n); end
            if n.leaf   push!(leaf,  n); end
        end
        new(size(leaf,1),size(node,1),size(edge,1),node,edge,1,[],hybrid,size(hybrid,1), #numTaxa,...,numHybrids
            [],[],[],[],leaf,[],[], #cladewiseorder,...,numht
            0,false,[],0,[],[],false,false) #numBad...
    end
    HybridNetwork() = new(0,0,0,[],[],0,[],[],0, # numTaxa ... numHybrid
                          [],[],[],[],[],[],[], # cladewiseorder...
                          0,false,[],0,[],[],false,false); # numBad ...
end


