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
- `ischild1` (boolean): `true` if `node[1]` is the child node of the edge,
  false if `node[1]` is the parent node of the edge
- `length`: branch length
- `hybrid` (boolean): whether the edge is a tree edge or a hybrid edge
  (in which case `ischild1` is important, even if the network is semi-directed)
- `gamma`: proportion of genetic material inherited by the child node via the edge;
  1.0 for a tree edge
- `ismajor` (boolean): whether the edge is the major path to the child node;
  `true` for tree edges, since a tree edge is the only path to its child node;
  normally true if `gamma>0.5`.

and other fields, used very internally
"""
mutable struct EdgeT{T<:ANode}
    number::Int
    length::Float64 # default 1.0
    hybrid::Bool
    y::Float64 # exp(-t), cannot set in constructor for congruence
    z::Float64 # 1-y , cannot set in constructor for congruence
    gamma::Float64 # set to 1.0 for tree edges, hybrid?gamma:1.0
    node::Vector{T}
    ischild1::Bool # used for hybrid edges to set the direction (default true)
    ismajor::Bool  # major edge treated as tree edge for network traversal
                   # true if gamma>.5, or if it is the original tree edge
    inte1::Int     # sometimes used for cycle number if level-1
    containroot::Bool # true if this edge may contain the root given the direction of hybrid edges
    boole1::Bool # reusable booleans, for an edge
    boole2::Bool
end
# outer constructors: ensure congruence among (length, y, z) and (gamma, hybrid, ismajor), and size(node)=2
function EdgeT{T}(number::Int, length::Float64=1.0) where {T<:ANode}
    y = exp(-length)
    EdgeT{T}(number,length,false,y,1.0-y,1.,T[],true,true,-1,true,true,false)
end
function EdgeT{T}(
  number::Int,
  length::Float64,
  hybrid::Bool,
  gamma::Float64,
  ismajor::Bool=(!hybrid || gamma>0.5)
) where {T<:ANode}
    y = exp(-length)
    EdgeT{T}(number,length,hybrid,y,1.0-y, hybrid ? gamma : 1.,T[],true,ismajor,-1,!hybrid,true,false)
end
function EdgeT{T}(number::Int, length::Float64,hybrid::Bool,gamma::Float64,node::Vector{T}) where {T<:ANode}
    size(node,1) == 2 || error("vector of nodes must have exactly 2 values")
    y = exp(-length)
    Edge{T}(number,length,hybrid,y,1.0-y,
        hybrid ? gamma : 1., node,true, !hybrid || gamma>0.5,
        -1,!hybrid,true,false)
end
function EdgeT{T}(
  number::Int,
  length::Float64,
  hybrid::Bool,
  gamma::Float64,
  node::Vector{T},
  ischild1::Bool,
  inte1::Int,
  containroot::Bool,
  boole1::Bool
) where {T<:ANode}
    size(node,1) == 2 || error("vector of nodes must have exactly 2 values")
    y = exp(-length)
    Edge{T}(number,length,hybrid,y,1.0-y,
        hybrid ? gamma : 1., node, ischild1, !hybrid || gamma>0.5,
        inte1, containroot, boole1, false)
end

# warning: fvalue, intn1, booln6/booln2 updated until the node is part of a network
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

- `isBadDiamondI` -> `booln2` and `isBadDiamondII` -> `booln3` (booleans): whether the node is a
  hybrid node where the reticulation forms a cycle of 4 nodes (diamond),
  and where both parents of the hybrid nodes are connected to a leaf.
  In a bad diamond of type I, the hybrid node itself is also connected
  to a leaf but the common neighbor of the 2 hybrid's parents is not connected
  to a leaf.
  In a bad diamond of type II, the hybrid node has an internal node as child,
  and the common neighbor of the 2 hybrid's parents is connected to a leaf.
- `isBadTriangle` -> `booln6`, `isVeryBadTriangle` -> `booln5` and `isExtBadTriangle` -> `booln4` (booleans):
  true if the reticulation forms a cycle of 3 nodes (triangle) and
  depending on the number of leaves attached these 3 nodes. The triangle means
  that the 2 parents of the hybrid node are directly related:
  one is the child of the other. `isBadTriangle` is true if the triangle is
  "good", as per Solís-Lemus & Ané (2016), that is, if all 3 nodes in the cycle
  are not connected to any leaves (the reticulation is detectable from quartet
  concordance factors, even though all branch lengths are not identifiable).
  `isVeryBadTriangle` -> `booln5` is true if 2 (or all) of the 3 nodes are connected to a
  leaf, in which case the reticulation is undetectable from unrooted gene tree
  topologies (thus it's best to exclude these reticulations from a search).
  `isBadTriangle` is true if exactly 1 of the 3 nodes is connected to a leaf.

For details see Solís-Lemus & Ané (2016, doi:10.1371/journal.pgen.1005896)
"""
mutable struct Node <: ANode
    number::Int
    leaf::Bool
    hybrid::Bool
    fvalue::Float64 # -1.0 by default. re-usable. used for gammaz or gamma2z in SNaQ
    edge::Vector{EdgeT{Node}}
    booln1::Bool # default: is incident to a hybrid edge?
    booln2::Bool # reusable booleans, for nodes. false by default
    booln3::Bool
    booln4::Bool
    booln5::Bool
    booln6::Bool
    intn1::Int   # default -1. sometimes used for cycle number if level-1
    prev::Union{Nothing,Node} # for traversal algorithms that track previous node
    intn2::Int   # default -1
    int8n3::Int8 # default -1
    name::AbstractString # default ""
end

const Edge = EdgeT{Node}

Node() = Node(-1,false,false,-1.,Edge[],false,false,false,false,false,false,-1,nothing,-1,-1,"")
Node(number::Int, leaf::Bool, hybrid::Bool=false) = Node(number,leaf,hybrid,-1.,Edge[],hybrid,false,false,false,false,false,-1.,nothing,-1,-1,"")
# set booln1 depending on edge:
Node(number::Int, leaf::Bool, hybrid::Bool, edge::Vector{Edge}) = Node(number,leaf,hybrid,-1.,edge,any(e->e.hybrid,edge),false,false,false,false,false,-1.,nothing,-1,-1,"")
Node(number::Int, leaf::Bool, hybrid::Bool,fvalue::Float64, edge::Vector{Edge}) = Node(number,leaf,hybrid,fvalue,edge,any(e->e.hybrid, edge),false,false,false,false,false,-1.,nothing,-1,-1,"")

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

- numtaxa (taxa are tips, i.e. nodes attached to a single edge)
- numnodes (total number of nodes: tips and internal nodes)
- numedges
- numhybrids (number of hybrid nodes)
- edge (array of Edges)
- node (array of Nodes)
- rooti (index of root in vector 'node'. May be artificial, for printing and traversal purposes only.)
- hybrid (array of Nodes: those are are hybrid nodes)
- leaf (array of Nodes: those that are leaves)
- fscore (score after fitting network to data, i.e. negative log pseudolik for SNaQ)
- isrooted (true or false)
"""
mutable struct HybridNetwork <: Network
    numtaxa::Int  # cannot set in constructor for congruence
    numnodes::Int
    numedges::Int
    node::Array{Node,1}
    edge::Array{Edge,1}
    rooti::Int # node[rooti] is the root node, default 1
    names::Array{String,1} # translate table for taxon names --but also includes hybrid names...
    hybrid::Array{Node,1} # array of hybrid nodes in network
    numhybrids::Int # number of hybrid nodes
    vec_int1::Vector{Int}  # vector of integers, e.g. to get "cladewise" preorder in main tree
    vec_bool::Array{Bool,1}
    vec_edge::Array{Edge,1}
    vec_node::Array{Node,1} # reusable, but used for preorder traversal
    leaf::Array{Node,1} # array of leaves
    vec_float::Vector{Float64}
    vec_int2::Vector{Int}
    intg1::Int
    boolg1::Bool
    vec_int3::Vector{Int}
    fscore::Float64 # in SNaQ: min -loglik that was found
    vec_int4::Vector{Int}
    partition::Vector{Partition} # to choose edges from a partition only to avoid intersecting cycles
    boolg2::Bool # default false
    isrooted::Bool # is network rooted? otherwise: semidirected fixit: check how it's used
    # inner constructor
    function HybridNetwork(node::Array{Node,1},edge::Array{Edge,1})
        hybrid=Node[];
        leaf=Node[];
        for n in node
            if n.hybrid push!(hybrid,n); end
            if n.leaf   push!(leaf,  n); end
        end
        new(size(leaf,1),size(node,1),size(edge,1),node,edge,1,[],hybrid,size(hybrid,1), #numtaxa,...,numhybrids
            [],[],[],[],leaf,[],[], #cladewiseorder,...,vec_int2
            0,false,[],0,[],[],false,false) #intg1...
    end
    HybridNetwork() = new(0,0,0,[],[],0,[],[],0, # numtaxa ... numHybrid
                          [],[],[],[],[],[],[], # cladewiseorder...
                          0,false,[],0,[],[],false,false); # intg1 ...
end

struct RootMismatch <: Exception
  msg::String
end
RootMismatch() = RootMismatch("")
Base.showerror(io::IO, e::RootMismatch) = print(io, "RootMismatch: ", e.msg);

