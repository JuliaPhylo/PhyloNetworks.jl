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

# type created from a HybridNetwork only to extract a given quartet
"""
    QuartetNetwork(net::HybridNetwork)

Subtype of `Network` abstract type.
A `QuartetNetwork` object is an internal type used to calculate the
expected CFs of quartets on a given network.
Attributes of the `QuartetNetwork` objects need not be updated at a given time (see below).

The procedure to calculate expected CFs for a given network is as follows:
1. A `QuartetNetwork` object is created for each `Quartet` using
   `extractQuartet!(net,d)` for `net::HybridNetwork` and `d::DataCF`
2. The vector `d.quartet` has all the `Quartet` objects, each with a `QuartetNetwork`
   object (`q.qnet`). Attibutes in `QuartetNetwork` are not updated at this point
3. Attributes in `QuartetNetwork` are partially updated when calculating the
   expected CF (`calculateExpCFAll!`). To calculate the expected CF for this quartet,
   we need to update the attributes: `which`, `typeHyb`, `t1`, `split`, `formula`, `expCF`.
   To do this, we need to modify the `QuartetNetwork` object (i.e. merge edges,...).
   But we do not want to modify it directly because it is connected to the original
   `net` via a map of the edges and nodes, so we use a deep copy:
   `qnet=deepcopy(q.qnet)` and then `calculateExpCFAll!(qnet)`.
   Attributes that are updated on the original `QuartetNetwork` object `q.qnet` are:
    - `q.qnet.hasEdge`: array of booleans of length equal to `net.edge` that shows which identifiable edges and gammas of `net` (`net.ht`) are in `qnet` (and still identifiable). Note that the first elements of the vector correspond to the gammas.
    - `q.qnet.index`: length should match the number of trues in `qnet.hasEdge`. It has the indexes in `qnet.edge` from the edges in `qnet.hasEdge`. Note that the first elements of the vector correspond to the gammas.
    - `q.qnet.edge`: list of edges in `QuartetNetwork`. Note that external edges in `net` are collapsed when they appear in `QuartetNetwork`, so only internal edges map directly to edges in `net`
    - `q.qnet.expCF`: expected CF for this `Quartet`


Why not modify the original `QuartetNetwork`? We wanted to keep the original
`QuartetNetwork` stored in `DataCF` with all the identifiable edges, to be able
to determine if this object had been changed or not after a certain optimization.

The process is:

1. Deep copy of full network to create `q.qnet` for `Quartet q`.
   This `QuartetNetwork` object has only 4 leaves now, but does not have merged edges
   (the identifiable ones) so that we can correspond to the edges in net.
   This `QuartetNetwork` does not have other attributes updated.
2. For the current set of branch lengths and gammas, we can update the attributes
   in `q.qnet` to compute the expected CF. The functions that do this will "destroy"
   the `QuartetNetwork` object by merging edges, removing nodes, etc... So, we do
   this process in `qnet=deepcopy(q.qnet)`, and at the end, only update `q.qnet.expCF`.
3. After we optimize branch lengths in the full network, we want to update the
   branch lengths in `q.qnet`. The edges need to be there (which is why we do
   not want to modify this `QuartetNetwork` object by merging edges), and
   we do not do a deep-copy of the full network again. We only change the values
   of branch lengths and gammas in `q.qnet`, and we can re-calculate the expCF
   by creating a deep copy `qnet=deepcopy(q.qnet)` and run the other functions
   (which merge edges, etc) to get the `expCF`.

Future work: there are definitely more efficient ways to do this (without the deep copies).
In addition, currently edges that are no longer identifiable in `QuartetNetwork`
do not appear in `hasEdge` nor `index`. Need to study this.

```jldoctest
julia> net0 = readTopology("(s17:13.76,(((s3:10.98,(s4:8.99,s5:8.99)I1:1.99)I2:0.47,(((s6:2.31,s7:2.31)I3:4.02,(s8:4.97,#H24:0.0::0.279)I4:1.36)I5:3.64,((s9:8.29,((s10:2.37,s11:2.37)I6:3.02,(s12:2.67,s13:2.67)I7:2.72)I8:2.89)I9:0.21,((s14:2.83,(s15:1.06,s16:1.06)I10:1.78)I11:2.14)#H24:3.52::0.72)I12:1.47)I13:1.48)I14:1.26,(((s18:5.46,s19:5.46)I15:0.59,(s20:4.72,(s21:2.40,s22:2.40)I16:2.32)I17:1.32)I18:2.68,(s23:8.56,(s1:4.64,s2:4.64)I19:3.92)I20:0.16)I21:3.98)I22:1.05);");

julia> net = readTopologyLevel1(writeTopology(net0)) ## need level1 attributes for functions below
HybridNetwork, Un-rooted Network
46 edges
46 nodes: 23 tips, 1 hybrid nodes, 22 internal tree nodes.
tip labels: s17, s3, s4, s5, ...
(s4:8.99,s5:8.99,(s3:10.0,((((s6:2.31,s7:2.31)I3:4.02,(s8:4.97,#H24:0.0::0.279)I4:1.36)I5:3.64,((s9:8.29,((s10:2.37,s11:2.37)I6:3.02,(s12:2.67,s13:2.67)I7:2.72)I8:2.89)I9:0.21,((s14:2.83,(s15:1.06,s16:1.06)I10:1.78)I11:2.14)#H24:3.52::0.721)I12:1.47)I13:1.48,((((s18:5.46,s19:5.46)I15:0.59,(s20:4.72,(s21:2.4,s22:2.4)I16:2.32)I17:1.32)I18:2.68,(s23:8.56,(s1:4.64,s2:4.64)I19:3.92)I20:0.16)I21:3.98,s17:10.0)I22:1.26)I14:0.47)I2:1.99)I1;

julia> q1 = Quartet(1,["s1", "s16", "s18", "s23"],[0.296,0.306,0.398])
number: 1
taxon names: ["s1", "s16", "s18", "s23"]
observed CF: [0.296, 0.306, 0.398]
pseudo-deviance under last used network: 0.0 (meaningless before estimation)
expected CF under last used network: Float64[] (meaningless before estimation)

julia> qnet = PhyloNetworks.extractQuartet!(net,q1)
taxa: ["s1", "s16", "s18", "s23"]
number of hybrid nodes: 1

julia> sum([e.istIdentifiable for e in net.edge]) ## 23 identifiable edges in net
23

julia> idedges = [ee.number for ee in net.edge[[e.istIdentifiable for e in net.edge]]];

julia> print(idedges)
[5, 6, 9, 11, 12, 13, 17, 20, 21, 22, 26, 27, 28, 29, 30, 31, 34, 38, 39, 40, 44, 45, 46]

julia> length(qnet.hasEdge) ## 24 = 1 gamma + 23 identifiable edges
24

julia> sum(qnet.hasEdge) ## 8 = 1 gamma + 7 identifiable edges in qnet
8

julia> print(idedges[qnet.hasEdge[2:end]]) ## 7 id. edges: [12, 13, 29, 30, 31, 45, 46]
[12, 13, 29, 30, 31, 45, 46]

julia> qnet.edge[qnet.index[1]].number ## 11 = minor hybrid edge
11
```
"""
mutable struct QuartetNetwork <: Network
    numTaxa::Int
    numNodes::Int
    numEdges::Int
    node::Array{Node,1}
    edge::Array{Edge,1}
    hybrid::Array{Node,1} # array of hybrid nodes in network
    leaf::Array{Node,1} # array of leaves
    numHybrids::Int # number of hybrid nodes
    hasEdge::Array{Bool,1} # array of boolean with all the original identifiable edges of HybridNetwork and gammas (net.ht)
    quartetTaxon::Array{String,1} # the quartet taxa in the order it represents. Points to same array as its Quartet.taxon
    which::Int8 # 0 it tree quartet, 1 is equivalent to tree quartet and 2 if two minor CF different, default -1
    typeHyb::Array{Int8,1} #array with the type of hybridization of each hybrid node in the quartet
    t1::Float64 # length of internal edge, used when qnet.which=1, default = -1
    names::Array{String,1} # taxon and node names, same order as in network.node
    split::Array{Int8,1} # split that denotes to which side each leaf is from the split, i.e. [1,2,2,1] means that leaf1 and 4 are on the same side of the split, default -1,-1,-1,-1
    formula::Array{Int8,1} # array for qnet.which=1 that indicates if the expCf is major (1) or minor (2) at qnet.expCF[i] depending on qnet.formula[i], default -1,-1,-1
    expCF::Array{Float64,1} # three expected CF in order 12|34, 13|24, 14|23 (matching obsCF from qnet.quartet), default [0,0,0]
    indexht::Vector{Int} # index in net.ht for each edge in qnet.ht
    changed::Bool # true if the expCF would be changed with the current parameters in the optimization, to recalculate, default true
    index::Vector{Int} # index in qnet.edge (qnet.node for gammaz) of the members in qnet.indexht to know how to find quickly in qnet
    # inner constructor
    function QuartetNetwork(net::HybridNetwork)
        net2 = deepcopy(net); #fixit: maybe we dont need deepcopy of all, maybe only arrays
        new(net2.numTaxa,net2.numNodes,net2.numEdges,net2.node,net2.edge,net2.hybrid,net2.leaf,net2.numHybrids, [true for e in net2.edge],[],-1,[], -1.,net2.names,Int8[-1,-1,-1,-1],Int8[-1,-1,-1],[0,0,0],[],true,[])
    end
    QuartetNetwork() = new(0,0,0,[],[],[],[],0,[],[],-1,[],-1.0,[],[],[],[],[],true,[])
end

abstract type AQuartet end

"""
    Quartet

type that saves the information on a given 4-taxon subset. It contains the following attributes:

- number: integer
- taxon: vector of taxon names, like t1 t2 t3 t4
- obsCF: vector of observed CF, in order 12|34, 13|24, 14|23
- logPseudoLik
- ngenes: number of gene trees used to compute the observed CF; -1.0 if unknown
- qnet: [`QuartetNetwork`](@ref), which saves the expCF after snaq estimation to
  emphasize that the expCF depend on a specific network, not the data

see also: [`QuartetT`](@ref) for quartet with data of user-defined type `T`,
using a mapping between quartet indices and quartet taxa.
"""
mutable struct Quartet <: AQuartet
    number::Int
    taxon::Array{String,1} # taxa 1234. qnet.quartetTaxon points to the same array.
    obsCF::Array{Float64,1} # three observed CF in order 12|34, 13|24, 14|23
    qnet::QuartetNetwork # quartet network for the current network (want to keep as if private attribute)
    logPseudoLik::Float64 # log pseudolik value for the quartet. 0.0 by default
    ngenes::Float64 # number of gene trees used to compute the obsCV, default -1.; Float in case ngenes is average
    # inner constructor: to guarantee obsCF are only three and add up to 1
    function Quartet(number::Integer,t1::AbstractString,t2::AbstractString,t3::AbstractString,t4::AbstractString,obsCF::Array{Float64,1})
        size(obsCF,1) != 3 ? error("observed CF vector should have size 3, not $(size(obsCF,1))") : nothing
        0.99 < sum(obsCF) < 1.02 || @warn "observed CF should add up to 1, not $(sum(obsCF))"
        new(number,[t1,t2,t3,t4],obsCF,QuartetNetwork(),0.0,-1.0);
    end
    function Quartet(number::Integer,t1::Array{String,1},obsCF::Array{Float64,1})
        size(obsCF,1) != 3 ? error("observed CF vector should have size 3, not $(size(obsCF,1))") : nothing
        0.99< sum(obsCF) < 1.02 || @warn "observed CF should add up to 1, not $(sum(obsCF))"
        size(t1,1) != 4 ? error("array of taxa should have size 4, not $(size(t1,1))") : nothing
        0.0 <= obsCF[1] <= 1.0 || error("obsCF must be between (0,1), but it is $(obsCF[1]) for $(t1)")
        0.0 <= obsCF[2] <= 1.0 || error("obsCF must be between (0,1), but it is $(obsCF[2]) for $(t1)")
        0.0 <= obsCF[3] <= 1.0 || error("obsCF must be between (0,1), but it is $(obsCF[3]) for $(t1)")
        new(number,t1,obsCF,QuartetNetwork(),0.0,-1.0);
    end
    Quartet() = new(0,[],[],QuartetNetwork(),0.0,-1.0)
end

"""
QuartetT{T}

Generic type for 4-taxon sets. Fields:
- `number`: rank of the 4-taxon set
- `taxonnumber`: static vector of 4 integers, assumed to be distinct and sorted
- `data`: object of type `T`

For easier look-up, a unique mapping is used between the rank (`number`) of a
4-taxon set and its 4 taxa (see [`quartetrank`](@ref) and [`nchoose1234`](@ref)):

rank-1 = (t1-1) choose 1 + (t2-1) choose 2 + (t3-1) choose 3 + (t4-1) choose 4

# examples

```jldoctest
julia> nCk = PhyloNetworks.nchoose1234(5)
6×4 Matrix{Int64}:
 0   0   0  0
 1   0   0  0
 2   1   0  0
 3   3   1  0
 4   6   4  1
 5  10  10  5

julia> PhyloNetworks.QuartetT(1,3,4,6, [.92,.04,.04, 100], nCk)
4-taxon set number 8; taxon numbers: 1,3,4,6
data: [0.92, 0.04, 0.04, 100.0]
```
"""
struct QuartetT{T} <: AQuartet where T
    number::Int
    taxonnumber::StaticArrays.SVector{4,Int}
    data::T
end
function Base.show(io::IO, obj::QuartetT{T}) where T
    disp = "4-taxon set number $(obj.number); taxon numbers: "
    disp *= join(obj.taxonnumber,",")
    disp *= "\ndata: "
    print(io, disp)
    print(io, obj.data)
end
function QuartetT(tn1::Int,tn2::Int,tn3::Int,tn4::Int, data::T, nCk::Matrix, checksorted=true::Bool) where T
    if checksorted
        (tn1<tn2 && tn2<tn3 && tn3<tn4) || error("taxon numbers must be sorted")
    end
    QuartetT{T}(quartetrank(tn1,tn2,tn3,tn4,nCk), SVector(tn1,tn2,tn3,tn4), data)
end

"""
    quartetrank(t1,t2,t3,t4, nCk::Matrix)
    quartetrank([t1,t2,t3,t4], nCk)

Return the rank of a four-taxon set with taxon numbers `t1,t2,t3,t4`,
assuming that `ti`s are positive integers such that t1<t2, t2<t3 and t3<t4
(assumptions not checked!).
`nCk` should be a matrix of "n choose k" binomial coefficients:
see [`nchoose1234`](@ref).

# examples

```jldoctest
julia> nCk = PhyloNetworks.nchoose1234(5)
6×4 Matrix{Int64}:
 0   0   0  0
 1   0   0  0
 2   1   0  0
 3   3   1  0
 4   6   4  1
 5  10  10  5

julia> PhyloNetworks.quartetrank([1,2,3,4], nCk)
1

julia> PhyloNetworks.quartetrank([3,4,5,6], nCk)
15
```
"""
@inline function quartetrank(tnum::AbstractVector, nCk::Matrix)
    quartetrank(tnum..., nCk)
end
@inline function quartetrank(t1::Int, t2::Int, t3::Int, t4::Int, nCk::Matrix)
    # rank-1 = t1-1 choose 1 + t2-1 choose 2 + t3-1 choose 3 + t4-1 choose 4
    return nCk[t1,1] + nCk[t2,2] + nCk[t3,3] + nCk[t4,4] + 1
end

"""
    nchoose1234(nmax)

`nmax+1 x 4` matrix containing the binomial coefficient
"n choose k" in row `n+1` and column `k`. In other words,
`M[i,k]` gives "i-1 choose k". It is useful to store these
values and look them up to rank (a large number of) 4-taxon sets:
see [`quartetrank`](@ref).
"""
function nchoose1234(nmax::Int)
    # compute nC1, nC2, nC3, nC4 for n in [0, nmax]: used for ranking quartets
    M = Matrix{Int}(undef, nmax+1, 4)
    for i in 1:(nmax+1)
        M[i,1] = i-1 # n choose 1 = n. row i is for n=i-1
    end
    M[1,2:4] .= 0 # 0 choose 2,3,4 = 0
    for i in 2:(nmax+1)
        for k in 2:4 # to choose k items in 1..n: the largest could be n, else <= n-1
            M[i,k] = M[i-1,k-1] + M[i-1,k]
        end
    end
    return M
end

# Data on quartet concordance factors -------

"""
    DataCF

type that contains the following attributes:

- quartet (vector of Quartets)
- numQuartets
- tree (vector of trees: empty if a table of CF was input instead of list of trees)
- numTrees (-1 if a table CF was input instead of list of trees)
- repSpecies (taxon names that were repeated in table of CF or input gene trees: used inside snaq for multiple alleles case)

The list of Quartet may be accessed with the attribute .quartet.
If the input was a list of trees, the HybridNetwork's can be accessed with the attribute .tree.
For example, if the DataCF object is named d, d.quartet[1] will show the first quartet
and d.tree[1] will print the first input tree.
"""
mutable struct DataCF # fixit
    quartet::Array{Quartet,1} # array of quartets read from CF output table or list of quartets in file
    numQuartets::Integer # number of quartets
    tree::Vector{HybridNetwork} #array of input gene trees
    numTrees::Integer # number of gene trees
    repSpecies::Vector{String} #repeated species in the case of multiple alleles
    DataCF(quartet::Array{Quartet,1}) = new(quartet,length(quartet),[],-1,[])
    DataCF(quartet::Array{Quartet,1},trees::Vector{HybridNetwork}) = new(quartet,length(quartet),trees,length(trees),[])
    DataCF() = new([],0,[],-1,[])
end

# aux type for the updateBL function
mutable struct EdgeParts
    edgenum::Int
    part1::Vector{Node}
    part2::Vector{Node}
    part3::Vector{Node}
    part4::Vector{Node}
end

struct RootMismatch <: Exception
    msg::String
end
RootMismatch() = RootMismatch("")
Base.showerror(io::IO, e::RootMismatch) = print(io, "RootMismatch: ", e.msg);
