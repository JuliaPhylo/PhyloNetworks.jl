
# Types for Julia implementation of pseudolikelihood estimation (Stage2)
# Claudia August 2014
# in julia: include("types.jl")
#
# Types: Edge, Node, HybridNetwork
################################################################

# procedure to create hybrid network:
# 1) create edges defined as hybrid or not
# 2) create nodes with such edges
# 3) setNode! to add nodes into edges
# 4) create hybrid network
# 5) updateGammaz!

# -------------- EDGE -------------------------#

abstract ANode

# warning: inCycle and containRoot are updated after edge is part of a network
type Edge
    number::Int64
    length::Float64 #default 1.0
    hybrid::Bool
    y::Float64 # exp(-t), cannot set in constructor for congruence
    z::Float64 # 1-y , cannot set in constructor for congruence
    gamma::Float64 # set to 1.0 for tree edges, hybrid?gamma:1.0
    node::Array{ANode,1} # we can also leave blank: node (see issues.jl)
    isChild1::Bool # used for hybrid edges to set the direction (default true)
    isMajor::Bool  # major edge treated as tree edge for network traversal
                   # true if gamma>.5, or if it is the original tree edge
    inCycle::Int64 # = Hybrid node number if this edge is part of a cycle created by such hybrid node
    # -1 if not part of cycle. used to add new hybrid edge. updated after edge is part of a network
    containRoot::Bool # true if this edge can contain a root given the direction of hybrid edges
                      # used to add new hybrid edge. updated after edge is part of a network
    istIdentifiable::Bool # true if the parameter t (length) for this edge is identifiable as part of a network
                          # updated after part of a network
    fromBadDiamondI::Bool # true if the edge came from deleting a bad diamond I hybridization
    # inner constructors: ensure congruence among (length, y, z) and (gamma, hybrid, isMajor), and size(node)=2
    Edge(number::Int64) = new(number,1.0,false,exp(-1.0),1-exp(-1.0),1.,[],true,true,-1,true,true,false)
    Edge(number::Int64, length::Float64) = new(number,length,false,exp(-length),1-exp(-length),1.,[],true,true,-1,true,true,false)
    Edge(number::Int64, length::Float64,hybrid::Bool,gamma::Float64)= new(number,length,hybrid,exp(-length),1-exp(-length),hybrid?gamma:1.,[],true,(!hybrid || gamma>0.5)?true:false,-1,!hybrid,true,false)
    Edge(number::Int64, length::Float64,hybrid::Bool,gamma::Float64,isMajor::Bool)= new(number,length,hybrid,exp(-length),1-exp(-length),hybrid?gamma:1.,[],true,isMajor,-1,!hybrid,true,false)
    function Edge(number::Int64, length::Float64,hybrid::Bool,gamma::Float64,node::Array{ANode,1})
        size(node,1) != 2 ?
        error("vector of nodes must have exactly 2 values") :
        new(number,length,hybrid,exp(-length),1-exp(-length),hybrid?gamma:1.,node,true,(!hybrid || gamma>0.5)?true:false,-1,!hybrid,true,false)
    end
    function Edge(number::Int64, length::Float64,hybrid::Bool,gamma::Float64,node::Array{ANode,1},isChild1::Bool, inCycle::Int32, containRoot::Bool, istIdentifiable::Bool)
        size(node,1) != 2 ?
        error("vector of nodes must have exactly 2 values") :
        new(number,length,hybrid,exp(-length),1-exp(-length),hybrid?gamma:1.,node,isChild1,(!hybrid || gamma>0.5)?true:false,inCycle,containRoot,istIdentifiable,false)
    end
end

# warning: gammaz, inCycle, isBadTriangle/Diamond updated until the node is part of a network
type Node <: ANode
    number::Int64
    leaf::Bool
    hybrid::Bool
    gammaz::Float64  # notes file for explanation. gammaz if tree node, gamma2z if hybrid node.
                     # updated after node is part of network with updateGammaz!
    edge::Array{Edge,1}
    hasHybEdge::Bool #is there a hybrid edge in edge? only needed when hybrid=false (tree node)
    isBadDiamondI::Bool # for hybrid node, is it bad diamond case I, update in updateGammaz!
    isBadDiamondII::Bool # for hybrid node, is it bad diamond case II, update in updateGammaz!
    isExtBadTriangle::Bool # for hybrid node, is it extremely bad triangle, udpate in updateGammaz!
    isVeryBadTriangle::Bool # for hybrid node, is it very bad triangle, udpate in updateGammaz!
    isBadTriangle::Bool # for hybrid node, is it very bad triangle, udpate in updateGammaz!
    inCycle::Int64 # = hybrid node if this node is part of a cycle created by such hybrid node, -1 if not part of cycle
    prev # previous node in cycle, used in updateInCycle. defined as "Any", set as "nothing" to begin with
    k::Int64 # number of nodes in cycle, only stored in hybrid node and updated after node becomes part of network
             # default -1
    typeHyb::Int64 # type of hybridization, needed for quartet network only, default -1
    name::AbstractString
    # inner constructor: set hasHybEdge depending on edge
    Node() = new(-1.,false,false,-1.,[],false,false,false,false,false,false,-1.,nothing,-1,-1,"")
    Node(number::Int64, leaf::Bool) = new(number,leaf,false,-1.,[],false,false,false,false,false,false,-1.,nothing,-1,-1,"")
    Node(number::Int64, leaf::Bool, hybrid::Bool) = new(number,leaf,hybrid,-1.,[],hybrid,false,false,false,false,false,-1.,nothing,-1,-1,"")
    Node(number::Int64, leaf::Bool, hybrid::Bool, edge::Array{Edge,1})=new(number,leaf,hybrid,-1.,edge,!all((e->!e.hybrid),edge),false,false,false,false,false,-1.,nothing,-1,-1,"")
    Node(number::Int64, leaf::Bool, hybrid::Bool,gammaz::Float64, edge::Array{Edge,1}) = new(number,leaf,hybrid,gammaz,edge,!all((e->!e.hybrid), edge),false,false,false,false,false,-1.,nothing,-1,-1,"")
end

# partition type
type Partition
    cycle::Vector{Int64} #hybrid node number for cycle (or cycles)
    edges::Vector{Edge} #edges in partition
end

abstract Network

# warning: no attempt to make sure the direction of edges matches with the root
# warning: no check if it is network or tree, node array can have no hybrids
# warning: nodes and edges need to be defined and linked before adding to a network
"""
`HybridNetwork type`
Explicit network or tree with the following attributes:

- numTaxa
- numNodes (total number of nodes)
- numEdges
- numHybrids (number of hybrid nodes)
- edge (array of Edges)
- node (array of Nodes)
- root (index of root in vector 'node'. May be artificial, for printing and traversal purposes only.)
- hybrids (array of Nodes: those are are hybrid nodes)
- leaf (array of Nodes: those that are leaves)
- loglik (negative log pseudolik after estimation)
- isRooted (true or false)
"""
type HybridNetwork <: Network
    numTaxa::Int64 # cannot set in constructor for congruence
    numNodes::Int64 # cannot set in constructor for congruence
    numEdges::Int64 # cannot set in constructor for congruence
    node::Array{Node,1}
    edge::Array{Edge,1}
    root::Int64 # node[root] is the root node, default 1
    names::Array{ASCIIString,1} # translate table for taxon names --but also includes hybrid names...
    hybrid::Array{Node,1} # array of hybrid nodes in network
    numHybrids::Int64 # number of hybrid nodes
    visited::Array{Bool,1} # reusable array of booleans
    edges_changed::Array{Edge,1} # reusable array of edges
    nodes_changed::Array{Node,1} # reusable array of nodes
    leaf::Array{Node,1} # array of leaves
    ht::Vector{Float64} # vector of parameters to optimize
    numht::Vector{Int64} # vector of number of the hybrid nodes and edges in ht e.g. [3,6,8,...], 2 hybrid nodes 3,6, and edge 8 is the 1st identifiable
    numBad::Int64 # number of bad diamond I hybrid nodes, set as 0
    hasVeryBadTriangle::Bool # true if the network has extremely/very bad triangles that should be ignored
    index::Vector{Int64} #index in net.edge, net.node of elements in net.ht to make updating easy
    loglik::Float64 # value of the min -loglik after optBL
    blacklist::Vector{Int64} # reusable array of integers, used in afterOptBL
    partition::Vector{Partition} # to choose edges from a partition only to avoid intersecting cycles
    cleaned::Bool # attribute to know if the network has been cleaned after readm default false
    isRooted::Bool #attribute to know if network is rooted (after root! or directEdges)
    # maxTaxNumber::Int32 --in case it's needed later when we prune taxa
    # inner constructor
    function HybridNetwork(node::Array{Node,1},edge::Array{Edge,1})
        hybrid=Node[];
        leaf=Node[];
        [n.hybrid?push!(hybrid,n):nothing for n in node];
        [n.leaf?push!(leaf,n):nothing for n in node];
        new(size(leaf,1),size(node,1),size(edge,1),node,edge,1,[],hybrid,size(hybrid,1),[],[],[],leaf,[],[],0,false,[],0,[],[],false,false)
    end
    function HybridNetwork(node::Array{Node,1},edge::Array{Edge,1},root::Int64)
        hybrid=Node[];
        leaf=Node[];
        [n.hybrid?push!(hybrid,n):nothing for n in node];
        [n.leaf?push!(leaf,n):nothing for n in node];
        new(size(leaf,1),size(node,1),size(edge,1),node,edge,root,[],hybrid,size(hybrid,1),[],[],[],leaf,[],[],0,false,[],0,[],[],false,false)
    end
    HybridNetwork() = new(0,0,0,[],[],0,[],[],0,[],[],[],[],[],[],0,false,[],0,[],[],false,false);
end

# type created from a HybridNetwork only to extract a given quartet
type QuartetNetwork <: Network
    numTaxa::Int64
    numNodes::Int64
    numEdges::Int64
    node::Array{Node,1}
    edge::Array{Edge,1}
    hybrid::Array{Node,1} # array of hybrid nodes in network
    leaf::Array{Node,1} # array of leaves
    numHybrids::Int64 # number of hybrid nodes
    hasEdge::Array{Bool,1} # array of boolean with all the original identifiable edges of HybridNetwork and gammas (net.ht)
    quartetTaxon::Array{ASCIIString,1} # the quartet taxa in the order it represents
    which::Int64 # 0 it tree quartet, 1 is equivalent to tree quartet and 2 if two minor CF different, default -1
    typeHyb::Array{Int64,1} #array with the type of hybridization of each hybrid node in the quartet
    t1::Float64 # length of internal edge, used when qnet.which=1, default = -1
    names::Array{ASCIIString,1} # translate table for taxon names
    split::Array{Int64,1} # split that denotes to which side each leaf is from the split, i.e. [1,2,2,1] means that leaf1 and 4 are on the same side of the split, default -1,-1,-1,-1
    formula::Array{Int64,1} # array for qnet.which=1 that indicates if the expCf is major (2) or minor (1) at qnet.expCF[i] depending on qnet.formula[i], default -1,-1,-1
    expCF::Array{Float64,1} # three expected CF in order 12|34, 13|24, 14|23 (matching obsCF from qnet.quartet), default [0,0,0]
    indexht::Vector{Int64} # index in net.ht for each edge in qnet.ht
    changed::Bool # true if the expCF would be changed with the current parameters in the optimization, to recalculate, default true
    index::Vector{Int64} # index in qnet.edge (qnet.node for gammaz) of the members in qnet.indexht to know how to find quickly in qnet
    # inner constructor
    function QuartetNetwork(net::HybridNetwork)
        net2 = deepcopy(net); #fixit: maybe we dont need deepcopy of all, maybe only arrays
        new(net2.numTaxa,net2.numNodes,net2.numEdges,net2.node,net2.edge,net2.hybrid,net2.leaf,net2.numHybrids, [true for e in net2.edge],[],-1,[], -1.,net2.names,[-1,-1,-1,-1],[-1,-1,-1],[0,0,0],[],true,[])
        #new(sum([n.leaf?1:0 for n in net.node]),size(net.node,1),size(net.edge,1),copy(net.node),copy(net.edge),copy(net.hybrid),size(net.hybrid,1), [true for e in net2.edge],[],-1,[],-1.,net2.names,[-1,-1,-1,-1],[-1,-1,-1],[],true,[])
    end
    function QuartetNetwork(net::HybridNetwork,quartet::Array{ASCIIString,1})
        net2 = deepcopy(net);
        new(net2.numTaxa,net2.numNodes,net2.numEdges,net2.node,net2.edge,net2.hybrid,net2.leaf,net2.numHybrids, [true for e in net2.edge],quartet,-1,[],-1.,net2.names,[-1,-1,-1,-1],[-1,-1,-1],[0,0,0],[],true,[])
        #new(sum([n.leaf?1:0 for n in net.node]),size(net.node,1),size(net.edge,1),copy(net.node),copy(net.edge),copy(net.hybrid),size(net.hybrid,1), [true for e in net2.edge],quartet,-1,[],-1.,net2.names,[-1,-1,-1,-1],[-1,-1,-1],[],true,[])
    end
    QuartetNetwork() = new(0,0,0,[],[],[],[],0,[],[],0,[],0.0,[],[],[],[],[],true,[])
end

"""
`Quartet type`

type that saves the information on a given 4-taxon subset. It contains the following attributes:

- number
- taxon (vector of taxon names)
- obsCF (vector of observed CF)
- logPseudoLik
- ngenes (number of gene trees used to compute the observed CF: -1 if unknown)
- qnet (internal topological structure that saves the expCF after snaq estimation to emphasize that the expCF depend on a specific network, not the data)
"""
type Quartet
    number::Int64
    taxon::Array{ASCIIString,1}
    obsCF::Array{Float64,1} # three observed CF in order 12|34, 13|24, 14|23
    qnet::QuartetNetwork # quartet network for the current network (want to keep as if private attribute)
    logPseudoLik::Float64 # log pseudolik value for the quartet
    ngenes::Number # number of gene trees used to compute the obsCV, default -1; changed to Number in case ngenes is average, so float
    # inner constructor: to guarantee obsCF are only three and add up to 1
    function Quartet(number::Int64,t1::AbstractString,t2::AbstractString,t3::AbstractString,t4::AbstractString,obsCF::Array{Float64,1})
        size(obsCF,1) != 3 ? error("observed CF vector should have size 3, not $(size(obsCF,1))") : nothing
        0.99 < sum(obsCF) < 1.02 || warn("observed CF should add up to 1, not $(sum(obsCF))")
        new(number,[t1,t2,t3,t4],obsCF,QuartetNetwork(),0,-1);
    end
    function Quartet(number::Int64,t1::Array{ASCIIString,1},obsCF::Array{Float64,1})
        size(obsCF,1) != 3 ? error("observed CF vector should have size 3, not $(size(obsCF,1))") : nothing
        0.99< sum(obsCF) < 1.02 || warn("observed CF should add up to 1, not $(sum(obsCF))")
        size(t1,1) != 4 ? error("array of taxa should have size 4, not $(size(t1,1))") : nothing
        0.0 <= obsCF[1] <= 1.0 || error("obsCF must be between (0,1), but it is $(obsCF[1]) for $(t1)")
        0.0 <= obsCF[2] <= 1.0 || error("obsCF must be between (0,1), but it is $(obsCF[2]) for $(t1)")
        0.0 <= obsCF[3] <= 1.0 || error("obsCF must be between (0,1), but it is $(obsCF[3]) for $(t1)")
        new(number,t1,obsCF,QuartetNetwork(),0,-1);
    end
    Quartet() = new(0,[],[],QuartetNetwork(),0,-1)
end

# Data -------

"""
`DataCF type`

type that contains the following attributes:

- quartet (vector of Quartets)
- numQuartets
- tree (vector of trees: empty if a table of CF was input instead of list of trees)
- numTrees (-1 if a table CF was input instead of list of trees)
- repSpecies (taxon names that were repeated in table of CF or input gene trees: used inside snaq for multiple alleles case)
"""
type DataCF
    quartet::Array{Quartet,1} # array of quartets read from CF output table or list of quartets in file
    numQuartets::Int64 # number of quartets
    tree::Vector{HybridNetwork} #array of input gene trees
    numTrees::Int64 # number of gene trees
    repSpecies::Vector{ASCIIString} #repeated species in the case of multiple alleles
    DataCF(quartet::Array{Quartet,1}) = new(quartet,length(quartet),[],-1,[])
    DataCF(quartet::Array{Quartet,1},trees::Vector{HybridNetwork}) = new(quartet,length(quartet),trees,length(trees),[])
    DataCF() = new([],0,[],-1,[])
end

# aux type for the updateBL function
type EdgeParts
    edgenum::Int64
    part1::Vector{Node}
    part2::Vector{Node}
    part3::Vector{Node}
    part4::Vector{Node}
end

