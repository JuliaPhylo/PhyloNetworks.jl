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
# 5) updateGammaz! updateGamma2z!

# -------------- EDGE -------------------------#

abstract ANode

# warning: inCycle and containRoot are updated after edge is part of a network
type Edge
    number::Int64
    length::Float64
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
    # inner constructors: ensure congruence among (length, y, z) and (gamma, hybrid, isMajor), and size(node)=2
    Edge(number::Int64) = new(number,1.0,false,exp(-1.0),1-exp(-1.0),1.,[],true,true,-1,true,true)
    Edge(number::Int64, length::Float64) = new(number,length,false,exp(-length),1-exp(-length),1.,[],true,true,-1,true,true)
    Edge(number::Int64, length::Float64,hybrid::Bool,gamma::Float64)= new(number,length,hybrid,exp(-length),1-exp(-length),hybrid?gamma:1.,[],true,(!hybrid || gamma>0.5)?true:false,-1,true,true)
        Edge(number::Int64, length::Float64,hybrid::Bool,gamma::Float64,isMajor::Bool)= new(number,length,hybrid,exp(-length),1-exp(-length),hybrid?gamma:1.,[],true,isMajor,-1,true,true)
    function Edge(number::Int64, length::Float64,hybrid::Bool,gamma::Float64,node::Array{ANode,1})
        size(node,1) != 2 ?
        error("vector of nodes must have exactly 2 values") :
        new(number,length,hybrid,exp(-length),1-exp(-length),hybrid?gamma:1.,node,true,(!hybrid || gamma>0.5)?true:false,-1,true,true)
    end
    function Edge(number::Int64, length::Float64,hybrid::Bool,gamma::Float64,node::Array{ANode,1},isChild1::Bool, inCycle::Int32, containRoot::Bool, istIdentifiable::Bool)
        size(node,1) != 2 ?
        error("vector of nodes must have exactly 2 values") :
        new(number,length,hybrid,exp(-length),1-exp(-length),hybrid?gamma:1.,node,isChild1,(!hybrid || gamma>0.5)?true:false,inCycle,containRoot,istIdentifiable)
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
    isBadTriangleI::Bool # for hybrid node, is it bad triangle case I, udpate in updateGammaz!
    isBadTriangleII::Bool # for hybrid node, is it bad triangle case II, udpate in updateGammaz!
    inCycle::Int64 # = hybrid node if this node is part of a cycle created by such hybrid node, -1 if not part of cycle
    prev # previous node in cycle, used in updateInCycle. defined as "Any", set as "nothing" to begin with
    k::Int64 # number of nodes in cycle, only stored in hybrid node and updated after node becomes part of network
             # default -1
    # inner constructor: set hasHybEdge depending on edge
    Node() = new(-1.,false,false,-1.,[],false,false,false,false,false,-1.,nothing,-1.)
    Node(number::Int64, leaf::Bool) = new(number,leaf,false,-1.,[],false,false,false,false,false,-1.,nothing,-1.)
    Node(number::Int64, leaf::Bool, hybrid::Bool) = new(number,leaf,hybrid,-1.,[],hybrid,false,false,false,false,-1.,nothing,-1.)
    Node(number::Int64, leaf::Bool, hybrid::Bool, edge::Array{Edge,1})=new(number,leaf,hybrid,-1.,edge,!all([!edge[i].hybrid for i=1:size(edge,1)]),false,false,false,false,-1.,nothing,-1.)
    Node(number::Int64, leaf::Bool, hybrid::Bool,gammaz::Float64, edge::Array{Edge,1}) = new(number,leaf,hybrid,gammaz,edge,!all([!e.hybrid for e in edge]),false,false,false,false,-1.,nothing,-1.)
end

# warning: no attempt to make sure the direction of edges matches with the root
# warning: no check if it is network or tree, node array can have no hybrids
# warning: nodes and edges need to be defined and linked before adding to a network
type HybridNetwork
    numTaxa::Int64 # cannot set in constructor for congruence
    numNodes::Int64 # cannot set in constructor for congruence
    numEdges::Int64 # cannot set in constructor for congruence
    node::Array{Node,1}
    edge::Array{Edge,1}
    root::Int64 # node[root] is the root node, default 1
    names::Array{ASCIIString,1} # translate table for taxon names
    hybrid::Array{Node,1} # array of hybrid nodes in network
    numHybrids::Int64 # number of hybrid nodes
    visited::Array{Bool,1} # reusable array of booleans
    edges_changed::Array{Edge,1} # reusable array of edges
    nodes_changed::Array{Node,1} # reusable array of nodes
    # maxTaxNumber::Int32 --in case it's needed later when we prune taxa
    # inner constructor
    function HybridNetwork(node::Array{Node,1},edge::Array{Edge,1})
        hybrid=Node[];
        [n.hybrid?push!(hybrid,n):nothing for n in node];
        new(sum([node[i].leaf?1:0 for i=1:size(node,1)]),size(node,1),size(edge,1),node,edge,1,[],hybrid,size(hybrid,1),[],[],[])
    end
    function HybridNetwork(node::Array{Node,1},edge::Array{Edge,1},root::Int64)
        hybrid=Node[];
        [n.hybrid?push!(hybrid,n):nothing for n in node];
        new(sum([node[i].leaf?1:0 for i=1:size(node,1)]),size(node,1),size(edge,1),node,edge,root,[],hybrid,size(hybrid,1),[],[],[])
    end
    HybridNetwork() = new(0,0,0,[],[],0,[],[],0,[],[],[]);
end

type Quartet
    number::Int64
    taxon1::ASCIIString
    taxon2::ASCIIString
    taxon3::ASCIIString
    taxon4::ASCIIString
    obsCF::Array{Float64,1} # three observed CF in order 12|34, 13|24, 14|23
    # inner constructor: to guarantee obsCF are only three and add up to 1
    function Quartet(number::Int64,t1::ASCIIString,t2::ASCIIString,t3::ASCIIString,t4::ASCIIString,obsCF::Array{Float64,1})
        size(obsCF,1) != 3 ? error("observed CF vector should have size 3, not $(size(obsCF,1))") : nothing
        sum(obsCF) != 1 ? error("observed CF should add up to 1, not $(sum(obsCF))") : nothing
        new(number,t1,t2,t3,t4,obsCF);
    end
end

# type created from a HybridNetwork only to extract a given quartet
type QuartetNetwork
    node::Array{Node,1}
    edge::Array{Edge,1}
    hybrid::Array{Node,1} # array of hybrid nodes in network
    numHybrids::Int64 # number of hybrid nodes
    hasEdge::Array{Bool,1} # array of boolean with all the original edges of HybridNetwork
    quartet::Quartet # the quartet it represent
    visited::Array{Bool,1} # reusable array of booleans
    edges_changed::Array{Edge,1} # reusable array of edges
    nodes_changed::Array{Node,1} # reusable array of nodes
    # inner constructor
    QuartetNetwork(net::HybridNetwork) = new(net.node,net.edge,net.hybrid,net.numHybrids,[true for e in net.edge],nothing,[],[],[])
    QuartetNetwork(net::HybridNetwork,quartet::Quartet) = new(net.node,net.edge,net.hybrid,net.numHybrids,[true for e in net.edge],quartet,[],[],[])
end
