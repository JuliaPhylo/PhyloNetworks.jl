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
                   # true if gamma>.5, cannot be set in constructor
    inCycle::Int64 # = Hybrid node number if this edge is part of a cycle created by such hybrid node
    # -1 if not part of cycle. used to add new hybrid edge. updated after edge is part of a network
    containRoot::Bool # true if this edge can contain a root given the direction of hybrid edges
                      # used to add new hybrid edge. updated after edge is part of a network
    # inner constructors: ensure congruence among (length, y, z) and (gamma, hybrid, isMajor), and size(node)=2
    Edge(number::Int64, length::Float64) = new(number,length,false,exp(-length),1-exp(-length),1.,[],true,true,-1,true)
    Edge(number::Int64, length::Float64,hybrid::Bool,gamma::Float64)= new(number,length,hybrid,exp(-length),1-exp(-length),hybrid?gamma:1.,[],true,(!hybrid || gamma>0.5)?true:false,-1,true)
    function Edge(number::Int64, length::Float64,hybrid::Bool,gamma::Float64,node::Array{ANode,1})
        size(node,1) != 2 ?
        error("vector of nodes must have exactly 2 values") :
        new(number,length,hybrid,exp(-length),1-exp(-length),hybrid?gamma:1.,node,true,(!hybrid || gamma>0.5)?true:false,-1,true)
    end
    function Edge(number::Int64, length::Float64,hybrid::Bool,gamma::Float64,node::Array{ANode,1},isChild1::Bool, inCycle::Int32, containRoot::Bool)
        size(node,1) != 2 ?
        error("vector of nodes must have exactly 2 values") :
        new(number,length,hybrid,exp(-length),1-exp(-length),hybrid?gamma:1.,node,isChild1,(!hybrid || gamma>0.5)?true:false,inCycle,containRoot)
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
    isBadDiamond::Bool # for hybrid node, is it bad diamond case, update in updateGammaz!
    isBadTriangle::Bool # for hybrid node, is it bad triangle case, udpate in updateGamma2z!
    inCycle::Int64 # = hybrid node if this node is part of a cycle created by such hybrid node, -1 if not part of cycle
    prev # previous node in cycle, used in updateInCycle. defined as "Any", set as "nothing" to begin with
    # inner constructor: set hasHybEdge depending on edge
    Node(number::Int64, leaf::Bool) = new(number,leaf,false,-1.,[],false,false,false,-1.,nothing)
    Node(number::Int64, leaf::Bool, hybrid::Bool) = new(number,leaf,hybrid,-1.,[],hybrid,false,false,-1.,nothing)
    Node(number::Int64, leaf::Bool, hybrid::Bool, edge::Array{Edge,1})=new(number,leaf,hybrid,-1.,edge,!all([!edge[i].hybrid for i=1:size(edge,1)]),false,false,-1.,nothing)
    Node(number::Int64, leaf::Bool, hybrid::Bool,gammaz::Float64, edge::Array{Edge,1}) = new(number,leaf,hybrid,gammaz,edge,!all([!edge[i].hybrid for i=1:size(edge,1)]),false,false,-1.,nothing)
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
     # maxTaxNumber::Int32 --in case it's needed later when we prune taxa
     # inner constructor
     HybridNetwork(node::Array{Node,1},edge::Array{Edge,1})=new(sum([node[i].leaf?1:0 for i=1:size(node,1)]),size(node,1),size(edge,1),node,edge,1)
     HybridNetwork(node::Array{Node,1},edge::Array{Edge,1},root::Int64)=new(sum([node[i].leaf?1:0 for i=1:size(node,1)]),size(node,1),size(edge,1),node,edge,root)
end
