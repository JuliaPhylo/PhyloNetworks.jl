# Stage2 pseudolikelihood implementation: sketch of classes (types) in Julia
# Claudia (June 2014)
# Classes based on "ane/public/quartetNetwork/notes" and "mbsumtree.h"
# Methods (functions) for each class are coded outside the class and they can
# be used by any class
#
# Ctrl S: fixit, todo, check
######################################################################

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

# warning: gammaz is updated until the node is part of a network
type Node <: ANode
    number::Int64
    leaf::Bool
    hybrid::Bool
    gammaz::Float64  # notes file for explanation. gammaz if tree node, gamma2z if hybrid node. updated after node is part of network
    edge::Array{Edge,1}
    hasHybEdge::Bool #is there a hybrid edge in edge? only needed when hybrid=false (tree node)
    # inner constructor: set hasHybEdge depending on edge
    Node(number::Int64, leaf::Bool) = new(number,leaf,false,-1.,[],false)
    Node(number::Int64, leaf::Bool, hybrid::Bool) = new(number,leaf,hybrid,-1.,[],hybrid)
    Node(number::Int64, leaf::Bool, hybrid::Bool, edge::Array{Edge,1})=new(number,leaf,hybrid,-1.,edge,!all([!edge[i].hybrid for i=1:size(edge,1)]))
    Node(number::Int64, leaf::Bool, hybrid::Bool,gammaz::Float64, edge::Array{Edge,1}) = new(number,leaf,hybrid,gammaz,edge,!all([!edge[i].hybrid for i=1:size(edge,1)]))
end

# warning: no attempt to make sure the direction of edges matches with the root
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


#------------- EDGE functions --------------------#

#Node needs to be defined as hybrid before adding to a hybrid edge
#First, an edge is defined as hybrid, and then the nodes are added to it
function setNode!(edge::Edge, node::Node)
  if(size(edge.node,1) == 2)
    error("vector of nodes already has 2 values");
  else
    push!(edge.node,node);
    n=size(edge.node,1);
    if(edge.hybrid)
	if(n==1)
            if(node.hybrid)
               edge.isChild1=true;
            else
               edge.isChild1=false;
	    end
	else
	    if(node.hybrid)
               if(edge.node[1].hybrid)
                  error("hybrid edge has two hybrid nodes");
               else
                  edge.isChild1=false;
	       end
	    else
	       if(!edge.node[1].hybrid)
	          error("hybrid edge has no hybrid nodes");
	       else
	          edge.isChild1=true;
	       end
	    end
	end
    end
  end
end

function setNode!(edge::Edge,node::Array{Node,1})
  size(node,1) != 2 ?
  error("vector of nodes must have exactly 2 values") :
  edge.node=node;
  if(edge.hybrid)
      if(node[1].hybrid)
          edge.isChild1=true;
      else
          if(node[2].hybrid)
              edge.isChild1=false;
          else
              error("hybrid edge without hybrid node");
          end
      end
   end
end



# -------------- NODE -------------------------#

function setEdge!(node::Node,edge::Edge)
   push!(node.edge,edge);
   edge.hybrid ? node.hasHybEdge=true : node.hasHybEdge=false;
end

function getOtherNode(edge::Edge,node::Node)
  isequal(edge.node[1],node) ? edge.node[2] : edge.node[1]
end
# -------------- NETWORK ----------------------- #

function getIndex(node::Node, net::HybridNetwork)
    i=1;
    while(i<=size(net.node,1) && !isequal(node,net.node[i]))
        i=i+1;
    end
    i>size(net.node,1)?error("node not in network"):return i;
end

function getIndex(edge::Edge, net::HybridNetwork)
    i=1;
    while(i<=size(net.edge,1) && !isequal(edge,net.edge[i]))
        i=i+1;
    end
    i>size(net.edge,1)?error("edge not in network"):return i;
end

function getIndex(bool::Bool, array::Array{Bool,1})
    i=1;
    while(i<=size(array,1) && !isequal(bool,array[i]))
        i=i+1;
    end
    i>size(array,1)?error("$(bool) not in array"):return i;
end

function getIndex(bool::Bool, array::Array{Any,1})
    i=1;
    while(i<=size(array,1) && !isequal(bool,array[i]))
        i=i+1;
    end
    i>size(array,1)?error("$(bool) not in array"):return i;
end


# search the hybrid node(s) in network
# check: do more tests
function searchHybridNode(net::HybridNetwork)
    suma=sum([net.node[i].hybrid?1:0 for i=1:size(net.node,1)]);
    k=getIndex(true,[net.node[i].hybrid for i=1:size(net.node,1)]);
    if(suma>1)
        a=[k];
        count=suma-1;
        index=k;
        vect=[net.node[i].hybrid for i=1:size(net.node,1)];
        while(count>0 && count<size(net.node,1))
            index==1 ? index=getIndex(true,[false,vect[2:size(net.node,1)]]) : index=getIndex(true,[vect[1:(index-1)],false,vect[(index+1):size(net.node,1)]])
            push!(a,index);
            count=count-1;
        end
        return a
    else
        return k
    end
end

# search the hybrid edge in network. assumes only one hybrid edge, or it returns the first one
# fixit: change to return the minor edge
function searchHybridEdge(net::HybridNetwork)
     getIndex(true,[net.edge[i].hybrid for i=1:size(net.edge,1)])
end

# old searchHybrid that traverses the graph looking for the hybrid node (assumes only one hybrid node)
# (assumes only one hybrid, stops after finding one hybrid node)
# visited=Array{Bool,1} created with same length as network.node and keeps track of which nodes have been visited
# hybrid=Array{Bool,1} created with same length as network.node and keeps track of which node is the hybrid
# (check: hybrid array needed to use because return function never worked!)
# check: this is better or having attribute "visit::Bool" in Node and Edge for traversals?
# maybe better to have attribute instead of creating new arrays everytime
# cecile: ideas of preorder, postorder, inorder. in tree travsersals, you don't need to keep track of visited or not
    # for network: isMajor edge treat as tree edge, then same traversal as tree
function searchHybrid(node::Node,visited::Array{Bool,1},net::HybridNetwork,hybrid::Array{Bool,1})
    visited[getIndex(node,net)]=true;
    if(node.hybrid)
        #println("this is the hybrid: $(getIndex(node,net))")
	hybrid[getIndex(node,net)]=true
    elseif(node.hasHybEdge)
         for(i=1:size(node.edge,1))
             node.edge[i].hybrid ? searchHybrid(getOtherNode(node.edge[i],node),visited,net,hybrid): nothing;
         end
    else
         i=1;
         other=getOtherNode(node.edge[i],node);
         while(i<=size(node.edge,1) && other.leaf && visited[getIndex(other,net)])
              i=i+1;
              other=getOtherNode(node.edge[i],node);
         end
         #println("vamos a ir a $(other.number)");
         searchHybrid(other,visited,net,hybrid);
    end
end

# todo: function to update edges (inCycle, containRoot) information after becoming part of a network,
# assume only one hybrid node to start there, another function can run this one for every hybrid node


# function to update gammaz in a network. assumes: only one hybrid node
# needs to have inCycle attributes updated already
# check: assume any tree node that has hybrid Edge has only one tree edge in cycle (true?)
# fixit: add possibility of more hybrid nodes!
# cecile: add hybrid node as parameter to update gammaz only for that hybrid event
# cecile: update gammaz for bad triangle only inside the bad triangle case,
    # note that the gammaz for hybrid node depends on the topology, not value of gamma
function updateGammaz!(net::HybridNetwork)
    index=searchHybridNode(net);
    node=net.node[index];
    if(node.hybrid)
        edge_maj=nothing;
        edge_min=nothing;
        edge_maj2=nothing;
        edge_min2=nothing;
        for(i=1:size(node.edge,1))
            if(isa(edge_maj,Nothing))
	       edge_maj=(node.edge[i].hybrid && node.edge[i].isMajor)? node.edge[i]:nothing;
            end
            if(isa(edge_min,Nothing))
               edge_min=(node.edge[i].hybrid && !node.edge[i].isMajor)? node.edge[i]:nothing;
            end
        end
        other_maj=getOtherNode(edge_maj,node);
        other_min=getOtherNode(edge_min,node);
        for(j=1:size(other_min.edge,1))
            if(isa(edge_min2,Nothing))
               edge_min2=(!other_min.edge[j].hybrid && other_min.edge[j].inCycle != -1)? other_min.edge[j]:nothing;
            end
        end
        for(j=1:size(other_maj.edge,1))
            if(isa(edge_maj2,Nothing))
              edge_maj2=(!other_maj.edge[j].hybrid && other_maj.edge[j].inCycle != -1)? other_maj.edge[j]:nothing;
            end
        end
        other_min.gammaz=edge_min.gamma*edge_min2.z;
        other_maj.gammaz=edge_maj.gamma*edge_maj2.z;
        node.gammaz=edge_maj.gamma^2*edge_maj.z+edge_min.gamma^2*edge_min2.z;
   end
end

# cecile: print more information
function printEdges(net::HybridNetwork)
    println("Edge#\tNode1\tNode2")
    for i in (1:net.numEdges)
        println("$(net.edge[i].number)\t$(net.edge[i].node[1].number)\t$(net.edge[i].node[2].number)")
    end;
end;

function printNodes(net::HybridNetwork)
    println("Node#\tEdges numbers")
    for i in (1:net.numNodes)
        print(net.node[i].number)
        for j in (1:length(net.node[i].edge))
            print("\t$(net.node[i].edge[j].number)")
        end;
        print("\n")
    end;
end;

# find the hybrid edges for a given hybrid node
function hybridEdges(node::Node)
   if(node.hybrid)
      n=size(node.edge,1);
      hybedges=Edge[];
      for(i=1:n)
        if(node.edge[i].hybrid)
	   push!(hybedges,node.edge[i]);
	end
      end
      if(size(hybedges,1)==2)
        return hybedges
      else
        println("node has more or less than 2 hybrid edges");
      end
   else
      println("node is not hybrid");
   end
end


# setLength using updateGammaz
# check: do we want edge as parameter? or its number?
function setLength!(edge::Edge, new_length::Float64,net::HybridNetwork)
  if(new_length<0)
      error("length has to be nonnegative");
  else
      edge.length = new_length;
      edge.y = exp(-new_length);
      edge.z = 1 - edge.y;
      updateGammaz!(net);
  end
end


# setGamma using updateGammaz
# check: we need to put edge as parameter still, maybe we can have a function to search for the hybrid edge?
#        maybe we don't have the actual edge as parameter, only its number? the number of the hybrid node?
function setGamma!(edge::Edge, new_gamma::Float64, net::HybridNetwork)
 if(edge.hybrid)
	if(0 < new_gamma < 1)
	     edge.gamma = new_gamma;
	     new_gamma<0.5 ? edge.isMajor=false : edge.isMajor=true;
	     edge.isChild1 ? ind=1 : ind=2 ; # hybrid edge pointing at node 1 or 2
             hybedges=hybridEdges(edge.node[ind]); # hybrid edges for hybrid node
	     if(!isempty(hybedges) && !isa(hybedges,Nothing)) # change 1-gamma in other hybrid edge
	          isequal(hybedges[1],edge) ? other=2 : other=1;
		  hybedges[other].gamma=1.0-new_gamma;
       	          new_gamma<0.5 ? hybedges[other].isMajor=true : hybedges[other].isMajor=false;
	    end
        updateGammaz!(net);
	else
	     error("gamma has to be between 0 and 1");
        end
  else
	error("cannot change gamma in a tree edge");
  end
end


# todo: function to create an hybrid edge:
# - make sure the hybridization is "identifiable": not between the same edge, or in a cherry
# - detect whether the new cycle would overlap with another cycle already in the network.
#   just check that the 2 edges to be connected are not already marked as
#   being on a cycle.
# - detect where the cycle is
# - mark edges along the cycle with the number of the hybrid edge/node
# - create the new nodes and edges, with correct hybrid labels
# - mark which edges can contain the root, check that the set of edges that
#   can contain the root is non-empty
# - check cycle configuration (value of k, and clade sizes ni)
#   if bad triangle: set gammaz and gamma2z for appropriate nodes
#   if bad diamond: gammaz for the "other" node (one just created) of each hybrid edge
#   if some parameters need to be set to 0:
# - identify the second hybrid edge, mark it as hybrid
# - depending on gamma, mark one of the 2 edges as the major "tree" edge

# todo: functions to propose a new network
# example: pick 2 edges and add a hybrid edge to link the 2
# todo: function to change direction of hybrid edge (hybrid edge=hybrid&&!isMajor),
#                    source or recipient of either hybrid edge, to propose new network


# todo: function readNetwork!(network::Network, string) # check string as parameter
# C function to read in tree (recursive) in mbsum*,
# maybe start reading a tree, and then add the hybrid edge(s)
# string will contain the parenthetical format, maybe not needed as parameter, but as return


# todo: function printTopology!(string, network::Network) # parameters

# todo: function network2Tree(network::Network) function to remove a hybrid edge and transform the network in tree?

# todo: function to reduce network to quartet: think of rules of how to remove hybrid edges, and when do we need to keep them and when not.

# todo: function to check that everything in network makes sense (gamma, t, gammaz, hybrid edges pointing at hybrid nodes, 2 hybrid edges: one major, one minor)

# todo: function to identify bad diamond/triangle in a network?
