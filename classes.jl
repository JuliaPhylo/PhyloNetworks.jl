# Stage2 pseudolikelihood implementation: sketch of classes (types) in Julia
# Claudia (June 2014)
# Classes based on "ane/public/quartetNetwork/notes" and "mbsumtree.h"
# Methods (functions) for each class are coded outside the class and they can
# be used by any class
#
# Working code: after tests, code goes to types.jl or functions.jl
#
# Ctrl S: fixit, todo, check
######################################################################

# types in "types.jl"
include("types.jl")

# functions in "functions.jl"
include("functions.jl")


# -------------- NETWORK ----------------------- #

# old searchHybrid that traverses the graph looking for the hybrid node (assumes only one hybrid node)
# (assumes only one hybrid, stops after finding one hybrid node)
# visited=Array{Bool,1} created with same length as network.node and keeps track of which nodes have been visited
# hybrid=Array{Bool,1} created with same length as network.node and keeps track of which node is the hybrid
# (check: hybrid array needed to use because return function never worked!)
# check: this is better or having attribute "visit::Bool" in Node and Edge for traversals?
# maybe better to have attribute instead of creating new arrays everytime
# cecile: ideas of preorder, postorder, inorder. in tree travsersals, you don't need to keep track of visited or not
    # for network: isMajor edge treat as tree edge, then same traversal as tree
# check: do we need traversals? we have node and edge arrays for the network.
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
# better to include node as parameter, like in updateGammaz
# check cecile notes for traversal and inCycle updating

# todo: function to identify bad diamond, only updateGammaz in bad diamond cases?

# todo: function to identify hybridization event as bad triangle, identify the network or the hybrid node?

# todo: function to update gammaz in the bad triangle case
# cecile: update gammaz for bad triangle only inside the bad triangle case,
# note that the gammaz for hybrid node depends on the topology, not value of gamma

# cecile: check updategammaz function, maybe we need two functions, one to update when changing length
# one to update when changing gamma? what i like about updategammaz is that you use that directly at the beginning
# of network, so maybe we should consider doing things ourselves inside setLength and setGamma, instead of calling
# update gamma

# cecile: attribute for bad diamond, bad triangle for node type, use in setLength, setGamma and updateGammaz


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

# setLength using updateGammaz
# check: do we want edge as parameter? or its number?
# fixit: need to update gamma2z as well for bad triangle case
function setLength!(edge::Edge, new_length::Float64,net::HybridNetwork)
  if(new_length<0)
      error("length has to be nonnegative");
  else
      edge.length = new_length;
      edge.y = exp(-new_length);
      edge.z = 1 - edge.y;
      if(edge.hybrid && edge.isMajor)
          #updateGamma2z for bad triangle
      elseif (!edge.hybrid && edge.inCycle!=-1 && !all([!edge.node[i].hasHybEdge for i=1:size(edge.node,1)])) # tree edge with node.hasHybEdge true
          updateGammaz!(net,getIndexNumNode(edge.inCycle,net));
      end
  end
end


# setGamma using updateGammaz
# check: we need to put edge as parameter still, maybe we can have a function to search for the hybrid edge?
#        maybe we don't have the actual edge as parameter, only its number? the number of the hybrid node?
# fixit: need to update gamma2z as well for bad triangle case
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
        updateGammaz!(net, getIndex(edge.node[ind],net));
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
