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

# examples
include("case_f_example.jl");
include("bad_triangle_example.jl")

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

# function that given a hybrid node, it gives you the minor hybrid edge
function getHybridEdge(node::Node)
    if(node.hybrid)
        a=nothing;
        for(i in 1:size(node.edge,1))
            (node.edge[i].hybrid && !node.edge[i].isMajor) ? a=node.edge[i] : nothing;
        end
        isa(a,Nothing) ? error("hybrid node does not have minor hybrid edge") : return a
    else
        error("node is not hybrid node")
    end
end

# todo: function to update edges (inCycle, containRoot) information after becoming part of a network,
# check cecile notes for traversal and inCycle updating
# outside:
# if(node.hybrid)
#    hybedge=getHybridEdge(node);
#    node1=getOtherNode(hybedge,node);
visited=[false for i=1:size(net.node,1)];
# fixit: still not working, two main problems:
# one, how to define how is child? so that it is not in infinite loop child<->parent
# and two, node.inCycle propagates outside the cycle, how to contain?
function updateInCycle!(net::HybridNetwork, node::Node, node1::Node,visited::Array{Bool})
    println("start in $(node.number)");
    visited[getIndex(node,net)]=true;
    for(i in 1:size(node.edge,1))
        if(node.edge[i].isMajor || !node.edge[i].hybrid)
            println("we go to edge: $(node.edge[i].number)");
            child=getOtherNode(node.edge[i],node);
            if(isequal(node1,child) || child.inCycle!=-1)
                println("condition true: $(node1.number)=$(child.number) OR $(child.inCycle) NOT -1");
                println("we will make this edge in cycle: $(node.edge[i].number)");
                node.edge[i].inCycle=1;
                println("we will make the current node in cycle: $(node.number)");
                node.inCycle=1;
            elseif(!visited[getIndex(child,net)])
                println("condition NOT true: $(node1.number) NOT $(child.number) OR $(child.inCycle) = -1");
                println("we will go now to child: $(child.number)");
                updateInCycle!(net,child,node1,visited);
            end
        end
    end
end



# cecile: check updategammaz function, maybe we need two functions, one to update when changing length
# one to update when changing gamma? what i like about updategammaz is that you use that directly at the beginning
# of network, so maybe we should consider doing things ourselves inside setLength and setGamma, instead of calling
# update gamma


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
