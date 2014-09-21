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

# needed modules:
using DataStructures # for updateInCycle with queue
using Base.Collections # for updateInCycle with priority queue

# examples
include("case_f_example.jl");
include("bad_triangle_example.jl");

include("tree_example.jl");
# -------------- NETWORK ----------------------- #

# function to delete a hybridization event
# input: hybrid node and network
function deleteHybrid!(node::Node,net::HybridNetwork)
    if(node.hybrid)
        hybedges = hybridEdges(node);
        if(hybedges[1].isMajor)
            hybedge1 = hybedges[1];
            hybedge2 = hybedges[2];
        else
            hybedge1 = hybedges[2];
            hybedge2 = hybedges[1];
        end
        other1 = getOtherNode(hybedge1,node);
        other2 = getOtherNode(hybedge2,node);
        for(e in node.edge)
            if(!e.hybrid)
                treeedge1 = e;
                other3 = getOtherNode(e,node);
            end
        end
        max_node = maximum([e.number for e in net.node]);
        max_edge = maximum([e.number for e in net.edge]);
        edge1 = Edge(max_edge+1,hybedge1.length+treeedge1.length);
        setNode!(edge1,[other1,other3]);
        removeEdge!(other1,hybedge1);
        setEdge!(other1,edge1);
        removeEdge!(other3,treeedge1);
        setEdge!(other3,edge1);
        index = getIndex(hybedge1,net);
        deleteat!(net.edge,index);
        index = getIndex(hybedge2,net);
        deleteat!(net.edge,index);
        index = getIndex(treeedge1,net);
        deleteat!(net.edge,index);
        index = getIndex(other2,net);
        deleteat!(net.edge,index);
        index = getIndex(node,net);
        deleteat!(net.edge,index);
        push!(net.edge,edge1);
        tree_edge = Edge[];
        for(e in other2.edge)
            if(!e.hybrid)
                push!(tree_edge,e);
            end
        end
        #assume node has only 2 tree edges
        treenode1 = getOtherNode(tree_edge[1],other2);
        treenode2 = getOtherNode(tree_edge[2],other2);
        edge2 = Edge(max_edge+2,tree_edge[1].length+tree_edge[2].length);
        setNode!(edge2,[treenode1,treenode2]);
        removeEdge!(treenode1,tree_edge[1]);
        removeEdge!(treenode2,tree_edge[2]);
        setEdge!(treenode1,edge2);
        setEdge!(treenode2,edge2);
    else
        error("node has to be hybrid")
    end
end



# cecile: check updategammaz function, maybe we need two functions, one to update when changing length
# one to update when changing gamma? what i like about updategammaz is that you use that directly at the beginning
# of network, so maybe we should consider doing things ourselves inside setLength and setGamma, instead of calling
# update gamma



# function to traverse the network
# simply prints the traversal path, can be modified to do other things
# needs:
visited  =  [false for i  =  1:size(net.node,1)];

function traverse(net::HybridNetwork, node::Node, visited::Array{Bool,1})
    println("estamos en $(node.number)");
    visited[getIndex(node,net)]  =  true;
    if(node.leaf)
        println("llegamos a leaf $(node.number)");
    else
        for(i in 1:size(node.edge,1))
            other  =  getOtherNode(node.edge[i],node);
            if(!visited[getIndex(other,net)])
                println("vamos a ir a $(other.number)");
                traverse(net,other,visited);
            end
        end
    end
end

# need function to check if after updateContainRoot! there is no place for the root
# careful because updateContainRoot changes things, so maybe we want to be careful and only change
# if the new hybridization is going to stay

# think of the process of adding a hybrid edge:
# updateInCycle: what happens if cycle intersects, can we go back?
# updateContainRoot: what happens if containRoot is empty, can we go back?

# todo: function to create an hybrid edge:
# - make sure the hybridization is "identifiable": not between the same edge, or in a cherry
# - detect whether the new cycle would overlap with another cycle already in the network.
#   just check that the 2 edges to be connected are not already marked as
#   being on a cycle: updateInCycle! returns false
# - detect where the cycle is: i think it always starts in the hybrid node, so simply use searchHybridNode, or use the hybrid node just created
# - mark edges along the cycle with the number of the hybrid edge/node: updateInCycle!
# - create the new nodes and edges, with correct hybrid labels
# - mark which edges can contain the root, check that the set of edges that
#   can contain the root is non-empty: updateContainRoot, still need function to check if empty
# - check cycle configuration (value of k, and clade sizes ni)
#   if bad triangle: set gammaz and gamma2z for appropriate nodes
#   if bad diamond: gammaz for the "other" node (one just created) of each hybrid edge
#   if some parameters need to be set to 0:
# - identify the second hybrid edge, mark it as hybrid
# - depending on gamma, mark one of the 2 edges as the major "tree" edge

# todo: functions to propose a new network
# example: pick 2 edges and add a hybrid edge to link the 2
# todo: function to change direction of hybrid edge (hybrid edge  =  hybrid&&!isMajor),
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
