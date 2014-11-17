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

# aux function to make a hybrid node a tree node
# used in deleteLeaf
# input: hybrid node
function makeNodeTree!(hybrid::Node)
    if(hybrid.hybrid)
        warn("we make node $(hybrid.number) a tree node, but it can still have hybrid edges pointing at it")
        hybrid.hybrid = true
        hybrid.gammaz = -1
        hybrid.isBadDiamondI = false
        hybrid.isBadDiamondII = false
        hybrid.isBadTriangleI = false
        hybrid.isBadTriangleII = false
        hybrid.k = -1
    else
        error("cannot make node $(hybrid.number) tree node because it already is")
    end
end

# aux function to make a hybrid edge tree edge
# used in deleteLeaf
# input: edge and hybrid node it pointed to
function makeEdgeTree!(edge::Edge, node::Node)
    if(edge.hybrid)
        if(node.hybrid)
            warn("we make edge $(edge.number) a tree edge, but it will still point to hybrid node $(node.number)")
            edge.hybrid = false
            edge.isMajor = true
            edge.gamma = 1.0
            getOtherNode(edge,node).hasHybEdge = false
        else
            error("need the hybrid node at which edge $(edge.number) is pointing to, node $(node.number) is tree node")
        end
    else
        error("cannot make edge $(edge.number) tree because it is tree already")
    end
end

# function to delete an internal node in an external edge
# input: network, internal node and leaf
# returns the new middle
function deleteIntLeaf!(net::HybridNetwork, middle::Node, leaf::Node)
    if(size(middle.edge,1) == 2)
        if(isequal(getOtherNode(middle.edge[1],middle),leaf))
            leafedge = middle.edge[1]
            otheredge = middle.edge[2]
        elseif(isequal(getOtherNode(middle.edge[2],middle),leaf))
            leafedge = middle.edge[2]
            otheredge = middle.edge[1]
        else
            error("leaf $(leaf.number) is not attached to internal node $(middle.number) by any edge")
        end
        othernode = getOtherNode(otheredge,middle)
        setLength!(leafedge,leafedge.length + otheredge.length)
        removeNode!(middle,leafedge)
        removeEdge!(othernode,otheredge)
        setNode!(leafedge,othernode)
        setEdge!(othernode,leafedge)
        deleteNode!(net,middle)
        deleteEdge!(net,otheredge)
        return othernode
    else
        error("internal node $(middle.number) does not have two edges only, it has $(size(middle.edge,1))")
    end
end


# function to delete a leaf from a network
# input: network, leaf node
# warning: it will delete from the actual network
#          need to create a copy before calling this
#          function
# fixit: still missing to test hybrid-> bad triangle II and
#        normal case (not bad triangle/diamond) of hasHybEdge
function deleteLeaf!(net::HybridNetwork, leaf::Node)
    if(leaf.leaf)
        if(size(leaf.edge,1) == 1)
            other = getOtherNode(leaf.edge[1],leaf);
            if(other.hybrid)
                edge1,edge2 = hybridEdges(other,leaf.edge[1]);
                if(!edge1.hybrid || !edge2.hybrid)
                    error("hybrid node $(other.node) does not have two hybrid edges, they are tree edges: $(edge1.number), $(edge2.number)")
                end
                other1 = getOtherNode(edge1,other);
                other2 = getOtherNode(edge2,other);
                removeEdge!(other1,edge1)
                removeEdge!(other2,edge2)
                deleteEdge!(net,edge1)
                deleteEdge!(net,edge2)
                deleteEdge!(net,leaf.edge[1])
                deleteNode!(net,other)
                deleteNode!(net,leaf)
                if(other.isBadTriangleII)
                    edge11 = isequal(getOtherNode(other1.edge[1],other1),other2) ? other1.edge[2] : other1.edge[1]
                    if(isequal(getOtherNode(other2.edge[1],other2),other1))
                        edge12 = other2.edge[2]
                        other3 = getOtherNode(other2.edge[2],other2)
                    else
                        edge12 = other2.edge[1]
                        other3 = getOtherNode(other2.edge[1],other2)
                    end
                    removeNode!(other1,edge11)
                    removeEdge!(other1,edge11)
                    setEdge!(other3,edge11)
                    setNode!(edge11,other3)
                    setLength!(edge11,other.gammaz+other1,gammaz*other2.gammaz)
                else
                    if(size(other1.edge,1) == 2)
                        node = getOtherNode(other1.edge[1],other1)
                        leaf1 = node.leaf ? node : getOtherNode(other1.edge[2],other1)
                        node = deleteIntLeaf!(net,other1,leaf1);
                    end
                    if(size(other2.edge,1) == 2)
                        node = getOtherNode(other2.edge[1],other1)
                        leaf1 = node.leaf ? node : getOtherNode(other2.edge[2],other2)
                        node = deleteIntLeaf!(net,other2,leaf1);
                    end
                end
            else
                if(other.hasHybEdge)
                    println("entro a caso has hyb edge")
                    edge1,edge2 = hybridEdges(other,leaf.edge[1]);
                    if(!edge1.hybrid && !edge2.hybrid)
                        error("node $(other.number) has hybrid edge attribute true, but the edges $(edge1.number), $(edge2.number) are not hybrid (and the third edge has a leaf $(leaf.number)")
                    end
                    other1 = getOtherNode(edge1,other);
                    other2 = getOtherNode(edge2,other);
                    if(other1.hybrid)
                        if(other1.isBadDiamondI)
                            edgemaj,edgemin,edgebla = hybridEdges(other1)
                            edge4 = isequal(edge1,edgemaj) ? edgemin : edgemaj
                            other3 = getOtherNode(edge4,other1)
                            edgebla,edge3,edge5 = hybridEdges(other3)
                            leaf5 = getOtherNode(edge5,other3)
                            removeNode!(other,edge2)
                            removeEdge!(other1,edge1)
                            setNode!(edge2,other1)
                            setEdge!(other1,edge2)
                            if(other3.gammaz != -1)
                                println("entro a cambiar length en edge $(edge2.number) con gammaz $(other3.gammaz)")
                                setLength!(edge2,-log(1-other3.gammaz))
                            else
                                error("hybrid node $(other1.number) is bad diamond, but for node $(other3.number), gammaz is not well updated, it is $(other3.gammaz)")
                            end
                            makeEdgeTree!(edge4,other1)
                            makeNodeTree!(other1)
                            removeEdge!(other2,edge3)
                            removeEdge!(other3,edge3)
                            deleteNode!(net,other)
                            deleteEdge!(net,edge1)
                            deleteEdge!(net,edge3)
                            removeNode!(other3,edge4)
                            removeEdge!(leaf5,edge5)
                            setNode!(edge4,leaf5)
                            setEdge!(leaf5,edge4)
                            deleteNode!(net,other3)
                            deleteEdge!(net,edge5)
                        elseif(other1.isBadTriangleI)
                            println("identify node $(other1.number) as bad triangle")
                            edgemaj,edgemin,treeedge = hybridEdges(other1)
                            println("treeedge is edge $(treeedge.number)")
                            edge3 = isequal(edge1,edgemaj) ? edgemin : edgemaj
                            removeNode!(other1,treeedge)
                            removeEdge!(other2,edge3)
                            removeEdge!(other1,edge1)
                            setEdge!(other2,treeedge)
                            setNode!(treeedge,other2)
                            if(other1.gammaz != -1)
                                setLength!(treeedge,-log(1-other1.gammaz))
                            else
                                error("hybrid node $(other1.number) is bad triangle I, but gammaz is not updated, it is $(other1.gammaz)")
                            end
                            deleteEdge!(net,edge2)
                            deleteEdge!(net,edge1)
                            deleteEdge!(net,edge3)
                            deleteNode!(net,other1)
                            deleteNode!(net,other)
                        else
                            removeEdge!(other,leaf.edge[1])
                            edgebla,edgebla,treeedge = hybridEdges(other1)
                            if(getOtherNode(treeedge,other1).leaf)
                                setLength!(edge1,edge2.length+edge1.length)
                                removeEdge!(other2,edge2)
                                removeNode!(other,edge1)
                                setNode!(edge1,other2)
                                setEdge!(other2,edge1)
                                deleteEdge!(net,edge2)
                                deleteNode!(net,other)
                            end
                        end
                        deleteNode!(net,leaf)
                        deleteEdge!(net,leaf.edge[1])
                    elseif(other2.hybrid)
                        if(other2.isBadDiamondI)
                            edgemaj,edgemin,edgebla = hybridEdges(other2)
                            edge4 = isequal(edge2,edgemaj) ? edgemin : edgemaj
                            other3 = getOtherNode(edge4,other2)
                            edgebla,edge3,edge5 = hybridEdges(other3)
                            leaf5 = getOtherNode(edge5,other3)
                            removeNode!(other,edge1)
                            removeEdge!(other2,edge2)
                            setNode!(edge1,other2)
                            setEdge!(other2,edge1)
                            if(other3.gammaz != -1)
                                setLength!(edge1,-log(1-other3.gammaz))
                            else
                                error("hybrid node $(other2.number) is bad diamond, but for node $(other3.number), gammaz is not well updated, it is $(other3.gammaz)")
                            end
                            makeEdgeTree!(edge4,other2)
                            makeNodeTree!(other2)
                            removeEdge!(other1,edge3)
                            removeEdge!(other3,edge3)
                            deleteNode!(net,other)
                            deleteEdge!(net,edge2)
                            deleteEdge!(net,edge3)
                            removeNode!(other3,edge4)
                            removeEdge!(leaf5,edge5)
                            setNode!(edge4,leaf5)
                            setEdge!(leaf5,edge4)
                            deleteNode!(net,other3)
                            deleteEdge!(net,edge5)
                        elseif(other2.isBadTriangleI)
                            edgemaj,edgemin,treeedge = hybridEdges(other2)
                            edge3 = isequal(edge2,edgemaj) ? edgemin : edgemaj
                            removeNode!(other2,treeedge)
                            removeEdge!(other1,edge3)
                            removeEdge!(other1,edge1)
                            setEdge!(other1,treeedge)
                            setNode!(treeedge,other1)
                            if(other2.gammaz != -1)
                                setLength!(treeedge,-log(1-other2.gammaz))
                            else
                                error("hybrid node $(other2.number) is bad triangle I, but gammaz is not updated, it is $(other2.gammaz)")
                            end
                            deleteEdge!(net,edge1)
                            deleteEdge!(net,edge2)
                            deleteEdge!(net,edge3)
                            deleteNode!(net,other2)
                            deleteNode!(net,other)
                        else
                            removeEdge!(other,leaf.edge[1])
                            edgebla,edgebla,treeedge = hybridEdges(other2)
                            if(getOtherNode(treeedge,other2).leaf)
                                setLength!(edge2,edge2.length+edge1.length)
                                removeEdge!(other1,edge1)
                                removeNode!(other,edge2)
                                setNode!(edge2,other1)
                                setEdge!(other1,edge2)
                                deleteEdge!(net,edge1)
                                deleteNode!(net,other)
                            end
                        end
                        deleteNode!(net,leaf)
                        deleteEdge!(net,leaf.edge[1])
                    else
                        error("node $(other.number) has hybrid edge, but neither of the other nodes $(other1.number), $(other2.number )are hybrid")
                    end
                else
                    edge1,edge2 = hybridEdges(other,leaf.edge[1]);
                    other1 = getOtherNode(edge1,other);
                    other2 = getOtherNode(edge2,other);
                    removeEdge!(other,leaf.edge[1])
                    deleteNode!(net,leaf)
                    deleteEdge!(net,leaf.edge[1])
                    if(other1.leaf || other2.leaf)
                        if(other1.leaf && other2.leaf)
                            error("just deleted a leaf $(leaf.number) and its two attached nodes are leaves also $(other1.number), $(other2.number)")
                        end
                        newleaf = other1.leaf ? other1 : other2
                        middle = other
                        while(size(middle.edge,1) == 2)
                            middle = deleteIntLeaf!(net,middle,newleaf)
                        end
                    end
                end
            end
        else
            error("strange leaf with $(size(leaf.edge,1)) edges instead of 1")
        end
    else
        error("node $(leaf.number) is not a leaf, cannot delete it")
    end
end


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
