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
include("bad_triangle_example.jl");

include("tree_example.jl");

# -------------- NETWORK ----------------------- #





# ---------------------- branch length optimization ---------------------------------


# function to update a QuartetNetwork for a given
# vector of parameters based on a boolean vector "changed"
# which shows which parameters have changed
function update!(qnet::QuartetNetwork,x::Vector{Float64}, net::HybridNetwork)
    ch = changed(net,x)
    length(x) == length(ch) || error("x (length $(length(x))) and changed $(length(changed)) should have the same length")
    length(ch) == length(qnet.hasEdge) || error("changed (length $(length(changed))) and qnet.hasEdge (length $(length(qnet.hasEdge))) should have same length")
    qnet.changed = false
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    for(i in 1:length(ch))
        qnet.changed |= (ch[i] & qnet.hasEdge[i])
    end
    if(qnet.changed)
        if(qnet.numHybrids == 1 && qnet.hybrid[1].isBadDiamondI) # qnet.indexht is only two values: gammaz1,gammaz2
            length(qnet.indexht) == 2 || error("strange qnet from bad diamond I with hybrid node, it should have only 2 elements: gammaz1,gammaz2, not $(length(qnet.indexht))")
            for(i in 1:2)
                0 <= x[qnet.indexht[i]] <= 1 || error("new gammaz value should be between 0,1: $(x[qnet.indexht[i]]).")
                qnet.node[qnet.index[i]].gammaz = x[qnet.indexht[i]]
            end
        else
            for i in 1:length(qnet.indexht)
                if(qnet.indexht[i] <= net.numHybrids - net.numBad)
                    0 <= x[qnet.indexht[i]] <= 1 || error("new gamma value should be between 0,1: $(x[qnet.indexht[i]]).")
                    qnet.edge[qnet.index[i]].hybrid || error("something odd here, optimizing gamma for tree edge $(qnet.edge[qnet.index[i]].number)")
                    setGamma!(qnet.edge[qnet.index[i]],x[qnet.indexht[i]])
                    node = qnet.edge[qnet.index[i]].node[qnet.edge[qnet.index[i]].isChild1 ? 1 : 2]
                    node.hybrid || error("hybrid edge $(qnet.edge[qnet.index[i]].number) pointing at tree node $(node.number)")
                    edges = hybridEdges(node,qnet.edge[qnet.index[i]])
                    length(edges) == 2 || error("strange here: node $(node.number) should have 3 edges and it has $(length(edges)+1).")
                    if(edges[1].hybrid && !edges[2].hybrid)
                        setGamma!(edges[1],1-x[qnet.indexht[i]])
                    elseif(edges[2].hybrid && !edges[1].hybrid)
                        setGamma!(edges[2],1-x[qnet.indexht[i]])
                    else
                        error("strange hybrid node $(node.number) with only one hybrid edge or with three hybrid edges")
                    end
                elseif(qnet.indexht[i] <= net.numHybrids - net.numBad + k)
                    setLength!(qnet.edge[qnet.index[i]],x[qnet.indexht[i]])
                else
                    0 <= x[qnet.indexht[i]] <= 1 || error("new gammaz value should be between 0,1: $(x[qnet.indexht[i]]).")
                    setLength!(qnet.edge[qnet.index[i]],-log(1-x[qnet.indexht[i]]))
                end
            end
        end
    end
end




# -------------------------------------------------------------------------------------------------
# ORIGINAL
# function to identify the QuartetNetwork as
# 1 (equivalent to tree), 2 (minor CF different)
# around a given hybrid node
# it also cleans the hybridizations of type 1
# returns 0,1,2
function identifyQuartet!(qnet::QuartetNetwork, node::Node)
    if(node.hybrid)
        k = sum([(n.inCycle == node.number && size(n.edge,1) == 3) ? 1 : 0 for n in qnet.node])
        if(k < 2)
            error("strange quartet network with a hybrid node $(node.number) but no cycle")
        elseif(k == 2)
            other = qnet.node[getIndex(true, [(n.inCycle == node.number && size(n.edge,1) == 3) for n in qnet.node])]
            edgemaj,edgemin,edge1 = hybridEdges(node)
            edgemin2,edgebla,edge2 = hybridEdges(other)
            if(getOtherNode(edge1,node).leaf || getOtherNode(edge2,other).leaf) # k=2, unidentifiable
                leaf = getOtherNode(edge1,node)
                middle = node
                if(!leaf.leaf)
                    leaf = getOtherNode(edge2,node)
                    middle = other
                end
                if(isequal(getOtherNode(edgemaj,node),other))
                    removeEdge!(node,edgemaj)
                    removeEdge!(other,edgemaj)
                    deleteEdge!(qnet,edgemaj)
                    makeNodeTree!(qnet,node)
                    deleteIntLeaf!(qnet,middle,leaf)
                elseif(isequal(getOtherNode(edgemin,node),other))
                    removeEdge!(node,edgemin)
                    removeEdge!(other,edgemin)
                    deleteEdge!(qnet,edgemin)
                    makeNodeTree!(qnet,node)
                    deleteIntLeaf!(qnet,middle,leaf)
                else
                    error("nodes $(node.number) and $(other.number) should be united by a hybrid edge but are not")
                end
                qnet.which = 0
            else

        elseif(k == 3)
            f
        elseif(k == 4)
            f
        else
            error("strange quartet network with $(k) nodes in cycle, maximum should be 4")
        end
    else
        error("cannot identify the hybridization around node $(node.number) because it is not hybrid node.")
    end
end




# function to identify the QuartetNetwork as one of the
# 6 possibilities
function identifyQuartet!(qnet::QuartetNetwork)
    if(qnet.which == -1)
        if(qnet.numHybrids == 0)
            qnet.which = 0
        elseif(qnet.numHybrids == 1)
            qnet.which = identifyQuartet!(qnet,qnet.hybrid[1])
        elseif(qnet.numHybrids > 1)
            for(n in qnet.hybrid)
                identifyQuartet!()
        else
            error("strange quartet network with negative number of hybrids: $(qnet.numHybrids).")
        end
    else
        error("Quartet has already been identified as $(qnet.which)")
    end
end

# -------------



# ------------------

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
