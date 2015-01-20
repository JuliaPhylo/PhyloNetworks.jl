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



# working on new in classes.jl
# function to update hasEdge attribute in a
# QuartetNetwork after leaves deleted
# with deleteLeaf!
function updateHasEdge!(qnet::QuartetNetwork, net::HybridNetwork)
    warn("function to compare edges depends on edges number being unique")
    warn("assumes no bad scenario, that is, all gammas and internal t are identifiable")
    edges = Bool[]
    for e in net.edge
        if e.istIdentifiable
            push!(edges,isEdgeNumIn(e,qnet.edge))
        end
    end
    qnet.hasEdge = vcat([isNodeNumIn(n,qnet.hybrid) for n in net.hybrid],edges)
end


# function to extract a quartet from a network
# input: QuartetNetwork (already created from HybridNetwork)
#        quartet: array with the 4 leaf nodes to keep
# return: QuartetNetwork with only 4 tips
# it updates qnet.hasEdge and qnet.indexht
function extractQuartet(net::HybridNetwork,quartet::Array{Node,1})
    if(size(quartet,1) != 4)
        error("quartet array should have 4 nodes, it has $(size(quartet,1))")
    else
        if(!quartet[1].leaf || !quartet[2].leaf || !quartet[3].leaf || !quartet[4].leaf)
            error("all four nodes to keep when extracting the quartet should be leaves: $([q.number for q in quartet])")
        else
            qnet = QuartetNetwork(net)
            leaves = copy(qnet.leaf)
            for(n in leaves)
                if(!isNodeNumIn(n,quartet))
                    println("delete leaf $(n.number)")
                    deleteLeaf!(qnet,n)
                end
            end
            # fixit: do something for bad cases when there are no identifiable parameters, but you could not detect them when deleting
            updateHasEdge!(qnet,net)
            parameters!(qnet,net)
            return qnet
        end
    end
end


# ---------------------- branch length optimization ---------------------------------

# function to update qnet.indexht based on net.numht
# warning: assumes net.numht is updated already with parameters!(net)
function parameters!(qnet::QuartetNetwork, net::HybridNetwork)
    if(size(net.numht,1) > 0)
        size(qnet.indexht,1) > 0 ? warn("deleting qnet.indexht to replace with info in net") : nothing
        n2 = net.numht[1:net.numHybrids]
        n = net.numht[net.numHybrids + 1 : length(net.numht)]
        qn2 = Int64[]
        qn = Int64[]
        for(e in qnet.edge)
            if(e.istIdentifiable)
                push!(qn, getIndex(e.number,n)+net.numHybrids)
            end
            if(e.hybrid && !e.isMajor)
                push!(qn2, getIndex(e.node[e.isChild1 ? 1 : 2].number,n2))
            end
        end
        qnet.indexht = vcat(qn2,qn)
    else
        error("net.numht not correctly updated, need to run parameters first")
    end
end


# function to update a QuartetNetwork for a given
# vector of parameters based on a boolean vector "changed"
# which shows which parameters have changed
function update!(qnet::QuartetNetwork,x::Vector{Float64}, ch::Vector{Bool})
    if(length(x) == length(ch))
        if(length(ch) == length(qnet.hasEdge))
            qnet.changed = false
            for(i in 1:length(ch))
                qnet.changed |= (ch[i] & qnet.hasEdge[i])
            end
            if(qnet.changed)
                i = 1
                j = 1
                for(e in qnet.edge)
                    if(e.istIdentifiable)
                        setLength!(e,x[qnet.indexht[i+qnet.numHybrids]])
                        i += 1
                    end
                    if(e.hybrid && !e.isMajor)
                        0 <= x[qnet.indexht[j]] <= 1 || error("new gamma value should be between 0,1: $(x[qnet.indexht[j]]).")
                        setGamma!(e,x[qnet.indexht[j]])
                        e.node[e.isChild1?1:2].hybrid || error("hybrid edge $(e.number) points at tree node.")
                        edges = hybridEdges(e.node[e.isChild1 ? 1 : 2],e)
                        length(edges) == 2 || error("strange here: node $(e.node[e.isChild1?1:2].number) should have 3 edges and it has $(length(edges)-1).")
                        if(edges[1].hybrid && edges[1].isMajor)
                            setGamma!(edges[1],1-x[qnet.indexht[j]])
                        elseif(edges[2].hybrid && edges[2].isMajor)
                            setGamma!(edges[2],1-x[qnet.indexht[j]])
                        else
                            error("strange hybrid node with only one hybrid edge $(e.number)")
                        end
                        j += 1
                    end
                end
            end
        else
            error("changed (length $(length(changed))) and qnet.hasEdge (length $(length(qnet.hasEdge))) should have same length")
        end
    else
        error("x (length $(length(x))) and changed $(length(changed)) should have the same length")
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
