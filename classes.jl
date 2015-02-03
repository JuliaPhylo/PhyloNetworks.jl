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

# -------------- heuristic search for topology -----------------------

move2int = Dict{Symbol,Int64}([:add=>1,:MVorigin=>2,:MVtarget=>3,:CHdir=>4,:delete=>5, :nni=>6])

function isTree(net::HybridNetwork)
    net.numHybrids == length(net.hybrid) || error("numHybrids does not match to length of net.hybrid")
    net.numHybrids != 0 || return true
    return false
end

# function to decide what next move to do when searching
# for topology that maximizes the P-loglik within the space of
# topologies with the same number of hybridizations
# possible moves: move origin/target, change direction hybrid edge, tree nni
# needs the network to know if it is a tree
function whichMove(net::HybridNetwork)
    if(isTree(net))
        return :nni
    else
        r = rand()
        if(r < 1/3)
            return :MVorigin
        elseif(r < 2/3)
            return :MVtarget
        else
            return :CHdir
        end
    end
end

#function to choose a hybrid node for the given moves
function chooseHybrid(net::HybridNetwork)
    !isTree(net) || error("net is a tree, cannot choose hybrid node")
    net.numHybrids > 1 || return net.hybrid[1]
    index1 = 0
    while(index1 == 0 || index1 > size(net.hybrid,1))
        index1 = iround(rand()*size(net.hybrid,1));
    end
    return net.hybrid[index1]
end

# function to propose a new topology given a move
# random = false uses the minor hybrid edge always
function proposedTop(move::Integer, T::HybridNetwork, random::Bool)
    1 <= move <= 6 || error("invalid move $(move)")
    newT = deepcopy(T)
    if(move == 1)
        success,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(newT)
    elseif(move == 2)
        node = chooseHybrid(newT)
        success,flag,nocycle,flag2 = moveOriginUpdate!(newT,node,random)
        flag3 = true
    elseif(move == 3)
        node = chooseHybrid(newT)
        success,flag,nocycle,flag2 = moveTargetUpdate!(newT,node,random) #fixit: add flag3
    elseif(move == 4)
        node = chooseHybrid(newT)
        success,flag2,flag3 = changeDirectionUpdate!(node,newT)
        flag = true
        nocycle = false
    elseif(move == 5)
        node = chooseHybrid(newT)
        deleteHybridizationUpdate!(newT,node,random)
        success = true
        flag = true
        flag2 = true
        flag3 = true
        nocycle = false
    elseif(move == 6)
        NNI!(newT)
    end
    !success || return newT # fixit: what if one of this fails? need to check that the code undoes what the move did
    error("new proposed topology failed for move $(move): success $(success), incycle $(flag), nocycle $(nocycle), updategammaz $(flag2), containroot $(flag3)")
end

proposedTop(move::Symbol, T::HybridNetwork, random::Bool) = proposedTop(try move2int[move] catch error("invalid move $(string(move))") end,T,random)

function optTopLevel(net0::HybridNetwork, epsilon::Float64, d::DataCF)
    epsilon > 0 || error("epsilon must be greater than zero: $(epsilon)")
    delta = epsilon - 1
    currT = net0
    currloglik,currxmin = optBL(net0,d) # do we want to updateParameters in net0?
    while(delta < epsilon)
        success = false
        while(!success) # will propose a new move after failure, not the same move
            move = whichMove(currT)
            try
                newT = proposedTop(move,currT,true)
            catch
                success = false
            end
            newT = proposedTop(move,currT,true)
            success = true
        end
        newloglik, newxmin = optBL(newT,d)
        if(newloglik > currloglik)
            delta = abs(newloglik - currloglik)
            currT = newT
            currloglik = newloglik
            currxmin = newxmin
        end
    end
    updateParameters!(newT)
    updateLik!(newT,newloglik)
    return newT
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
