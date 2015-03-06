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

# function that will add a hybridization with addHybridizationUpdate,
# if success=false, it will try to move the hybridization before
# declaring failure
function addHybridizationUpdateSmart!(net::HybridNetwork)
    success, hybrid, flag, nocycle, flag2, flag3 = addHybridizationUpdate!(net)
    if(!success)
        while(nocycle || !flag) #incycle failed
            deleteHybrid!(hybrid,net,true)
            success, hybrid, flag, nocycle, flag2, flag3 = addHybridizationUpdate!(net)
        end
        if(!flag2) #gammaz failed
            flag3, edgesRoot = updateContainRoot!(net,hybrid);
            if(flag3)
                success = moveOriginUpdateRepeat!(net,hybrid,true)
            else
                success = moveTargetUpdateRepeat!(net,hybrid,true)
            end
        else
            if(!flag3) #containRoot failed
                success,hybrid = changeDirectionUpdate!(net,hybrid)
            end
        end
        if(!success)
            deleteHybridizationUpdate!(net,hybrid)
        end
    end
    return success
end

# function to delete a hybrid, and then add a new hybrid:
# deleteHybridizationUpdate and addHybridizationUpdate,
# will keep trying to add until success or after N times
# returns success
function moveHybrid!(net::HybridNetwork, edge::Edge)
    edge.hybrid || error("edge $(edge.number) cannot be deleted because it is not hybrid")
    warn("moving hybrid for edge $(edge.number)")
    node = edge.node[edge.isChild1 ? 1 : 2]
    node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
    deleteHybridizationUpdate!(net, node)
    success = addHybridizationUpdateSmart!(net)
    N = 100
    i = 1
    while(!success && i<N)
        success = addHybridizationUpdateSmart!(net)
        i += 1
    end
    if(success)
        newf,newx = optBL!(net,d)
    end
    return success
end


# function to deal with h=0,1 case by:
# - change direction of hybrid edge, do optBL again
# - if failed, call moveHybrid
# input: net (network to be altered)
# returns success
function gammaZero!(net::HybridNetwork, edge::Edge)
    currTloglik = net.loglik
    edge.hybrid || error("edge $(edge.number) should be hybrid edge because it corresponds to a gamma (or gammaz) in net.ht")
    warn("gamma zero situation found for hybrid edge $(edge.number) with gamma $(edge.gamma)")
    node = edge.node[edge.isChild1 ? 1 : 2]
    node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
    success, newnode = changeDirectionUpdate!(net,node)
    if(success)
        newf,newx = optBL!(net,d)
        if(net.loglik <= currTloglik && !approxEq(edge.gamma,0.0) && !approxEq(edge.gamma,1.0))
            println("changing direction fixed the gamma zero situation")
            success2 = true
        else
            warn("changing direction does not fix the gamma zero situation, need to move hybrid")
            success,node = changeDirectionUpdate!(net,newnode)
            success || error("strange thing, changed direction and success, but lower loglik; want to undo changeDirection, and success=false! Hybrid node is $(node.number)")
            success2 = moveHybrid!(net,edge)
        end
    else
        warn("changing direction was not possible to fix the gamma zero situation (success=false), need to move hybrid")
        success2 = moveHybrid!(net,edge)
    end
    return success2
end

# function to check if h or t (in currT.ht) are 0 (or 1 for h)
# returns false if found hz=1.0, or if could not add new hybrid; true ow
function afterOptBL!(currT::HybridNetwork, d::DataCF)
    nh = currT.ht[1 : currT.numHybrids - currT.numBad]
    k = sum([e.istIdentifiable ? 1 : 0 for e in currT.edge])
    nt = currT.ht[currT.numHybrids - currT.numBad + 1 : currT.numHybrids - currT.numBad + k]
    nhz = currT.ht[currT.numHybrids - currT.numBad + k + 1 : length(currT.ht)]
    indh = currT.index[1 : currT.numHybrids - currT.numBad]
    indt = currT.index[currT.numHybrids - currT.numBad + 1 : currT.numHybrids - currT.numBad + k]
    indhz = currT.index[currT.numHybrids - currT.numBad + k + 1 : length(currT.ht)]
    flagh,flagt,flaghz = isValid(nh,nt,nhz)
    samenet = true
    successchange = true
    if(!flagh)
        for(i in 1:length(nh))
            if(approxEq(nh[i],0.0) || approxEq(nh[i],1.0))
                edge = currT.edge[indh[i]]
                approxEq(edge.gamma,nh[i]) || error("edge $(edge.number) gamma $(edge.gamma) should match the gamma in net.ht $(nh[i]) and it does not")
                successchange = gammaZero!(currT, edge)
                samenet = false
                break
            end
        end
    elseif(!flagt)
        for(i in 1:length(nt))
            if(approxEq(nt[i],0.0))
                edge = currT.edge[indt[i]]
                approxEq(edge.length,nt[i]) || error("edge $(edge.number) length $(edge.length) should match the length in net.ht $(nt[i]) and it does not")
                if(all([!n.hasHybEdge for n in edge.node]))
                    successchange = NNI!(currT,edge)
                    if(successchange)
                        newf,newx = optBL!(currT,d)
                        samenet = false
                        break
                    end
                end
            end
        end
    elseif(!flaghz)
        i = 1
        while(i <= length(nhz))
            if(approxEq(nhz[i],0.0))
                nodehz = currT.node[indhz[i]]
                approxEq(nodehz.gammaz,nhz[i]) || error("nodehz $(nodehz.number) gammaz $(nodehz.gammaz) should match the gammaz in net.ht $(nhz[i]) and it does not")
                edges = hybridEdges(nodehz)
                edges[1].hybrid || error("bad diamond I situation, node $(nodehz.number) has gammaz $(nodehz.gammaz) so should be linked to hybrid edge, but it is not")
                successchange = gammaZero!(currT,edges[1])
                samenet = false
                break
            elseif(approxEq(nhz[i],1.0))
                approxEq(nhz[i+1],0.0) || error("gammaz for node $(currT.node[indhz[i]].number) is $(nhz[i]) but the other gammaz is $(nhz[i+1]), the sum should be less than 1.0")
                nodehz = currT.node[indhz[i+1]]
                approxEq(nodehz.gammaz,nhz[i+1]) || error("nodehz $(nodehz.number) gammaz $(nodehz.gammaz) should match the gammaz in net.ht $(nhz[i+1]) and it does not")
                edges = hybridEdges(nodehz)
                edges[1].hybrid || error("bad diamond I situation, node $(nodehz.number) has gammaz $(nodehz.gammaz) so should be linked to hybrid edge, but it is not")
                successchange = gammaZero!(currT,edges[1])
                samenet = false
                break
            else
                if(approxEq(nhz[i+1],0.0))
                    nodehz = currT.node[indhz[i+1]]
                    approxEq(nodehz.gammaz,nhz[i+1]) || error("nodehz $(nodehz.number) gammaz $(nodehz.gammaz) should match the gammaz in net.ht $(nhz[i+1]) and it does not")
                    edges = hybridEdges(nodehz)
                    edges[1].hybrid || error("bad diamond I situation, node $(nodehz.number) has gammaz $(nodehz.gammaz) so should be linked to hybrid edge, but it is not")
                    successchange = gammaZero!(currT,edges[1])
                    samenet = false
                    break
                elseif(approxEq(nhz[i+1],1.0))
                    error("gammaz for node $(currT.node[indhz[i]].number) is $(nhz[i]) but the other gammaz is $(nhz[i+1]), the sum should be less than 1.0")
                end
            end
            i += 2
        end
    end
    return successchange,samenet,flagh,flagt,flaghz
end


# todo: recheck afterOptBL, debug in test_optBL, test_optTopLevelParts (do optBL and then afteroptBL)
# add afterOptBL into optBL, also add other minor changes (lower tol, etc)
# redo tests in optBL and optTopLevel


# numerical optimization of branch lengths given a network (or tree)
# and data (set of quartets with obsCF)
# using BOBYQA from NLopt package
function optBL!(net::HybridNetwork, d::DataCF)
    ht = parameters!(net); # branches/gammas to optimize: net.ht, net.numht
    extractQuartet!(net,d) # quartets are all updated: hasEdge, expCF, indexht
    k = length(net.ht)
    net.numBad >= 0 || error("network has negative number of bad hybrids")
    opt = NLopt.Opt(net.numBad == 0 ? :LN_BOBYQA : :LN_COBYLA,k) # :LD_MMA if use gradient, :LN_COBYLA for nonlinear/linear constrained optimization derivative-free, :LN_BOBYQA for bound constrained derivative-free
    # criterion based on prof Bates code
    NLopt.ftol_rel!(opt,1e-12) # relative criterion -12
    NLopt.ftol_abs!(opt,1e-12) # absolute critetion -8
    NLopt.xtol_abs!(opt,1e-10) # criterion on parameter value changes -10
    NLopt.lower_bounds!(opt, zeros(k))
    NLopt.upper_bounds!(opt,upper(net))
    count = 0
    function obj(x::Vector{Float64},g::Vector{Float64}) # added g::Vector{Float64} for gradient, ow error
        count += 1
        calculateExpCFAll!(d,x,net) # update qnet branches and calculate expCF
        update!(net,x) # update net.ht
        val = logPseudoLik(d)
        println("f_$count: $(round(val,5)), x: $(x)")
        return val
    end
    NLopt.min_objective!(opt,obj)
    if(net.numBad == 1)
        function inequalityGammaz(x::Vector{Float64},g::Vector{Float64})
            val = calculateIneqGammaz(x,net,1)
            return val
        end
        NLopt.inequality_constraint!(opt,inequalityGammaz,1e-8)
    elseif(net.numBad > 1)
        function inequalityGammaz(result::Vector{Float64},x::Vector{Float64},g::Matrix{Float64})
            for(i in 1:net.numBad)
                result[i] = calculateIneqGammaz(x,net,i)
            end
        end
        NLopt.inequality_constraint!(opt,inequalityGammaz)
    end
    fmin, xmin, ret = NLopt.optimize(opt,ht)
    println("got $(round(fmin,5)) at $(round(xmin,5)) after $(count) iterations (returned $(ret))")
    #println("net.ht is $(round(net.ht,5))")
    updateParameters!(net)
    updateLik!(net,fmin)
    return fmin,xmin
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
