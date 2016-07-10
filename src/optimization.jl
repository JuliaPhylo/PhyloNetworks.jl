# functions for numerical/heuristic optimization
# originally in functions.jl
# Claudia March 2015


const move2int = Dict{Symbol,Int}(:add=>1,:MVorigin=>2,:MVtarget=>3,:CHdir=>4,:delete=>5, :nni=>6)
const int2move = Dict{Int,Symbol}([move2int[k]=>k for k in keys(move2int)])
# if changes are made here, make the same in the docstring for snaq! below.
const fAbs = 1e-6 #1e-10 prof Bates, 1e-6
const fRel = 1e-5 # 1e-12 prof Bates, 1e-5
const xAbs = 1e-4 # 0.001 in phylonet, 1e-10 prof Bates, 1e-4
const xRel = 1e-3 # 0.01 in phylonet, 1e-10 prof Bates, 1e-3
const numFails = 100 # number of failed proposals allowed before stopping the procedure (like phylonet)
const numMoves = Int[] #empty to be calculated inside based on coupon's collector
const multiplier = 10000 # for loglik absolute tol (multiplier*fAbs)

# ---------------------- branch length optimization ---------------------------------

# function to get the branch lengths/gammas to optimize for a given network
# warning: order of parameters (h,t,gammaz)
# updates net.numht also with the number of hybrid nodes and number of identifiable edges (n2,n,hzn)
function parameters(net::Network)
    t = Float64[]
    h = Float64[]
    n = Int[]
    n2 = Int[]
    hz = Float64[]
    hzn = Int[]
    indxt = Int[]
    indxh = Int[]
    indxhz = Int[]
    for(e in net.edge)
        if(e.istIdentifiable)
            push!(t,e.length)
            push!(n,e.number)
            push!(indxt, getIndex(e,net))
        end
        if(e.hybrid && !e.isMajor)
            node = e.node[e.isChild1 ? 1 : 2]
            node.hybrid || error("strange thing, hybrid edge $(e.number) pointing at tree node $(node.number)")
            if(!node.isBadDiamondI)
                push!(h,e.gamma)
                push!(n2,e.number)
                push!(indxh, getIndex(e,net))
            else
                if(node.isBadDiamondI)
                    edges = hybridEdges(node)
                    push!(hz,getOtherNode(edges[1],node).gammaz)
                    push!(hz,getOtherNode(edges[2],node).gammaz)
                    push!(hzn,parse(Int,string(string(node.number),"1")))
                    push!(hzn,parse(Int,string(string(node.number),"2")))
                    push!(indxhz,getIndex(getOtherNode(edges[1],node),net))
                    push!(indxhz,getIndex(getOtherNode(edges[2],node),net))
                end
            end
        end
    end
    size(t,1) == 0 ? warn("net does not have identifiable branch lengths") : nothing
    return vcat(h,t,hz),vcat(n2,n,hzn),vcat(indxh,indxt,indxhz)
end

function parameters!(net::Network)
    #warn("deleting net.ht,net.numht and updating with current edge lengths (numbers)")
    net.ht,net.numht,net.index = parameters(net)
    return net.ht
end


# function to update qnet.indexht,qnet.index based on net.numht
# warning: assumes net.numht is updated already with parameters!(net)
function parameters!(qnet::QuartetNetwork, net::HybridNetwork)
    global DEBUG
    size(net.numht,1) > 0 || error("net.numht not correctly updated, need to run parameters first")
    if(DEBUG)
        size(qnet.indexht,1) == 0 ||  println("deleting qnet.indexht to replace with info in net")
    end
    nh = net.numht[1 : net.numHybrids - net.numBad]
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    nt = net.numht[net.numHybrids - net.numBad + 1 : net.numHybrids - net.numBad + k]
    nhz = net.numht[net.numHybrids - net.numBad + k + 1 : length(net.numht)]
    qnh = Int[]
    qnt = Int[]
    qnhz = Int[]
    qindxh = Int[]
    qindxt = Int[]
    qindxhz = Int[]
    if(qnet.numHybrids == 1 && qnet.hybrid[1].isBadDiamondI)
        ind1 = parse(Int,string(string(qnet.hybrid[1].number),"1"))
        ind2 = parse(Int,string(string(qnet.hybrid[1].number),"2"))
        i = getIndex(ind1,nhz)
        edges = hybridEdges(qnet.hybrid[1])
        push!(qnhz,i+net.numHybrids-net.numBad+k)
        push!(qnhz,i+1+net.numHybrids-net.numBad+k)
        push!(qindxhz,getIndex(getOtherNode(edges[1],qnet.hybrid[1]),qnet))
        push!(qindxhz,getIndex(getOtherNode(edges[2],qnet.hybrid[1]),qnet))
    else
        found = false
        for(n in qnet.hybrid)
            if(n.isBadDiamondI)
                ind1 = parse(Int,string(string(n.number),"1"))
                ind2 = parse(Int,string(string(n.number),"2"))
                i = getIndex(ind1,nhz)
                edges = hybridEdges(n)
                push!(qnhz,i+net.numHybrids-net.numBad+k)
                push!(qnhz,i+1+net.numHybrids-net.numBad+k)
                push!(qindxhz,getIndex(getOtherNode(edges[1],n),qnet))
                push!(qindxhz,getIndex(getOtherNode(edges[2],n),qnet))
                found = true
                break
            end
        end
        if(!found)
            all((n -> !n.isBadDiamondI),qnet.hybrid) || error("cannot have bad diamond I hybrid nodes in this qnet, case dealt separately before")
            for(e in qnet.edge)
                if(e.istIdentifiable)
                    try
                        getIndex(e.number,nt)
                    catch
                        error("identifiable edge $(e.number) in qnet not found in net")
                    end
                    push!(qnt, getIndex(e.number,nt) + net.numHybrids - net.numBad)
                    push!(qindxt, getIndex(e,qnet))
                end
                if(!e.istIdentifiable && all((n->!n.leaf),e.node) && !e.hybrid && e.fromBadDiamondI) # tree edge not identifiable but internal with length!=0 (not bad diamII nor bad triangle)
                    try
                        getIndex(e.number,nhz)
                    catch
                        error("internal edge $(e.number) corresponding to gammaz in qnet not found in net.ht")
                    end
                    push!(qnhz, getIndex(e.number,nhz) + net.numHybrids - net.numBad + k)
                    push!(qindxhz, getIndex(e,qnet))
                end
                if(e.hybrid && !e.isMajor)
                    node = e.node[e.isChild1 ? 1 : 2]
                    node.hybrid || error("strange hybrid edge $(e.number) poiting to tree node $(node.number)")
                    found = true
                    try
                        getIndex(e.number,nh)
                    catch
                        found = false
                    end
                    found  ? push!(qnh, getIndex(e.number,nh)) : nothing
                    found ? push!(qindxh, getIndex(e,qnet)) : nothing
                end
            end # for qnet.edge
        end # not found
    end
    qnet.indexht = vcat(qnh,qnt,qnhz)
    qnet.index = vcat(qindxh,qindxt,qindxhz)
    length(qnet.indexht) == length(qnet.index) || error("strange in setting qnet.indexht and qnet.index, they do not have same length")
end


# function to compare a vector of parameters with the current vector in net.ht
# to know which parameters were changed
function changed(net::HybridNetwork, x::Vector{Float64})
    if(length(net.ht) == length(x))
        #DEBUG && println("inside changed with net.ht $(net.ht) and x $(x)")
        return [!approxEq(net.ht[i],x[i]) for i in 1:length(x)]
    else
        error("net.ht (length $(length(net.ht))) and vector x (length $(length(x))) need to have same length")
    end
end


# function to update a QuartetNetwork for a given
# vector of parameters based on a boolean vector "changed"
# which shows which parameters have changed
function update!(qnet::QuartetNetwork,x::Vector{Float64}, net::HybridNetwork)
    global DEBUG
    ch = changed(net,x)
    length(x) == length(ch) || error("x (length $(length(x))) and changed $(length(changed)) should have the same length")
    length(ch) == length(qnet.hasEdge) || error("changed (length $(length(changed))) and qnet.hasEdge (length $(length(qnet.hasEdge))) should have same length")
    qnet.changed = false
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    for(i in 1:length(ch))
        qnet.changed |= (ch[i] & qnet.hasEdge[i])
    end
    #DEBUGC && println("inside update!, qnet.changed is $(qnet.changed), ch $(ch) and qnet.hasEdge $(qnet.hasEdge), $(qnet.quartetTaxon), numHyb $(qnet.numHybrids)")
    if(qnet.changed)
        if(any([n.isBadDiamondI for n in qnet.hybrid])) # qnet.indexht is only two values: gammaz1,gammaz2 #FIXIT: this could crash if hybrid for bad diamond should disappear after cleaning qnet
            DEBUG && println("it is inside update! and identifies that ht changed and it is inside the bad diamond I case")
            length(qnet.indexht) == 2 || error("strange qnet from bad diamond I with hybrid node, it should have only 2 elements: gammaz1,gammaz2, not $(length(qnet.indexht))")
            for(i in 1:2)
                0 <= x[qnet.indexht[i]] <= 1 || error("new gammaz value should be between 0,1: $(x[qnet.indexht[i]]).")
                if(DEBUG)
                    x[qnet.indexht[1]] + x[qnet.indexht[2]] <= 1 || warn("new gammaz should add to less than 1: $(x[qnet.indexht[1]] + x[qnet.indexht[2]])")
                end
                qnet.node[qnet.index[i]].gammaz = x[qnet.indexht[i]]
            end
        else
            for(i in 1:length(qnet.indexht))
                if(qnet.indexht[i] <= net.numHybrids - net.numBad)
                    0 <= x[qnet.indexht[i]] <= 1 || error("new gamma value should be between 0,1: $(x[qnet.indexht[i]]).")
                    qnet.edge[qnet.index[i]].hybrid || error("something odd here, optimizing gamma for tree edge $(qnet.edge[qnet.index[i]].number)")
                    setGamma!(qnet.edge[qnet.index[i]],x[qnet.indexht[i]], true)
                elseif(qnet.indexht[i] <= net.numHybrids - net.numBad + k)
                    setLength!(qnet.edge[qnet.index[i]],x[qnet.indexht[i]])
                else
                    DEBUGC && println("updating qnet parameters, found gammaz case when hybridization has been removed")
                    0 <= x[qnet.indexht[i]] <= 1 || error("new gammaz value should be between 0,1: $(x[qnet.indexht[i]]).")
                    #x[qnet.indexht[i]] + x[qnet.indexht[i]+1] <= 1 || warn("new gammaz value should add to less than 1: $(x[qnet.indexht[i]])  $(x[qnet.indexht[i]+1]).")
                    if(approxEq(x[qnet.indexht[i]],1.0))
                        setLength!(qnet.edge[qnet.index[i]],10.0)
                    else
                        setLength!(qnet.edge[qnet.index[i]],-log(1-x[qnet.indexht[i]]))
                    end
                end
            end
        end
    end
end

# function to update the branch lengths/gammas for a network
# warning: order of parameters (h,t)
function update!(net::HybridNetwork, x::Vector{Float64})
    if(length(x) == length(net.ht))
        net.ht = deepcopy(x) # to avoid linking them
    else
        error("net.ht (length $(length(net.ht))) and x (length $(length(x))) must have the same length")
    end
end

# function to update the branch lengths and gammas in a network after
# the optimization
# warning: optBL need to be run before, note that
# xmin will not be in net.ht (net.ht is one step before)
function updateParameters!(net::HybridNetwork, xmin::Vector{Float64})
    length(xmin) == length(net.ht) || error("xmin vector should have same length as net.ht $(length(net.ht)), not $(length(xmin))")
    net.ht = xmin
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    for(i in 1:length(net.ht))
        if(i <= net.numHybrids - net.numBad)
            0 <= net.ht[i] <= 1 || error("new gamma value should be between 0,1: $(net.ht[i]).")
            net.edge[net.index[i]].hybrid || error("something odd here, optimizing gamma for tree edge $(net.edge[net.index[i]].number)")
            setGamma!(net.edge[net.index[i]],net.ht[i], true)
        elseif(i <= net.numHybrids - net.numBad + k)
            setLength!(net.edge[net.index[i]],net.ht[i])
        else
            0 <= net.ht[i] <= 1 || error("new gammaz value should be between 0,1: $(net.ht[i]).")
            net.node[net.index[i]].gammaz = net.ht[i]
        end
    end
end

# function to update the attribute net.loglik
function updateLik!(net::HybridNetwork, l::Float64)
    net.loglik = l
end

# function for the upper bound of ht
function upper(net::HybridNetwork)
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    return vcat(ones(net.numHybrids-net.numBad), repeat([10],inner=[k]),
                ones(length(net.ht)-k-net.numHybrids+net.numBad))
end

# function to calculate the inequality gammaz1+gammaz2 <= 1
function calculateIneqGammaz(x::Vector{Float64}, net::HybridNetwork, ind::Integer, verbose::Bool)
    global DEBUG
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    hz = x[net.numHybrids - net.numBad + k + 1 : length(x)]
    if(verbose || DEBUG)
        println("enters calculateIneqGammaz with hz $(hz), and hz[ind*2] + hz[ind*2-1] - 1 = $(hz[ind*2] + hz[ind*2-1] - 1)")
    end
    hz[ind*2] + hz[ind*2-1] - 1
end

# numerical optimization of branch lengths given a network (or tree)
# and data (set of quartets with obsCF)
# using BOBYQA from NLopt package
# warning: this function assumes that the network has all the good attributes set. It will not be efficient to re-read the network inside
# to guarantee all the correct attributes, because this is done over and over inside snaq
# also, net is modified inside to set its attribute net.loglik equal to the min
function optBL!(net::HybridNetwork, d::DataCF, verbose::Bool, ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64)
    global DEBUG
    (ftolRel > 0 && ftolAbs > 0 && xtolAbs > 0 && xtolRel > 0) || error("tolerances have to be positive, ftol (rel,abs), xtol (rel,abs): $([ftolRel, ftolAbs, xtolRel, xtolAbs])")
    (DEBUG || verbose) && println("OPTBL: begin branch lengths and gammas optimization, ftolAbs $(ftolAbs), ftolRel $(ftolRel), xtolAbs $(xtolAbs), xtolRel $(xtolRel)")
    ht = parameters!(net); # branches/gammas to optimize: net.ht, net.numht
    extractQuartet!(net,d) # quartets are all updated: hasEdge, expCF, indexht
    k = length(net.ht)
    net.numBad >= 0 || error("network has negative number of bad hybrids")
    #opt = NLopt.Opt(net.numBad == 0 ? :LN_BOBYQA : :LN_COBYLA,k) # :LD_MMA if use gradient, :LN_COBYLA for nonlinear/linear constrained optimization derivative-free, :LN_BOBYQA for bound constrained derivative-free
    opt = NLopt.Opt(:LN_BOBYQA,k) # :LD_MMA if use gradient, :LN_COBYLA for nonlinear/linear constrained optimization derivative-free, :LN_BOBYQA for bound constrained derivative-free
    # criterion based on prof Bates code
    NLopt.ftol_rel!(opt,ftolRel) # relative criterion -12
    NLopt.ftol_abs!(opt,ftolAbs) # absolute critetion -8, later changed to -10
    NLopt.xtol_rel!(opt,xtolRel) # criterion on parameter value changes -10
    NLopt.xtol_abs!(opt,xtolAbs) # criterion on parameter value changes -10
    NLopt.maxeval!(opt,1000) # max number of iterations
    NLopt.lower_bounds!(opt, zeros(k))
    NLopt.upper_bounds!(opt,upper(net))
    count = 0
    function obj(x::Vector{Float64},g::Vector{Float64}) # added g::Vector{Float64} for gradient, ow error
        if(DEBUG || verbose) #|| net.numBad > 0) #we want to see what happens with bad diamond I
            println("inside obj with x $(x)")
        end
        count += 1
        calculateExpCFAll!(d,x,net) # update qnet branches and calculate expCF
        update!(net,x) # update net.ht
        val = logPseudoLik(d)
        if(DEBUG || verbose) #|| net.numBad > 0)#we want to see what happens with bad diamond I
            println("f_$count: $(round(val,5)), x: $(x)")
        end
        return val
    end
    NLopt.min_objective!(opt,obj)
    ## if(net.numBad == 1)
    ##     function inequalityGammaz(x::Vector{Float64},g::Vector{Float64})
    ##         val = calculateIneqGammaz(x,net,1,verbose)
    ##         return val
    ##     end
    ##     NLopt.inequality_constraint!(opt,inequalityGammaz)
    ## elseif(net.numBad > 1)
    ##     function inequalityGammaz(result::Vector{Float64},x::Vector{Float64},g::Matrix{Float64})
    ##         i = 1
    ##         while(i < net.numBad)
    ##             result[i] = calculateIneqGammaz(x,net,i,verbose)
    ##             i += 2
    ##         end
    ##     end
    ##     NLopt.inequality_constraint!(opt,inequalityGammaz)
    ## end
    (DEBUG || verbose) && println("OPTBL: starting point $(ht)")
    fmin, xmin, ret = NLopt.optimize(opt,ht)
    (DEBUG || verbose) && println("got $(round(fmin,5)) at $(round(xmin,5)) after $(count) iterations (returned $(ret))")
    updateParameters!(net,xmin)
    net.loglik = fmin
    #return fmin,xmin
end

optBL!(net::HybridNetwork, d::DataCF) = optBL!(net, d, false, fRel, fAbs, xRel, xAbs)
optBL!(net::HybridNetwork, d::DataCF, verbose::Bool) = optBL!(net, d,verbose, fRel, fAbs, xRel, xAbs)
optBL!(net::HybridNetwork, d::DataCF, ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64) = optBL!(net, d, false, ftolRel, ftolAbs, xtolRel, xtolAbs)

# rename optBL for a more user-friendly name
"""
`topologyMaxQPseudolik!(net::HybridNetwork, d::DataCF)`

Estimate the branch lengths and inheritance probabilities (γ's) for a given network topology.
The network is *not* modified, only the object `d` is, with updated expected concordance factors.

Ouput: new network, with optimized parameters (branch lengths and gammas).
The maximized quartet pseudo-deviance is the negative log pseudo-likelihood,
up to an additive constant, such that a perfect fit corresponds to a deviance of 0.0.
This is also an attribute of the network, which can be accessed with `net.loglik`.

Optional arguments (default value):

- verbose (false): if true, information on the numerical optimization is printed to screen
- ftolRel (1e-5), ftolAbs (1e-6), xtolRel (1e-3), xtolAbs (1e-4):
  absolute and relative tolerance values for the pseudo-deviance function
  and the parameters
"""
function topologyMaxQPseudolik!(net::HybridNetwork, d::DataCF; verbose=false::Bool, ftolRel=fRel::Float64, ftolAbs=fAbs::Float64, xtolRel=xRel::Float64, xtolAbs=xAbs::Float64)
    net = readTopologyUpdate(writeTopologyLevel1(net)) # update everything for level 1
    try
        checkNet(net)
    catch err
        err.msg = "starting topology not a level 1 network:\n" * err.msg
        rethrow(err)
    end
    optBL!(net, d, verbose, ftolRel, ftolAbs, xtolRel,xtolAbs)
    if(net.numBad > 0) # to keep gammaz info in parenthetical description of bad diamond I
        for(n in net.hybrid)
            setGammaBLfromGammaz!(n,net) # get t and γ that are compatible with estimated gammaz values
        end
    end
    setNonIdBL!(net)
    return net
end

# function to delete a hybrid, and then add a new hybrid:
# deleteHybridizationUpdate and addHybridizationUpdate,
# closeN=true will try move origin/target, if false, will delete/add new hybrid
# default is closeN =true
# origin=true, moves origin, if false, moves target. option added to control not keep coming
# to the same network over and over
# returns success
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
function moveHybrid!(net::HybridNetwork, edge::Edge, closeN ::Bool, origin::Bool,N::Integer, movesgamma::Vector{Int})
    global DEBUG
    edge.hybrid || error("edge $(edge.number) cannot be deleted because it is not hybrid")
    node = edge.node[edge.isChild1 ? 1 : 2];
    node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
    DEBUG && println("MOVE: moving hybrid for edge $(edge.number)")
    if(closeN )
        if(origin)
            movesgamma[2] += 1
            success = moveOriginUpdateRepeat!(net,node,true)
            movesgamma[8] += success ? 1 : 0
        else
            movesgamma[3] += 1
            success = moveTargetUpdateRepeat!(net,node,true)
            movesgamma[9] += success ? 1 : 0
        end
    else
        movesgamma[1] += 1
        deleteHybridizationUpdate!(net, node, true, true)
        success = addHybridizationUpdateSmart!(net,N)
        movesgamma[7] += success ? 1 : 0
    end
    return success
end

# function to deal with h=0,1 case by:
# - change direction of hybrid edge, do optBL again
# - if failed, call moveHybrid
# input: net (network to be altered)
# closeN =true will try move origin/target, if false, will delete/add new hybrid
# origin=true, moves origin, if false, moves target. option added to control not keep coming
# to the same network over and over
# returns success
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
function gammaZero!(net::HybridNetwork, d::DataCF, edge::Edge, closeN ::Bool, origin::Bool, N::Integer, movesgamma::Vector{Int})
    global CHECKNET, DEBUG
    currTloglik = net.loglik
    edge.hybrid || error("edge $(edge.number) should be hybrid edge because it corresponds to a gamma (or gammaz) in net.ht")
    DEBUG && println("gamma zero situation found for hybrid edge $(edge.number) with gamma $(edge.gamma)")
    node = edge.node[edge.isChild1 ? 1 : 2];
    node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
    success = changeDirectionUpdate!(net,node) #changes dir of minor
    movesgamma[4] += 1
    if(success)
        DEBUG && printEverything(net)
        CHECKNET && checkNet(net)
        movesgamma[10] += 1
        optBL!(net,d,false)
        flags = isValid(net)
        if(net.loglik <= currTloglik && flags[1] && flags[3])
            DEBUG && println("changing direction fixed the gamma zero situation")
            success2 = true
        else
            DEBUG && println("changing direction does not fix the gamma zero situation, need to undo change direction and move hybrid")
            success = changeDirectionUpdate!(net,node)
            success || error("strange thing, changed direction and success, but lower loglik; want to undo changeDirection, and success=false! Hybrid node is $(node.number)")
            DEBUG && printEverything(net)
            CHECKNET && checkNet(net)
            success2 = moveHybrid!(net,edge,closeN ,origin,N, movesgamma)
        end
    else
        DEBUG && println("changing direction was not possible to fix the gamma zero situation (success=false), need to move hybrid")
        DEBUG && printEverything(net)
        CHECKNET && checkNet(net)
        success2 = moveHybrid!(net,edge,closeN ,origin,N,movesgamma)
    end
    return success2
end

# function to check if h or t (in currT.ht) are 0 (or 1 for h)
# closeN =true will try move origin/target, if false, will delete/add new hybrid
# origin=true, moves origin, if false, moves target. option added to control not keep coming
# to the same network over and over
# returns successchange=false if could not add new hybrid; true ow
# returns successchange,flagh,flagt,flaghz
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
function afterOptBL!(currT::HybridNetwork, d::DataCF,closeN ::Bool, origin::Bool,verbose::Bool, N::Integer, movesgamma::Vector{Int})
    global CHECKNET, DEBUG
    !isTree(currT) || return false,true,true,true
    nh = currT.ht[1 : currT.numHybrids - currT.numBad]
    k = sum([e.istIdentifiable ? 1 : 0 for e in currT.edge])
    nt = currT.ht[currT.numHybrids - currT.numBad + 1 : currT.numHybrids - currT.numBad + k]
    nhz = currT.ht[currT.numHybrids - currT.numBad + k + 1 : length(currT.ht)]
    indh = currT.index[1 : currT.numHybrids - currT.numBad]
    indt = currT.index[currT.numHybrids - currT.numBad + 1 : currT.numHybrids - currT.numBad + k]
    indhz = currT.index[currT.numHybrids - currT.numBad + k + 1 : length(currT.ht)]
    flagh,flagt,flaghz = isValid(nh,nt,nhz)
    !reduce(&,[flagh,flagt,flaghz]) || return false,true,true,true
    DEBUG && println("begins afterOptBL because of conflicts: flagh,flagt,flaghz=$([flagh,flagt,flaghz])")
    successchange = true
    if(!flagh)
        for(i in 1:length(nh))
            if(approxEq(nh[i],0.0) || approxEq(nh[i],1.0))
                edge = currT.edge[indh[i]]
                approxEq(edge.gamma,nh[i]) || error("edge $(edge.number) gamma $(edge.gamma) should match the gamma in net.ht $(nh[i]) and it does not")
                successchange = gammaZero!(currT, d,edge,closeN ,origin,N,movesgamma)
                !successchange || break
            end
        end
    elseif(!flagt)
        for(i in 1:length(nt))
            if(approxEq(nt[i],0.0))
                edge = currT.edge[indt[i]]
                approxEq(edge.length,nt[i]) || error("edge $(edge.number) length $(edge.length) should match the length in net.ht $(nt[i]) and it does not")
                if(all((n->!n.hasHybEdge), edge.node))
                    movesgamma[6] += 1
                    successchange = NNI!(currT,edge)
                    movesgamma[12] += successchange ? 1 : 0
                    !successchange || break
                else
                    DEBUG && println("MOVE: need NNI on edge $(edge.number) because length is $(edge.length), but edge is attached to nodes that have hybrid edges")
                    successchange = false
                end
            end
        end
    elseif(!flaghz)
        currT.numBad > 0 || error("not a bad diamond I and flaghz is $(flaghz), should be true")
        i = 1
        while(i <= length(nhz))
            if(approxEq(nhz[i],0.0))
                nodehz = currT.node[indhz[i]]
                approxEq(nodehz.gammaz,nhz[i]) || error("nodehz $(nodehz.number) gammaz $(nodehz.gammaz) should match the gammaz in net.ht $(nhz[i]) and it does not")
                edges = hybridEdges(nodehz)
                edges[1].hybrid || error("bad diamond I situation, node $(nodehz.number) has gammaz $(nodehz.gammaz) so should be linked to hybrid edge, but it is not")
                successchange = gammaZero!(currT,d,edges[1],closeN ,origin,N,movesgamma)
                break
            elseif(approxEq(nhz[i],1.0))
                approxEq(nhz[i+1],0.0) || error("gammaz for node $(currT.node[indhz[i]].number) is $(nhz[i]) but the other gammaz is $(nhz[i+1]), the sum should be less than 1.0 (afteroptbl)")
                nodehz = currT.node[indhz[i+1]]
                approxEq(nodehz.gammaz,nhz[i+1]) || error("nodehz $(nodehz.number) gammaz $(nodehz.gammaz) should match the gammaz in net.ht $(nhz[i+1]) and it does not")
                edges = hybridEdges(nodehz)
                edges[1].hybrid || error("bad diamond I situation, node $(nodehz.number) has gammaz $(nodehz.gammaz) so should be linked to hybrid edge, but it is not")
                successchange = gammaZero!(currT,d,edges[1],closeN ,origin,N,movesgamma)
                break
            else
                if(approxEq(nhz[i+1],0.0))
                    nodehz = currT.node[indhz[i+1]];
                    approxEq(nodehz.gammaz,nhz[i+1]) || error("nodehz $(nodehz.number) gammaz $(nodehz.gammaz) should match the gammaz in net.ht $(nhz[i+1]) and it does not")
                    edges = hybridEdges(nodehz);
                    edges[1].hybrid || error("bad diamond I situation, node $(nodehz.number) has gammaz $(nodehz.gammaz) so should be linked to hybrid edge, but it is not")
                    successchange = gammaZero!(currT,d,edges[1],closeN ,origin,N,movesgamma)
                    break
                elseif(approxEq(nhz[i+1],1.0))
                    error("gammaz for node $(currT.node[indhz[i]].number) is $(nhz[i]) but the other gammaz is $(nhz[i+1]), the sum should be less than 1.0 (afteroptbl)")
                end
            end
            i += 2
        end
    end
    if(successchange)
        DEBUG && println("afterOptBL SUCCESSFUL change, need to run again to see if new topology is valid")
        DEBUG && printEverything(currT)
        CHECKNET && checkNet(currT)
        optBL!(currT,d,verbose)
        flagh,flagt,flaghz = isValid(currT)
        DEBUG && println("new flags: flagh,flagt,flaghz $([flagh,flagt,flaghz])")
    end
    return successchange,flagh,flagt,flaghz
end

# function to repeat afterOptBL every time after changing something
# N: number of times it will delete/add hybrid if closeN =false
# origin=true, moves origin, if false, moves target. option added to control not keep coming
# to the same network over and over
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
function afterOptBLRepeat!(currT::HybridNetwork, d::DataCF, N::Integer,closeN ::Bool, origin::Bool,
                           verbose::Bool, movesgamma::Vector{Int})
    global DEBUG
    success,flagh,flagt,flaghz = afterOptBL!(currT,d,closeN ,origin,verbose,N, movesgamma)
    DEBUG && println("inside afterOptBLRepeat, after afterOptBL once, we get: success, flags: $([success,flagh,flagt,flaghz])")
    if(!closeN )
        DEBUG && println("closeN  is $(closeN ), and gets inside the loop for repeating afterOptBL")
        i = 1
        while(success && !reduce(&,[flagh,flagt,flaghz]) && i<N) #fixit: do we want to check loglik in this step also?
            success,flagh,flagt,flaghz = afterOptBL!(currT,d,closeN ,origin,verbose,N,movesgamma)
            i += 1
        end
        if(DEBUG)
            i < N || println("tried afterOptBL $(i) times")
        end
    end
    return success,flagh,flagt,flaghz
end

# function to check if we have to keep checking loglik
# returns true if we have to stop because loglik not changing much anymore, or it is close to 0.0 already
function stopLoglik(currloglik::Float64, newloglik::Float64, ftolAbs::Float64, M::Number)
    return (abs(currloglik-newloglik) <= M*ftolAbs) || (newloglik <= M*ftolAbs)
end

# function to repeat afterOptBL every time after changing something
# closeN =true will try move origin/target, if false, will delete/add new hybrid
# default is closeN =true
# returns new approved currT (no gammas=0.0)
# N: number of times failures of accepting loglik is allowed, M: multiplier for loglik tolerance (M*ftolAbs)
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
function afterOptBLAll!(currT::HybridNetwork, d::DataCF, N::Integer,closeN ::Bool, M::Number, ftolAbs::Float64, verbose::Bool, movesgamma::Vector{Int},ftolRel::Float64, xtolRel::Float64, xtolAbs::Float64)
    global DEBUG
    DEBUG && println("afterOptBLAll: checking if currT has gamma(z) = 0.0(1.0): currT.ht $(currT.ht)")
    currloglik = currT.loglik
    currT.blacklist = Int[];
    origin = (rand() > 0.5) #true=moveOrigin, false=moveTarget
    startover = true
    tries = 0
    N2 = N > 10 ? N/10 : 1 #num of failures of badlik around a gamma=0.0, t=0.0
    while(startover && tries < N)
        tries += 1
        DEBUG && println("inside afterOptBLALL: number of tries $(tries) out of $(N) possible")
        badliks = 0
        if(currT.loglik < M*ftolAbs) #curr loglik already close to 0.0
            startover = false
        else
            backCurrT0 = false
            while(badliks < N2) #will try a few options around currT
                DEBUG && println("tried $(badliks) bad likelihood options so far out of $(N2)")
                currT0 = deepcopy(currT)
                origin = !origin #to guarantee not going back to previous topology
                success,flagh,flagt,flaghz = afterOptBLRepeat!(currT,d,N,closeN ,origin,verbose,movesgamma)
                all((e->!(e.hybrid && e.inCycle == -1)), currT.edge) || error("found hybrid edge with inCycle == -1")
                DEBUG && println("inside afterOptBLAll, after afterOptBLRepeat once we get: success, flags: $([success,flagh,flagt,flaghz])")
                if(!success) #tried to change something but failed
                    DEBUG && println("did not change anything inside afterOptBL: could be nothing needed change or tried but couldn't anymore. flagh, flagt, flaghz = $([flagh,flagt,flaghz])")
                    if(reduce(&,[flagh,flagt,flaghz])) #currT was ok
                        startover = false
                    elseif(!flagh || !flaghz) #currT was bad but could not change it, need to go down a level
                        !isTree(currT) || error("afterOptBL should not give reject=true for a tree")
                        DEBUG && println("current topology has numerical parameters that are not valid: gamma=0(1), gammaz=0(1); need to move down a level h-1")
                        moveDownLevel!(currT)
                        optBL!(currT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
                        startover = true
                    elseif(!flagt)
                        startover = false
                    end
                else #changed something
                    DEBUG && println("changed something inside afterOptBL: flagh, flagt, flaghz = $([flagh,flagt,flaghz]). oldloglik $(currloglik), newloglik $(currT.loglik)")
                    if(currT.loglik > currloglik) #|| abs(currT.loglik-currloglik) <= M*ftolAbs) #fixit: allowed like this because of changeDir that does not change much the lik but can fix h=0
                        DEBUG && println("worse likelihood, back to currT")
                        startover = true
                        backCurrT0 = true
                    else
                        DEBUG && println("better likelihood, jump to new topology and startover")
                        backCurrT0 = false
                        movesgamma[13] += 1
                        if(reduce(&,[flagh,flagt,flaghz]))
                            startover = false
                        else
                            currloglik = currT.loglik
                            startover = true
                        end
                    end
                end
                if(backCurrT0)
                    currT = currT0
                    startover = true
                    badliks += 1
                else
                    badliks = N+1
                end
            end
            if(backCurrT0) # leaves while for failed loglik
                DEBUG && println("tried to fix gamma zero situation for $(badliks) times and could not")
                flagh,flagt,flaghz = isValid(currT)
                if(!flagh || !flaghz)
                    !isTree(currT) || error("afterOptBL should not give reject=true for a tree")
                    DEBUG && println("current topology has numerical parameters that are not valid: gamma=0(1), t=0, gammaz=0(1); need to move down a level h-1")
                    movesgamma[5] += 1
                    movesgamma[11] += 1
                    moveDownLevel!(currT)
                    optBL!(currT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
                    startover = true
                else
                    DEBUG && println("the only problem were lengths equal to zero, so we will keep them")
                    startover = false
                end
            end
        end
    end
    if(tries >= N)
        DEBUG && println("afterOptBLAll ended because it tried $(tries) times with startover $(startover)")
        DEBUG && writeTopologyLevel1(currT,true)
        flagh,flagt,flaghz = isValid(currT)
        if(!flagh || !flaghz)
            DEBUG && println("gammaz zero situation still in currT, need to move down one level to h-1")
            moveDownLevel!(currT)
            DEBUG && printEdges(currT)
            DEBUG && printPartitions(currT)
            #printNodes(currT)
            DEBUG && println(writeTopologyLevel1(currT,true))
            optBL!(currT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
        end
    end
    currT.blacklist = Int[];
    return currT
end


# -------------- heuristic search for topology -----------------------

function isTree(net::HybridNetwork)
    net.numHybrids == length(net.hybrid) || error("numHybrids does not match to length of net.hybrid")
    net.numHybrids != 0 || return true
    return false
end

# function to adjust the weight of addHybrid if net is in a much lower layer
# net.numHybrids<<hmax
# takes as input the vector of weights for each move (add,mvorigin,mvtarget,chdir,delete,nni)
function adjustWeight(net::HybridNetwork,hmax::Integer,w::Vector{Float64})
    if(hmax - net.numHybrids > 0)
        hmax >= 0 || error("hmax must be non negative: $(hmax)")
        length(w) == 6 || error("length of w should be 6 as there are only 6 moves: $(w)")
        approxEq(sum(w),1.0) || error("vector of move weights should add up to 1: $(w),$(sum(w))")
        all((i->(0<=i<=1)), w) || error("weights must be nonnegative and less than one $(w)")
        suma = w[5]+w[2]+w[3]+w[4]+w[6]
        v = zeros(6)
        k = hmax - net.numHybrids
        for(i in 1:6)
            if(i == 1)
                v[i] = w[1]*k/(suma + w[1]*k)
            else
                v[i] = w[i]/(suma + w[1]*k)
            end
        end
        return v
    end
    return w
end

# function to adjust weights (v) based on movesfail and Nmov, to
# avoid proposing over and over moves that cannot work: (add,mvorigin, mvtarget, chdir, delete,nni)
#returns false if sum(v)=0, no more moves available
# needs the net and hmax to decide if there are available moves
function adjustWeightMovesfail!(v::Vector{Float64}, movesfail::Vector{Int}, Nmov::Vector{Int}, net::HybridNetwork, hmax::Integer)
    length(v) ==length(movesfail) || error("v and movesfail must have same length")
    length(Nmov) ==length(movesfail) || error("Nmov and movesfail must have same length")
    for(i in 1:length(v))
        v[i] = v[i]*(movesfail[i]<Nmov[i] ? 1 : 0)
    end
    if(hmax == 0)
        isTree(net) || error("hmax is $(hmax) but net is not a tree")
        v[6] == 0 && return false #nni
    else
        if(0 < net.numHybrids < hmax)
            sum(v) != 0 || return false #all moves
        elseif(net.numHybrids == 0)
            v[1] == 0 && v[6] == 0 && return false #nni or add
        elseif(net.numHybrids == hmax)
            sum(v[2:4]) + v[6] != 0 || return false #all moves except add/delete
        end
    end
    suma = sum(v)
    for(i in 1:length(v))
        v[i] = v[i]/suma
    end
    return true
end

# function to decide what next move to do when searching
# for topology that maximizes the P-loglik within the space of
# topologies with the same number of hybridizations
# possible moves: move origin/target, change direction hybrid edge, tree nni
# needs the network to know the current numHybrids
# takes as input the vector of weights for each move (add,mvorigin, mvtarget, chdir, delete,nni)
# and dynamic=true, adjusts the weight for addHybrid if net is in a lower layer (net.numHybrids<<hmax)
# movesfail and Nmov are to count number of fails in each move
function whichMove(net::HybridNetwork,hmax::Integer,w::Vector{Float64}, dynamic::Bool, movesfail::Vector{Int}, Nmov::Vector{Int})
    global DEBUG
    hmax >= 0 || error("hmax must be non negative: $(hmax)")
    length(w) == 6 || error("length of w should be 6 as there are only 6 moves: $(w)")
    approxEq(sum(w),1.0) || error("vector of move weights should add up to 1: $(w),$(sum(w))")
    all((i->(0<=i<=1)),w) || error("weights must be nonnegative and less than one $(w)")
    if(hmax == 0)
        isTree(net) || error("hmax is $(hmax) but net is not a tree")
        flag = adjustWeightMovesfail!(w,movesfail,Nmov,net,hmax)
        flag || return :none
        return :nni
    else
        r = rand()
        if(dynamic)
            v = adjustWeight(net,hmax,w)
        else
            v = w
        end
        DEBUG && println("weights before adjusting by movesfail $(v)")
        flag = adjustWeightMovesfail!(v,movesfail,Nmov,net,hmax)
        DEBUG && println("weights after adjusting by movesfail $(v)")
        flag || return :none
        if(0 < net.numHybrids < hmax)
            if(r < v[1])
                return :add
            elseif(r < v[1]+v[2])
                return :MVorigin
            elseif(r < v[1]+v[2]+v[3])
                return :MVtarget
            elseif(r < v[1]+v[2]+v[3]+v[4])
                return :CHdir
            elseif(r < v[1]+v[2]+v[3]+v[4]+v[5])
                return :delete
            else
                return :nni
            end
        elseif(net.numHybrids == 0)
            suma = v[1]+v[6]
            if(r < (v[1])/suma)
                return :add
            else
                return :nni
            end
        else # net.numHybrids == hmax
            suma = v[5]+v[2]+v[3]+v[4]+v[6]
            if(r < v[2]/suma)
                return :MVorigin
            elseif(r < (v[3]+v[2])/suma)
                return :MVtarget
            elseif(r < (v[4]+v[2]+v[3])/suma)
                return :CHdir
            elseif(r < (v[5]+v[2]+v[3]+v[4])/suma)
                return :delete
            else
                return :nni
            end
        end
    end
end

whichMove(net::HybridNetwork,hmax::Integer,movesfail::Vector{Int}, Nmov::Vector{Int}) = whichMove(net,hmax,[1/5,1/5,1/5,1/5,0.0,1/5], true,movesfail, Nmov)
whichMove(net::HybridNetwork,hmax::Integer,w::Vector{Float64},movesfail::Vector{Int}, Nmov::Vector{Int}) = whichMove(net,hmax,w, true,movesfail, Nmov)

#function to choose a hybrid node for the given moves
function chooseHybrid(net::HybridNetwork)
    global DEBUG
    !isTree(net) || error("net is a tree, cannot choose hybrid node")
    net.numHybrids > 1 || return net.hybrid[1]
    index1 = 0
    while(index1 == 0 || index1 > size(net.hybrid,1))
        index1 = round(Integer,rand()*size(net.hybrid,1));
    end
    DEBUG && println("chosen hybrid node for network move: $(net.hybrid[index1].number)")
    return net.hybrid[index1]
end

# function to propose a new topology given a move
# random = false uses the minor hybrid edge always
# count to know in which step we are, N for NNI trials
# order in movescount as in IF here (add,mvorigin,mvtarget,chdir,delete,nni)
function proposedTop!(move::Integer, newT::HybridNetwork,random::Bool, count::Integer, N::Integer, movescount::Vector{Int}, movesfail::Vector{Int})
    global CHECKNET, DEBUG
    1 <= move <= 6 || error("invalid move $(move)") #fixit: if previous move rejected, do not redo it!
    DEBUG && println("current move: $(int2move[move])")
    if(move == 1)
        success = addHybridizationUpdateSmart!(newT,N)
    elseif(move == 2)
        node = chooseHybrid(newT)
        success = moveOriginUpdateRepeat!(newT,node,random)
        CHECKNET && isBadTriangle(node) && success && error("success is $(success) in proposedTop, but node $(node.number) is very bad triangle")
    elseif(move == 3)
        node = chooseHybrid(newT)
        success = moveTargetUpdateRepeat!(newT,node,random)
        CHECKNET && isBadTriangle(node) && success && error("success is $(success) in proposedTop, but node $(node.number) is very bad triangle")
    elseif(move == 4)
        node = chooseHybrid(newT)
        success = changeDirectionUpdate!(newT,node, random)
        CHECKNET && isBadTriangle(node) && success && error("success is $(success) in proposedTop, but node $(node.number) is very bad triangle")
    elseif(move == 5)
        node = chooseHybrid(newT)
        deleteHybridizationUpdate!(newT,node)
        success = true
    elseif(move == 6)
        success = NNIRepeat!(newT,N)
    end
    movescount[move] += 1
    movescount[move+6] += success ? 1 : 0
    movesfail[move] += success ? 0 : 1
    DEBUG && println("success $(success), movescount (add,mvorigin,mvtarget,chdir,delete,nni) proposed: $(movescount[1:6]); successful: $(movescount[7:12]); movesfail: $(movesfail)")
    !success || return true
    DEBUG && println("new proposed topology failed in step $(count) for move $(int2move[move])")
    DEBUG && printEverything(newT)
    CHECKNET && checkNet(newT)
    return false
end


proposedTop!(move::Symbol, newT::HybridNetwork, random::Bool, count::Integer,N::Integer, movescount::Vector{Int},movesfail::Vector{Int}) = proposedTop!(try move2int[move] catch error("invalid move $(string(move))") end,newT, random,count,N, movescount,movesfail)

# function to calculate Nmov, number max of tries per move
# order: (add,mvorigin,mvtarget,chdir,delete,nni)
function calculateNmov!(net::HybridNetwork, N::Vector{Int})
    if(isempty(N))
        N = zeros(Int, 6)
    else
        length(N) == 6 || error("vector Nmov should have length 6: $(N)")
    end
    if(isTree(net))
        N[1] = ceil(coupon(binom(numTreeEdges(net),2))) #add
        N[2] = 1
        N[3] = 1
        N[4] = 1
        N[5] = 1 #delete
        N[6] = ceil(coupon(4*numIntTreeEdges(net))) #nni
    else
        N[1] = ceil(coupon(binom(numTreeEdges(net),2))) #add
        N[2] = ceil(coupon(2*4*net.numHybrids)) #mvorigin
        N[3] = ceil(coupon(2*4*net.numHybrids)) #mtarget
        N[4] = ceil(coupon(2*net.numHybrids)) #chdir
        N[5] = 10000 #delete
        N[6] = ceil(coupon(4*numIntTreeEdges(net))) #nni
    end
end


# function to optimize on the space of networks with the same (or fewer) numHyb
# currT, the starting network will be modified inside
# Nmov: vector with max number of tries per move (add,mvorigin,mvtarget,chdir,delete,nni)
# Nfail: number of failure networks with lower loglik before aborting
# M: multiplier to stop the search if loglik close to M*ftolAbs, or if absDiff less than M*ftolAbs
# hmax: max number of hybrids allowed
# closeN =true if gamma=0.0 fixed only around neighbors with move origin/target
# sout=IOStream to capture the output from the functions called, only used when DEBUG=true and REDIRECT=true, default STDOUT
# logfile=IOStream to capture the information on the heurisitc optimization, default STDOUT
function optTopLevel!(currT::HybridNetwork, M::Number, Nfail::Integer, d::DataCF, hmax::Integer,
                      ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
                      verbose::Bool, closeN ::Bool, Nmov0::Vector{Int}, sout::IO, logfile::IO, writelog::Bool)
    global CHECKNET, DEBUG, REDIRECT
    DEBUG && println("OPT: begins optTopLevel with hmax $(hmax)")
    M > 0 || error("M must be greater than zero: $(M)")
    Nfail > 0 || error("Nfail must be greater than zero: $(Nfail)")
    isempty(Nmov0) || all((n-> (n > 0)), Nmov0) || error("Nmov must be greater than zero: $(Nmov0)")
    if(DEBUG && REDIRECT) #for debugging
        originalSTDOUT = STDOUT
        originalSTDERR = STDERR
        info("redirecting STDOUT and STDERR") # fixit: move this to snaqDebug?
        redirect_stdout(sout)
        redirect_stderr(sout)
    end
    DEBUG && printEverything(currT)
    CHECKNET && checkNet(currT)
    count = 0
    movescount = zeros(Int,18) #1:6 number of times moved proposed, 7:12 number of times success move (no intersecting cycles, etc.), 13:18 accepted by loglik
    movesgamma = zeros(Int,13) #number of moves to fix gamma zero: proposed, successful, movesgamma[13]: total accepted by loglik
    movesfail = zeros(Int,6) #count of failed moves for current topology
    failures = 0
    stillmoves = true
    if(isempty(Nmov0))
        Nmov = zeros(Int,6)
    else
        Nmov = deepcopy(Nmov0)
    end
    all((e->!(e.hybrid && e.inCycle == -1)), currT.edge) || error("found hybrid edge with inCycle == -1")
    optBL!(currT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
    currT = afterOptBLAll!(currT, d, Nfail,closeN , M, ftolAbs, verbose,movesgamma,ftolRel,xtolRel,xtolAbs)
    absDiff = M*ftolAbs + 1
    newT = deepcopy(currT)
    DEBUG && printEdges(newT)
    DEBUG && printPartitions(newT)
    #DEBUG && printNodes(newT)
    DEBUG && println(writeTopologyLevel1(newT,true))
    writelog && write(logfile, "\nBegins heuristic optimization of network------\n")
    while(absDiff > M*ftolAbs && failures < Nfail && currT.loglik > M*ftolAbs && stillmoves) #stops if close to zero because of new deviance form of the pseudolik
        count += 1
        DEBUG && println("--------- loglik_$(count) = $(round(currT.loglik,6)) -----------")
        if(isempty(Nmov0)) #if empty, not set by user
            calculateNmov!(newT,Nmov)
        end
        DEBUG && println("will propose move with movesfail $(movesfail), Nmov $(Nmov)")
        move = whichMove(newT,hmax,movesfail,Nmov)
        if(move != :none)
            flag = proposedTop!(move,newT,true, count,10, movescount,movesfail) #N=10 because with 1 it never finds an edge for nni
            if(flag) #no need else because newT always undone if failed
                accepted = false
                all((e->!(e.hybrid && e.inCycle == -1)), newT.edge) || error("found hybrid edge with inCycle == -1")
                DEBUG && println("proposed new topology in step $(count) is ok to start optBL")
                DEBUG && printEverything(newT)
                CHECKNET && checkNet(newT)
                optBL!(newT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
                DEBUG && println("OPT: comparing newT.loglik $(newT.loglik), currT.loglik $(currT.loglik)")
                if(newT.loglik < currT.loglik && abs(newT.loglik-currT.loglik) > M*ftolAbs) #newT better loglik: need to check for error or keeps jumping back and forth
                    newloglik = newT.loglik
                    newT = afterOptBLAll!(newT, d, Nfail,closeN , M, ftolAbs,verbose,movesgamma,ftolRel, xtolRel,xtolAbs) #needed return because of deepcopy inside
                    DEBUG && println("loglik before afterOptBL $(newloglik), newT.loglik now $(newT.loglik), loss in loglik by fixing gamma(z)=0.0(1.0): $(newloglik>newT.loglik ? 0 : abs(newloglik-newT.loglik))")
                    accepted = true
                else
                    accepted = false
                end
                if(accepted)
                    absDiff = abs(newT.loglik - currT.loglik)
                    DEBUG && println("proposed new topology with better loglik in step $(count): oldloglik=$(round(currT.loglik,3)), newloglik=$(round(newT.loglik,3)), after $(failures) failures")
                    currT = deepcopy(newT)
                    failures = 0
                    movescount[move2int[move]+12] += 1
                    movesfail = zeros(Int,6) #count of failed moves for current topology
                else
                    DEBUG && println("rejected new topology with worse loglik in step $(count): currloglik=$(round(currT.loglik,3)), newloglik=$(round(newT.loglik,3)), with $(failures) failures")
                    failures += 1
                    movesfail[move2int[move]] += 1
                    newT = deepcopy(currT)
                end
                DEBUG && printEdges(newT)
                DEBUG && printPartitions(newT)
                #printNodes(newT)
                DEBUG && println(writeTopologyLevel1(newT,true))
                DEBUG && println("ends step $(count) with absDiff $(accepted? absDiff : 0.0) and failures $(failures)")
            end
        else
            stillmoves = false
        end
        # here would set the debug file up to beginning
        if(DEBUG && REDIRECT)
            seekstart(sout)
        end
    end
    if(ftolAbs > 1e-7 || ftolRel > 1e-7 || xtolAbs > 1e-7 || xtolRel > 1e-7)
        writelog && write(logfile,"\nfound best network, now we re-optimize branch lengths and gamma more precisely")
        optBL!(newT,d,verbose,1e-12,1e-10,1e-10,1e-10)
    end
    assignhybridnames!(newT)
    if(absDiff <= M*ftolAbs)
        writelog && write(logfile,"\nSTOPPED by absolute difference criteria")
    elseif(currT.loglik <= M*ftolAbs)
        writelog && write(logfile,"\nSTOPPED by loglik close to zero criteria")
    elseif(!stillmoves)
        writelog && write(logfile,"\nSTOPPED for not having more moves to propose: movesfail $(movesfail), Nmov $(Nmov)")
    else
        writelog && write(logfile,"\nSTOPPED by number of failures criteria")
    end
    ## if(newT.loglik > M*ftolAbs) #not really close to 0.0, based on absTol also
    ##     write(logfile,"\nnewT.loglik $(newT.loglik) not really close to 0.0 based on loglik abs. tol. $(M*ftolAbs), you might need to redo with another starting point")
    ## end
    if(newT.numBad > 0) # if bad diamond I, need to keep gammaz info in parenthetical description
        for(n in newT.hybrid)
            setGammaBLfromGammaz!(n,newT) # get t and γ that are compatible with estimated gammaz values
        end
    end
    writelog && write(logfile,"\nEND optTopLevel: found minimizer topology at step $(count) (failures: $(failures)) with -loglik=$(round(newT.loglik,5)) and ht_min=$(round(newT.ht,5))")
    writelog && printCounts(movescount,movesgamma,logfile)
    DEBUG && printEdges(newT)
    DEBUG && printPartitions(newT)
    DEBUGC && printNodes(newT)
    DEBUG && println(writeTopologyLevel1(newT,true))
    if(DEBUG && REDIRECT)
        redirect_stdout(originalSTDOUT)
        redirect_stderr(originalSTDERR)
        info("end of redirection")
    end
    return newT
end

optTopLevel!(currT::HybridNetwork, d::DataCF, hmax::Integer) = optTopLevel!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false,true,numMoves,STDOUT, STDOUT,true)
optTopLevel!(currT::HybridNetwork, d::DataCF, hmax::Integer, verbose::Bool) = optTopLevel!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, verbose,true,numMoves,STDOUT,STDOUT,true)
optTopLevel!(currT::HybridNetwork, M::Number, Nfail::Integer, d::DataCF, hmax::Integer,ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64, verbose::Bool, closeN ::Bool, Nmov0::Vector{Int}) = optTopLevel!(currT, M, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN , Nmov0,STDOUT,STDOUT,true)



# function to print count of number of moves after optTopLevel
# s: IOStream to decide if you want to print to a file or to screen, by default to screen
function printCounts(movescount::Vector{Int}, movesgamma::Vector{Int},s::Union{IOStream,Base.TTY})
    length(movescount) == 18 || error("movescount should have length 18, not $(length(movescount))")
    length(movesgamma) == 13 || error("movesgamma should have length 13, not $(length(movescount))")
    print(s,"\nPERFORMANCE: total number of moves (proposed, successful, accepted) in general, and to fix gamma=0.0,t=0.0 cases\n")
    print(s,"\t--------moves general--------\t\t\t\t --------moves gamma,t------\t\n")
    print(s,"move\t Num.Proposed\t Num.Successful\t Num.Accepted\t | Num.Proposed\t Num.Successful\t Num.Accepted\n")
    names = ["add","mvorigin","mvtarget","chdir","delete","nni"]
    for(i in 1:6)
        if(i == 2 || i == 3)
            print(s,"$(names[i])\t $(movescount[i])\t\t $(movescount[i+6])\t\t $(movescount[i+12])\t |\t $(movesgamma[i])\t\t $(movesgamma[i+6])\t\t --\n")
        elseif(i == 1)
            print(s,"$(names[i])\t\t $(movescount[i])\t\t $(movescount[i+6])\t\t $(movescount[i+12])\t |\t NA \t\t NA \t\t NA \n")
        else
            print(s,"$(names[i])\t\t $(movescount[i])\t\t $(movescount[i+6])\t\t $(movescount[i+12])\t |\t $(movesgamma[i])\t\t $(movesgamma[i+6])\t\t --\n")
        end
    end
    suma = sum(movescount[7:12]);
    suma2 = sum(movesgamma[7:12]) == 0 ? 1 : sum(movesgamma[7:12])
    print(s,"Total\t\t $(sum(movescount[1:6]))\t\t $(sum(movescount[7:12]))\t\t $(sum(movescount[13:18]))\t |\t $(sum(movesgamma[1:6]))\t\t $(sum(movesgamma[7:12]))\t\t $(movesgamma[13])\n")
    print(s,"Proportion\t -- \t\t $(round(sum(movescount[7:12])/suma,1))\t\t $(round(sum(movescount[13:18])/suma,1))\t |\t -- \t\t $(round(sum(movesgamma[7:12])/suma2,1))\t\t $(round(movesgamma[13]/suma2,1))\n")
end

printCounts(movescount::Vector{Int}, movesgamma::Vector{Int}) = printCounts(movescount, movesgamma,STDOUT)

# function to print count of number of moves to a file
# file=name of file where to save
function printCounts(movescount::Vector{Int}, movesgamma::Vector{Int},file::AbstractString)
    s = open(file,"w")
    printCounts(movescount,movesgamma,s)
    close(s)
end


# function to move down onw level to h-1
# caused by gamma=0,1 or gammaz=0,1
function moveDownLevel!(net::HybridNetwork)
    global CHECKNET, DEBUG
    !isTree(net) ||error("cannot delete hybridization in a tree")
    DEBUG && println("MOVE: need to go down one level to h-1=$(net.numHybrids-1) hybrids because of conflicts with gamma=0,1")
    DEBUG && printEverything(net)
    CHECKNET && checkNet(net)
    nh = net.ht[1 : net.numHybrids - net.numBad]
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    nt = net.ht[net.numHybrids - net.numBad + 1 : net.numHybrids - net.numBad + k]
    nhz = net.ht[net.numHybrids - net.numBad + k + 1 : length(net.ht)]
    indh = net.index[1 : net.numHybrids - net.numBad]
    indhz = net.index[net.numHybrids - net.numBad + k + 1 : length(net.ht)]
    flagh,flagt,flaghz = isValid(nh,nt,nhz)
    if(!flagh)
        for(i in 1:length(nh))
            if(approxEq(nh[i],0.0) || approxEq(nh[i],1.0))
                edge = net.edge[indh[i]]
                node = edge.node[edge.isChild1 ? 1 : 2];
                node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
                deleteHybridizationUpdate!(net,node)
            end
            break
        end
    elseif(!flaghz)
        net.numBad > 0 || error("not a bad diamond I and flaghz is $(flaghz), should be true")
        i = 1
        while(i <= length(nhz))
            if(approxEq(nhz[i],0.0))
                nodehz = net.node[indhz[i]]
                approxEq(nodehz.gammaz,nhz[i]) || error("nodehz $(nodehz.number) gammaz $(nodehz.gammaz) should match the gammaz in net.ht $(nhz[i]) and it does not")
                edges = hybridEdges(nodehz)
                edges[1].hybrid || error("bad diamond I situation, node $(nodehz.number) has gammaz $(nodehz.gammaz) so should be linked to hybrid edge, but it is not")
                node = edges[1].node[edges[1].isChild1 ? 1 : 2];
                node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
                deleteHybridizationUpdate!(net,node)
                break
            elseif(approxEq(nhz[i],1.0))
                approxEq(nhz[i+1],0.0) || error("gammaz for node $(net.node[indhz[i]].number) is $(nhz[i]) but the other gammaz is $(nhz[i+1]), the sum should be less than 1.0 (movedownlevel)")
                nodehz = net.node[indhz[i+1]]
                approxEq(nodehz.gammaz,nhz[i+1]) || error("nodehz $(nodehz.number) gammaz $(nodehz.gammaz) should match the gammaz in net.ht $(nhz[i+1]) and it does not")
                edges = hybridEdges(nodehz)
                edges[1].hybrid || error("bad diamond I situation, node $(nodehz.number) has gammaz $(nodehz.gammaz) so should be linked to hybrid edge, but it is not")
                node = edges[1].node[edges[1].isChild1 ? 1 : 2];
                node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
                deleteHybridizationUpdate!(net,node)
                break
            else
                if(approxEq(nhz[i+1],0.0))
                    nodehz = net.node[indhz[i+1]];
                    approxEq(nodehz.gammaz,nhz[i+1]) || error("nodehz $(nodehz.number) gammaz $(nodehz.gammaz) should match the gammaz in net.ht $(nhz[i+1]) and it does not")
                    edges = hybridEdges(nodehz);
                    edges[1].hybrid || error("bad diamond I situation, node $(nodehz.number) has gammaz $(nodehz.gammaz) so should be linked to hybrid edge, but it is not")
                    node = edges[1].node[edges[1].isChild1 ? 1 : 2];
                    node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
                    deleteHybridizationUpdate!(net,node)
                    break
                elseif(approxEq(nhz[i+1],1.0))
                    error("gammaz for node $(net.node[indhz[i]].number) is $(nhz[i]) but the other gammaz is $(nhz[i+1]), the sum should be less than 1.0 (movedownlevel)")
                end
            end
            i += 2
        end
    end
    DEBUG && printEverything(net)
    CHECKNET && checkNet(net)
end

# checks if there are problems in estimated net.ht:
# returns flag for h, flag for t, flag for hz
function isValid(net::HybridNetwork)
    nh = net.ht[1 : net.numHybrids - net.numBad]
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    nt = net.ht[net.numHybrids - net.numBad + 1 : net.numHybrids - net.numBad + k]
    nhz = net.ht[net.numHybrids - net.numBad + k + 1 : length(net.ht)]
    #println("isValid on nh $(nh), nt $(nt), nhz $(nhz)")
    return all((n->(0<n<1 && !approxEq(n,0.0) && !approxEq(n,1.0))), nh), all((n->(n>0 && !approxEq(n,0.0))), nt), all((n->(0<n<1 && !approxEq(n,0.0) && !approxEq(n,1.0))), nhz)
end

# checks if there are problems in estimated net.ht:
# returns flag for h, flag for t, flag for hz
function isValid(nh::Vector{Float64},nt::Vector{Float64},nhz::Vector{Float64})
    #println("isValid on nh $(nh), nt $(nt), nhz $(nhz)")
    return all((n->(0<n<1 && !approxEq(n,0.0) && !approxEq(n,1.0))), nh), all((n->(n>0 && !approxEq(n,0.0))), nt), all((n->(0<n<1 && !approxEq(n,0.0) && !approxEq(n,1.0))), nhz)
end

# function to do optTopLevel for each level below hmax
# returns the best topology overall
# uses the best topology in h-1 as starting point for the search in h
# no test in between! fixit: add a test to decide to move from h-1 to h?
function optTop!(currT::HybridNetwork, M::Number, Nfail::Integer, d::DataCF, hmax::Integer,ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64, verbose::Bool, closeN ::Bool, Nmov0::Vector{Int},s::IO)
    global DEBUG
    hmax >= 0 || error("hmax cannot be negative $(hmax)")
    currT.numHybrids == 0 || write(s,"\ncurrT has already $(currT.numHybrids), so search will not start in tree, but in the space of networks with $(currT.numHybrids) hybrids")
    currT.numHybrids <= hmax || error("currT has more hybrids: $(currT.numHybrids) than hmax $(hmax)")
    bestT = HybridNetwork[]

    maxNet = HybridNetwork();
    maxNet.loglik = 1.e15;

    push!(bestT,currT)
    i = 1
    for(h in currT.numHybrids:hmax)
        startT = deepcopy(bestT[i]);
        write(s,"\n")
        write(s,"\nStart search on space of $(h) hybridizations with previous opt topology $(writeTopologyLevel1(startT,true))")
        j = 0
        if(h > 0)
            success = false
            while(!success && j < Nfail) #try to add new hybrid
                success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(startT)
                if(!success)
                    write(s,"\n$(j)th failed new added hybridization, will delete and try again until $(Nfail) failures")
                    j += 1
                    deleteHybridizationUpdate!(startT,hybrid)
                end
            end
            if(!success)
                write(s,"\ncould not find a place for new hybridization after $(Nfail) attempts, will stop search here with $(h) hybridizations, instead of hmax= $(hmax)")
                maxNet.loglik < 1.e15 || error("never updated maxNet")
                return maxNet
            end
        end
        newT = optTopLevel!(startT, M, Nfail, d, h,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN , Nmov0)
        write(s,"\nBEST NETWORK at $(h) hybridizations: $(writeTopologyLevel1(newT,true)) with pseudologlik $(round(newT.loglik,5))")
        write(s,"\nfailed attempts to place a new hybridization to begin with: $(j)")
        push!(bestT,deepcopy(newT))
        DEBUG && println("maxNet.loglik $(maxNet.loglik), newT.loglik $(newT.loglik)")
        if(newT.loglik < maxNet.loglik)
            DEBUG && println("entra a q newT tiene mejor loglik: maxNet.loglik $(maxNet.loglik), newT.loglik $(newT.loglik)")
            maxNet = deepcopy(newT)
        end
        DEBUG && println("maxNet.loglik $(maxNet.loglik)")
        i += 1
    end
    maxNet.loglik < 1.e15 || error("never updated maxNet")
    DEBUG && println("typeof maxNet: $(typeof(maxNet))")
    return maxNet
end



optTop!(currT::HybridNetwork, d::DataCF, hmax::Integer) = optTop!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false,true,numMoves,STDOUT)
optTop!(currT::HybridNetwork, d::DataCF, hmax::Integer, verbose::Bool) = optTop!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, verbose,true,numMoves,STDOUT)
optTop!(currT::HybridNetwork, M::Number, Nfail::Integer, d::DataCF, hmax::Integer,ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64, verbose::Bool, closeN ::Bool, Nmov0::Vector{Int}) = optTop!(currT, M, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN , Nmov0, STDOUT)

"""
Road map for various functions behind snaq!

    snaq!
    optTopRuns!
    optTopRun1!
    optTopLevel!
    optBL!
    snaqDebug

All return their optimized network.

- snaq! calls optTopRuns! once, after a deep copy of the starting network.
- optTopRuns! calls optTopRun1! several (nrun) times.
  assumes level-1 network with >0 branch lengths.
  each call gets the same starting network.
- optTopRun1! calls optTopLevel! once, after deep copying + changing the starting network slightly.
- optTopLevel! calls optBL! various times and proposes new network with various moves.
- snaqDebug calls optTopRun1! once

Functions that are no longer used, and not tested as much:

    optTop! : optimizes at h=0, then h=1 etc.
"""
# function to repeat optTopLevel for a certain number of runs
# used to run optTop but now we want to compare to phylonet
# runs will add one extra to desired because the first run is not counted towards time (because of compilation time), default runs=10
# outgroup is needed to root appropriately to use java distance function
# rootname is for output files: log, err, out, default "optTopRuns"
# seed= integer to set the seed, default 0, so clocktime used
# probST is the probability of starting in the starting topology currT per run (1-probST is the probability of doing a NNI move) default 0.3
# updateBL=true if we want to update the missing BL to average CF
# does *not* modify currT0
function optTopRuns!(currT0::HybridNetwork, M::Number, Nfail::Integer, d::DataCF, hmax::Integer,
                     ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
                     verbose::Bool, closeN ::Bool, Nmov0::Vector{Int}, runs::Integer,
                     outgroup::AbstractString, rootname::AbstractString, seed::Integer, probST::Float64)
    writelog = true
    if (rootname != "")
        juliaerr = string(rootname,".err")
        errfile = open(juliaerr,"w")
        julialog = string(rootname,".log")
        logfile = open(julialog,"w")
        juliaout = string(rootname,".out")
    else
      writelog = false
      logfile = STDOUT # used in call to optTopRun1!
    end

    str = """optimization of topology, BL and inheritance probabilities using:
              hmax = $(hmax),
              tolerance parameters: ftolRel=$(ftolRel), ftolAbs=$(ftolAbs),
                                    xtolAbs=$(xtolAbs), xtolRel=$(xtolRel).
              max number of failed proposals = $(Nfail), multiplier M = $(M).
             """
    if outgroup != "none"
        str *= "Outgroup: $(outgroup) (for rooting at the final step)\n"
    end
    str *= (writelog ? "rootname for files: $(rootname)\n" : "no output files\n")
    str *= "BEGIN: $(runs) runs on starting tree $(writeTopologyLevel1(currT0,true))\n"
    if (writelog)
      write(logfile,str)
      flush(logfile)
    end
    print(STDOUT,str)
    print(STDOUT, "$(Libc.strftime(time()))\n")
    # time printed to logfile at start of every run, not here.

    maxNet = HybridNetwork();
    maxNet.loglik = 1.e15;

    failed = Int[] #seeds with bugs
    bestnet = HybridNetwork[];
    #runs += 1 #extra run for compilation in v0.3

    if(seed == 0)
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) #better seed based on clock
    end
    if (writelog)
      write(logfile,"\nmain seed $(seed)\n")
      flush(logfile)
    else print(STDOUT,"\nmain seed $(seed)\n"); end
    srand(seed)
    seeds = [seed;round(Integer,floor(rand(runs-1)*100000))]

    tic();
    for(i in 1:runs)
        writelog && write(logfile,"seed: $(seeds[i]) for run $(i)\n$(Libc.strftime(time()))\n")
        writelog && flush(logfile)
        print(STDOUT,"seed: $(seeds[i]) for run $(i)\n")
        gc();
        try
            writelog && write(logfile,"\n BEGIN SNaQ for run $(i), seed $(seeds[i]) and hmax $(hmax)")
            verbose && print(STDOUT,"\n BEGIN SNaQ for run $(i), seed $(seeds[i]) and hmax $(hmax)")
            best = optTopRun1!(currT0, M, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN , Nmov0,seeds[i],logfile,writelog,probST);
            writelog && write(logfile,"\n FINISHED SNaQ for run $(i), -loglik of best $(best.loglik)\n")
            verbose && print(STDOUT,"\n FINISHED SNaQ for run $(i), -loglik of best $(best.loglik)\n")
            if (writelog)
              if(outgroup == "none")
                write(logfile,writeTopologyLevel1(best,true)) #no outgroup
              else
                write(logfile,writeTopologyLevel1(best,outgroup)) #outgroup
              end
            flush(logfile)
            end
            push!(bestnet, best)
            if(best.loglik < maxNet.loglik)
                maxNet = best
            end
            DEBUG && println("typeof maxNet $(typeof(maxNet)), -loglik of maxNet $(maxNet.loglik)")
        catch(err)
            if (writelog)
            write(logfile,"\n ERROR found on SNaQ for run $(i) seed $(seeds[i]): $(err)")
            flush(logfile)
            write(errfile,"\n ERROR found on SNaQ for run $(i) seed $(seeds[i]): $(err)")
            flush(errfile)
            push!(failed,seeds[i])
            else
                warn("\n ERROR found on SNaQ for run $(i) seed $(seeds[i]): $(err)")
            end
        end
        writelog && write(logfile,"\n---------------------\n")
        writelog && flush(logfile)
    end
    t=toc();
    if (writelog)
    write(errfile, "\n Total errors: $(length(failed)) in seeds $(failed)")
    write(logfile,"\n$(Libc.strftime(time()))")
    close(errfile)
    else
        verbose && print(STDOUT,"""\n Total errors: $(length(failed)) in seeds $(failed)
                                     $(Libc.strftime(time()))""")
    end

    maxNet.numTaxa > 0 || error("the best network is empty!")

    ## need to do this before setting BL to -1
    if (writelog && !isTree(maxNet)) ## only do networks file if maxNet is not tree
        println("best network and networks with different hybrid/gene flow directions printed to .networks file")
        julianet = string(rootname,".networks")
        s = open(julianet,"w")
        otherNet = []
        try
            otherNet = undirectedOtherNetworks(maxNet, outgroup=outgroup, insideSnaq=true) # do not use rootMaxNet
        catch
            write(s,"""Bug found when trying to obtain networks with modified hybrid/gene flow direction.
                       To help debug these cases and get other similar estimated networks for your analysis,
                       please send the estimated network in parenthetical format to claudia@stat.wisc.edu
                       with the subject BUG IN NETWORKS FILE. You can get this network from the .out file.
                       You can also post this problem to the google group, or github issues. Thank you!\n""")
        end
        write(s,"$(writeTopologyLevel1(maxNet,true)), with -loglik $(maxNet.loglik) (best network found, remaining sorted by log-pseudolik; the smaller, the better)\n")
        # best network is included first: for score comparison with other networks
        foundBad = false
        for(n in otherNet)
            try
                optBL!(n,d) ##optBL MUST have network with all the attributes, and undirectedOtherNetworks will return "good" networks that way
                if(n.numBad > 0) # to keep gammaz info in parenthetical description of bad diamond I
                    for(nod in n.hybrid)
                        setGammaBLfromGammaz!(nod,n) # get t and γ that are compatible with estimated gammaz values
                    end
                end
                setNonIdBL!(n)
            catch
                n.loglik = -1
                foundBad = true
            end
        end
        ## to sort otherNet by loglik value:
        ind = sortperm([n.loglik for n in otherNet])
        otherNet = otherNet[ind]
        for(n in otherNet)
            write(s,"$(writeTopologyLevel1(n,true)), with -loglik $(n.loglik)\n")
        end
        foundBad && write(s,"Problem found when optimizing branch lengths for some networks, left loglik as -1. Please report this issue to claudia@stat.wisc.edu, google group or github issues. Thank you!")
        close(s)
    end

    setNonIdBL!(maxNet)
    writelog &&
    write(logfile,"\nMaxNet is $(writeTopologyLevel1(maxNet,true)) \nwith -loglik $(maxNet.loglik)\n")
    print(STDOUT,"\nMaxNet is $(writeTopologyLevel1(maxNet,true)) \nwith -loglik $(maxNet.loglik)\n")

    # new output file with (nearly) non-identifiable networks: (run before maxNet gets rooted)

    if outgroup != "none"
        try
            checkRootPlace!(maxNet,outgroup=outgroup) ## keeps all attributes
        catch err
            if isa(err, RootMismatch)
                 println("RootMismatch: ", err.msg,
                 """\nThe estimated network has hybrid edges that are incompatible with the desired outgroup.
                    Reverting to an admissible root position.
                    """)
            else
                println("error trying to reroot: ", err.msg);
            end
            checkRootPlace!(maxNet,verbose=false) # message about problem already printed above
        end
    else
        checkRootPlace!(maxNet,verbose=false) #leave root in good place after snaq
    end
    s = writelog ? open(juliaout,"w") : STDOUT
    str = writeTopologyLevel1(maxNet) * """
     -Ploglik = $(maxNet.loglik)
     Dendroscope: $(writeTopologyLevel1(maxNet,di=true))
     Elapsed time: $(t) seconds in $(runs-length(failed)) successful runs
    -------
    List of estimated networks for all runs (sorted by log-pseudolik; the smaller, the better):
    """
    ## to sort bestnet by loglik value:
    ind = sortperm([n.loglik for n in bestnet])
    bestnet = bestnet[ind]
    for(n in bestnet)
      str *= " "
      str *= (outgroup == "none" ? writeTopologyLevel1(n,true) :
                                   writeTopologyLevel1(n,outgroup))
      str *= ", with -loglik $(n.loglik)\n"
    end
    str *= "-------\n"
    write(s,str);
    writelog && close(s) # to close juliaout file (but not STDOUT!)
    writelog && close(logfile)

    return maxNet
end

optTopRuns!(currT::HybridNetwork, d::DataCF, hmax::Integer, runs::Integer, outgroup::AbstractString, rootname::AbstractString) = optTopRuns!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false, true, numMoves, runs, outgroup,rootname,0,0.3)
optTopRuns!(currT::HybridNetwork, d::DataCF, hmax::Integer, runs::Integer, outgroup::AbstractString) = optTopRuns!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false, true, numMoves, runs, outgroup,"optTopRuns",0,0.3)
optTopRuns!(currT::HybridNetwork, d::DataCF, hmax::Integer, runs::Integer) = optTopRuns!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false, true, numMoves, runs, "none", "optTopRuns",0,0.3)
optTopRuns!(currT::HybridNetwork, d::DataCF, hmax::Integer) = optTopRuns!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false, true, numMoves, 10, "none", "optTopRuns",0,0.3)

# function that runs one run inside optTopRuns: it changes the starting topology currT0 and calls optTopLevel
# it has an extra parameter to optTopRuns: seed to set the seed inside and change the starting topology
# currT0
# does *not* modify currT0. Modifies data d only.
function optTopRun1!(currT0::HybridNetwork, M::Number, Nfail::Integer, d::DataCF, hmax::Integer,
                     ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
                     verbose::Bool, closeN ::Bool, Nmov0::Vector{Int},seed::Integer,
                     logfile::IO, writelog::Bool, probST::Float64, sout::IO)
    srand(seed)
    currT = deepcopy(currT0);
    if(probST<1.0 && rand() < 1-probST) # modify starting tree by a nni move
        suc = NNIRepeat!(currT,10); #will try 10 attempts to do an nni move, if set to 1, hard to find it depending on currT
        writelog && suc && write(logfile," changed starting topology by NNI move\n")
        if(!isTree(currT))
            if(rand() < 1-probST) # modify starting network by mvorigin, mvtarget with equal prob
                if(currT.numHybrids == 1)
                    ind = 1
                else
                    ind = 0
                    while(ind == 0 || ind > length(currT.hybrid))
                        ind = round(Integer,rand()*length(currT.hybrid));
                    end
                end
                if(rand()<0.5)
                    suc = moveOriginUpdateRepeat!(currT,currT.hybrid[ind],true)
                    writelog && suc && write(logfile,"\n changed starting network by move origin")
                else
                    suc = moveTargetUpdateRepeat!(currT,currT.hybrid[ind],true)
                    writelog && suc && write(logfile,"\n changed starting network by move target")
                end

            end
        end
    end
    gc();
    optTopLevel!(currT, M, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN , Nmov0,sout,logfile,writelog)
end

optTopRun1!(currT0::HybridNetwork, M::Number, Nfail::Integer, d::DataCF, hmax::Integer,ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64, verbose::Bool, closeN ::Bool, Nmov0::Vector{Int},seed::Integer,logfile::IO,writelog::Bool,probST::Float64) = optTopRun1!(currT0, M, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0,seed,logfile,writelog,probST, STDOUT)
optTopRun1!(currT::HybridNetwork, d::DataCF, hmax::Integer) = optTopRun1!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false, true, numMoves, 0,STDOUT,true,0.3, STDOUT)
optTopRun1!(currT::HybridNetwork, d::DataCF, hmax::Integer, seed::Integer) = optTopRun1!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false, true, numMoves,seed,STDOUT,true,0.3, STDOUT)


# function SNaQ: it calls directly optTopRuns but has a prettier name
# it differs from optTopRuns in that it creates a deepcopy of the starting topology,
#    check that it's of level 1 and updates its BL.
# only currT and d are necessary, all others are optional and have default values
"""
`snaq!(T::HybridNetwork, d::DataCF)`

Estimate the network (or tree) to fit observed quartet concordance factors (CFs)
stored in a DataCF object, using maximum pseudo-likelihood. A level-1 network is assumed.
The search starts from topology `T`,
which can be a tree or a network with no more than `hmax` hybrid nodes.
The function name ends with ! because it modifies the CF data `d` by updating its
attributes `expCF`: CFs expected under the network model.
It does *not* modify `T`.
The quartet pseudo-deviance is the negative log pseudo-likelihood,
up to an additive constant, such that a perfect fit corresponds to a deviance of 0.0.

Output:

- estimated network in file `.out` (also in `.log`): best network overall and list of
  networks from each individual run.
- the best network and modifications of it, in file `.networks`.
  All networks in this file have the same undirected topology as the best network,
  but have different hybrid/gene flow directions. These other networks are reported with
  their pseudo-likelihood scores, because
   non-identifiability issues can cause them to have very similar scores, and because
   SNaQ was shown to estimate the undirected topology accurately but not the direction of
   hybridization in cases of near non-identifiability.
- if any error occurred, file `.err` provides information (seed) to reproduce the error.

There are many optional arguments, including

- hmax: maximum number of hybridizations allowed (default 1)
- verbose: if true, print information about the numerical optimization
- runs: number of independent starting points for the search (default 10)
- outgroup: outgroup taxon to root the estimated topology at the very end
- filename: root name for the output files. Default is "snaq". If empty (""),
  files are *not* created, progress log goes to the screen only (standard out).
- seed: seed to replicate a given search

The following optional arguments control when to stop the optimization of branch lengths
and γ's on each individual candidate network. Defaults are in parentheses:

- ftolRel (1e-5) and ftolAbs (1e-6): relative and absolute differences of the network score
  between the current and proposed parameters,
- xtolRel (1e-3) and xtolAbs (1e-4): relative and absolute differences between the current
  and proposed parameters.

Greater values will result in a less thorough but faster search. These parameters are used
when evaluating candidate networks only. Branch lengths and γ's are optimized on the last
"best" network with different and very thorough tolerance parameters (1e-12 for ftolRel
and 1e-10 for the others).

The following optional arguments control when to stop proposing new network topologies:

- Nfail (100): maximum number of times that new topologies are proposed and rejected (in a row).
- M (10000): the proposed network is accepted if its score is better than the current score by
  at least M*ftolAbs.

Lower values of Nfail and greater values of M would result in a less thorough but faster search.

See also: `topologyMaxQPseudolik!` to optimize parameters on a fixed topology,
and `topologyQPseudolik!` to get the deviance (pseudo log-likelihood up to a constant)
of a fixed topology with fixed parameters.

Reference:
Claudia Solís-Lemus and Cécile Ané (2016).
Inferring phylogenetic networks with maximum pseudolikelihood under incomplete lineage sorting.
[PLoS Genetics](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896)
12(3):e1005896
"""
function snaq!(currT0::HybridNetwork, d::DataCF; hmax=1::Integer, M=multiplier::Number, Nfail=numFails::Integer,ftolRel=fRel::Float64, ftolAbs=fAbs::Float64, xtolRel=xRel::Float64, xtolAbs=xAbs::Float64, verbose=false::Bool, closeN=true::Bool, Nmov0=numMoves::Vector{Int}, runs=10::Integer, outgroup="none"::AbstractString, filename="snaq"::AbstractString, seed=0::Integer, probST=0.3::Float64, updateBL=true::Bool)
    0.0<=probST<=1.0 || error("probability to keep the same starting topology should be between 0 and 1: $(probST)")
    sameTaxa(d,currT0) || error("some taxon names in quartets do not appear on the starting topology")
    currT0.numTaxa >= 5 || error("cannot estimate hybridizations in topologies with fewer than 5 taxa, this topology has $(currT0.numTaxa) taxa")
    # need a clean starting net. fixit: maybe we need to be more thorough here
    # yes, need to check that everything is ok because it could have been cleaned and then modified
    startnet = readTopologyUpdate(writeTopologyLevel1(currT0)) # update all level-1 things
    flag = checkNet(startnet,true) # light checking only
    flag && error("starting topology suspected not level-1")
    try
        checkNet(startnet)
    catch err
        err.msg = "starting topology not a level 1 network:\n" * err.msg
        rethrow(err)
    end
    if(updateBL && isTree(startnet))
        updateBL!(startnet,d)
    end
    # for the case of multiple alleles: expand into two leaves quartets like sp1 sp1 sp2 sp3.
    if(!isempty(d.repSpecies))
        expandLeaves!(d.repSpecies,startnet)
    end
    optTopRuns!(startnet, M, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0, runs, outgroup, filename,seed,probST)
end


# function to run SNaQ for a given seed (a given run) to debug
# the parameter runs is ignored (put here to make it equivalent to snaq)
"""
`snaqDebug(currT::HybridNetwork,d::DataCF)`

function to replicate a given run that produces error according to the .err file generated by snaq.
The same settings used in that run should be used in this function, specially the seed. See the readme file online for more details.
"""
function snaqDebug(currT0::HybridNetwork, d::DataCF; hmax=1::Integer, M=multiplier::Number, Nfail=numFails::Integer,
                   ftolRel=fRel::Float64, ftolAbs=fAbs::Float64, xtolRel=xRel::Float64, xtolAbs=xAbs::Float64,
                   verbose=false::Bool, closeN=true::Bool, Nmov0=numMoves::Vector{Int}, runs=10::Integer,
                   outgroup="none"::AbstractString, rootname="snaqDebug"::AbstractString, seed=0::Integer, probST=0.3::Float64,
                   debug=true::Bool, redirect=true::Bool)
    currT0 = readTopologyUpdate(writeTopologyLevel1(currT0)) # update all level-1 things
    flag = checkNet(currT0,true)
    flag && error("starting topology suspected not level-1")

    global DEBUG, REDIRECT
    debuglog = string(rootname,".debug.log")
    if (debug) DEBUG = true; end
    if (redirect)
        REDIRECT = true
        println("file for redirection: $debuglog")
    end
    sout = open(debuglog,"w")
    write(sout,"DEBUG for seed $(seed)")
    write(sout,"\n$(Libc.strftime(time()))")
    flush(sout)
    julialog = string(rootname,".log")
    logfile = open(julialog,"w")
    write(logfile,"optimization of topology and BL for hmax = $(hmax), and tolerance parameters: ftolRel=$(ftolRel), ftolAbs= $(ftolAbs), xtolAbs= $(xtolAbs), xtolRel= $(xtolRel). Max number of failed proposals is $(Nfail), and multiplier M is $(M).")
    write(logfile,"\n Outgroup defined as $(outgroup) and rootname for files $(rootname)")
    write(logfile,"\nBEGIN: 1 run on starting tree $(writeTopologyLevel1(currT0,true))")
    write(logfile,"\n$(Libc.strftime(time()))")
    flush(logfile)
    res = nothing
    try
        res = optTopRun1!(currT0, M, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0,seed,logfile,true,probST,sout)
    catch(err)
        flush(sout)
        write(sout,"\n ERROR found on optTopRun1 for run with seed $(seed): $(err)")
        write(sout,"\n------------------------------------------------------------\n")
        write(logfile,"\n ERROR found on optTopRun1 for run with seed $(seed): $(err)")
        flush(logfile)
        flush(sout)
    end
    close(sout)
    close(logfile)
    if (debug) DEBUG = false; end  # change back to default
    if (redirect) REDIRECT = false; end
    return res
end


