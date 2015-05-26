# functions for numerical/heuristic optimization
# originally in functions.jl
# Claudia March 2015


const move2int = Dict{Symbol,Int64}([:add=>1,:MVorigin=>2,:MVtarget=>3,:CHdir=>4,:delete=>5, :nni=>6])
const int2move = (Int64=>Symbol)[move2int[k]=>k for k in keys(move2int)]
const fAbs = 1e-6 #1e-10 prof Bates, 1e-6
const fRel = 1e-5 # 1e-12 prof Bates, 1e-5
const xAbs = 1e-4 # 0.001 in phylonet, 1e-10 prof Bates, 1e-4
const xRel = 1e-3 # 0.01 in phylonet, 1e-10 prof Bates, 1e-3
const numFails = 100 # like phylonet
const numMoves = Int64[] #empty to be calculated inside based on coupon's collector
const multiplier = 10000 # for loglik absolute tol (multiplier*fAbs)

# ---------------------- branch length optimization ---------------------------------

# function to get the branch lengths/gammas to optimize for a given network
# warning: order of parameters (h,t,gammaz)
# updates net.numht also with the number of hybrid nodes and number of identifiable edges (n2,n,hzn)
function parameters(net::Network)
    t = Float64[]
    h = Float64[]
    n = Int64[]
    n2 = Int64[]
    hz = Float64[]
    hzn = Int64[]
    indxt = Int64[]
    indxh = Int64[]
    indxhz = Int64[]
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
                    push!(hzn,int(string(string(node.number),"1")))
                    push!(hzn,int(string(string(node.number),"2")))
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
    size(net.numht,1) > 0 || error("net.numht not correctly updated, need to run parameters first")
    size(qnet.indexht,1) == 0 ||  println("deleting qnet.indexht to replace with info in net")
    nh = net.numht[1 : net.numHybrids - net.numBad]
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    nt = net.numht[net.numHybrids - net.numBad + 1 : net.numHybrids - net.numBad + k]
    nhz = net.numht[net.numHybrids - net.numBad + k + 1 : length(net.numht)]
    qnh = Int64[]
    qnt = Int64[]
    qnhz = Int64[]
    qindxh = Int64[]
    qindxt = Int64[]
    qindxhz = Int64[]
    if(qnet.numHybrids == 1 && qnet.hybrid[1].isBadDiamondI)
        ind1 = int(string(string(qnet.hybrid[1].number),"1"))
        ind2 = int(string(string(qnet.hybrid[1].number),"2"))
        i = getIndex(ind1,nhz)
        edges = hybridEdges(qnet.hybrid[1])
        push!(qnhz,i+net.numHybrids-net.numBad+k)
        push!(qnhz,i+1+net.numHybrids-net.numBad+k)
        push!(qindxhz,getIndex(getOtherNode(edges[1],qnet.hybrid[1]),qnet))
        push!(qindxhz,getIndex(getOtherNode(edges[2],qnet.hybrid[1]),qnet))
    else
        all([!n.isBadDiamondI for n in qnet.hybrid]) || error("cannot have bad diamond I hybrid nodes in this qnet, case dealt separately before")
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
            if(!e.istIdentifiable && all([!n.leaf for n in e.node]) && !e.hybrid && e.fromBadDiamondI) # tree edge not identifiable but internal with length!=0 (not bad diamII nor bad triangle)
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
    end
    qnet.indexht = vcat(qnh,qnt,qnhz)
    qnet.index = vcat(qindxh,qindxt,qindxhz)
    length(qnet.indexht) == length(qnet.index) || error("strange in setting qnet.indexht and qnet.index, they do not have same length")
end


# function to compare a vector of parameters with the current vector in net.ht
# to know which parameters were changed
function changed(net::HybridNetwork, x::Vector{Float64})
    if(length(net.ht) == length(x))
        return [!approxEq(net.ht[i],x[i]) for i in 1:length(x)]
    else
        error("net.ht (length $(length(net.ht))) and vector x (length $(length(x))) need to have same length")
    end
end


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
                x[qnet.indexht[1]] + x[qnet.indexht[2]] <= 1 || warn("new gammaz should add to less than 1: $(x[qnet.indexht[1]] + x[qnet.indexht[2]])")
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
                    #println("updating qnet parameters, found gammaz case when hybridization has been removed")
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
    return vcat(ones(net.numHybrids-net.numBad),DataFrames.rep(10,k),ones(length(net.ht)-k-net.numHybrids+net.numBad))
end

# function to calculate the inequality gammaz1+gammaz2 <= 1
function calculateIneqGammaz(x::Vector{Float64}, net::HybridNetwork, ind::Int64, verbose::Bool)
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    hz = x[net.numHybrids - net.numBad + k + 1 : length(x)]
    if(verbose)
        println("enters calculateIneqGammaz with hz $(hz), and hz[ind*2] + hz[ind*2-1] - 1 = $(hz[ind*2] + hz[ind*2-1] - 1)")
    end
    hz[ind*2] + hz[ind*2-1] - 1
end

# numerical optimization of branch lengths given a network (or tree)
# and data (set of quartets with obsCF)
# using BOBYQA from NLopt package
function optBL!(net::HybridNetwork, d::DataCF, verbose::Bool, ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64)
    (ftolRel > 0 && ftolAbs > 0 && xtolAbs > 0 && xtolRel > 0) || error("tolerances have to be positive, ftol (rel,abs), xtol (rel,abs): $([ftolRel, ftolAbs, xtolRel, xtolAbs])")
    println("OPTBL: begin branch lengths and gammas optimization, ftolAbs $(ftolAbs), ftolRel $(ftolRel), xtolAbs $(xtolAbs), xtolRel $(xtolRel)")
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
        if(verbose) #|| net.numBad > 0) #we want to see what happens with bad diamond I
            println("inside obj with x $(x)")
        end
        count += 1
        calculateExpCFAll!(d,x,net) # update qnet branches and calculate expCF
        update!(net,x) # update net.ht
        val = logPseudoLik(d)
        if(verbose) #|| net.numBad > 0)#we want to see what happens with bad diamond I
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
    println("OPTBL: starting point $(ht)")
    fmin, xmin, ret = NLopt.optimize(opt,ht)
    println("got $(round(fmin,5)) at $(round(xmin,5)) after $(count) iterations (returned $(ret))")
    updateParameters!(net,xmin)
    net.loglik = fmin
    #return fmin,xmin
end

optBL!(net::HybridNetwork, d::DataCF) = optBL!(net, d, false, fRel, fAbs, xRel, xAbs)
optBL!(net::HybridNetwork, d::DataCF, verbose::Bool) = optBL!(net, d,verbose, fRel, fAbs, xRel, xAbs)
optBL!(net::HybridNetwork, d::DataCF, ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64) = optBL!(net, d, false, ftolRel, ftolAbs, xtolRel, xtolAbs)


# function that will add a hybridization with addHybridizationUpdate,
# if success=false, it will try to move the hybridization before
# declaring failure
# blacklist used in afterOptBLAll
function addHybridizationUpdateSmart!(net::HybridNetwork, blacklist::Bool, N::Int64)
    println("MOVE: addHybridizationUpdateSmart")
    success, hybrid, flag, nocycle, flag2, flag3 = addHybridizationUpdate!(net, blacklist)
    i = 0
    if(!success)
        while((nocycle || !flag) && i < N) #incycle failed
            warn("MOVE: added hybrid causes conflict with previous cycle, need to delete and add another")
            deleteHybrid!(hybrid,net,true)
            success, hybrid, flag, nocycle, flag2, flag3 = addHybridizationUpdate!(net, blacklist)
        end
        if(nocycle || !flag)
            warn("MOVE: added hybridization $(i) times trying to avoid incycle conflicts, but failed")
        else
            if(!flag2) #gammaz failed
                warn("MOVE: added hybrid has problem with gammaz (not identifiable bad triangle)")
                flag3, edgesRoot = updateContainRoot!(net,hybrid);
                if(flag3)
                    warn("MOVE: we will move origin to fix the gammaz situation")
                    success = moveOriginUpdateRepeat!(net,hybrid,true)
                else
                    warn("MOVE: we will move target to fix the gammaz situation")
                    success = moveTargetUpdateRepeat!(net,hybrid,true)
                end
            else
                if(!flag3) #containRoot failed
                    warn("MOVE: added hybrid causes problems with containRoot, will change the direction to fix it")
                    success = changeDirectionUpdate!(net,hybrid) #change dir of minor
                end
            end
        end
        if(!success)
            warn("MOVE: could not fix the added hybrid by any means, we will delete it now")
            deleteHybridizationUpdate!(net,hybrid)
        end
    end
    !success || println("MOVE: added hybridization SUCCESSFUL: new hybrid $(hybrid.number)")
    return success
end

addHybridizationUpdateSmart!(net::HybridNetwork, N::Int64) = addHybridizationUpdateSmart!(net, false,N)

# function to delete a hybrid, and then add a new hybrid:
# deleteHybridizationUpdate and addHybridizationUpdate,
# close=true will try move origin/target, if false, will delete/add new hybrid
# default is close=true
# origin=true, moves origin, if false, moves target. option added to control not keep coming
# to the same network over and over
# returns success
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
function moveHybrid!(net::HybridNetwork, edge::Edge, close::Bool, origin::Bool,N::Int64, movesgamma::Vector{Int64})
    edge.hybrid || error("edge $(edge.number) cannot be deleted because it is not hybrid")
    node = edge.node[edge.isChild1 ? 1 : 2];
    node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
    println("MOVE: moving hybrid for edge $(edge.number)")
    if(close)
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
# close=true will try move origin/target, if false, will delete/add new hybrid
# origin=true, moves origin, if false, moves target. option added to control not keep coming
# to the same network over and over
# returns success
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
function gammaZero!(net::HybridNetwork, d::DataCF, edge::Edge, close::Bool, origin::Bool, N::Int64, movesgamma::Vector{Int64})
    currTloglik = net.loglik
    edge.hybrid || error("edge $(edge.number) should be hybrid edge because it corresponds to a gamma (or gammaz) in net.ht")
    println("gamma zero situation found for hybrid edge $(edge.number) with gamma $(edge.gamma)")
    node = edge.node[edge.isChild1 ? 1 : 2];
    node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
    success = changeDirectionUpdate!(net,node) #changes dir of minor
    movesgamma[4] += 1
    if(success)
        movesgamma[10] += 1
        optBL!(net,d,false)
        flags = isValid(net)
        if(net.loglik <= currTloglik && flags[1] && flags[3])
            println("changing direction fixed the gamma zero situation")
            success2 = true
        else
            println("changing direction does not fix the gamma zero situation, need to undo change direction and move hybrid")
            success = changeDirectionUpdate!(net,node)
            success || error("strange thing, changed direction and success, but lower loglik; want to undo changeDirection, and success=false! Hybrid node is $(node.number)")
            success2 = moveHybrid!(net,edge,close,origin,N, movesgamma)
        end
    else
        println("changing direction was not possible to fix the gamma zero situation (success=false), need to move hybrid")
        success2 = moveHybrid!(net,edge,close,origin,N,movesgamma)
    end
    return success2
end

# function to check if h or t (in currT.ht) are 0 (or 1 for h)
# close=true will try move origin/target, if false, will delete/add new hybrid
# origin=true, moves origin, if false, moves target. option added to control not keep coming
# to the same network over and over
# returns successchange=false if could not add new hybrid; true ow
# returns successchange,flagh,flagt,flaghz
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
function afterOptBL!(currT::HybridNetwork, d::DataCF,close::Bool, origin::Bool,verbose::Bool, N::Int64, movesgamma::Vector{Int64})
    !isTree(currT) || return false,true,true,true
    nh = currT.ht[1 : currT.numHybrids - currT.numBad]
    k = sum([e.istIdentifiable ? 1 : 0 for e in currT.edge])
    nt = currT.ht[currT.numHybrids - currT.numBad + 1 : currT.numHybrids - currT.numBad + k]
    nhz = currT.ht[currT.numHybrids - currT.numBad + k + 1 : length(currT.ht)]
    indh = currT.index[1 : currT.numHybrids - currT.numBad]
    indt = currT.index[currT.numHybrids - currT.numBad + 1 : currT.numHybrids - currT.numBad + k]
    indhz = currT.index[currT.numHybrids - currT.numBad + k + 1 : length(currT.ht)]
    flagh,flagt,flaghz = isValid(nh,nt,nhz)
    !all([flagh,flagt,flaghz]) || return false,true,true,true
    println("begins afterOptBL because of conflicts: flagh,flagt,flaghz=$([flagh,flagt,flaghz])")
    successchange = true
    if(!flagh)
        for(i in 1:length(nh))
            if(approxEq(nh[i],0.0) || approxEq(nh[i],1.0))
                edge = currT.edge[indh[i]]
                approxEq(edge.gamma,nh[i]) || error("edge $(edge.number) gamma $(edge.gamma) should match the gamma in net.ht $(nh[i]) and it does not")
                successchange = gammaZero!(currT, d,edge,close,origin,N,movesgamma)
                !successchange || break
            end
        end
    elseif(!flagt)
        for(i in 1:length(nt))
            if(approxEq(nt[i],0.0))
                edge = currT.edge[indt[i]]
                approxEq(edge.length,nt[i]) || error("edge $(edge.number) length $(edge.length) should match the length in net.ht $(nt[i]) and it does not")
                if(all([!n.hasHybEdge for n in edge.node]))
                    movesgamma[6] += 1
                    successchange = NNI!(currT,edge)
                    movesgamma[12] += successchange ? 1 : 0
                    !successchange || break
                else
                    println("MOVE: need NNI on edge $(edge.number) because length is $(edge.length), but edge is attached to nodes that have hybrid edges")
                    successchange = false
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
                successchange = gammaZero!(currT,d,edges[1],close,origin,N,movesgamma)
                break
            elseif(approxEq(nhz[i],1.0))
                approxEq(nhz[i+1],0.0) || error("gammaz for node $(currT.node[indhz[i]].number) is $(nhz[i]) but the other gammaz is $(nhz[i+1]), the sum should be less than 1.0")
                nodehz = currT.node[indhz[i+1]]
                approxEq(nodehz.gammaz,nhz[i+1]) || error("nodehz $(nodehz.number) gammaz $(nodehz.gammaz) should match the gammaz in net.ht $(nhz[i+1]) and it does not")
                edges = hybridEdges(nodehz)
                edges[1].hybrid || error("bad diamond I situation, node $(nodehz.number) has gammaz $(nodehz.gammaz) so should be linked to hybrid edge, but it is not")
                successchange = gammaZero!(currT,d,edges[1],close,origin,N,movesgamma)
                break
            else
                if(approxEq(nhz[i+1],0.0))
                    nodehz = currT.node[indhz[i+1]];
                    approxEq(nodehz.gammaz,nhz[i+1]) || error("nodehz $(nodehz.number) gammaz $(nodehz.gammaz) should match the gammaz in net.ht $(nhz[i+1]) and it does not")
                    edges = hybridEdges(nodehz);
                    edges[1].hybrid || error("bad diamond I situation, node $(nodehz.number) has gammaz $(nodehz.gammaz) so should be linked to hybrid edge, but it is not")
                    successchange = gammaZero!(currT,d,edges[1],close,origin,N,movesgamma)
                    break
                elseif(approxEq(nhz[i+1],1.0))
                    error("gammaz for node $(currT.node[indhz[i]].number) is $(nhz[i]) but the other gammaz is $(nhz[i+1]), the sum should be less than 1.0")
                end
            end
            i += 2
        end
    end
    if(successchange)
        #println(writeTopology(currT))
        #printEdges(currT)
        optBL!(currT,d,verbose)
    end
    !successchange || println("afterOptBL SUCCESSFUL change, need to run again to see if new topology is valid")
    if(successchange)
        flagh,flagt,flaghz = isValid(currT)
        println("new flags: flagh,flagt,flaghz $([flagh,flagt,flaghz])")
    end
    return successchange,flagh,flagt,flaghz
end

# function to repeat afterOptBL every time after changing something
# N: number of times it will delete/add hybrid if close=false
# origin=true, moves origin, if false, moves target. option added to control not keep coming
# to the same network over and over
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
function afterOptBLRepeat!(currT::HybridNetwork, d::DataCF, N::Int64,close::Bool, origin::Bool, verbose::Bool, movesgamma::Vector{Int64})
    success,flagh,flagt,flaghz = afterOptBL!(currT,d,close,origin,verbose,N, movesgamma)
    println("inside afterOptBLRepeat, after afterOptBL once, we get: success, flags: $([success,flagh,flagt,flaghz])")
    if(!close)
        println("close is $(close), and gets inside the loop for repeating afterOptBL")
        i = 1
        while(success && !all([flagh,flagt,flaghz]) && i<N) #fixit: do we want to check loglik in this step also?
            success,flagh,flagt,flaghz = afterOptBL!(currT,d,close,origin,verbose,N,movesgamma)
            i += 1
        end
        i < N || println("tried afterOptBL $(i) times")
    end
    return success,flagh,flagt,flaghz
end

# function to check if we have to keep checking loglik
# returns true if we have to stop because loglik not changing much anymore, or it is close to 0.0 already
function stopLoglik(currloglik::Float64, newloglik::Float64, ftolAbs::Float64, M::Number)
    return (abs(currloglik-newloglik) <= M*ftolAbs) || (newloglik <= M*ftolAbs)
end

# function to repeat afterOptBL every time after changing something
# close=true will try move origin/target, if false, will delete/add new hybrid
# default is close=true
# returns new approved currT (no gammas=0.0)
# N: number of times failures of accepting loglik is allowed, M: multiplier for loglik tolerance (M*ftolAbs)
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
function afterOptBLAll!(currT::HybridNetwork, d::DataCF, N::Int64,close::Bool, M::Number, ftolAbs::Float64, verbose::Bool, movesgamma::Vector{Int64},ftolRel::Float64, xtolRel::Float64, xtolAbs::Float64)
    println("afterOptBLAll: checking if currT has gamma(z) = 0.0(1.0): currT.ht $(currT.ht)")
    currloglik = currT.loglik
    currT.blacklist = Int64[];
    origin = (rand() > 0.5) #true=moveOrigin, false=moveTarget
    startover = true
    tries = 0
    N2 = N > 5 ? N/5 : 1 #num of failures of badlik around a gamma=0.0
    while(startover && tries < N)
        tries += 1
        println("inside afterOptBLALL: number of tries $(tries) out of $(N) possible")
        badliks = 0
        if(currT.loglik < M*ftolAbs) #curr loglik already close to 0.0
            startover = false
        else
            backCurrT0 = false
            while(badliks < N2) #will try a few options around currT
                println("tried $(badliks) bad likelihood options so far out of $(N2)")
                currT0 = deepcopy(currT)
                origin = !origin #to guarantee not going back to previous topology
                success,flagh,flagt,flaghz = afterOptBLRepeat!(currT,d,N,close,origin,verbose,movesgamma)
                println("inside afterOptBLAll, after afterOptBLRepeat once we get: success, flags: $([success,flagh,flagt,flaghz])")
                if(!success) #tried to change something but failed
                    println("did not change anything inside afterOptBL: could be nothing needed change or tried but couldn't anymore. flagh, flagt, flaghz = $([flagh,flagt,flaghz])")
                    if(all([flagh,flagt,flaghz])) #currT was ok
                        startover = false
                    elseif(!flagh || !flaghz) #currT was bad but could not change it, need to go down a level
                        !isTree(currT) || error("afterOptBL should not give reject=true for a tree")
                        println("current topology has numerical parameters that are not valid: gamma=0(1), gammaz=0(1); need to move down a level h-1")
                        moveDownLevel!(currT)
                        optBL!(currT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
                        startover = true
                    elseif(!flagt)
                        startover = false
                    end
                else #changed something
                    println("changed something inside afterOptBL: flagh, flagt, flaghz = $([flagh,flagt,flaghz]). oldloglik $(currloglik), newloglik $(currT.loglik)")
                    #printEdges(currT)
                    #printNodes(currT)
                    println(writeTopology(currT))
                    if(currT.loglik > currloglik) #|| abs(currT.loglik-currloglik) <= M*ftolAbs) #fixit: allowed like this because of changeDir that does not change much the lik but can fix h=0
                        println("worse likelihood, back to currT")
                        startover = true
                        backCurrT0 = true
                    else
                        println("better likelihood, jump to new topology and startover")
                        movesgamma[13] += 1
                        if(all([flagh,flagt,flaghz]))
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
                println("tried to fix gamma zero situation for $(badliks) times and could not")
                flagh,flagt,flaghz = isValid(currT)
                if(!flagh || !flaghz)
                    !isTree(currT) || error("afterOptBL should not give reject=true for a tree")
                    println("current topology has numerical parameters that are not valid: gamma=0(1), t=0, gammaz=0(1); need to move down a level h-1")
                    movesgamma[5] += 1
                    movesgamma[11] += 1
                    moveDownLevel!(currT)
                    #printEdges(currT)
                    #printNodes(currT)
                    println(writeTopology(currT))
                    optBL!(currT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
                    startover = true
                else
                    println("the only problem were lengths equal to zero, so we will keep them")
                    startover = false
                end
            end
        end
    end
    if(tries >= N)
        println("afterOptBLAll ended because it tried $(tries) times with startover $(startover)")
        writeTopology(currT)
        flagh,flagt,flaghz = isValid(currT)
        if(!flagh || !flaghz)
            println("gammaz zero situation still in currT, need to move down one level to h-1")
            moveDownLevel!(currT)
            #printEdges(currT)
            #printNodes(currT)
            println(writeTopology(currT))
            optBL!(currT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
        end
    end
    currT.blacklist = Int64[];
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
function adjustWeight(net::HybridNetwork,hmax::Int64,w::Vector{Float64})
    if(hmax - net.numHybrids > 0)
        hmax >= 0 || error("hmax must be non negative: $(hmax)")
        length(w) == 6 || error("length of w should be 6 as there are only 6 moves: $(w)")
        approxEq(sum(w),1.0) || error("vector of move weights should add up to 1: $(w),$(sum(w))")
        all([0<=i<=1 for i in w]) || error("weights must be nonnegative and less than one $(w)")
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
function adjustWeightMovesfail!(v::Vector{Float64}, movesfail::Vector{Int64}, Nmov::Vector{Int64}, net::HybridNetwork, hmax::Int64)
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
function whichMove(net::HybridNetwork,hmax::Int64,w::Vector{Float64}, dynamic::Bool, movesfail::Vector{Int64}, Nmov::Vector{Int64})
    hmax >= 0 || error("hmax must be non negative: $(hmax)")
    length(w) == 6 || error("length of w should be 6 as there are only 6 moves: $(w)")
    approxEq(sum(w),1.0) || error("vector of move weights should add up to 1: $(w),$(sum(w))")
    all([0<=i<=1 for i in w]) || error("weights must be nonnegative and less than one $(w)")
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
        println("weights before adjusting by movesfail $(v)")
        flag = adjustWeightMovesfail!(v,movesfail,Nmov,net,hmax)
        println("weights after adjusting by movesfail $(v)")
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

whichMove(net::HybridNetwork,hmax::Int64,movesfail::Vector{Int64}, Nmov::Vector{Int64}) = whichMove(net,hmax,[1/5,1/5,1/5,1/5,0.0,1/5], true,movesfail, Nmov)
whichMove(net::HybridNetwork,hmax::Int64,w::Vector{Float64},movesfail::Vector{Int64}, Nmov::Vector{Int64}) = whichMove(net,hmax,w, true,movesfail, Nmov)

#function to choose a hybrid node for the given moves
function chooseHybrid(net::HybridNetwork)
    !isTree(net) || error("net is a tree, cannot choose hybrid node")
    net.numHybrids > 1 || return net.hybrid[1]
    index1 = 0
    while(index1 == 0 || index1 > size(net.hybrid,1))
        index1 = iround(rand()*size(net.hybrid,1));
    end
    #println("chosen hybrid node for network move: $(net.hybrid[index1].number)")
    return net.hybrid[index1]
end

# function to propose a new topology given a move
# random = false uses the minor hybrid edge always
# count to know in which step we are, N for NNI trials
# order in movescount as in IF here (add,mvorigin,mvtarget,chdir,delete,nni)
function proposedTop!(move::Integer, newT::HybridNetwork,random::Bool, count::Int64, N::Int64, movescount::Vector{Int64}, movesfail::Vector{Int64})
    1 <= move <= 6 || error("invalid move $(move)") #fixit: if previous move rejected, do not redo it!
    println("current move: $(int2move[move])")
    if(move == 1)
        success = addHybridizationUpdateSmart!(newT,N)
    elseif(move == 2)
        node = chooseHybrid(newT)
        success = moveOriginUpdateRepeat!(newT,node,random)
    elseif(move == 3)
        node = chooseHybrid(newT)
        success = moveTargetUpdateRepeat!(newT,node,random)
    elseif(move == 4)
        node = chooseHybrid(newT)
        success = changeDirectionUpdate!(newT,node, random)
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
    println("success $(success), movescount (add,mvorigin,mvtarget,chdir,delete,nni) proposed: $(movescount[1:6]); successful: $(movescount[7:12]); movesfail: $(movesfail)")
    !success || return true
    println("new proposed topology failed in step $(count) for move $(int2move[move])")
    return false
end


proposedTop!(move::Symbol, newT::HybridNetwork, random::Bool, count::Int64,N::Int64, movescount::Vector{Int64},movesfail::Vector{Int64}) = proposedTop!(try move2int[move] catch error("invalid move $(string(move))") end,newT, random,count,N, movescount,movesfail)

# function to calculate Nmov, number max of tries per move
# order: (add,mvorigin,mvtarget,chdir,delete,nni)
function calculateNmov!(net::HybridNetwork, N::Vector{Int64})
    if(isempty(N))
        N = rep(0,6)
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
# close=true if gamma=0.0 fixed only around neighbors with move origin/target
# ret=true: return the network, default is false
function optTopLevel!(currT::HybridNetwork, M::Number, Nfail::Int64, d::DataCF, hmax::Int64,ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64, verbose::Bool, close::Bool, Nmov0::Vector{Int64},ret::Bool)
    M > 0 || error("M must be greater than zero: $(M)")
    Nfail > 0 || error("Nfail must be greater than zero: $(Nfail)")
    isempty(Nmov0) || all([n > 0 for n in Nmov0]) || error("Nmov must be greater than zero: $(Nmov0)")
    count = 0
    movescount = rep(0,18) #1:6 number of times moved proposed, 7:12 number of times success move (no intersecting cycles, etc.), 13:18 accepted by loglik
    movesgamma = rep(0,13) #number of moves to fix gamma zero: proposed, successful, movesgamma[13]: total accepted by loglik
    movesfail = rep(0,6) #count of failed moves for current topology
    failures = 0
    stillmoves = true
    if(isempty(Nmov0))
        Nmov = rep(0,6)
    else
        Nmov = deepcopy(Nmov0)
    end
    optBL!(currT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
    currT = afterOptBLAll!(currT, d, Nfail,close, M, ftolAbs, verbose,movesgamma,ftolRel,xtolRel,xtolAbs)
    absDiff = M*ftolAbs + 1
    newT = deepcopy(currT)
    #printEdges(newT)
    #printNodes(newT)
    println(writeTopology(newT))
    while(absDiff > M*ftolAbs && failures < Nfail && currT.loglik > M*ftolAbs && stillmoves) #stops if close to zero because of new deviance form of the pseudolik
        count += 1
        println("--------- loglik_$(count) = $(round(currT.loglik,6)) -----------")
        if(isempty(Nmov0)) #if empty, not set by user
            calculateNmov!(newT,Nmov)
        end
        println("will propose move with movesfail $(movesfail), Nmov $(Nmov)")
        move = whichMove(newT,hmax,movesfail,Nmov)
        if(move != :none)
            flag = proposedTop!(move,newT,true, count,10, movescount,movesfail) #N=10 because with 1 it never finds an edge for nni
            if(flag)
                accepted = false
                println("accepted proposed new topology in step $(count)")
                #printEdges(newT)
                #printNodes(newT)
                println(writeTopology(newT))
                optBL!(newT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
                println("OPT: comparing newT.loglik $(newT.loglik), currT.loglik $(currT.loglik)")
                if(newT.loglik < currT.loglik && abs(newT.loglik-currT.loglik) > M*ftolAbs) #newT better loglik: need to check for error or keeps jumping back and forth
                    newloglik = newT.loglik
                    newT = afterOptBLAll!(newT, d, Nfail,close, M, ftolAbs,verbose,movesgamma,ftolRel, xtolRel,xtolAbs)
                    println("loglik before afterOptBL $(newloglik), newT.loglik now $(newT.loglik), loss in loglik by fixing gamma(z)=0.0(1.0): $(newloglik>newT.loglik ? 0 : abs(newloglik-newT.loglik))")
                    accepted = true
                else
                    accepted = false
                end
                if(accepted)
                    absDiff = abs(newT.loglik - currT.loglik)
                    println("proposed new topology with better loglik in step $(count): oldloglik=$(round(currT.loglik,3)), newloglik=$(round(newT.loglik,3)), after $(failures) failures")
                    currT = deepcopy(newT)
                    failures = 0
                    movescount[move2int[move]+12] += 1
                    movesfail = rep(0,6) #count of failed moves for current topology
                else
                    println("rejected new topology with worse loglik in step $(count): currloglik=$(round(currT.loglik,3)), newloglik=$(round(newT.loglik,3)), with $(failures) failures")
                    failures += 1
                    movesfail[move2int[move]] += 1
                    newT = deepcopy(currT)
                end
                #printEdges(newT)
                #printNodes(newT)
                println(writeTopology(newT))
                println("ends step $(count) with absDiff $(accepted? absDiff : 0.0) and failures $(failures)")
            end
        else
            stillmoves = false
        end
    end
    if(ftolAbs > 1e-7 || ftolRel > 1e-7 || xtolAbs > 1e-7 || xtolRel > 1e-7)
        println("found best network, now we re-optimize branch lengths and gamma more precisely")
        optBL!(newT,d,verbose,1e-12,1e-10,1e-10,1e-10)
    end
    if(absDiff <= M*ftolAbs)
        println("STOPPED by absolute difference criteria")
    elseif(currT.loglik <= M*ftolAbs)
        println("STOPPED by loglik close to zero criteria")
    elseif(!stillmoves)
        println("STOPPED for not having more moves to propose: movesfail $(movesfail), Nmov $(Nmov)")
    else
        println("STOPPED by number of failures criteria")
    end
    if(newT.loglik > M*ftolAbs) #not really close to 0.0, based on absTol also
        println("newT.loglik $(newT.loglik) not really close to 0.0 based on loglik abs. tol. $(M*ftolAbs), you might need to redo with another starting point")
    end
    if(newT.numBad > 0) #need to undogammaz if newT has bad diamond I to use gammaz as proxy of gamma for writeTopology
        for(n in newT.hybrid)
            if(n.isBadDiamondI)
                undoGammaz!(n,newT)
            end
        end
    end
    println("END optTopLevel: found minimizer topology at step $(count) (failures: $(failures)) with -loglik=$(round(newT.loglik,5)) and ht_min=$(round(newT.ht,5))")
    printCounts(movescount,movesgamma)
    #printEdges(newT)
    #printNodes(newT)
    println(writeTopology(newT))
    !ret || return newT
end

optTopLevel!(currT::HybridNetwork, d::DataCF, hmax::Int64) = optTopLevel!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false,true,numMoves,false)
optTopLevel!(currT::HybridNetwork, d::DataCF, hmax::Int64,ret::Bool) = optTopLevel!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false,true,numMoves,ret)
optTopLevel!(currT::HybridNetwork, d::DataCF, hmax::Int64, verbose::Bool) = optTopLevel!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, verbose,true,numMoves,false)
optTopLevel!(currT::HybridNetwork, d::DataCF, hmax::Int64, verbose::Bool,ret::Bool) = optTopLevel!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, verbose,true,numMoves,ret)
optTopLevel!(currT::HybridNetwork, M::Number, Nfail::Int64, d::DataCF, hmax::Int64,ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64, verbose::Bool, close::Bool, Nmov0::Vector{Int64}) = optTopLevel!(currT, M, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, close, Nmov0, false)



# function to print count of number of moves after optTopLeve
function printCounts(movescount::Vector{Int64}, movesgamma::Vector{Int64})
    length(movescount) == 18 || error("movescount should have length 18, not $(length(movescount))")
    length(movesgamma) == 13 || error("movesgamma should have length 13, not $(length(movescount))")
    println("PERFORMANCE: total number of moves (proposed, successful, accepted) in general, and to fix gamma=0.0,t=0.0 cases")
    println("\t--------moves general--------\t\t\t\t --------moves gamma--------\t")
    println("move\t Num.Proposed\t Num.Successful\t Num.Accepted\t | Num.Proposed\t Num.Successful\t Num.Accepted")
    names = ["add","mvorigin","mvtarget","chdir","delete","nni"]
    for(i in 1:6)
        if(i == 2 || i == 3)
            println("$(names[i])\t $(movescount[i])\t\t $(movescount[i+6])\t\t $(movescount[i+12])\t |\t $(movesgamma[i])\t\t $(movesgamma[i+6])\t\t --")
        else
            println("$(names[i])\t\t $(movescount[i])\t\t $(movescount[i+6])\t\t $(movescount[i+12])\t |\t $(movesgamma[i])\t\t $(movesgamma[i+6])\t\t --")
        end
    end
    suma = sum(movescount[1:6]);
    suma2 = sum(movesgamma[1:6]) == 0 ? 1 : sum(movesgamma[1:6])
    println("Total\t\t $(sum(movescount[1:6]))\t\t $(sum(movescount[7:12]))\t\t $(sum(movescount[13:18]))\t |\t $(sum(movesgamma[1:6]))\t\t $(sum(movesgamma[7:12]))\t\t $(movesgamma[13])")
    println("Proportion\t $(round(sum(movescount[1:6])/suma,1))\t\t $(round(sum(movescount[7:12])/suma,1))\t\t $(round(sum(movescount[13:18])/suma,1))\t |\t $(round(sum(movesgamma[1:6])/suma2,1))\t\t $(round(sum(movesgamma[7:12])/suma2,1))\t\t $(round(movesgamma[13]/suma2,1))")
end

# function to move down onw level to h-1
# caused by gamma=0,1 or gammaz=0,1
function moveDownLevel!(net::HybridNetwork)
    !isTree(net) ||error("cannot delete hybridization in a tree")
    println("MOVE: need to go down one level to h-1=$(net.numHybrids-1) hybrids because of conflicts with gamma=0,1")
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
                approxEq(nhz[i+1],0.0) || error("gammaz for node $(net.node[indhz[i]].number) is $(nhz[i]) but the other gammaz is $(nhz[i+1]), the sum should be less than 1.0")
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
                    error("gammaz for node $(net.node[indhz[i]].number) is $(nhz[i]) but the other gammaz is $(nhz[i+1]), the sum should be less than 1.0")
                end
            end
            i += 2
        end
    end
end

# checks if there are problems in estimated net.ht:
# returns flag for h, flag for t, flag for hz
function isValid(net::HybridNetwork)
    nh = net.ht[1 : net.numHybrids - net.numBad]
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    nt = net.ht[net.numHybrids - net.numBad + 1 : net.numHybrids - net.numBad + k]
    nhz = net.ht[net.numHybrids - net.numBad + k + 1 : length(net.ht)]
    #println("isValid on nh $(nh), nt $(nt), nhz $(nhz)")
    return all([(0<n<1 && !approxEq(n,0.0) && !approxEq(n,1.0)) for n in nh]), all([(n>0 && !approxEq(n,0.0)) for n in nt]), all([(0<n<1 && !approxEq(n,0.0) && !approxEq(n,1.0)) for n in nhz])
end

# checks if there are problems in estimated net.ht:
# returns flag for h, flag for t, flag for hz
function isValid(nh::Vector{Float64},nt::Vector{Float64},nhz::Vector{Float64})
    #println("isValid on nh $(nh), nt $(nt), nhz $(nhz)")
    return all([(0<n<1 && !approxEq(n,0.0) && !approxEq(n,1.0)) for n in nh]), all([(n>0 && !approxEq(n,0.0)) for n in nt]), all([(0<n<1 && !approxEq(n,0.0) && !approxEq(n,1.0)) for n in nhz])
end

# function to do optTopLevel for each level below hmax
# returns array of HybridNetworks with currT in first position, and best topology for each level later
# uses the best topology in h-1 as starting point for the search in h
# no test in between! fixit: add a test to decide to move from h-1 to h?
function optTop!(currT::HybridNetwork, M::Number, Nfail::Int64, d::DataCF, hmax::Int64,ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64, verbose::Bool, close::Bool, Nmov0::Vector{Int64},ret::Bool)
    hmax >= 0 || error("hmax cannot be negative $(hmax)")
    currT.numHybrids == 0 || warn("currT has already $(currT.numHybrids), so search will not start in tree, but in the space of networks with $(currT.numHybrids) hybrids")
    currT.numHybrids <= hmax || error("currT has more hybrids: $(currT.numHybrids) than hmax $(hmax)")
    bestT = HybridNetwork[]
    push!(bestT,currT)
    i = 1
    for(h in currT.numHybrids:hmax)
        startT = deepcopy(bestT[i])
        newT = optTopLevel!(startT, M, Nfail, d, h,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, close, Nmov0,true)
        println("BEST NETWORK at $(h) hybridizations: $(writeTopology(newT))")
        push!(bestT,deepcopy(newT))
        i += 1
    end
    !ret || return bestT
end



optTop!(currT::HybridNetwork, d::DataCF, hmax::Int64) = optTop!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false,true,numMoves,false)
optTop!(currT::HybridNetwork, d::DataCF, hmax::Int64,ret::Bool) = optTop!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false,true,numMoves,ret)
optTop!(currT::HybridNetwork, d::DataCF, hmax::Int64, verbose::Bool) = optTop!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, verbose,true,numMoves,false)
optTop!(currT::HybridNetwork, d::DataCF, hmax::Int64, verbose::Bool,ret::Bool) = optTop!(currT, multiplier, numFails, d, hmax,fRel, fAbs, xRel, xAbs, verbose,true,numMoves,ret)
optTop!(currT::HybridNetwork, M::Number, Nfail::Int64, d::DataCF, hmax::Int64,ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64, verbose::Bool, close::Bool, Nmov0::Vector{Int64}) = optTop!(currT, M, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, close, Nmov0, false)
