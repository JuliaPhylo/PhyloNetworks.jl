# functions for numerical/heuristic optimization
# originally in functions.jl
# Claudia March 2015


const move2int = Dict{Symbol,Int}(:add=>1,:MVorigin=>2,:MVtarget=>3,:CHdir=>4,:delete=>5, :nni=>6)
const int2move = Dict{Int,Symbol}(move2int[k]=>k for k in keys(move2int))
"""
default values for tolerance parameters,
used in the optimization of branch lengths (fAbs, fRel, xAbs, xRel)
and in the acceptance of topologies (likAbs, numFails).

if changes are made here, **make the same** in the docstring for snaq! below

version | fAbs | fRel | xAbs | xRel | numFails | likAbs | multiplier
--------|------|------|------|------|----------|--------|-----------
v0.5.1  | 1e-6 | 1e-6 | 1e-3 | 1e-2 |     75   |  1e-6  |
v0.3.0  | 1e-6 | 1e-5 | 1e-4 | 1e-3 |    100   |  0.01  |
v0.0.1  | 1e-6 | 1e-5 | 1e-4 | 1e-3 |    100   |        | 10000
older   | 1e-10| 1e-12| 1e-10| 1e-10|

v0.5.1: based on Nan Ji's work. same xAbs and xRel as in phylonet (as of 2015).
earlier: multiplier was used; later: likAbs = multiplier*fAbs)
"older": values from GLM.jl, Prof Bates

default values used on a single topology, to optimize branch lengths
and gammas, at the very end of snaq!, and by
topologyMaxQPseudolik! since v0.5.1.

version | fAbsBL | fRelBL | xAbsBL | xRelBL
--------|--------|--------|--------|-------
v0.0.1  | 1e-10  | 1e-12  | 1e-10  | 1e-10
"""
const fAbs = 1e-6
const fRel = 1e-6
const xAbs = 1e-3
const xRel = 1e-2
const numFails = 75 # number of failed proposals allowed before stopping the procedure (like phylonet)
const numMoves = Int[] #empty to be calculated inside based on coupon's collector
const likAbs = 1e-6 # loglik absolute tolerance to accept new topology

const fAbsBL = 1e-10
const fRelBL = 1e-12
const xAbsBL = 1e-10
const xRelBL = 1e-10
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
    for e in net.edge
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
    # size(t,1) > 0 || @warn "net does not have identifiable branch lengths"
    return vcat(h,t,hz),vcat(n2,n,hzn),vcat(indxh,indxt,indxhz)
end

function parameters!(net::Network)
    #@warn "deleting net.ht,net.numht and updating with current edge lengths (numbers)"
    net.ht,net.numht,net.index = parameters(net)
    return net.ht
end


# function to update qnet.indexht,qnet.index based on net.numht
# warning: assumes net.numht is updated already with parameters!(net)
function parameters!(qnet::QuartetNetwork, net::HybridNetwork)
    size(net.numht,1) > 0 || error("net.numht not correctly updated, need to run parameters first")
    @debug (size(qnet.indexht,1) == 0 ? "" :
        "deleting qnet.indexht to replace with info in net")
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
        i = findfirst(isequal(ind1), nhz)
        i != nothing || error("ind1 not found in nhz")
        edges = hybridEdges(qnet.hybrid[1])
        push!(qnhz,i+net.numHybrids-net.numBad+k)
        push!(qnhz,i+1+net.numHybrids-net.numBad+k)
        push!(qindxhz,getIndex(getOtherNode(edges[1],qnet.hybrid[1]),qnet))
        push!(qindxhz,getIndex(getOtherNode(edges[2],qnet.hybrid[1]),qnet))
    else
        found = false
        for n in qnet.hybrid
            if(n.isBadDiamondI)
                ind1 = parse(Int,string(string(n.number),"1"))
                ind2 = parse(Int,string(string(n.number),"2"))
                i = findfirst(isequal(ind1), nhz)
                i != nothing || error("ind1 not found in nhz")
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
            for e in qnet.edge
                if(e.istIdentifiable)
                    enum_in_nt = findfirst(isequal(e.number), nt)
                    if isnothing(enum_in_nt)
                        error("identifiable edge $(e.number) in qnet not found in net")
                    end
                    push!(qnt, enum_in_nt + net.numHybrids - net.numBad)
                    push!(qindxt, getIndex(e,qnet))
                end
                if(!e.istIdentifiable && all((n->!n.leaf),e.node) && !e.hybrid && e.fromBadDiamondI) # tree edge not identifiable but internal with length!=0 (not bad diamII nor bad triangle)
                    enum_in_nhz = findfirst(isequal(e.number), nhz)
                    if isnothing(enum_in_nhz)
                        error("internal edge $(e.number) corresponding to gammaz in qnet not found in net.ht")
                    end
                    push!(qnhz, enum_in_nhz + net.numHybrids - net.numBad + k)
                    push!(qindxhz, getIndex(e,qnet))
                end
                if(e.hybrid && !e.isMajor)
                    node = e.node[e.isChild1 ? 1 : 2]
                    node.hybrid || error("strange hybrid edge $(e.number) poiting to tree node $(node.number)")
                    enum_in_nh = findfirst(isequal(e.number), nh)
                    found = (enum_in_nh != nothing)
                    found ? push!(qnh, enum_in_nh) : nothing
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
        #@debug "inside changed with net.ht $(net.ht) and x $(x)"
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
    for i in 1:length(ch)
        qnet.changed |= (ch[i] & qnet.hasEdge[i])
    end
    #DEBUGC && @debug "inside update!, qnet.changed is $(qnet.changed), ch $(ch) and qnet.hasEdge $(qnet.hasEdge), $(qnet.quartetTaxon), numHyb $(qnet.numHybrids)"
    if(qnet.changed)
        if(any([n.isBadDiamondI for n in qnet.hybrid])) # qnet.indexht is only two values: gammaz1,gammaz2 #FIXIT: this could crash if hybrid for bad diamond should disappear after cleaning qnet
            @debug "it is inside update! and identifies that ht changed and it is inside the bad diamond I case"
            length(qnet.indexht) == 2 || error("strange qnet from bad diamond I with hybrid node, it should have only 2 elements: gammaz1,gammaz2, not $(length(qnet.indexht))")
            for i in 1:2
                0 <= x[qnet.indexht[i]] <= 1 || error("new gammaz value should be between 0,1: $(x[qnet.indexht[i]]).")
                @debug (x[qnet.indexht[1]] + x[qnet.indexht[2]] <= 1 ? "" :
                        "warning: new gammaz should add to less than 1: $(x[qnet.indexht[1]] + x[qnet.indexht[2]])")
                qnet.node[qnet.index[i]].gammaz = x[qnet.indexht[i]]
            end
        else
            for i in 1:length(qnet.indexht)
                if(qnet.indexht[i] <= net.numHybrids - net.numBad)
                    0 <= x[qnet.indexht[i]] <= 1 || error("new gamma value should be between 0,1: $(x[qnet.indexht[i]]).")
                    qnet.edge[qnet.index[i]].hybrid || error("something odd here, optimizing gamma for tree edge $(qnet.edge[qnet.index[i]].number)")
                    setGamma!(qnet.edge[qnet.index[i]],x[qnet.indexht[i]], true)
                elseif(qnet.indexht[i] <= net.numHybrids - net.numBad + k)
                    setLength!(qnet.edge[qnet.index[i]],x[qnet.indexht[i]])
                else
                    DEBUGC && @debug "updating qnet parameters, found gammaz case when hybridization has been removed"
                    0 <= x[qnet.indexht[i]] <= 1 || error("new gammaz value should be between 0,1: $(x[qnet.indexht[i]]).")
                    #x[qnet.indexht[i]] + x[qnet.indexht[i]+1] <= 1 || @warn "new gammaz value should add to less than 1: $(x[qnet.indexht[i]])  $(x[qnet.indexht[i]+1])."
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
    for i in 1:length(net.ht)
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
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    hz = x[net.numHybrids - net.numBad + k + 1 : length(x)]
    if verbose # goes to stdout
        println("enters calculateIneqGammaz with hz $(hz), and hz[ind*2] + hz[ind*2-1] - 1 = $(hz[ind*2] + hz[ind*2-1] - 1)")
    else # goes to logger (if debug messages are turned on by user)
        @debug "enters calculateIneqGammaz with hz $(hz), and hz[ind*2] + hz[ind*2-1] - 1 = $(hz[ind*2] + hz[ind*2-1] - 1)"
    end
    hz[ind*2] + hz[ind*2-1] - 1
end

# numerical optimization of branch lengths given a network (or tree)
# and data (set of quartets with obsCF)
# using BOBYQA from NLopt package
# warning: this function assumes that the network has all the good attributes set. It will not be efficient to re-read the network inside
# to guarantee all the correct attributes, because this is done over and over inside snaq
# also, net is modified inside to set its attribute net.loglik equal to the min
"""
`optBL` road map

Function that optimizes the numerical parameters (branch lengths and inheritance probabilities) for a given network. This function is called multiple times inside `optTopLevel!`.

- Input: network `net`, data `d`
- Numerical tolerances: `ftolAbs, ftolRel, xtolAbs, xtolRel`
- Function based on `MixedModels` `fit` function
- The function assumes `net` has all the right attributes, and cannot check this inside because it would be inefficient

Procedure:

- `ht = parameters!(net)` extracts the vector of parameters to estimate `(h,t,gammaz)`, and sets as `net.ht`; identifies a bad diamond I, sets `net.numht` (vector of hybrid node numbers for h, edge numbers for t, hybrid node numbers for gammaz), and `net.index` to keep track of the vector of parameters to estimate
- `extractQuartet!(net,d)` does the following for all quartets in `d.quartet`:
   - Extract quartet by deleting all leaves not in q -> create `QuartetNetwork` object saved in `q.qnet`
   - This network is ugly and does not have edges collapsed. This is done to keep a one-to-one correspondence between the edges in `q.qnet` and the edges in `net` (if we remove nodes with only two edges, we will lose this correspondence)
   - Calculate expected CF with `calculateExpCFAll` for a copy of `q.qnet`. We do this copy because we want to keep `q.qnet` as it is (without collapsed edges into one). The function will then save the `expCF` in `q.qnet.expCF`
- `calculateExpCFAll!(qnet)` will
   - identify the type of quartet as type 1 (equivalent to a tree) or type 2 (minor CF different).
     Here the code will first clean up any hybrid node by removing nodes with only two edges before
     identifying the `qnet` (because identification depends on neighbor nodes to hybrid node);
     later, set `qnet.which` (1 or 2), `node.prev` (neighbor node to hybrid node),
     updates `node.k` (number of nodes in hybridization cycle, this can change after deleting the nodes with only two edges),
     `node.typeHyb` (1,2,3,4,5 depending on the number of nodes in the hybridization cycle
     and the origin/target of the minor hybrid edge; this attribute is never used).
   - eliminate hybridization: this will remove type 1 hybridizations first.
     If `qnet.which=1`, then the `qnet` is similar to a tree quartet,
     so it will calculate the internal length of the tree quartet: `qnet.t1`.
   - update split for `qnet.which=1`, to determine which taxa are together.
     For example, for the quartet 12|34, the split is [1,1,2,2] or [2,2,1,1],
     that is, taxon 1 and 2 are on the same side of the split. This will update `qnet.split`
   - update formula for `qnet.which=1` to know the order of minorCF and majorCF in the vector `qnet.expCF`. That is, if the quartet is 1342 (order in `qnet.quartet.taxon`), then the expected CF should match the observed CF in 13|42, 14|32, 12|34 and the `qnet` is 12|34 (given by `qnet.split`), `qnet.formula` will be [2,2,1] minor, minor, major
   - `calculateExpCF!(qnet)` for `qnet.which=1`, it will do `1-2/3exp(-qnet.t1)` if `qnet.formula[i]==1`, and `1/3exp(qnet.t1)` if `qnet.formula[i]==2`. For `qnet.which=2`, we need to make sure that there is only one hybrid node, and compute the major, minor1,minor2 expected CF in the order 12|34, 13|24, 14|23 of the taxa in `qnet.quartet.taxon`

Then we create a `NLopt` object with algorithm BOBYQA and k parameters (length of ht).
We define upper and lower bounds and define the objective function that should only depend on `x=(h,t,gz)` and g (gradient, which we do not have, but still need to put as argument).

The objective function `obj(x,g)` calls

- `calculateExpCFAll!(d,x,net)` needs to be run after `extractQuartet(net,d)` that will update `q.qnet` for all quartet.
   Assumes that `qnet.indexht` is updated already: we only need to do this at the beginning of `optBL!` because the topology is fixed at this point)
   - First it will update the edge lengths according to x
   - If the `q.qnet.changed=true` (that is, any of `qnet` branches changed value), we need to call `calculateExpCFAll!(qnet)` on a copy of `q.qnet` (again because we want to leave `q.qnet` with the edge correspondence to `net`)
- `update!(net,x)` simply saves the new x in `net.ht`

Finally, we call `NLopt.optimize`, and we update the `net.loglik` and `net.ht` at the end.
After `optBL`, we want to call `afterOptBLAll` (or `afterOptBLAllMultipleAlleles`) to check if there are `h==0,1`; `t==0`; `hz==0,1`.
"""
function optBL!(net::HybridNetwork, d::DataCF, verbose::Bool, ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64)
    (ftolRel > 0 && ftolAbs > 0 && xtolAbs > 0 && xtolRel > 0) || error("tolerances have to be positive, ftol (rel,abs), xtol (rel,abs): $([ftolRel, ftolAbs, xtolRel, xtolAbs])")
    if verbose println("OPTBL: begin branch lengths and gammas optimization, ftolAbs $(ftolAbs), ftolRel $(ftolRel), xtolAbs $(xtolAbs), xtolRel $(xtolRel)");
    else @debug        "OPTBL: begin branch lengths and gammas optimization, ftolAbs $(ftolAbs), ftolRel $(ftolRel), xtolAbs $(xtolAbs), xtolRel $(xtolRel)"; end
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
        if verbose #|| net.numBad > 0)#we want to see what happens with bad diamond I
            println("f_$count: $(round(val, digits=5)), x: $(x)")
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
    if verbose println("OPTBL: starting point $(ht)")     # to stdout
    else @debug        "OPTBL: starting point $(ht)"; end # to logger if debug turned on by user
    fmin, xmin, ret = NLopt.optimize(opt,ht)
    if verbose println("got $(round(fmin, digits=5)) at $(round.(xmin, digits=5)) after $(count) iterations (returned $(ret))")
    else @debug        "got $(round(fmin, digits=5)) at $(round.(xmin, digits=5)) after $(count) iterations (returned $(ret))"; end
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
    tmp1, tmp2 = taxadiff(d,net)
    length(tmp1)==0 || error("these taxa appear in one or more quartets, but not in the starting topology: $tmp1")
    if length(tmp2)>0
        s = "these taxa will be deleted from the starting topology, they have no quartet CF data:\n"
        for tax in tmp2 s *= " $tax"; end
        @warn s
        net = deepcopy(net)
        for tax in tmp2
            deleteleaf!(net, tax)
        end
    end
    net = readTopologyUpdate(writeTopologyLevel1(net)) # update everything for level 1
    try
        checkNet(net)
    catch err
        err.msg = "starting topology not a level 1 network:\n" * err.msg
        rethrow(err)
    end
    if(!isempty(d.repSpecies))
      expandLeaves!(d.repSpecies, net)
      net = readTopologyLevel1(writeTopologyLevel1(net)) # dirty fix to multiple alleles problem with expandLeaves
    end
    optBL!(net, d, verbose, ftolRel, ftolAbs, xtolRel,xtolAbs)
    if(net.numBad > 0) # to keep gammaz info in parenthetical description of bad diamond I
        for n in net.hybrid
            setGammaBLfromGammaz!(n,net) # get t and γ that are compatible with estimated gammaz values
        end
    end
    setNonIdBL!(net)
    if(!isempty(d.repSpecies))
      mergeLeaves!(net)
    end
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
"""
`moveHybrid` road map

Function that tries to fix a gamma zero problem (`h==0,1; t==0; hz==0,1`) after changing direction of hybrid edge failed.
This function is called in `gammaZero`.

Arguments:
- `closeN=true` will try move origin/target on all neighbors (first choose minor/major edge at random, then make list of all neighbor edges and tries to put the hybrid node in all the neighbors until successful move); if false, will delete and add hybrid until successful move up to N times (this is never tested)

Returns true if change was successful (not testing `optBL` again), and false if we could not move anything
"""
function moveHybrid!(net::HybridNetwork, edge::Edge, closeN ::Bool, origin::Bool,N::Integer, movesgamma::Vector{Int})
    edge.hybrid || error("edge $(edge.number) cannot be deleted because it is not hybrid")
    node = edge.node[edge.isChild1 ? 1 : 2];
    node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
    @debug "MOVE: moving hybrid for edge $(edge.number)"
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
"""
`gammaZero` road map

Function that tries to fix a gamma zero problem (`h==0,1; t==0; hz==0,1`)
1) First tries to do `changeDirection`
2) If not successful from start, we call `moveHybrid`
3) If successful move (change direction), we call `optBL` and check if we fixed the problem
4) If problem fixed and we do not have worse pseudolik, we return `success=true`
5) If still problem or worse pseudolik, we call `moveHybrid`

** Important: ** Any function (`afterOptBL`) calling `gammaZero` is assuming that it only made a change, so if the returned value is true, then a change was made, and the other function needs to run `optBL` and check that all parameters are 'valid'. If the returned value is false, then no change was possible and we need to remove a hybridization if the problem is h==0,1; hz==0,1. If the problem is t==0, we ignore this problem.
"""
function gammaZero!(net::HybridNetwork, d::DataCF, edge::Edge, closeN ::Bool, origin::Bool, N::Integer, movesgamma::Vector{Int})
    global CHECKNET
    currTloglik = net.loglik
    edge.hybrid || error("edge $(edge.number) should be hybrid edge because it corresponds to a gamma (or gammaz) in net.ht")
    @debug "gamma zero situation found for hybrid edge $(edge.number) with gamma $(edge.gamma)"
    node = edge.node[edge.isChild1 ? 1 : 2];
    node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
    success = changeDirectionUpdate!(net,node) #changes dir of minor
    movesgamma[4] += 1
    if(success)
        @debug begin printEverything(net); "printed everything" end
        CHECKNET && checkNet(net)
        movesgamma[10] += 1
        optBL!(net,d,false)
        flags = isValid(net)
        if(net.loglik <= currTloglik && flags[1] && flags[3])
            @debug "changing direction fixed the gamma zero situation"
            success2 = true
        else
            @debug "changing direction does not fix the gamma zero situation, need to undo change direction and move hybrid"
            success = changeDirectionUpdate!(net,node)
            success || error("strange thing, changed direction and success, but lower loglik; want to undo changeDirection, and success=false! Hybrid node is $(node.number)")
            @debug begin printEverything(net); "printed everything" end
            CHECKNET && checkNet(net)
            success2 = moveHybrid!(net,edge,closeN ,origin,N, movesgamma)
        end
    else
        @debug "changing direction was not possible to fix the gamma zero situation (success=false), need to move hybrid"
        @debug begin printEverything(net); "printed everything" end
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
# returns successchange,flagh,flagt,flaghz (flag_=false if problem with gamma, t=0 or gammaz)
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
"""
`afterOptBL` road map

Function that will check if there are `h==0,1;t==0,hz==0,1` cases in a network after calling `optBL!`.

Arguments:
- `closeN=true` will move origin/target, if false, add/delete N times before giving up (we have only tested `closeN=true`)
- `origin=true` will move origin, false will move target. We added this to avoid going back and forth between the same networks
- `movesgamma` vector of counts of number of times each move is proposed to fix a gamma zero problem: `(add,mvorigin,mvtarget,chdir,delete,nni)`

Procedure:

- First we split the `ht` vector in `nh,nt,nhz` (gammas, lengths, gammaz)
- If we find a `h==0,1`, we loop through `nh` to find a hybrid edge with h==0 or 1 and want to try to fix this by doing:
  - `gammaZero!(currT,d,edge,closeN,origin,N,movesgamma)` which returns true if there was a successful change, and we stop the loop
- If we find a `t==0`, we loop through all `nt` to find such edge, and do NNI move on this edge; return true if change successful and we stop the loop
- If we find a `hz==0,1`, we loop through `nhz` to find such hybrid edge and call `gammaZero` again

- If we did a successful change, we run `optBL` again, and recheck if there are no more problems.
- Returns successchange, flagh, flagt,flaghz (flag=true means no problems)

- If it is the multiple alleles case, it will not try to fix `h==0,1;hz==0,1` because it can reach a case that violates the multiple alleles condition. If we add a check here, things become horribly slow and inefficient, so we just delete a hybridization that has `h==0,1;hz==0,1`

** Important: ** `afterOptBL` is doing only one change, but we need to repeat multiple times to be sure that we fix all the gamma zero problems, which is why we call `afterOptBLRepeat`
"""
function afterOptBL!(currT::HybridNetwork, d::DataCF,closeN ::Bool, origin::Bool,verbose::Bool, N::Integer, movesgamma::Vector{Int})
    global CHECKNET
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
    @debug "begins afterOptBL because of conflicts: flagh,flagt,flaghz=$([flagh,flagt,flaghz])"
    successchange = true
    if(!flagh)
        for i in 1:length(nh)
            if(approxEq(nh[i],0.0) || approxEq(nh[i],1.0))
                edge = currT.edge[indh[i]]
                approxEq(edge.gamma,nh[i]) || error("edge $(edge.number) gamma $(edge.gamma) should match the gamma in net.ht $(nh[i]) and it does not")
                successchange = gammaZero!(currT, d,edge,closeN ,origin,N,movesgamma)
                !successchange || break
            end
        end
    elseif(!flagt)
        for i in 1:length(nt)
            if(approxEq(nt[i],0.0))
                edge = currT.edge[indt[i]]
                approxEq(edge.length,nt[i]) || error("edge $(edge.number) length $(edge.length) should match the length in net.ht $(nt[i]) and it does not")
                if(all((n->!n.hasHybEdge), edge.node))
                    movesgamma[6] += 1
                    successchange = NNI!(currT,edge)
                    movesgamma[12] += successchange ? 1 : 0
                    !successchange || break
                else
                    @debug "MOVE: need NNI on edge $(edge.number) because length is $(edge.length), but edge is attached to nodes that have hybrid edges"
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
    if successchange
        @debug begin
            printEverything(currT);
            "afterOptBL SUCCESSFUL change, need to run again to see if new topology is valid"
        end
        CHECKNET && checkNet(currT)
        optBL!(currT,d,verbose)
        flagh,flagt,flaghz = isValid(currT)
        @debug "new flags: flagh,flagt,flaghz $([flagh,flagt,flaghz])"
    end
    return successchange,flagh,flagt,flaghz
end

# function to repeat afterOptBL every time after changing something
# N: number of times it will delete/add hybrid if closeN =false
# origin=true, moves origin, if false, moves target. option added to control not keep coming
# to the same network over and over
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
"""
`afterOptBLRepeat` road map

`afterOptBL` is doing only one change, but we need to repeat multiple times to be sure that we fix all the gamma zero problems, which is why we call `afterOptBLRepeat`.
This function will repeat `afterOptBL` every time a successful change happened; this is done only
if `closeN=false`, because we would delete/add hybridizations and need to stop after tried N times. If `closeN=true` (default), then `afterOptBLRepeat` only does one `afterOptBL`, because in this case, only the neighbor edges need to be tested, and this would have been done already in `gammaZero`.
"""
function afterOptBLRepeat!(currT::HybridNetwork, d::DataCF, N::Integer,closeN ::Bool, origin::Bool,
                           verbose::Bool, movesgamma::Vector{Int})
    success,flagh,flagt,flaghz = afterOptBL!(currT,d,closeN ,origin,verbose,N, movesgamma)
    @debug "inside afterOptBLRepeat, after afterOptBL once, we get: success, flags: $([success,flagh,flagt,flaghz])"
    if !closeN
        @debug "closeN  is $(closeN ), and gets inside the loop for repeating afterOptBL"
        i = 1
        while(success && !reduce(&,[flagh,flagt,flaghz]) && i<N) #fixit: do we want to check loglik in this step also?
            success,flagh,flagt,flaghz = afterOptBL!(currT,d,closeN ,origin,verbose,N,movesgamma)
            i += 1
        end
        @debug (i < N ? "" : "tried afterOptBL $(i) times")
    end
    return success,flagh,flagt,flaghz
end

# function to check if we have to keep checking loglik
# returns true if we have to stop because loglik not changing much anymore, or it is close to 0.0 already
## function stopLoglik(currloglik::Float64, newloglik::Float64, ftolAbs::Float64, M::Number)
##     return (abs(currloglik-newloglik) <= M*ftolAbs) || (newloglik <= M*ftolAbs)
## end

# function similar to afterOptBLALl but for multiple alleles
# it will not try as hard to fix gamma,t=0 problem
# it will only call moveDownLevel if problem with gamma or gammaz
# fixit: later we can modify this function so that it
# will call the original afterOptBLAll if it is safe
# (i.e. there is no chance the alleles will be involved in any change)
function afterOptBLAllMultipleAlleles!(currT::HybridNetwork, d::DataCF, N::Integer,closeN ::Bool, ftolAbs::Float64, verbose::Bool, movesgamma::Vector{Int},ftolRel::Float64, xtolRel::Float64, xtolAbs::Float64)
    global CHECKNET
    !isempty(d.repSpecies) || error("calling afterOptBLAllMultipleAlleles but this is not a case with multple alleles")
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
    @debug "begins afterOptBL because of conflicts: flagh,flagt,flaghz=$([flagh,flagt,flaghz])"
    ## fixit: we could check here currT and see if we can call the original afterOptBLAll
    if(!flagh || !flaghz)
        moveDownLevel!(currT)
        optBL!(currT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
    end
    if(CHECKNET)
        checkNet(currT)
        checkTop4multAllele(currT) || error("network after moveDownLevel does not satisfy multiple alleles condition")
    end
end


# function to repeat afterOptBL every time after changing something
# closeN =true will try move origin/target, if false, will delete/add new hybrid
# default is closeN =true
# returns new approved currT (no gammas=0.0)
# N: number of times failures of accepting loglik is allowed, liktolAbs: tolerance to stop the search for lik improvement
# movesgama: vector of count of number of times each move is proposed to fix gamma zero situation:(add,mvorigin,mvtarget,chdir,delete,nni)
# movesgamma[13]: total number of accepted moves by loglik
"""
`afterOptBLAll` road map

After `optBL`, we want to call `afterOptBLAll` (or `afterOptBLAllMultipleAlleles`) to check if there are `h==0,1`; `t==0`; `hz==0,1`. This function will try to fix the gamma zero problem, but if it cannot, it will call `moveDownLevel`, to delete the hybridization from the network.

Procedure:

While `startover=true` and `tries<N`
- While `badliks < N2` (number of bad pseudolikelihoods are less than `N2`)
  - Run `success = afterOptBLRepeat`
  - If `success = true` (it changed something):
    - If worse pseudolik, then go back to original topology `currT`, set `startover=true` and `badliks++`
    - If better pseudolik, then check flags. If all good, then `startover=false`; otherwise `startover = true`
  - If `success = false` (nothing changed), then set `badliks=N2+1` (to end the while on `currT`)
    - If all flags are ok, then `startover = false`
    - If bad h or hz, then call `moveDownLevel` (delete one hybridization), and set `startover = true` (maybe deleting that hybridization did not fix other gamma zero problems)
    - If bad t, then set `startover = false`
- If left second while by back to original `currT`, and still bad h/hz, then move down one level, and `startover=true`; otherwise `startover=false`
If first while ends by `tries>N`, then it checks one last time the flags, if bad h/hz will move down one level, and exit
"""
function afterOptBLAll!(currT::HybridNetwork, d::DataCF, N::Integer,closeN ::Bool, liktolAbs::Float64, ftolAbs::Float64, verbose::Bool, movesgamma::Vector{Int},ftolRel::Float64, xtolRel::Float64, xtolAbs::Float64)
    @debug "afterOptBLAll: checking if currT has gamma (gammaz) = 0.0(1.0): currT.ht $(currT.ht)"
    currloglik = currT.loglik
    currT.blacklist = Int[];
    origin = (rand() > 0.5) #true=moveOrigin, false=moveTarget
    startover = true
    tries = 0
    N2 = N > 10 ? N/10 : 1 #num of failures of badlik around a gamma=0.0, t=0.0
    while(startover && tries < N)
        tries += 1
        @debug "inside afterOptBLALL: number of tries $(tries) out of $(N) possible"
        badliks = 0
        if(currT.loglik < liktolAbs) #curr loglik already close to 0.0
            startover = false
        else
            backCurrT0 = false
            while(badliks < N2) #will try a few options around currT
                @debug "tried $(badliks) bad likelihood options so far out of $(N2)"
                currT0 = deepcopy(currT)
                origin = !origin #to guarantee not going back to previous topology
                success,flagh,flagt,flaghz = afterOptBLRepeat!(currT,d,N,closeN ,origin,verbose,movesgamma)
                all((e->!(e.hybrid && e.inCycle == -1)), currT.edge) || error("found hybrid edge with inCycle == -1")
                @debug "inside afterOptBLAll, after afterOptBLRepeat once we get: success, flags: $([success,flagh,flagt,flaghz])"
                if !success #tried to change something but failed
                    @debug "did not change anything inside afterOptBL: could be nothing needed change or tried but couldn't anymore. flagh, flagt, flaghz = $([flagh,flagt,flaghz])"
                    if reduce(&,[flagh,flagt,flaghz]) #currT was ok
                        startover = false
                    elseif !flagh || !flaghz #currT was bad but could not change it, need to go down a level
                        !isTree(currT) || error("afterOptBL should not give reject=true for a tree")
                        @debug "current topology has numerical parameters that are not valid: gamma=0(1), gammaz=0(1); need to move down a level h-1"
                        moveDownLevel!(currT)
                        optBL!(currT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
                        startover = true
                    elseif !flagt
                        startover = false
                    end
                else #changed something
                    @debug "changed something inside afterOptBL: flagh, flagt, flaghz = $([flagh,flagt,flaghz]). oldloglik $(currloglik), newloglik $(currT.loglik)"
                    if currT.loglik > currloglik #|| abs(currT.loglik-currloglik) <= liktolAbs) #fixit: allowed like this because of changeDir that does not change much the lik but can fix h=0
                        @debug "worse likelihood, back to currT"
                        startover = true
                        backCurrT0 = true
                    else
                        @debug "better likelihood, jump to new topology and startover"
                        backCurrT0 = false
                        movesgamma[13] += 1
                        if reduce(&,[flagh,flagt,flaghz])
                            startover = false
                        else
                            currloglik = currT.loglik
                            startover = true
                        end
                    end
                end
                if backCurrT0
                    currT = currT0
                    startover = true
                    badliks += 1
                else
                    badliks = N+1 ## exit second while
                end
            end
            if backCurrT0 # leaves while for failed loglik
                @debug "tried to fix gamma zero situation for $(badliks) times and could not"
                flagh,flagt,flaghz = isValid(currT)
                if(!flagh || !flaghz)
                    !isTree(currT) || error("afterOptBL should not give reject=true for a tree")
                    @debug "current topology has numerical parameters that are not valid: gamma=0(1), t=0, gammaz=0(1); need to move down a level h-1"
                    movesgamma[5] += 1
                    movesgamma[11] += 1
                    moveDownLevel!(currT)
                    optBL!(currT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
                    startover = true
                else
                    @debug "the only problem were lengths equal to zero, so we will keep them"
                    startover = false
                end
            end
        end
    end
    if tries >= N
        @debug "afterOptBLAll ended because it tried $(tries) times with startover $(startover)"
        @debug writeTopologyLevel1(currT,true)
        flagh,flagt,flaghz = isValid(currT)
        if(!flagh || !flaghz)
            @debug "gammaz zero situation still in currT, need to move down one level to h-1"
            moveDownLevel!(currT)
            @debug begin
                printEdges(currT)
                printPartitions(currT)
                #printNodes(currT)
                writeTopologyLevel1(currT,true)
            end
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
        for i in 1:6
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
    for i in 1:length(v)
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
    for i in 1:length(v)
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
        @debug "weights before adjusting by movesfail $(v)"
        flag = adjustWeightMovesfail!(v,movesfail,Nmov,net,hmax)
        @debug "weights after adjusting by movesfail $(v)"
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
    !isTree(net) || error("net is a tree, cannot choose hybrid node")
    net.numHybrids > 1 || return net.hybrid[1]
    index1 = 0
    while(index1 == 0 || index1 > size(net.hybrid,1))
        index1 = round(Integer,rand()*size(net.hybrid,1));
    end
    @debug "chosen hybrid node for network move: $(net.hybrid[index1].number)"
    return net.hybrid[index1]
end

# function to propose a new topology given a move
# random = false uses the minor hybrid edge always
# count to know in which step we are, N for NNI trials
# order in movescount as in IF here (add,mvorigin,mvtarget,chdir,delete,nni)
# multAll = true if d.repSpecies is not empty, checked outside
"""
`proposedTop!(move,newT,random,count,N,movescount,movesfail,multall)` road map

Function to change the current network `newT` by a given `move`, and checks that the move was successful (correct attributes). If not successful, `newT` is changed back to its original state, except for the case of multiple alleles.

**Note** that the update of attributes by each move is not done in all the network, but only in the local edges that were changed by the move. This is efficient (and makes a move easy to undo), but makes the code of each move function very clunky.

Arguments:

- move chosen from `whichMove` as described in `optTopLevel`
- `newT` is the topology that will be modified inside with the move
- `random=true`: chooses minor hybrid edge with prob 1-h, and major edge with prob h, if false, always chooses minor hybrid edge
- `count`: simply which likelihood step we are in in the optimization at `optTopLevel`
- `movescount` and `movesfail`: vector of counts of number of moves proposed
- `multall=true` if multiple alleles case: we need to check if the move did not violate the multiple alleles condition (sister alleles together and no gene flow into the alleles). This is inefficient because we are proposing moves that we can reject later, instead of being smart about the moves we propose: for example, move origin/target could rule out some neighbors that move gene flow into the alleles, the same for add hybridization; nni move can check if it is trying to separate the alleles)

Moves:

- `addHybridizationUpdate(newT,N)`:
will choose a partition first (to avoid choosing edges that will create a non level-1 network)
will choose two edges from this partition randomly, will not allow two edges in a cherry (non-identifiable), or sister edges that are not identifiable
(the blacklist was a way to keep track of "bad edges" were we should not waste time trying to put hybridizations, it has never been used nor tested). Also choose gamma from U(0,0.5). The "Update" in the function name means that it creates the new hybrid, and also updates all the attributes of `newT`

- `node = chooseHybrid(newT)` choose a hybrid randomly for the next moves:
- `moveOriginUpdateRepeat!(newT,node,random)`
will choose randomly the minor/major hybrid edge to move (if `random=true`); will get the list of all neighbor edges where to move the origin, will move the origin and update all the attributes and check if the move was successful (not conflicting attributes); if not, will undo the move, and try with a different neighbor until it runs out of neighbors. Return true if the move was successful.

- `moveTargetUpdateRepeat!(newT,node,random)`
same as move origin but moving the target

- `changeDirectionUpdate!(newT,node,random)`
chooses minor/major hybrid edge at random (if `random=true), and changes the direction, and updates all the attributes. Checks if the move was successful (returns true), or undoes the change and returns false.

- `deleteHybridizationUpdate!(newT,node)`
removes the hybrid node, updates the attributes, no need to check any attributes, always successful move

- NNIRepeat!(newT,N)
choose an edge for nni that does not have a neighbor hybrid. It will try to find such an edge N times, and if it fails, it will return false (unsuccessful move). N=10 by default. If N=1, it rarely finds such an edge if the network is small or complex. The function cannot choose an external edge. it will update locally the attributes.

** Important: ** All the moves undo what they did if the move was not successful, so at the end you either have a `newT` with a new move and with all good attributes, or the same `newT` that started. This is important to avoid having to do deepcopy of the network before doing the move.
Also, after each move, when we update the attributes, we do not update the attributes of the whole network, we only update the attributes of the edges that were affected by the move. This saves time, but makes the code quite clunky.
Only the case of multiple alleles the moves does not undo what it did, because it finds out that it failed after the function is over, so just need to treat this case special.
"""
function proposedTop!(move::Integer, newT::HybridNetwork,random::Bool, count::Integer, N::Integer, movescount::Vector{Int}, movesfail::Vector{Int}, multall::Bool)
    global CHECKNET
    1 <= move <= 6 || error("invalid move $(move)") #fixit: if previous move rejected, do not redo it!
    @debug "current move: $(int2move[move])"
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
    if(multall)
        success2 = checkTop4multAllele(newT)
        @debug "entered to check topology for mult allele: $(success2)"
        success &= success2
    end
    movescount[move] += 1
    movescount[move+6] += success ? 1 : 0
    movesfail[move] += success ? 0 : 1
    @debug "success $(success), movescount (add,mvorigin,mvtarget,chdir,delete,nni) proposed: $(movescount[1:6]); successful: $(movescount[7:12]); movesfail: $(movesfail)"
    !success || return true
    @debug "new proposed topology failed in step $(count) for move $(int2move[move])"
    @debug begin printEverything(newT); "printed everything" end
    CHECKNET && checkNet(newT)
    return false
end


proposedTop!(move::Symbol, newT::HybridNetwork, random::Bool, count::Integer,N::Integer, movescount::Vector{Int},movesfail::Vector{Int}, multall::Bool) =
    proposedTop!( try move2int[move] catch; error("invalid move $(string(move))") end,
        newT, random,count,N, movescount,movesfail, multall)

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
# liktolAbs: to stop the search if loglik close to liktolAbs, or if absDiff less than liktolAbs
# hmax: max number of hybrids allowed
# closeN =true if gamma=0.0 fixed only around neighbors with move origin/target
# logfile=IOStream to capture the information on the heurisitc optimization, default stdout
"""
`optTopLevel` road map

Function that does most of the heavy-lifting of `snaq`. It optimizes the pseudolikelihood for a given starting topology, and returns the best network.
Assumes that the starting topology is level-1 network, and has all the attributes correctly updated.

Input parameters:

- Starting topology `currT`, input data `DataCF` `d`, maximum number of hybridizations `hmax`
- Numerical optimization parameters: `liktolAbs, Nfail, ftolRel, ftolAbs, xtolRel, xtolAbs`
- Print parameters: `verbose, logfile, writelog`
- Parameters to tune the search in space of networks: `closeN=true` only propose move origin/target to neighbor edges (coded, but not tested with `closeN=false`), `Nmov0` vector with maximum number of trials allowed per type of move `(add, mvorigin, mvtarget, chdir, delete, nni)`, by default computed inside with coupon’s collector formulas

The optimization procedure keeps track of
- `movescount`: count of proposed moves,
- `movesgamma`: count of proposed moves to fix a gamma zero situation (see below for definition of this situation),
- `movesfail`: count of failed moves by violation of level-1 network (`inCycle` attribute) or worse pseudolikelihood than current,
- `failures`: number of failed proposals that had a worse pseudolikelihood

Optimization procedure:

While the difference between current loglik and proposed loglik is greater than `liktolAbs`,
or `failures<Nfail`, or `stillmoves=true`:

- `Nmov` is updated based on `newT`. The type of move proposed will depend on `newT` (which is the same as `currT` at this point). For example, if `currT` is a tree, we cannot propose move origin/target.

- `move = whichMove` selects randomly a type of move, depending on `Nmov,movesfail,hmax,newT` with weights 1/5 by default for all, and 0 for delete. These weights are adjusted depending on `newT.numHybrids` and `hmax`. If `newT.numHybrids` is far from `hmax`, we give higher probability to adding a new hybrid (we want to reach the `hmax` sooner, maybe not the best strategy, easy to change).
   Later, we adjust the weights by `movesfail` (first, give weight of 0 if `movesfail[i]>Nmov[i]`, that is, if we reached the maximum possible number of moves allowed for a certain type) and then increase the probability of the other moves.
   So, unless one move has `w=0`, nothing changes. This could be improved by using the outlier quartets to guide the proposal of moves.

- `whichMove` will choose a move randomly from the weights, it will return `none` if no more moves allowed, in which case, the optimization ends

- `flag=proposedTop!(move, newT)` will modify `newT` based on `move`.
  The function `proposedTop` will return `flag=true` if the move was successful (the move succeeded by `inCycle`, `containRoot`, available edge to make the move (more details in `proposedTop`)).
  If `flag=false`, then `newT` is cleaned, except for the case of multiple alleles.
  The function `proposedTop` keeps count of `movescount` (successful move), `movesfail` (unsuccessful move),

  Options:

  `random=true`: moves major/minor hybrid edge with prob h,1-h, respectively

  `N=10`: number of trials for NNI edge.

- if(flag)
  Optimize branch lengths with `optBL`

  If `newT.loglik` is better than `currT.loglik` by `liktolAbs`, jump to `newT` (`accepted=true`) and fix `gamma=0, t=0` problems (more info on `afterOptBL`)

  If(accepted)
    `failures=0`, `movesfail=zeros`, `movescount` for successful move +1

end while

After choosing the best network `newT`, we do one last more thorough optimization of branch lengths with `optBL`,
we change non identifiable branch lengths to -1 and return `newT`
"""
function optTopLevel!(currT::HybridNetwork, liktolAbs::Float64, Nfail::Integer, d::DataCF, hmax::Integer,
                      ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
                      verbose::Bool, closeN ::Bool, Nmov0::Vector{Int}, logfile::IO, writelog::Bool)
    global CHECKNET
    @debug "OPT: begins optTopLevel with hmax $(hmax)"
    liktolAbs > 0 || error("liktolAbs must be greater than zero: $(liktolAbs)")
    Nfail > 0 || error("Nfail must be greater than zero: $(Nfail)")
    isempty(Nmov0) || all((n-> (n > 0)), Nmov0) || error("Nmov must be greater than zero: $(Nmov0)")
    if(!isempty(d.repSpecies))
        checkTop4multAllele(currT) || error("starting topology does not fit multiple alleles condition")
    end
    @debug begin printEverything(currT); "printed everything" end
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
    if(!isempty(d.repSpecies)) ## muliple alleles case
        afterOptBLAllMultipleAlleles!(currT, d, Nfail,closeN , ftolAbs, verbose,movesgamma,ftolRel,xtolRel,xtolAbs)
    else
        currT = afterOptBLAll!(currT, d, Nfail,closeN , liktolAbs, ftolAbs, verbose,movesgamma,ftolRel,xtolRel,xtolAbs) #needed return because of deepcopy inside
    end
    absDiff = liktolAbs + 1
    newT = deepcopy(currT)
    @debug begin
        printEdges(newT)
        printPartitions(newT)
        println("++++")
        writeTopologyLevel1(newT,true)
    end
    writelog && write(logfile, "\nBegins heuristic optimization of network------\n")
    while(absDiff > liktolAbs && failures < Nfail && currT.loglik > liktolAbs && stillmoves) #stops if close to zero because of new deviance form of the pseudolik
        if CHECKNET && !isempty(d.repSpecies)
            checkTop4multAllele(currT) || error("currT is not good for multiple alleles")
        end
        count += 1
        @debug "--------- loglik_$(count) = $(round(currT.loglik, digits=6)) -----------"
        if isempty(Nmov0) #if empty, not set by user
            calculateNmov!(newT,Nmov)
        end
        @debug "will propose move with movesfail $(movesfail), Nmov $(Nmov)"
        move = whichMove(newT,hmax,movesfail,Nmov)
        @debug "++++"
        if move != :none
            if !isempty(d.repSpecies) # need the original newT in case the proposed top fails by multiple alleles condition
                newT0 = deepcopy(newT)
            end
            flag = proposedTop!(move,newT,true, count,10, movescount,movesfail,!isempty(d.repSpecies)) #N=10 because with 1 it never finds an edge for nni
            if(flag) #no need else in general because newT always undone if failed, but needed for multiple alleles
                accepted = false
                all((e->!(e.hybrid && e.inCycle == -1)), newT.edge) || error("found hybrid edge with inCycle == -1")
                @debug "proposed new topology in step $(count) is ok to start optBL"
                @debug begin printEverything(newT); "printed everything" end
                CHECKNET && checkNet(newT)
                optBL!(newT,d,verbose,ftolRel, ftolAbs, xtolRel, xtolAbs)
                @debug "OPT: comparing newT.loglik $(newT.loglik), currT.loglik $(currT.loglik)"
                if(newT.loglik < currT.loglik && abs(newT.loglik-currT.loglik) > liktolAbs) #newT better loglik: need to check for error or keeps jumping back and forth
                    newloglik = newT.loglik
                    if(!isempty(d.repSpecies)) ## multiple alleles
                        afterOptBLAllMultipleAlleles!(newT, d, Nfail,closeN , ftolAbs,verbose,movesgamma,ftolRel, xtolRel,xtolAbs)
                    else
                        newT = afterOptBLAll!(newT, d, Nfail,closeN , liktolAbs, ftolAbs,verbose,movesgamma,ftolRel, xtolRel,xtolAbs) #needed return because of deepcopy inside
                    end
                    @debug "loglik before afterOptBL $(newloglik), newT.loglik now $(newT.loglik), loss in loglik by fixing gamma (gammaz)=0.0(1.0): $(newloglik>newT.loglik ? 0 : abs(newloglik-newT.loglik))"
                    accepted = true
                else
                    accepted = false
                end
                if(accepted)
                    absDiff = abs(newT.loglik - currT.loglik)
                    @debug "proposed new topology with better loglik in step $(count): oldloglik=$(round(currT.loglik, digits=3)), newloglik=$(round(newT.loglik, digits=3)), after $(failures) failures"
                    currT = deepcopy(newT)
                    failures = 0
                    movescount[move2int[move]+12] += 1
                    movesfail = zeros(Int,6) #count of failed moves for current topology
                else
                    @debug "rejected new topology with worse loglik in step $(count): currloglik=$(round(currT.loglik, digits=3)), newloglik=$(round(newT.loglik, digits=3)), with $(failures) failures"
                    failures += 1
                    movesfail[move2int[move]] += 1
                    newT = deepcopy(currT)
                end
                @debug begin
                    printEdges(newT)
                    printPartitions(newT)
                    #printNodes(newT)
                    println("++++")
                    println(writeTopologyLevel1(newT,true))
                    "ends step $(count) with absDiff $(accepted ? absDiff : 0.0) and failures $(failures)"
                end
            else
                if(!isempty(d.repSpecies))
                    newT = newT0 ## only need to go back with multiple alleles, bc other functions in proposedTop undo what they did
                end
            end
        else
            stillmoves = false
        end
        @debug "--------- loglik_$(count) end: earlier log can be discarded ----"
    end
    if ftolAbs > 1e-7 || ftolRel > 1e-7 || xtolAbs > 1e-7 || xtolRel > 1e-7
        writelog && write(logfile,"\nfound best network, now we re-optimize branch lengths and gamma more precisely")
        optBL!(newT,d,verbose, fRelBL,fAbsBL,xRelBL,xAbsBL)
    end
    assignhybridnames!(newT)
    if(absDiff <= liktolAbs)
        writelog && write(logfile,"\nSTOPPED by absolute difference criteria")
    elseif(currT.loglik <= liktolAbs)
        writelog && write(logfile,"\nSTOPPED by loglik close to zero criteria")
    elseif(!stillmoves)
        writelog && write(logfile,"\nSTOPPED for not having more moves to propose: movesfail $(movesfail), Nmov $(Nmov)")
    else
        writelog && write(logfile,"\nSTOPPED by number of failures criteria")
    end
    ## if(newT.loglik > liktolAbs) #not really close to 0.0, based on absTol also
    ##     write(logfile,"\nnewT.loglik $(newT.loglik) not really close to 0.0 based on loglik abs. tol. $(liktolAbs), you might need to redo with another starting point")
    ## end
    if(newT.numBad > 0) # if bad diamond I, need to keep gammaz info in parenthetical description
        for n in newT.hybrid
            setGammaBLfromGammaz!(n,newT) # get t and γ that are compatible with estimated gammaz values
        end
    end
    writelog && write(logfile,"\nEND optTopLevel: found minimizer topology at step $(count) (failures: $(failures)) with -loglik=$(round(newT.loglik, digits=5)) and ht_min=$(round.(newT.ht, digits=5))")
    writelog && printCounts(movescount,movesgamma,logfile)
    @debug begin
        printEdges(newT)
        printPartitions(newT)
        printNodes(newT)
        writeTopologyLevel1(newT,true)
    end
    if CHECKNET && !isempty(d.repSpecies)
        checkTop4multAllele(newT) || error("newT not suitable for multiple alleles at the very end")
    end
    return newT
end

optTopLevel!(currT::HybridNetwork, d::DataCF, hmax::Integer) = optTopLevel!(currT, likAbs, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false,true,numMoves, stdout,true)
optTopLevel!(currT::HybridNetwork, d::DataCF, hmax::Integer, verbose::Bool) = optTopLevel!(currT, likAbs, numFails, d, hmax,fRel, fAbs, xRel, xAbs, verbose,true,numMoves,stdout,true)
optTopLevel!(currT::HybridNetwork, liktolAbs::Float64, Nfail::Integer, d::DataCF, hmax::Integer,ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64, verbose::Bool, closeN ::Bool, Nmov0::Vector{Int}) = optTopLevel!(currT, liktolAbs, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN , Nmov0,stdout,true)



# function to print count of number of moves after optTopLevel
# s: IOStream to decide if you want to print to a file or to screen, by default to screen
function printCounts(movescount::Vector{Int}, movesgamma::Vector{Int},s::Union{IOStream,Base.TTY})
    length(movescount) == 18 || error("movescount should have length 18, not $(length(movescount))")
    length(movesgamma) == 13 || error("movesgamma should have length 13, not $(length(movescount))")
    print(s,"\nPERFORMANCE: total number of moves (proposed, successful, accepted) in general, and to fix gamma=0.0,t=0.0 cases\n")
    print(s,"\t--------moves general--------\t\t\t\t --------moves gamma,t------\t\n")
    print(s,"move\t Num.Proposed\t Num.Successful\t Num.Accepted\t | Num.Proposed\t Num.Successful\t Num.Accepted\n")
    names = ["add","mvorigin","mvtarget","chdir","delete","nni"]
    for i in 1:6
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
    print(s,"Proportion\t -- \t\t $(round(sum(movescount[7:12])/suma, digits=1))\t\t $(round(sum(movescount[13:18])/suma, digits=1))\t |\t -- \t\t $(round(sum(movesgamma[7:12])/suma2, digits=1))\t\t $(round(movesgamma[13]/suma2, digits=1))\n")
end

printCounts(movescount::Vector{Int}, movesgamma::Vector{Int}) = printCounts(movescount, movesgamma,stdout)

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
    global CHECKNET
    !isTree(net) ||error("cannot delete hybridization in a tree")
    @debug "MOVE: need to go down one level to h-1=$(net.numHybrids-1) hybrids because of conflicts with gamma=0,1"
    @debug begin printEverything(net); "printed everything" end
    CHECKNET && checkNet(net)
    nh = net.ht[1 : net.numHybrids - net.numBad]
    k = sum([e.istIdentifiable ? 1 : 0 for e in net.edge])
    nt = net.ht[net.numHybrids - net.numBad + 1 : net.numHybrids - net.numBad + k]
    nhz = net.ht[net.numHybrids - net.numBad + k + 1 : length(net.ht)]
    indh = net.index[1 : net.numHybrids - net.numBad]
    indhz = net.index[net.numHybrids - net.numBad + k + 1 : length(net.ht)]
    flagh,flagt,flaghz = isValid(nh,nt,nhz)
    if(!flagh)
        for i in 1:length(nh)
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
    @debug begin printEverything(net); "printed everything" end
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


# optTopRuns! and optTopRun1! do *not* modify currT0
"""
Road map for various functions behind [`snaq!`](@ref)

    snaq!
    optTopRuns!
    optTopRun1!
    optTopLevel!
    optBL!

All return their optimized network.

- [`snaq!`](@ref) calls `optTopRuns!` once, after a deep copy of the starting network.
  If the data contain multiple alleles from a given species, `snaq!` first
  expands the leaf for that species into 2 separate leaves, and merges them
  back into a single leaf after calling `optTopRuns!`.
- `optTopRuns!` calls [`optTopRun1!`](@ref) several (`nrun`) times.
  assumes level-1 network with >0 branch lengths.
  assumes same tips in network as in data: i.e. 2 separate tips per species
                                           that has multiple alleles.
  each call to `optTopRun1!` gets the same starting network.
- `optTopRun1!` calls `optTopLevel!` once, after deep copying + changing the starting network slightly.
- `optTopLevel!` calls `optBL!` various times and proposes new network with various moves.

"""
function optTopRuns!(currT0::HybridNetwork, liktolAbs::Float64, Nfail::Integer, d::DataCF, hmax::Integer,
                     ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
                     verbose::Bool, closeN ::Bool, Nmov0::Vector{Int}, runs::Integer,
                     outgroup::AbstractString, rootname::AbstractString, seed::Integer, probST::Float64)
    writelog = true
    writelog_1proc = false
    if (rootname != "")
        julialog = string(rootname,".log")
        logfile = open(julialog,"w")
        juliaout = string(rootname,".out")
        if Distributed.nprocs() == 1
            writelog_1proc = true
            juliaerr = string(rootname,".err")
            errfile = open(juliaerr,"w")
        end
    else
      writelog = false
      logfile = stdout # used in call to optTopRun1!
    end
    str = """optimization of topology, BL and inheritance probabilities using:
              hmax = $(hmax),
              tolerance parameters: ftolRel=$(ftolRel), ftolAbs=$(ftolAbs),
                                    xtolAbs=$(xtolAbs), xtolRel=$(xtolRel).
              max number of failed proposals = $(Nfail), liktolAbs = $(liktolAbs).
             """
    if outgroup != "none"
        str *= "Outgroup: $(outgroup) (for rooting at the final step)\n"
    end
    str *= (writelog ? "rootname for files: $(rootname)\n" : "no output files\n")
    str *= "BEGIN: $(runs) runs on starting tree $(writeTopologyLevel1(currT0,true))\n"
    if Distributed.nprocs()>1
        str *= "       using $(Distributed.nprocs()) processors\n"
    end
    if (writelog)
      write(logfile,str)
      flush(logfile)
    end
    print(stdout,str)
    print(stdout, Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s") * "\n")
    # if 1 proc: time printed to logfile at start of every run, not here.

    if(seed == 0)
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) #better seed based on clock
    end
    if (writelog)
      write(logfile,"\nmain seed $(seed)\n")
      flush(logfile)
    else print(stdout,"\nmain seed $(seed)\n"); end
    Random.seed!(seed)
    seeds = [seed;round.(Integer,floor.(rand(runs-1)*100000))]
    if writelog && !writelog_1proc
        for i in 1:runs # workers won't write to logfile
            write(logfile, "seed: $(seeds[i]) for run $(i)\n")
        end
        flush(logfile)
    end

    tstart = time_ns()
    bestnet = Distributed.pmap(1:runs) do i # for i in 1:runs
        logstr = "seed: $(seeds[i]) for run $(i), $(Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s"))\n"
        print(stdout, logstr)
        msg = "\nBEGIN SNaQ for run $(i), seed $(seeds[i]) and hmax $(hmax)"
        if writelog_1proc # workers can't write on streams opened by master
            write(logfile, logstr * msg)
            flush(logfile)
        end
        verbose && print(stdout, msg)
        GC.gc();
        try
            best = optTopRun1!(currT0, liktolAbs, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs,
                       verbose, closeN , Nmov0,seeds[i],logfile,writelog_1proc,probST);
            logstr *= "\nFINISHED SNaQ for run $(i), -loglik of best $(best.loglik)\n"
            verbose && print(stdout, logstr)
            if writelog_1proc
              logstr = writeTopologyLevel1(best,outgroup=outgroup, printID=true, multall=!isempty(d.repSpecies))
              logstr *= "\n---------------------\n"
              write(logfile, logstr)
              flush(logfile)
            end
            return best
        catch(err)
            msg = "\nERROR found on SNaQ for run $(i) seed $(seeds[i]): $(err)\n"
            logstr = msg * "\n---------------------\n"
            if writelog_1proc
                write(logfile, logstr)
                flush(logfile)
                write(errfile, msg)
                flush(errfile)
            end
            @warn msg # returns: nothing
        end
    end
    tend = time_ns() # in nanoseconds
    telapsed = round(convert(Int64, tend-tstart) * 1e-9, digits=2) # in seconds
    writelog_1proc && close(errfile)
    msg = "\n" * Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s")
    if writelog
        write(logfile, msg)
    elseif verbose
        print(stdout, msg)
    end
    filter!(n -> n !== nothing, bestnet) # remove "nothing", failed runs
    if length(bestnet)>0
        ind = sortperm([n.loglik for n in bestnet])
        bestnet = bestnet[ind]
        maxNet = bestnet[1]::HybridNetwork # tell type to compiler
    else
        error("all runs failed")
    end

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
        write(s,"$(writeTopologyLevel1(maxNet,printID=true, multall=!isempty(d.repSpecies))), with -loglik $(maxNet.loglik) (best network found, remaining sorted by log-pseudolik; the smaller, the better)\n")
        # best network is included first: for score comparison with other networks
        foundBad = false
        for n in otherNet
            try
                optBL!(n,d) ##optBL MUST have network with all the attributes, and undirectedOtherNetworks will return "good" networks that way
                if(n.numBad > 0) # to keep gammaz info in parenthetical description of bad diamond I
                    for nod in n.hybrid
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
        for n in otherNet
            write(s,"$(writeTopologyLevel1(n,printID=true, multall=!isempty(d.repSpecies))), with -loglik $(n.loglik)\n")
        end
        foundBad && write(s,"Problem found when optimizing branch lengths for some networks, left loglik as -1. Please report this issue on github. Thank you!")
        close(s)
    end

    setNonIdBL!(maxNet)
    writelog &&
    write(logfile,"\nMaxNet is $(writeTopologyLevel1(maxNet,printID=true, multall=!isempty(d.repSpecies))) \nwith -loglik $(maxNet.loglik)\n")
    print(stdout,"\nMaxNet is $(writeTopologyLevel1(maxNet,printID=true, multall=!isempty(d.repSpecies))) \nwith -loglik $(maxNet.loglik)\n")

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
    s = writelog ? open(juliaout,"w") : stdout
    str = writeTopologyLevel1(maxNet, printID=true,multall=!isempty(d.repSpecies)) * """
     -Ploglik = $(maxNet.loglik)
     Dendroscope: $(writeTopologyLevel1(maxNet,di=true, multall=!isempty(d.repSpecies)))
     Elapsed time: $(telapsed) seconds, $(runs) attempted runs
    -------
    List of estimated networks for all runs (sorted by log-pseudolik; the smaller, the better):
    """
    for n in bestnet
      str *= " "
      str *= (outgroup == "none" ? writeTopologyLevel1(n,printID=true, multall=!isempty(d.repSpecies)) :
                                   writeTopologyLevel1(n,outgroup=outgroup, printID=true, multall=!isempty(d.repSpecies)))
      str *= ", with -loglik $(n.loglik)\n"
    end
    str *= "-------\n"
    write(s,str);
    writelog && close(s) # to close juliaout file (but not stdout!)
    writelog && close(logfile)

    return maxNet
end

optTopRuns!(currT::HybridNetwork, d::DataCF, hmax::Integer, runs::Integer, outgroup::AbstractString, rootname::AbstractString) = optTopRuns!(currT, likAbs, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false, true, numMoves, runs, outgroup,rootname,0,0.3)
optTopRuns!(currT::HybridNetwork, d::DataCF, hmax::Integer, runs::Integer, outgroup::AbstractString) = optTopRuns!(currT, likAbs, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false, true, numMoves, runs, outgroup,"optTopRuns",0,0.3)
optTopRuns!(currT::HybridNetwork, d::DataCF, hmax::Integer, runs::Integer) = optTopRuns!(currT, likAbs, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false, true, numMoves, runs, "none", "optTopRuns",0,0.3)
optTopRuns!(currT::HybridNetwork, d::DataCF, hmax::Integer) = optTopRuns!(currT, likAbs, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false, true, numMoves, 10, "none", "optTopRuns",0,0.3)

# picks a modification of starting topology and calls optTopLevel
# the seed is used as is
# does *not* modify currT0. Modifies data d only.
"""
    optTopRun1!(net, liktolAbs, Nfail, d::DataCF, hmax, etc.)

The function will run 1 run by modifying the starting topology and
calling `optTopLevel`. See [`optTopRuns!`](@ref) for a roadmap.

`probST` (default in snaq is 0.3) is the probability of starting one run
at the same input tree. So, with probability `1-probST`, we will change the
topology by a NNI move on a tree edge without neighbor hybrid.
If the starting topology is a network, then with probability `1-probST`
it will also modify one randomly chosen hybrid edge: with prob 0.5,
the function will move origin, with prob 0.5 will do move target.

If there are multiple alleles (`d.repSpecies` not empty),
then the function has to check that the starting topology
does not violate the multiple alleles condition.

After modifying the starting topology with NNI and/or move origin/target,
`optTopLevel` is called.
"""
function optTopRun1!(currT0::HybridNetwork, liktolAbs, Nfail::Integer, d::DataCF, hmax::Integer,
                     ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
                     verbose::Bool, closeN ::Bool, Nmov0::Vector{Int},seed::Integer,
                     logfile::IO, writelog::Bool, probST::Float64)
    Random.seed!(seed)
    currT = deepcopy(currT0);
    if(probST<1.0 && rand() < 1-probST) # modify starting tree by a nni move
        suc = NNIRepeat!(currT,10); #will try 10 attempts to do an nni move, if set to 1, hard to find it depending on currT
        if(!isempty(d.repSpecies))
            suc2 = checkTop4multAllele(currT)
            suc &= suc2
            if(!suc2)
                currT = deepcopy(currT0)
            end
        end
        writelog && suc && write(logfile," changed starting topology by NNI move\n")
        if(!isTree(currT))
            if(rand() < 1-probST) # modify starting network by mvorigin, mvtarget with equal prob
                currT0 = deepcopy(currT) # to go back if new topology does not work for mult alleles
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
                    if(!isempty(d.repSpecies))
                        suc2 = checkTop4multAllele(currT)
                        suc &= suc2
                        if(!suc2)
                            currT = deepcopy(currT0)
                        end
                    end
                    writelog && suc && write(logfile,"\n changed starting network by move origin")
                else
                    suc = moveTargetUpdateRepeat!(currT,currT.hybrid[ind],true)
                    if(!isempty(d.repSpecies))
                        suc2 = checkTop4multAllele(currT)
                        suc &= suc2
                        if(!suc2)
                            currT = deepcopy(currT0)
                        end
                    end
                    writelog && suc && write(logfile,"\n changed starting network by move target")
                end
            end
        end
    end
    GC.gc();
    optTopLevel!(currT, liktolAbs, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN , Nmov0,logfile,writelog)
end

optTopRun1!(currT::HybridNetwork, d::DataCF, hmax::Integer) = optTopRun1!(currT, likAbs, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false, true, numMoves, 0,stdout,true,0.3)
optTopRun1!(currT::HybridNetwork, d::DataCF, hmax::Integer, seed::Integer) = optTopRun1!(currT, likAbs, numFails, d, hmax,fRel, fAbs, xRel, xAbs, false, true, numMoves,seed,stdout,true,0.3)


# function SNaQ: it calls directly optTopRuns but has a prettier name
# it differs from optTopRuns in that it creates a deepcopy of the starting topology,
#    check that it's of level 1 and updates its BL.
# only currT and d are necessary, all others are optional and have default values
"""
    snaq!(T::HybridNetwork, d::DataCF)

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

- `hmax` (default 1): maximum number of hybridizations allowed
- `verbose` (default false): if true, print information about the numerical optimization
- `runs` (default 10): number of independent starting points for the search
- `outgroup` (default none): outgroup taxon to root the estimated topology at the very end
- `filename` (default "snaq"): root name for the output files (`.out`, `.err`). If empty (""),
  files are *not* created, progress log goes to the screen only (standard out).
- `seed` (default 0 to get it from the clock): seed to replicate a given search
- `probST` (default 0.3): probability to start from `T` at each given run.
  With problability 1-probST, the search is started from an NNI modification of `T`
  along a tree edge with no hybrid neighbor,
  with a possible modification of one reticulation if `T` has one.
- `updateBL` (default true): If true and if `T` is a tree, the branch lengths in `T`
  are first optimized roughly with [`updateBL!`](@ref) by using the average CF of
  all quartets defining each branch and back-calculating the coalescent units.

The following optional arguments control when to stop the optimization of branch lengths
and γ's on each individual candidate network. Defaults are in parentheses:

- `ftolRel` (1e-6) and `ftolAbs` (1e-6): relative and absolute differences of the network score
  between the current and proposed parameters,
- `xtolRel` (1e-2) and `xtolAbs` (1e-3): relative and absolute differences between the current
  and proposed parameters.

Greater values will result in a less thorough but faster search. These parameters are used
when evaluating candidate networks only.
The following optional arguments control when to stop proposing new network topologies:

- `Nfail` (75): maximum number of times that new topologies are proposed and rejected (in a row).
- `liktolAbs` (1e-6): the proposed network is accepted if its score is better than the current score by
  at least liktolAbs.

Lower values of `Nfail` and greater values of `liktolAbs` and `ftolAbs` would result in a less thorough but faster search.

At the end, branch lengths and γ's are optimized on the last "best" network
with different and very thorough tolerance parameters:
1e-12 for ftolRel, 1e-10 for ftolAbs, xtolRel, xtolAbs.

See also: [`topologyMaxQPseudolik!`](@ref) to optimize parameters on a fixed topology,
and [`topologyQPseudolik!`](@ref) to get the deviance (pseudo log-likelihood up to a constant)
of a fixed topology with fixed parameters.

Reference:  
Claudia Solís-Lemus and Cécile Ané (2016).
Inferring phylogenetic networks with maximum pseudolikelihood under incomplete lineage sorting.
[PLoS Genetics](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896)
12(3):e1005896
"""
function snaq!(currT0::HybridNetwork, d::DataCF;
      hmax=1::Integer, liktolAbs=likAbs::Float64, Nfail=numFails::Integer,
      ftolRel=fRel::Float64, ftolAbs=fAbs::Float64, xtolRel=xRel::Float64, xtolAbs=xAbs::Float64,
      verbose=false::Bool, closeN=true::Bool, Nmov0=numMoves::Vector{Int},
      runs=10::Integer, outgroup="none"::AbstractString, filename="snaq"::AbstractString,
      seed=0::Integer, probST=0.3::Float64, updateBL=true::Bool)
    0.0<=probST<=1.0 || error("probability to keep the same starting topology should be between 0 and 1: $(probST)")
    currT0.numTaxa >= 5 || error("cannot estimate hybridizations in topologies with fewer than 5 taxa, this topology has $(currT0.numTaxa) taxa")
    typemax(Int) > length(d.quartet) ||
    @warn "the number of rows / 4-taxon sets exceeds the max integer of type $Int ($(typemax(Int))). High risk of overflow errors..."
    # need a clean starting net. fixit: maybe we need to be more thorough here
    # yes, need to check that everything is ok because it could have been cleaned and then modified
    tmp1, tmp2 = taxadiff(d,currT0)
    length(tmp1)==0 || error("these taxa appear in one or more quartets, but not in the starting topology: $tmp1")
    if length(tmp2)>0
        s = "these taxa will be deleted from the starting topology, they have no quartet CF data:\n"
        for tax in tmp2 s *= " $tax"; end
        @warn s
        currT0 = deepcopy(currT0)
        for tax in tmp2
            deleteleaf!(currT0, tax)
        end
    end
    startnet = readTopologyUpdate(writeTopologyLevel1(currT0)) # update all level-1 things
    flag = checkNet(startnet,true) # light checking only
    flag && error("starting topology suspected not level-1")
    try
        checkNet(startnet)
    catch err
        err.msg = "starting topology not a level 1 network:\n" * err.msg
        rethrow(err)
    end
    if updateBL && isTree(startnet)
        updateBL!(startnet,d)
    end
    # for the case of multiple alleles: expand into two leaves quartets like sp1 sp1 sp2 sp3.
    if !isempty(d.repSpecies)
        expandLeaves!(d.repSpecies,startnet)
        startnet = readTopologyLevel1(writeTopologyLevel1(startnet)) # dirty fix to multiple alleles problem with expandLeaves
    end
    net = optTopRuns!(startnet, liktolAbs, Nfail, d, hmax, ftolRel,ftolAbs, xtolRel,xtolAbs,
                      verbose, closeN, Nmov0, runs, outgroup, filename,seed,probST)
    if(!isempty(d.repSpecies))
        mergeLeaves!(net)
    end
    return net
end


## function to modify the starting topology currT0 by an NNI move with probability 1-probST
## if starting topology is a network: also do move Origin/Target with prob 1-probST, each 50-50 chance
## if outgroup!="none" means that the root placement matters at the end (on outgroup), used for max parsimony
function findStartingTopology!(currT0::HybridNetwork, probST::Float64, multAll::Bool, writelog::Bool, logfile::IO;
                               outgroup="none"::Union{AbstractString,Integer})
    currT = deepcopy(currT0);
    if probST<1.0 && rand() < 1-probST # modify starting tree by a nni move
        suc = NNIRepeat!(currT,10); #will try 10 attempts to do an nni move, if set to 1, hard to find it depending on currT
        if multAll
            suc2 = checkTop4multAllele(currT)
            suc &= suc2
            if !suc2
                currT = deepcopy(currT0)
            end
        end
        writelog && suc && write(logfile," changed starting topology by NNI move\n")
        if !isTree(currT)
            if rand() < 1-probST # modify starting network by mvorigin, mvtarget with equal prob
                currT0 = deepcopy(currT) # to go back if new topology does not work for mult alleles
                ind = rand(1:currT.numHybrids)
                mymove = ( rand()<0.5 ? "origin" : "target" )
                mymove_fun! = (mymove=="origin" ? moveOriginUpdateRepeat! : moveTargetUpdateRepeat!)
                suc = mymove_fun!(currT,currT.hybrid[ind],true)
                if multAll
                    suc2 = checkTop4multAllele(currT)
                    suc &= suc2
                    if !suc2
                        currT = deepcopy(currT0)
                    end
                end
                writelog && suc && write(logfile,"\n changed starting network by move $(mymove)")
            end
        end
    end
    GC.gc();
    if outgroup != "none"
        currT1 = deepcopy(currT) ##just to check, but we do not want a rooted network at the end
        try
            rootatnode!(currT1,outgroup)
        catch err
            if isa(err, RootMismatch)
                 println("RootMismatch: ", err.msg,
                 """\nThe starting topology has hybrid edges that are incompatible with the desired outgroup.
                    Reverting to original starting topology.
                    """)
            else
                println("error trying to reroot: ", err.msg);
            end
            currT = deepcopy(currT0)
        end
    end
    return currT
end
