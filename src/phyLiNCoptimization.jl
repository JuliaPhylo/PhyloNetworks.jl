const moveweights = Distributions.aweights([0.5, 0.3, 0.2])
const likAbsAddHybLiNC = 0.5
const likAbsDelHybLiNC = -0.1
const amin = 0.05
const amax = 500.0

# optimization wrapper functions #
"""
    phyLiNC!(net::HybridNetwork, fastafile::String, modSymbol::Symbol,
        maxhybrid=1::Int64, no3cycle=true::Bool, unzip=true::Bool,
        nohybridladder=true::Bool, maxmoves=100::Int64, nreject=75::Int64,
        verbose=false::Bool, NLoptMethod=:LD_MMA::Symbol,
        ftolRel=fRelBL::Float64, ftolAbs=fAbsBL::Float64,
        xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64,
        constraints=TopologyConstraint[]::Vector{TopologyConstraint},
        alphamin=amin::Float64, alphamax=amax::Float64)

Estimate a phylogenetic network (or tree) from concatenated fasta data using
maximum likelihood. Any level network is accepted. The search starts from the
given `net` topology, using local search to optimize structure (nearest-neighbor
interchange moves, add hybridizations, and remove hybridizations) and non-linear
optimization to optimize evolutionary rates, rate variation across sites, branch
lengths and gammas.

Return a StatisticalSubstitutionModel object.

There are many optional arguments, including

`maxhybrid` (default = 1): maximum number of hybridizations allowed
`no3cycle` (default = true): prevents 3-cycles to preserve the identifiability of
the network structure
`unzip` (default = true): sets edges below hybrid nodes to zero to preserve the
identifiability of branch lengths.
`nohybridladder` (default = true): prevents hybrid ladder in network. If this is
true, input network must not have hybrid ladders.
`maxmoves` (default = 100): maximum number of moves before branch lengths,
hybrid gamma values, evolutionary rates, and rate variation parameters are
reestimated.
`verbose` (default = false): if true, print information about the numerical optimization
`constraints` (default none): topology constraints to follow during structure
optimization. Created using [`TopologyConstraint`] (@ref)
#? Should we create constraints in this function or ask users to do so beforehand?

The following optional arguments control when to stop the optimization of branch
lengths and gamma values on each individual candidate network. Defaults in
parentheses.
ftolRel (1e-6) and ftolAbs (1e-6): relative and absolute differences of the
network score between the current and proposed parameters
xtolRel (1e-2) and xtolAbs (1e-3): relative and absolute differences between the
current and proposed parameters.
Greater values will result in a less thorough but faster search. These parameters
are used when evaluating candidate networks only.

The following optional arguments control when to stop proposing new
network topologies:

`nreject` (default = 75): the maximum number of times that new topologies are
proposed and rejected in a row. As in snaq, optimize structure runs until
rejections = nreject.
liktolAbs (1e-6): the proposed network is accepted if its score is better than
the current score by at least liktolAbs.
Lower values of `nreject` and greater values of `liktolAbs` and `ftolAbs` would
result in a less thorough but faster search.

#? Cécile, should we allow users to set a seed to replicate a search?
seed (default 0 to get it from the clock): seed to replicate a given search

#TODO in future, run independently with multiple starting topologies?
"""
function phyLiNC!(net::HybridNetwork, fastafile::String, modSymbol::Symbol,
    maxhybrid=1::Int64, no3cycle=true::Bool, unzip=true::Bool,
    nohybridladder=true::Bool, maxmoves=100::Int64, nreject=75::Int64,
    verbose=false::Bool, NLoptMethod=:LD_MMA::Symbol, ftolRel=fRelBL::Float64,
    ftolAbs=fAbsBL::Float64, xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64,
    constraints=TopologyConstraint[]::Vector{TopologyConstraint},
    alphamin=amin::Float64, alphamax=amax::Float64)

    obj = StatisticalSubstitutionModel(net, fastafile, modSymbol, maxhybrid)

    checknetworkbeforeLiNC!(obj.net, maxhybrid, no3cycle, unzip, nohybridladder,
    constraints)

    startingBL!(obj.net, unzip, obj.trait, obj.siteweight)

    discrete_corelikelihood!(obj)

    fit!(obj; optimizeQ=true, optimizeRVAS=true, maxeval=20) #rough optim of rates and alpha

    done = false
    while !done # break out of this loop only is if rejections = nreject.
        done = optimizestructure!(obj, maxmoves, maxhybrid, no3cycle, unzip,
            nohybridladder, nreject, verbose, constraints)
        fit!(obj; optimizeQ=true, optimizeRVAS=true, ftolRel=1e-2, ftolAbs=1e-2,
             xtolRel=1e-1, xtolAbs=1e-2)
        optimizeBL!(obj, obj.net, obj.net.edge, unzip, verbose, 1000, NLoptMethod,
                    ftolRel, ftolAbs, xtolRel, xtolAbs)
        optimizeallgammas!(obj, obj.net, unzip, verbose, 1000, NLoptMethod,
                           ftolRel, ftolAbs, xtolRel, xtolAbs)
    end
    return obj
end

"""
    checknetworkbeforeLiNC!(net::HybridNetwork, maxhybrid::Int64, no3cycle::Bool,
        unzip::Bool, nohybridladder::Bool,
        constraints=TopologyConstraint[]::Vector{TopologyConstraint})

Remove nodes of degree 2, including root. According to user-given options, check
for 3-cycles, hybrid ladders, max number of hybrids, and `unzip` (set all edges
below hybrids to length zero).

```jldoctest
julia> maxhybrid = 3;

julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple.aln"));

julia> obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69, maxhybrid);

julia> PhyloNetworks.checknetworkbeforeLiNC!(obj.net, maxhybrid, true, true, true)
HybridNetwork, Rooted Network
8 edges
8 nodes: 4 tips, 1 hybrid nodes, 3 internal tree nodes.
tip labels: A, B, C, D
((A:2.0,(B:0.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0,D:2.5);


```
"""
function checknetworkbeforeLiNC!(net::HybridNetwork, maxhybrid::Int64, no3cycle::Bool,
    nohybridladder::Bool, unzip::Bool,
    constraints=TopologyConstraint[]::Vector{TopologyConstraint})

    checkspeciesnetwork!(net, constraints) || # checks for polytomies, constraint violations, nodes of degree 2
        error("individuals specified in species constraint are not grouped together in one polytomy.")
    if no3cycle
        !contains3cycles(net) || error("Options indicate there should be no 3-cycles
        in the returned network, but the input network contains
        one or more 3-cycles.")
    end
    if nohybridladder
        !hashybridladder(net) || error("Options indicate there should be no
        hybrid ladders in the returned network, but the input network contains
        one or more hybrid ladders.")
    end
    if length(net.hybrid) > maxhybrid
        error("Options indicate a maximum number of hybrids of $(maxhybrid), but
        the input network contains $(length(net.hybrid)) hybrids. Please
        increase maxhybrid to $(length(net.hybrid)) or provide input network
        with $(maxhybrid) (or fewer) hybrids.")
    end
    if unzip && length(net.hybrid) > 0
        constrainededges = [getChildEdge(h) for h in net.hybrid]
        setlengths!(constrainededges, zeros(length(constrainededges)))
    end
    return net
end

"""
    optimizestructure!(obj::SSM, maxmoves::Int64, maxhybrid::Int64,
        no3cycle::Bool, unzip::Bool, nohybridladder::Bool, nreject=75::Int64,
        verbose=false::Bool,
        constraints=TopologyConstraint[]::Vector{TopologyConstraint})

Alternate nni moves, hybrid moves, and root changes. Optimizes local branch
lengths and hybrid gammas after each move, then decides whether or not to accept
the move by comparing likelihoods. After adding or removing a hybrid, updates
`SSM` object's displayed `trees` and their attributes.
Return `done` boolean indicating if number of rejections was met.

The percent of nni moves, hybrid moves, and root changes to be performed is
0.5, 0.3, and 0.2 respectively.

For a description of optional arguments, see [`phyLiNC`](@ref).

Assumptions:
- `checknetworkbeforeLiNC` and `discrete_corelikelihood!` have been called on
obj.net and obj.
- starting with a network without 2- and 3- cycles
(checked by `checknetworkbeforeLiNC`)

Note: When removing a hybrid edge, always removes the minor edge.
"""
function optimizestructure!(obj::SSM, maxmoves::Int64, maxhybrid::Int64,
    no3cycle::Bool, unzip::Bool, nohybridladder::Bool, nreject=75::Int64,
    verbose=false::Bool,
    constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    nmoves = 0
    rejections = 0
    while nmoves < maxmoves && rejections < nreject # both should be true to continue
        currLik = obj.loglik
        movechoice = sample(["nni", "hybrid", "root"], moveweights)
        if movechoice == "nni"
            edgefound = false
            blacklist = Edge[]
            while !edgefound # randomly select interior edge
                if length(blacklist) == length(obj.net.edge)
                    verbose &&
                    println("There are no nni moves possible in this network.")
                end
                eindex = Random.randperm(length(obj.net.edge))[1]
                e1 = obj.net.edge[eindex]
                if !(e1 in blacklist) # else go back to top of nni while loop
                    undoinfo = nni!(obj.net, e1,nohybridladder,no3cycle,constraints)
                    if !isnothing(undoinfo)
                        nmoves += 1
                        edgefound = true
                        #optimize BL and gamma locally
                        discrete_corelikelihood!(obj) # update loglik
                        if obj.loglik - currLik < likAbs
                            nni!(undoinfo...) # undo move
                            rejections += 1
                        else
                            rejections = 0 # reset
                        end
                    else # if move unsuccessful, search for edge until successful
                        push!(blacklist, e1)
                    end
                end
            end
        elseif movechoice == "hybrid" # perform hybrid move
            if maxhybrid == 0
                @debug("The maximum number of hybrids allowed is $maxhybrid,
                so hybrid moves are not legal on this network.")
            elseif length(obj.net.hybrid) == 0
                add = true # add hybrid
            elseif length(obj.net.hybrid) == maxhybrid
                add = false # remove hybrid
            elseif length(obj.net.hybrid) > maxhybrid
                error("""The network has more hybrids than allowed. maxhybrid =
                 $maxhybrid, but network has $(obj.net.hybrid) hybrids.""")
            else  # meaning we can do either move
                add = (rand() > 0.3) # add hybrid with 70% probability
            end
            if add
                added = addhybridedge_LiNC!(obj, currLik, maxhybrid, no3cycle,
                        unzip, nohybridladder, verbose, constraints)
                nmoves +=1
                if isnothing(added)
                    verbose && println("Cannot add a hybrid to the network.")
                    movechoice = "add hybrid (unsuccessful attempt)"
                elseif added
                    movechoice = "add hybrid"
                    rejections = 0 # reset
                else
                    movechoice = "add hybrid (but deleted afterward)"
                    rejections += 1
                end
            else # delete hybrid
                deleted = deletehybridedge_LiNC!(obj, currLik, maxhybrid,
                        no3cycle, unzip, nohybridladder, verbose, constraints)
                nmoves += 1
                if isnothing(deleted)
                    verbose && println("""Cannot delete a hybrid to the network
                     without violating a topology constraint.""")
                    movechoice = "delete hybrid (unsuccessful attempt)"
                elseif deleted
                    movechoice = "delete hybrid"
                    rejections = 0 # reset
                else
                    movechoice = "delete hybrid (but added back)"
                    rejections += 1
                end
            end
        else # change root (doesn't affect likelihood)
            originalroot = obj.net.root
            changednet = moveroot!(obj.net, constraints)
            nmoves += 1
            if !changednet
                @debug("Cannot perform a root change move on current network.")
            end
        end
        verbose && println("""loglik = $(loglikelihood(obj)) after move of type
        $movechoice, $nmoves total moves, and $rejections rejected moves""")
    end
    return rejections >= nreject # done if rejections >= nreject
end

"""
    addhybridedge_LiNC!(obj::SSM, currLik::Float64, maxhybrid::Int64,
        no3cycle::Bool, unzip::Bool, nohybridladder::Bool, verbose::Bool,
        constraints::Vector{TopologyConstraint})

Completes checks, adds hybrid in a random location, updates SSM object, and
optimizes branch lengths and gammas locally as part of PhyLiNC optimization.

Return true if accepted add hybrid move. If move not accepted, return false.
If cannot add a hybrid, return nothing.

Assumptions:
- called by [`optimizestructure!`](@ref)
"""
function addhybridedge_LiNC!(obj::SSM, currLik::Float64, maxhybrid::Int64,
    no3cycle::Bool, unzip::Bool, nohybridladder::Bool, verbose::Bool,
    constraints::Vector{TopologyConstraint})
    origobj = deepcopy(obj) # hold old SSM in case we remove new hybrid
        # TODO can we use the same memory space for this every time?
    result = addhybridedge!(obj.net, nohybridladder, no3cycle, constraints)
    for h in obj.net.hybrid # TODO in future, remove these debugging statements
        @debug !any([n.leaf for n in getMajorParentEdge(h).node]) ||
            println("""edge $(getMajorParentEdge(h).number) is a major
                    hybrid edge and the stem edge of a leaf""")
        @debug !any([n.leaf for n in getMinorParentEdge(h).node]) ||
            println("""edge $(getMinorParentEdge(h).number) is a minor hybrid
                    edge and the stem edge of a leaf""")
        end
    if !isnothing(result)
        newhybridnode, newhybridedge = result
        updateSSM!(obj)
        optimizelocalBL!(obj, obj.net, newhybridedge, unzip, verbose) # updates obj.loglik
        optimizelocalgammas!(obj, obj.net, newhybridedge, unzip, verbose)
        if obj.loglik - currLik < likAbsAddHybLiNC # improvement too small or negative: undo
            obj = origobj
            return false
        else
            return true
        end
    else
        return nothing
    end
end

"""
    deletehybridedge_LiNC!(obj::SSM, currLik::Float64, maxhybrid::Int64,
        no3cycle::Bool, unzip::Bool, nohybridladder::Bool, verbose::Bool,
        constraints::Vector{TopologyConstraint})

Deletes a random hybrid edge, completes checks, and updates SSM object as part of
PhyLiNC optimization.

Return true if accepted delete hybrid move. If move not accepted, return false.

Assumptions:
- called by [`optimizestructure!`](@ref)
"""
function deletehybridedge_LiNC!(obj::SSM, currLik::Float64, maxhybrid::Int64,
    no3cycle::Bool, unzip::Bool, nohybridladder::Bool, verbose::Bool,
    constraints::Vector{TopologyConstraint})
    hybridnode = obj.net.hybrid[Random.rand(1:length(obj.net.hybrid))]
    minorhybridedge = getMinorParentEdge(hybridnode)
    if length(constraints) > 0
        blacklist = Node[]
        edgefound = false
        while length(blacklist) < obj.net.numHybrids && !edgefound
            hybridnode = obj.net.hybrid[Random.rand(1:length(obj.net.hybrid))]
            hybridnode ∉ blacklist || continue # if node already attempted, jump to top of while
            edgefound = true
            for c in constraints
                if minorhybridedge == c.edge # edge to remove is stem edge
                    push!(blacklist, hybridnode)
                    edgefound = false
                    break # out of constraint loop
                end
            end
        end
        if !edgefound # tried to remove all hybrids
            return nothing
        end
    end
    savededges, savedhybridedges = savelocalBLgamma(obj.net, minorhybridedge)
    setGamma!(minorhybridedge, 0.0)
    optimizelocalBL!(obj, obj.net, minorhybridedge, unzip, verbose) # updates loglik
    optimizelocalgammas!(obj, obj.net, minorhybridedge, unzip, verbose)
    if obj.loglik - currLik > likAbsDelHybLiNC # -0.1: loglik can decrease for parsimony
        e1 = getChildEdge(hybridnode)
        deletehybridedge!(obj.net, minorhybridedge, false, false)
        updateSSM!(obj, true)
        return true
    else # keep hybrid
        resetlocalBLgamma!(obj.net, savededges, savedhybridedges)
        obj.loglik = currLik
        return false
    end
end
"""
    updateSSM!(obj::SSM, renumber=false::Bool)

After adding or removing a hybrid, displayed trees will change. Updates
the displayed tree list. Return SSM object.

if `renumber`, reorder edge and internal node numbers. Only need
to renumber after deleting a hybrid (which could remove edges and nodes
from the middle of the edge and node lists).

Assumptions:
- The SSM object has cache arrays of size large enough, that is,
  the constructor [`StatisticalSubstitutionModel`](@ref) was previously
  called with maxhybrid equal or greater than in `obj.net`.
  `obj.priorltw` is not part of the "cache" arrays.

Warnings:
Does not update the likelihood.

```jldoctest
julia> maxhybrid = 3;

julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple.aln"));

julia> obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69, maxhybrid);

julia> PhyloNetworks.checknetworkbeforeLiNC!(obj.net, maxhybrid, true, true, true);

julia> PhyloNetworks.discrete_corelikelihood!(obj);

julia> PhyloNetworks.addhybridedge!(obj.net, obj.net.edge[8], obj.net.edge[1], true, 0.0, 0.4);

julia> PhyloNetworks.updateSSM!(obj);

julia> writeTopology(obj.net)
"(((B:0.0)#H1:0.1::0.9,(A:1.0)#H2:1.0::0.6):1.5,(C:0.6,#H1:1.0::0.1):1.0,(D:1.25,#H2:0.0::0.4):1.25);"
```
"""
function updateSSM!(obj::SSM, renumber=false::Bool)
    if renumber # traits are in leaf.number order, so leaf nodes not reordered
        resetNodeNumbers!(obj.net; checkPreorder=true, internalonly=true)
        resetEdgeNumbers!(obj.net)
    end
    # extract displayed trees
    obj.displayedtree = displayedTrees(obj.net, 0.0; keepNodes=true)
    nnodes = length(obj.net.node)
    for tree in obj.displayedtree
        preorder!(tree) # no need to call directEdges! before: already done on net
        length(tree.nodes_changed) == nnodes ||
            error("displayed tree with too few nodes: $(writeTopology(tree))")
        length(tree.edge) == length(obj.net.edge)-obj.net.numHybrids ||
            error("displayed tree with too few edges: $(writeTopology(tree))")
    end
    # log tree weights: sum log(γ) over edges, for each displayed tree
    obj.priorltw = inheritanceWeight.(obj.displayedtree)
    @debug begin
        all(!ismissing, obj.priorltw) ? "" :
        "one or more inheritance γ's are missing or negative. fix using setGamma!(network, edge)"
    end
    return obj
end

## Optimize Branch Lengths and Gammas ##
"""
    startingBL!(net::HybridNetwork, unzip::Bool,
                trait::AbstractVector{Vector{Union{Missings.Missing,Int}}},
                siteweight=ones(length(trait[1]))::AbstractVector{Float64})

Calibrate branch lengths in `net` by minimizing the mean squared error
between the JC-adjusted pairwise distance between taxa, and network-predicted
pairwise distances, using [`calibrateFromPairwiseDistances!`](@ref).
`siteweight[k]` gives the weight of site (or site pattern) `k` (default: all 1s).
`unzip` = true sets all edges below a hybrid node to length zero.

Assumptions:

- all species have the same number of traits (sites): `length(trait[i])` constant
- `trait[i]` is for leaf with `node.number = i` in `net`, and
  `trait[i][j] = k` means that leaf number `i` has state index `k` for trait `j`.
  These indices are those used in a substitution model:
  kth value of `getlabels(model)`.
- Hamming distances are < 0.75 with four states, or < (n-1)/n for n states.
  If not, all pairwise hamming distances are scaled by `.75/(m*1.01)` where `m`
  is the maximum observed hamming distance, to make them all < 0.75.
"""
function startingBL!(net::HybridNetwork, unzip::Bool,
        trait::AbstractVector{Vector{Union{Missings.Missing,Int}}},
        siteweight=ones(length(trait[1]))::AbstractVector{Float64})
    nspecies = net.numTaxa
    M = zeros(Float64, nspecies, nspecies) # pairwise distances initialized to 0
    # count pairwise differences, then multiply by pattern weight
    ncols = length(trait[1]) # assumption: all species have same # columns
    length(siteweight) == ncols ||
      error("$(length(siteweight)) site weights but $ncols columns in the data")
    for i in 2:nspecies
        species1 = trait[i]
        for j in 1:(i-1)
            species2 = trait[j]
            for col in 1:ncols
                if !(ismissing(species1[col]) || ismissing(species2[col])) &&
                    (species1[col] != species2[col])
                    M[i, j] += siteweight[col]
                end
            end
            M[j,i] = M[i,j]
        end
    end
    Mp = M ./ sum(siteweight) # to get proportion of sites, for each pair

    # estimate pairwise evolutionary distances using extended Jukes Cantor model
    nstates = mapreduce(x -> maximum(skipmissing(x)), max, trait)
    maxdist = (nstates-1)/nstates
    Mp[:] = Mp ./ max(maxdist, maximum(Mp*1.01)) # values in [0,0.9901]: log(1-Mp) well defined
    dhat = - maxdist .* log.( 1.0 .- Mp)

    taxonnames = [net.leaf[i].name for i in sortperm([n.number for n in net.leaf])]
    # taxon names: to tell the calibration that row i of dhat if for taxonnames[i]
    # ASSUMPTION: trait[i][j] = trait j for taxon at node number i: 'node.number' = i
    calibrateFromPairwiseDistances!(net, dhat, taxonnames,
        forceMinorLength0=false, ultrametric=false)

    if unzip && length(net.hybrid) > 0
        constrainededges = [getChildEdge(h) for h in net.hybrid]
        setlengths!(constrainededges, zeros(length(constrainededges)))
    end
    return net
end

"""
    optimizelocalBL!(obj::SSM, net::HybridNetwork, edge::Edge, unzip::Bool)

Optimize branch lengths in `net` locally around `edge`. Update all edges that
share a node with `edge` (including itself).
If `unzip = true`, constrain branch lengths to zero below hybrid edges.
Return vector of updated `edges`.

Used after `nni!` or `addhybridedge!` moves to update local branch lengths.

```jldoctest
julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
HybridNetwork, Rooted Network
9 edges
9 nodes: 4 tips, 1 hybrid nodes, 4 internal tree nodes.
tip labels: A, B, C, D
(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);


julia> fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple.aln"));

julia> obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69);

julia> obj.net.edge[4]
PhyloNetworks.Edge:
 number:4
 length:1.5
 attached to 2 node(s) (parent first): 6 8


julia> PhyloNetworks.optimizelocalBL!(obj, obj.net, obj.net.edge[4], true);

julia> obj.net.edge[4]
PhyloNetworks.Edge:
 number:4
 length:0.0
 attached to 2 node(s) (parent first): 6 8
```
"""
function optimizelocalBL!(obj::SSM, net::HybridNetwork, edge::Edge, unzip::Bool,
    verbose=false::Bool)
    edges = Edge[]
    for n in edge.node
        for e in n.edge # all edges sharing a node with `edge` (including self)
            if !(e in edges)
                push!(edges, e)
            end
        end
    end
    optimizeBL!(obj, net, edges, unzip, verbose, 10, :LD_MMA, fRelBL, fAbsBL,
    xRelBL, xAbsBL) #maxeval = 10
    return edges
end

"""
    optimizeBL!(obj::SSM, net::HybridNetwork, edges::Vector{Edge}, unzip::Bool,
                verbose=false::Bool, maxeval=1000::Int64,
                NLoptMethod=:LD_MMA::Symbol, ftolRel=fRelBL::Float64,
                ftolAbs=fAbsBL::Float64, xtolRel=xRelBL::Float64,
                xtolAbs=xAbsBL::Float64)

Optimize branch lengths for edges in vector `edges`.
If `unzip = true`, constrain branch lengths to zero below hybrid edges.
Return vector of updated `edges`.
"""
function optimizeBL!(obj::SSM, net::HybridNetwork, edges::Vector{Edge},
    unzip::Bool, verbose=false::Bool, maxeval=1000::Int64, NLoptMethod=:LD_MMA::Symbol,
    ftolRel=fRelBL::Float64, ftolAbs=fAbsBL::Float64, xtolRel=xRelBL::Float64,
    xtolAbs=xAbsBL::Float64)
    if unzip && length(net.hybrid) > 0
        constrainededges = [getChildEdge(h) for h in net.hybrid]
        setlengths!(constrainededges, zeros(length(constrainededges)))
        edges = setdiff(edges, constrainededges) # edges - constrainededges
    end
    if any([e.length < 0.0 for e in edges])
        for e in edges #? is there a faster way to do this?
            if e.length < 0.0
                e.length = 0.0
            end
        end
    end
    counter = [0]
    function loglikfunBL(lengths::Vector{Float64}, grad::Vector{Float64})
        counter[1] += 1
        setlengths!(edges, lengths)
        res = discrete_corelikelihood!(obj)
        # verbose && println("loglik: $res, branch lengths: $(lengths)")
        length(grad) == 0 || error("gradient not implemented")
        return res
    end
    # set-up optimization object for BL parameter
    NLoptMethod=:LN_COBYLA # no gradient
    # :LN_COBYLA for (non)linear constraits, :LN_BOBYQA for bound constraints
    nparBL = length(edges)
    optBL = NLopt.Opt(NLoptMethod, nparBL)
    NLopt.ftol_rel!(optBL,ftolRel) # relative criterion
    NLopt.ftol_abs!(optBL,ftolAbs) # absolute criterion
    NLopt.xtol_rel!(optBL,xtolRel)
    NLopt.xtol_abs!(optBL,xtolAbs)
    NLopt.maxeval!(optBL, maxeval) # max number of iterations
    # NLopt.maxtime!(optBL, t::Real)
    NLopt.lower_bounds!(optBL, zeros(length(edges)))
    counter[1] = 0
    NLopt.max_objective!(optBL, loglikfunBL)
    fmax, xmax, ret = NLopt.optimize(optBL, getlengths(edges))
    verbose && println("BL: got $(round(fmax, digits=5)) at
    BL = $(round.(xmax, digits=5)) after $(counter[1]) iterations
    (return code $(ret))")
    return edges
end

"""
    optimizeallgammas!(obj::SSM, net::HybridNetwork, unzip::Bool,
                        verbose=false::Bool, NLoptMethod=:LD_MMA::Symbol,
                        ftolRel=fRelBL::Float64, ftolAbs=fAbsBL::Float64,
                        xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64)

Optimize all gammas in a network. Creates a list of containing one parent edge
per hybrid then calls optimizegammas! on that list.
"""
function optimizeallgammas!(obj::SSM, net::HybridNetwork, unzip::Bool,
    verbose=false::Bool, maxeval=1000::Int64, NLoptMethod=:LD_MMA::Symbol, ftolRel=fRelBL::Float64,
    ftolAbs=fAbsBL::Float64, xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64)
    edges = Edge[]
    for h in net.hybrid
        push!(edges, getMajorParentEdge(h))
    end
    if length(edges) == 0
        @debug "no gammas in network to optimize because there are $(length(obj.net.hybrid)) hybrids"
    else
        optimizegammas!(obj, net, edges, unzip, verbose, maxeval, :LD_MMA,
        fRelBL, fAbsBL, xRelBL, xAbsBL)
    end
end

"""
    optimizelocalgammas!(obj::SSM, net::HybridNetwork, edge::Edge, unzip::Bool,
                            verbose=false::Bool)

Optimize gammas in `net` locally around `edge`. Update all edges that share a
node with `edge` (including itself). Does not include edges' partners because
the `setGamma!` updates partners automatically.

Return modified edges.

Used after `nni!` or `addhybridedge!` moves to update local gammas.

Assumptions:
- correct `isChild1` field for `edge` and for hybrid edges
- no in-coming polytomy: a node has 0, 1 or 2 parents, no more

```jldoctest
julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple.aln"));

julia> obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69);

julia> obj.net.hybrid[1].edge
3-element Array{PhyloNetworks.Edge,1}:
 PhyloNetworks.Edge:
 number:2
 length:1.0
 attached to 2 node(s) (parent first): 9 2

 PhyloNetworks.Edge:
 number:3
 length:0.1
 major hybrid edge with gamma=0.9
 attached to 2 node(s) (parent first): 8 9

 PhyloNetworks.Edge:
 number:6
 length:1.0
 minor hybrid edge with gamma=0.1
 attached to 2 node(s) (parent first): 7 9

julia> PhyloNetworks.optimizelocalgammas!(obj, net, obj.net.hybrid[1].edge[2], true);

julia> obj.net.hybrid[1].edge
3-element Array{PhyloNetworks.Edge,1}:
 PhyloNetworks.Edge:
 number:2
 length:1.0
 attached to 2 node(s) (parent first): 9 2

 PhyloNetworks.Edge:
 number:3
 length:0.1
 major hybrid edge with gamma=0.9000000000259154
 attached to 2 node(s) (parent first): 8 9

 PhyloNetworks.Edge:
 number:6
 length:1.0
 minor hybrid edge with gamma=0.0999999999740846
 attached to 2 node(s) (parent first): 7 9

````
"""
function optimizelocalgammas!(obj::SSM, net::HybridNetwork, edge::Edge,
    unzip::Bool, verbose=false::Bool)
    edges = Edge[]
    for n in edge.node
        for e in n.edge # edges that share a node with `edge` (including self)
            if e.hybrid && !(e in edges) && !(getPartner(e) in edges)
                push!(edges, e)
            end
        end
    end
    if length(edges) == 0
        @debug "no local gammas to optimize around edge $edge"
    else
        optimizegammas!(obj, net, edges, unzip, verbose, 10, :LD_MMA, fRelBL,
            fAbsBL, xRelBL, xAbsBL) # maxeval = 10
    end
end

"""
    optimizegammas!(obj::SSM, net::HybridNetwork, edges::Vector{Edge},
                    unzip::Bool, verbose=false::Bool, maxeval=1000::Int64,
                    NLoptMethod=:LD_MMA::Symbol, ftolRel=fRelBL::Float64,
                    ftolAbs=fAbsBL::Float64, xtolRel=xRelBL::Float64,
                    xtolAbs=xAbsBL::Float64)

Optimize gammas for hybrid edges in vector `edges`.
Return vector of updated `edges`.

Assumptions:
- Edges vector does not contain hybrid partners.

Warning: Do not call directly. Instead use optimizelocalgammas!
or optimizeallgamma! functions.
"""
function optimizegammas!(obj::SSM, net::HybridNetwork, edges::Vector{Edge},
    unzip::Bool, verbose=false::Bool, maxeval=1000::Int64,
    NLoptMethod=:LD_MMA::Symbol, ftolRel=fRelBL::Float64,
    ftolAbs=fAbsBL::Float64, xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64)
    if any([(e.gamma < 0.0 || e.gamma > 1.0) for e in edges])
        for e in edges #? is there a faster way to do this?
            if e.gamma < 0.0 || e.gamma > 1.0
                setGamma!(e, 0.5)
            end
        end
    end
    counter = [0]
    function loglikfungamma(gammas::Vector{Float64}, grad::Vector{Float64})
        counter[1] += 1
        setmultiplegammas!(edges, gammas)
        res = discrete_corelikelihood!(obj)
        # verbose && println("loglik: $res, gammas: $(gammas)")
        length(grad) == 0 || error("gradient not implemented")
        return res
    end
    # set-up optimization object for gamma parameter
    NLoptMethod=:LN_COBYLA # no gradient
    # :LN_COBYLA for (non)linear constraits, :LN_BOBYQA for bound constraints
    npargamma = length(edges)
    optgamma = NLopt.Opt(NLoptMethod, npargamma)
    NLopt.ftol_rel!(optgamma,ftolRel) # relative criterion
    NLopt.ftol_abs!(optgamma,ftolAbs) # absolute criterion
    NLopt.xtol_rel!(optgamma,xtolRel)
    NLopt.xtol_abs!(optgamma,xtolAbs)
    NLopt.maxeval!(optgamma, maxeval) # max number of iterations
    #NLopt.initial_step!(optgamma, 0.05) # step size
    # NLopt.maxtime!(optgamma, t::Real)
    NLopt.lower_bounds!(optgamma, zeros(Float64, npargamma))
    NLopt.upper_bounds!(optgamma, ones(Float64, npargamma))
    counter[1] = 0
    NLopt.max_objective!(optgamma, loglikfungamma)
    fmax, xmax, ret = NLopt.optimize(optgamma, [e.gamma for e in edges])
    verbose && println("gamma: got $(round(fmax, digits=5)) at
    $(round.(xmax, digits=5)) after $(counter[1]) iterations
    (return code $(ret))")
    return edges
end

## Prep Functions ##

# fixit: new type for two (dependent) binary traits
# need new algorithm for model at hybrid nodes: displayed trees aren't enough
"""
    learnlabels(model::Symbol, dat::DataFrame)

Return unique non-missing values in `dat`, and check that these labels
can be used to construct of substitution model of type `model`.

# examples:

```jldoctest
julia> using DataFrames

julia> dat = DataFrame(trait1 = ["A", "C", "A", missing]); # 4×1 DataFrame

julia> PhyloNetworks.learnlabels(:BTSM, dat)
2-element Array{String,1}:
 "A"
 "C"

julia> PhyloNetworks.learnlabels(:JC69, dat)
2-element Array{String,1}:
 "A"
 "C"
```
"""
function learnlabels(modSymbol::Symbol, dat::AbstractDataFrame)
    labels = mapreduce(x -> unique(skipmissing(x)), union, eachcol(dat))
    if modSymbol == :BTSM
        length(labels) == 2 || error("Binary Trait Substitution Model supports traits with two states. These data have do not have two states.")
    elseif modSymbol == :TBTSM
        unique(skipmissing(dat[!,1])) == 2 && unique(skipmissing(dat[!,2]) == 2) ||
          error("Two Binary Trait Substitution Model supports two traits with two states each.")
    elseif modSymbol in [:HKY85, :JC69]
        typeof(labels) == Array{DNA,1} ||
        (typeof(labels) == Array{Char,1}   &&       all(in.(uppercase.(labels), "-ABCDGHKMNRSTVWY"))) ||
        (typeof(labels) == Array{String,1} && all(occursin.(uppercase.(labels), "-ABCDGHKMNRSTVWY"))) ||
        # "ACGT" would dissallow ambiguous sites
          error("$modSymbol requires that trait data are dna bases A, C, G, and T")
    end
    return labels
end

"""
    startingrate(net)

Estimate an evolutionary rate appropriate for the branch lengths in the network,
which should be a good starting value before optimization in `fitdiscrete`,
assuming approximately 1 change across the entire tree.
If all edge lengths are missing, set starting rate to 1/(number of taxa).
"""
function startingrate(net::HybridNetwork)
    totaledgelength = 0.0
    for e in net.edge
        if e.length > 0.0
            totaledgelength += e.length
        end
    end
    if totaledgelength == 0.0 # as when all edge lengths are missing
        totaledgelength = net.numTaxa
    end
    return 1.0/totaledgelength
end

"""
    symboltomodel(network, modsymbol::Symbol, data::DataFrame,
                    siteweights::Vector)

Return a statistical substitution model (SSM) with appropriate state labels
and a rate appropriate for the branch lengths in `net`
(see [`startingrate`](@ref)).
The `data` frame must have the actual trait/site data in columns 2 and up,
as when the species names are in column 1.
For DNA data, the relative rate model is returned, with a
stationary distribution equal to empirical frequencies for DNA data.
"""
function symboltomodel(net::HybridNetwork, modsymbol::Symbol, data::DataFrame,
        siteweights=repeat([1.], inner=size(data,2))::AbstractVector)
    rate = startingrate(net)
    actualdat = view(data, :, 2:size(data,2))
    labels = learnlabels(modsymbol, actualdat)
    if modsymbol == :JC69
        return JC69([1.0], true) # 1.0 instead of rate because relative version
    elseif modsymbol == :HKY85 # transition/transversion rate ratio
        return HKY85([1.0], empiricalDNAfrequencies(actualdat, siteweights),
        true)
    elseif modsymbol == :ERSM
        return EqualRatesSubstitutionModel(length(labels), rate, labels);
    elseif modsymbol == :BTSM
        return BinaryTraitSubstitutionModel([rate, rate], labels)
    elseif modsymbol == :TBTSM
        return TwoBinaryTraitSubstitutionModel([rate, rate, rate, rate, rate,
        rate, rate, rate], labels)
    else
        error("model $modsymbol is unknown or not implemented yet")
    end
end

"""
    savelocalBLgamma(net::HybridNetwork, edge::Edge)

Saves local branch lengths and gammas before they're optimized so they can be
reset.
Return a tuple of two dictionaries holding edges and branch lengths
and edges and gammas, respectively. Each dictionary is of this format:
Dict{Edge, Float64}
"""
function savelocalBLgamma(net::HybridNetwork, edge::Edge)
    localedges = Dict{Edge, Float64}()
    for n in edge.node
        for e in n.edge # all edges sharing a node with `edge` (including self)
            if !(haskey(localedges, e))
                localedges[e] = e.length
            end
        end
    end
    localhybridedges = Dict{Edge, Float64}()
    for n in edge.node
        for e in n.edge # edges that share a node with `edge` (including self)
            if e.hybrid && !haskey(localhybridedges, e) && !haskey(localhybridedges, getPartner(e))
                localhybridedges[e] = e.gamma
            end
        end
    end
    return localedges, localhybridedges
end

"""
    resetlocalBLgamma!(net::HybridNetwork, localedges::Dict{Edge, Float64},
                        localhybridedges::Dict{Edge, Float64})

Reset local branch lengths and gammas to undo a local optimization.
Return net.
"""
function resetlocalBLgamma!(net::HybridNetwork, localedges::Dict{Edge, Float64},
    localhybridedges::Dict{Edge, Float64})
    for edge in keys(localedges) # set lengths
        edge.length = localedges[edge]
    end
    for hybridedge in keys(localhybridedges) # set gammas
        hybridedge.gamma = localhybridedges[hybridedge]
    end
    return net
end
