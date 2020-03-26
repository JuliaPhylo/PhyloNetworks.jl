# any change to these constants must be documented in phyLiNC!
const moveweights_LiNC = Distributions.aweights([0.4, 0.2, 0.2, 0.2])
const movelist_LiNC = ["nni", "addhybrid", "deletehybrid", "root"]
const likAbsAddHybLiNC = 0.5 #= loglik improvement required to retain a hybrid
  greater values raise the standard for newly-proposed hybrids,
  leading to fewer proposed hybrids accepted during the search =#
const likAbsDelHybLiNC = -0.1 #= loglik decrease allowed when removing a hybrid
  lower (more negative) values of lead to more hybrids removed during the search =#
const alphaRASmin = 0.02
const alphaRASmax = 50.0

"""
    CacheGammaLiNC

Type to store intermediate values during γ optimization,
to limit time spent on garbage collection.
"""
struct CacheGammaLiNC
    "conditional likelihood under focus hybrid edge, length nsites"
    clike::Vector{Float64}
    "conditional likelihood under partner edge, length nsites"
    clikp::Vector{Float64}
    """whether displayed trees have the focus hybrid edge (true)
       the partner edge (false) or none of them (missing),
       which can happen on non-tree-child networks"""
    hase::Vector{Union{Missing, Bool}}
end
function CacheGammaLiNC(obj::SSM)
    CacheGammaLiNC(
        Vector{Float64}(undef, obj.nsites),
        Vector{Float64}(undef, obj.nsites),
        Vector{Union{Missing, Bool}}(undef, length(obj.displayedtree))
    )
end

"""
    phyLiNC!(net::HybridNetwork, fastafile::String, modSymbol::Symbol,
             rateCategories=1::Int)

Estimate a phylogenetic network from concatenated DNA data using
maximum likelihood, ignoring incomplete lineage sorting
(phyLiNC: phylogenetic Likelihood Network from Concatenated data).
The network is constrained to have `maxhybrid` reticulations at most,
but can be of any level.
The search starts at (or near) the network `net`,
using a local hill-climbing search to optimize the topology
(nearest-neighbor interchange moves, add hybridizations,
and remove hybridizations). Also optimized are evolutionary rates,
amount of rate variation across sites, branch lengths and inheritance γs.
This search strategy is run `nruns` times, and the best of the `nruns`
networks is returned.

Return a [`StatisticalSubstitutionModel`](@ref) object, say `obj`, which
contains the estimated network in `obj.net`.

The length of the edge below a reticulation is not identifiable.
Therefore, phyLiNC estimates the canonical version of the network: with
reticulations **unzipped**: edges below reticulations are set to 0, and
hybrid edges (parental lineages) have estimated lengths that are
increased accordingly.

Optional arguments include (default value in parenthesis):
- `rateCategories` (1): number of categories to use in estimating variable evolutionary
  rates using a discretized gamma model. If `rateCategories` = 1, no rate variation
  is included. To allow for rate variation, four categories is typically used.
  See [`RateVariationAcrossSites`] (@ref)
- `nruns` (10): number of independent starting points for the search
- `filename` ("phyLiNC"): root name for the output files (`.out`, `.err`).
  If empty (""), files are *not* created, progress log goes to the screen only
  (standard out).
- `maxhybrid` (1): maximum number of hybridizations allowed
- `no3cycle` (true): prevents 3-cycles, which are (almost) not
  identifiable
- `nohybridladder` (true): prevents hybrid ladder in network. If true,
  the input network must not have hybrid ladders.
- `maxmoves` (100): maximum number of topology moves before branch lengths,
  hybrid γ values, evolutionary rates, and rate variation parameters are
  reestimated.
- `verbose` (true): set to false to turn off screen output
- `seed` (default 0 to get it from the clock): seed to replicate a given search
- `probST` (0.5): probability to use `net` as the starting topology
  for each given run. If probST < 1, the starting topology is k NNI moves
  away from `net`, where k is drawn from a geometric distribution: p (1-p)ᵏ,
  with success probability p = `probST`.
- `speciesfile` (""): path to a csv file containing two columns:
  species and individuals used to create one or more species topology
  constraints to meet during the search.
- `cladefile` (""): path to a csv file containing two columns:
  clades and individuals used to create one or more clade topology
  constraints to meet during the search.
  (NOTE: clade contraints not yet implemented.)
- `alphamin` (0.02): minimum value for shape parameter alpha in rate variation
  across sites model.
- `alphamax` (50.0): maximum value for shape parameter alpha in rate variation
  across sites model.

The following optional arguments control when to stop the optimization of branch
lengths and gamma values on each individual candidate network. Defaults in
parentheses.
- `ftolRel` (1e-6) and `ftolAbs` (1e-6): relative and absolute differences of the
  network score between the current and proposed parameters
- `xtolRel` (1e-5) and `xtolAbs` (1e-5): relative and absolute differences
  between the current and proposed parameters.
Greater values will result in a less thorough but faster search. These parameters
are used when evaluating candidate networks only.
Regardless of these arguments, once a final topology is chosen, branch lenghts
are optimized using stricter tolerances (1e-10, 1e-12, 1e-10, 1e-10) for better
estimates.

The following optional arguments control when to stop proposing new
network topologies:

- `nreject` (75): maximum number of times that new topologies are
  proposed and rejected in a row. Lower values of `nreject` result in a less
  thorough but faster search.
"""
function phyLiNC!(net::HybridNetwork, fastafile::String, modSymbol::Symbol,
                  rateCategories=1::Int;
                  maxhybrid=1::Int, no3cycle=true::Bool,
                  nohybridladder=true::Bool,
                  speciesfile=""::AbstractString,
                  cladefile=""::AbstractString, verbose=true::Bool,
                  kwargs...)
    # create constraints and new network if speciesfile
    if !isempty(speciesfile)
        net, constraints = mapindividuals(net, speciesfile)
    else
        constraints = TopologyConstraint[]
    end
    if !isempty(cladefile)
        error("Clade constraints not yet implemented.")
    end
    # create starting object for all runs
    obj = StatisticalSubstitutionModel(net, fastafile, modSymbol, rateCategories, maxhybrid)
    #= after SSM(), update constraint taxon names, taxonnums, edge, and node
       because some leaves may be pruned, and check_matchtaxonnames calls
       resetNodeNumbers! (changing leaf node numbers) and resetEdgeNumbers! =#
    for i in 1:length(constraints)
        constraints[i] = PhyloNetworks.TopologyConstraint(constraints[i].type,
                                            constraints[i].taxonnames, obj.net)
    end
    checknetwork_LiNC!(obj.net, maxhybrid, no3cycle, nohybridladder, constraints, verbose)
    #= checknetwork removes degree-2 nodes (including root) and 2- and 3-cycles
       and requires that the network is preordered.
       Warning: need to call updateSSM after using checknetwork_LiNC =#
    updateSSM!(obj, true; constraints=constraints)
    startingBL!(obj.net, true, obj.trait, obj.siteweight) # true: to unzip
    phyLiNC!(obj; maxhybrid=maxhybrid, no3cycle=no3cycle, nohybridladder=nohybridladder,
            constraints=constraints, verbose=verbose, kwargs...)
end
function phyLiNC!(obj::SSM;
                  maxhybrid=1::Int, no3cycle=true::Bool,
                  nohybridladder=true::Bool, maxmoves=100::Int, nreject=75::Int,
                  nruns=10::Int, filename="phyLiNC"::AbstractString,
                  verbose=true::Bool, seed=0::Int, probST=0.5::Float64,
                  ftolRel=1e-6::Float64, ftolAbs=1e-6::Float64,
                  xtolRel=1e-10::Float64, xtolAbs=1e-5::Float64,
                  constraints=TopologyConstraint[]::Vector{TopologyConstraint},
                  alphamin=alphaRASmin::Float64, alphamax=alphaRASmax::Float64)
    writelog = true
    writelog_1proc = false
    if filename != ""
        logfile = open(string(filename,".log"),"w")
        juliaout = string(filename,".out")
        if Distributed.nprocs() == 1
            writelog_1proc = true
            juliaerr = string(filename,".err")
            errfile = open(juliaerr,"w")
        end
    else
      writelog = false
      logfile = stdout
    end
    γcache = CacheGammaLiNC(obj)
    # rough optimization of rates and alpha, for better starting values used by all runs
    obj.loglik = -Inf
    fit!(obj; optimizeQ=(nparams(obj.model) > 0),
         optimizeRVAS=(nparams(obj.ratemodel) > 0),
         verbose=false, maxeval=20,
         ftolRel=ftolRel, ftolAbs=ftolAbs, xtolRel=xtolRel, xtolAbs=xtolAbs)
    @debug "loglik = $(loglikelihood(obj)) at the start"
    str = """
    PhyLiNC network estimation starting. Parameters:
       maxhybrid = $(maxhybrid)
       nohybridladder = $nohybridladder
       no3cycle = $no3cycle
       probST= $probST
       max number of moves per cycle = $maxmoves
       max number of consecutive failed proposals = $(nreject)
       optimization tolerance: ftolRel=$(ftolRel), ftolAbs=$(ftolAbs),
                               xtolAbs=$(xtolAbs), xtolRel=$(xtolRel)."""
    str *= (writelog ? "\n   filename for log and err files: $(filename)" :
                       "\n   no output files\n\n")
    str *= "\n$(nruns) run(s) starting near network topology:\n$(writeTopology(obj.net))\nstarting model:\n" *
            string(obj.model) * string(obj.ratemodel)
    # fixit: add info about constraints: type and tip names for each constraint
    if Distributed.nprocs()>1
        str *= "using $(Distributed.nworkers()) worker processors\n"
    end
    if writelog
      write(logfile,str)
      flush(logfile)
    end
    verbose && print(stdout,str)
    verbose && print(stdout, "Time: " * Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s") * "\n")
    # if 1 proc: time printed to logfile at start of every run, not here.

    if seed == 0
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) # seed based on clock
    end
    if writelog
      write(logfile,"main seed = $(seed)\n---------------------\n")
      flush(logfile)
    end
    verbose && print(stdout,"main seed = $(seed)\n---------------------\n")
    Random.seed!(seed)
    seeds = [seed; round.(Integer,floor.(rand(nruns-1)*100000))]

    if writelog && !writelog_1proc
        for i in 1:nruns # workers won't write to logfile
            write(logfile, "For run $(i), seed = $(seeds[i])\n")
        end
        flush(logfile)
    end

    tstart = time_ns()
    startingnet = obj.net
    startingconstraints = constraints
    netvector = Distributed.pmap(1:nruns) do i # for i in 1:nruns
        msg = "BEGIN PhyLiNC run $(i)\nseed = $(seeds[i])\ntime = $(Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s"))\n"
        if writelog_1proc # workers can't write on streams opened by master
            write(logfile, msg)
            flush(logfile)
        end
        verbose && print(stdout, msg)
        GC.gc()
        try
            obj.net = deepcopy(startingnet)
            obj.model = deepcopy(obj.model)
            obj.ratemodel = deepcopy(obj.ratemodel)
            constraints = deepcopy(startingconstraints) # problem: nodes were copied on previous line. when nodes copied again on this line, they are now different
            phyLiNCone!(obj, maxhybrid, no3cycle, nohybridladder, maxmoves,
                        nreject, verbose, writelog_1proc, logfile, seeds[i], probST,
                        constraints, ftolRel, ftolAbs, xtolRel, xtolAbs,
                        alphamin, alphamax, γcache)
            logstr = "\nFINISHED. loglik = $(obj.loglik)\n"
            verbose && print(stdout, logstr)
            if writelog_1proc
                logstr *= writeTopology(obj.net)
                logstr *= "\n---------------------\n"
                write(logfile, logstr)
                flush(logfile)
            end
            obj.net.loglik = obj.loglik
            return [obj.net, obj.model, obj.ratemodel]
        catch err
            msg = "\nERROR found on PhyLiNC for run $(i) seed $(seeds[i]):\n" *
                  sprint(showerror,err)
            logstr = msg * "\n---------------------\n"
            print(stacktrace(catch_backtrace()))
            println()
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
    telapsed = round(convert(Int, tend-tstart) * 1e-9, digits=2) # in seconds
    writelog_1proc && close(errfile)
    msg = "All runs complete.\nend time: " * Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s") *
            "\ntime elapsed: $telapsed seconds" *"\n---------------------\n"
    writelog && write(logfile, msg)
    verbose && print(stdout, msg)

    #= post processing of the networks
      type of netvector: Array{Union{Nothing, Array{Int64,1}},1}
      each item holds [net, model, ratemodel] if sucessful, nothing if failed =#
    filter!(n -> n !== nothing, netvector) # remove "nothing": failed runs
    !isempty(netvector) || error("all runs failed")
    sort!(netvector, by = x -> x[1].loglik, rev=true)
    # @debug "loglik from all runs:" [n[1].loglik for n in netvector]
    obj.net = netvector[1][1]::HybridNetwork # best network, tell type to compiler
    obj.model = netvector[1][2]
    obj.ratemodel = netvector[1][3]
    obj.loglik = obj.net.loglik
    updateSSM!(obj, true; constraints = constraints) # topology has changed, need to update displayedtree, priorltw
    logstr = "Best topology:\n$(writeTopology(obj.net))\n" *
              "with loglik $(obj.loglik) under:\n" * string(obj.model) *
              string(obj.ratemodel) * "---------------------\n" *
              "Final optimization of branch lengths and gammas on this network... "
    if writelog
        write(logfile, logstr)
        flush(logfile)
    end
    verbose && print(stdout,logstr)
    # warning: tolerance values from constants, not user-specified
    optimizeBL_LiNC!(obj, copy(obj.net.edge), fRelBL, fAbsBL, xRelBL, xAbsBL,
                     max(10*length(obj.net.edge), 1000))
    optimizeallgammas_LiNC!(obj, fAbsBL, γcache, 100) &&
        updateSSM!(obj, true; constraints=constraints)
    toptimend = time_ns() # in nanoseconds
    telapsed = round(convert(Int, toptimend-tstart) * 1e-9, digits=2) # in seconds
    logstr = "complete.\nFinal log-likelihood: $(obj.loglik)\n" *
             "Final network:\n" * "$(writeTopology(obj.net))\n" *
             "Total time elapsed: $telapsed seconds (includes final branch length and gamma optimization)\n" *
             "Final time: " * Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s") * "\n"
    if writelog
        write(logfile, logstr)
        flush(logfile)
    end
    verbose && print(stdout,logstr)
    return obj
end

"""
    phyLiNCone!(obj::SSM, maxhybrid::Int, no3cycle::Bool,
                nohybridladder::Bool, maxmoves::Int, nrejectmax::Int,
                verbose::Bool, writelog_1proc::Bool, logfile::IO,
                seed::Int, probST::Float64,
                constraints::Vector{TopologyConstraint},
                ftolRel::Float64, ftolAbs::Float64,
                xtolRel::Float64, xtolAbs::Float64,
                alphamin::Float64, alphamax::Float64,
                γcache::CacheGammaLiNC

Estimate one phylogenetic network (or tree) from concatenated DNA data,
like [`phyLiNC!`](@ref), but doing one run only, and taking as input an
StatisticalSubstitutionModel object `obj`. The starting network is `obj.net`
and is assumed to meet all the requirements.

`writelog_1proc` is passed by phyLiNC! an indicates if a log should be written.
If the number of processors is > 1, this will be false because workers can't
write on streams opened by master. `logfile` will be stdout if `writelog_1proc`
is false. Otherwise, it will be the log file created by `phyLiNC!`.

See [`phyLiNC!`](@ref) for other arguments.
"""
function phyLiNCone!(obj::SSM, maxhybrid::Int, no3cycle::Bool,
                    nohybridladder::Bool, maxmoves::Int, nrejectmax::Int,
                    verbose::Bool, writelog_1proc::Bool, logfile::IO,
                    seed::Int, probST::Float64,
                    constraints::Vector{TopologyConstraint},
                    ftolRel::Float64, ftolAbs::Float64,
                    xtolRel::Float64, xtolAbs::Float64,
                    alphamin::Float64, alphamax::Float64,
                    γcache::CacheGammaLiNC)

    Random.seed!(seed)
    # update files inside topology constraint to match the deepcopied net in obj.net created by pmap
    updateconstraintfields!(constraints, obj.net)
    if probST < 1.0 # modify starting tree by k nni moves (if possible), k=0 or more
        numNNI = rand(Geometric(probST)) # number of NNIs follows a geometric distribution: p (1 - p)^k
        for i in 1:numNNI
            nni_LiNC!(obj, no3cycle, nohybridladder, constraints, ftolRel,
                      ftolAbs, xtolRel, xtolAbs, γcache)
        end
        logstr = "changed starting topology by $numNNI attempted NNI move(s)\n"
        writelog_1proc && write(logfile, logstr)
        verbose && print(stdout, logstr)
    end
    logstr = "starting at $(writeTopology(obj.net))\n"
    if writelog_1proc
        write(logfile, logstr)
        flush(logfile)
    end
    verbose && print(stdout, logstr)
    optQ =  nparams(obj.model) > 0
    optRAS = nparams(obj.ratemodel) > 0
    nrejected = 0
    while nrejected < nrejectmax
        nrejected = optimizestructure!(obj, maxmoves, maxhybrid, no3cycle, nohybridladder,
                                  nrejected, nrejectmax, constraints,
                                  ftolRel,ftolAbs, xtolRel,xtolAbs, γcache)
        @debug "after optimizestructure returns, the likelihood is $(obj.loglik), nrejected = $nrejected"
        fit!(obj; optimizeQ=optQ, optimizeRVAS=optRAS, maxeval=20,
             ftolRel=ftolRel, ftolAbs=ftolAbs, xtolRel=xtolRel, xtolAbs=xtolAbs)
        @debug "after fit! runs, the likelihood is $(obj.loglik)"
        for i in Random.shuffle(1:obj.net.numEdges)
            e = obj.net.edge[i]
            optimizelocalBL_LiNC!(obj, e, ftolRel,ftolAbs,xtolRel,xtolAbs)
            e.hybrid || continue
            optimizelocalgammas_LiNC!(obj, e, ftolAbs, γcache)
        end
        @debug "after global BL and gamma optimization, the likelihood is $(obj.loglik)"
        for h in obj.net.hybrid # check for gammas close to zero
            minorhybridedge = getMinorParentEdge(h)
            if minorhybridedge.gamma == 0.0 # delete this edge, updateSSM!
                deletehybridedge!(obj.net, minorhybridedge, false, true) # don't keep nodes; unroot
                updateSSM!(obj, true; constraints=constraints) # renumber = true
                @debug "deleted hybrid edge with γ=0 at hybrid node number $(h.number)"
            end
        end
    end
    return obj
end

"""
    checknetwork_LiNC!(net::HybridNetwork, maxhybrid::Int, no3cycle::Bool,
        nohybridladder::Bool,
        constraints=TopologyConstraint[]::Vector{TopologyConstraint},
        verbose::Bool=false)

Check that `net` is an adequate starting network before phyLiNC:
remove nodes of degree 2 (possibly including the root);
check that `net` meets the topological `constraints`,
has no polytomies (except at species constraints),
and `maxhybrid` of fewer reticulations.
According to user-given options, also check for the absence of
3-cycles and/or hybrid ladders.

```jldoctest
julia> maxhybrid = 3;

julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> preorder!(net) # for correct unzipping in checknetwork_LiNC!

julia> PhyloNetworks.checknetwork_LiNC!(net, maxhybrid, true, true)
HybridNetwork, Rooted Network
8 edges
8 nodes: 4 tips, 1 hybrid nodes, 3 internal tree nodes.
tip labels: A, B, C, D
((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0,D:2.5);

```
"""
function checknetwork_LiNC!(net::HybridNetwork, maxhybrid::Int, no3cycle::Bool,
    nohybridladder::Bool,
    constraints=TopologyConstraint[]::Vector{TopologyConstraint},
    verbose::Bool=false)

    if maxhybrid > 0
        net.numTaxa >= 3 ||
            error("cannot estimate hybridizations in topologies with 3 or fewer tips: $(net.numTaxa) tips here.")
    end
    # checks for polytomies, constraint violations, nodes of degree 2
    checkspeciesnetwork!(net, constraints) ||
        error("The species or clade constraints are not satisfied in the starting network.")
        # checkspeciesnetwork removes nodes of degree 2, need to renumber with updateSSM
    if no3cycle
        !contain3cycles(net, no3cycle) || (verbose && @warn("Options indicate there should
        be no 3-cycles in the returned network, but the input network contains
        one or more 3-cycles (after removing any nodes of degree 2 including the
        root). These 3-cycles have been removed."))
        # if nodes removed, need to renumber with updateSSM
    end
    if nohybridladder
        !hashybridladder(net) || error("Options indicate there should be no
        hybrid ladders in the returned network, but the input network contains
        one or more hybrid ladders.")
    end
    if length(net.hybrid) > maxhybrid
        error("Options indicate a maximum of $(maxhybrid) reticulations, but
        the input network contains $(length(net.hybrid)) hybrid nodes. Please
        increase maxhybrid to $(length(net.hybrid)) or provide an input network
        with $(maxhybrid) or fewer reticulations.")
    end
    return net
end

"""
    optimizestructure!(obj::SSM, maxmoves::Integer, maxhybrid::Integer,
                       no3cycle::Bool, nohybridladder::Bool,
                       nreject::Integer, nrejectmax::Integer,
                       constraints::Vector{TopologyConstraint},
                       ftolRel::Float64, ftolAbs::Float64,
                       xtolRel::Float64, xtolAbs::Float64,
                       γcache::CacheGammaLiNC)

Alternate NNI moves, hybrid additions / deletions, and root changes based on their
respective weights in `moveweights_LiNC` ([0.4, 0.2, 0.2, 0.2]).
Branch lengths and hybrid γs around the NNI focal edge or around the
added / deleted hybrid edge are optimized (roughly) on the proposed
network. Each proposed network is accepted --or not-- if the likelihood
improves (or doesn't decrease much for hybrid deletion).
After adding or removing a hybrid, `obj` is updated, to have correct
displayed `trees` and node/edge numberings: done by [`nni_LiNC!`](@ref),
[`addhybridedgeLiNC!`](@ref) and [`deletehybridedgeLiNC!`](@ref).

Output: number of consecutive rejections so far.

The percent of nni moves, hybrid additions, hybrid deletions, and root changes
to be performed is in `PhyloNetworks.moveweights_LiNC`.

- `maxmoves`: maximum number of moves to be performed, including root changes,
  which don't actually change the semi-directed topology and likelihood.
- `nreject`: number of consecutive rejections (ignoring root changes), prior
  to starting the search (from a prior function call)
- `nrejectmax`: the search stops when there has been this number of moves that
  have been rejected in a row (ignoring root changes)

For a description of other arguments, see [`phyLiNC!`](@ref).

Assumptions:
- `checknetworkbeforeLiNC` and `discrete_corelikelihood!` have been called on
  `obj.net`.
- starting with a network without 2- and 3- cycles
  (checked by `checknetworkbeforeLiNC`)

Note: When removing a hybrid edge, always removes the minor edge.
"""
function optimizestructure!(obj::SSM, maxmoves::Integer, maxhybrid::Integer,
    no3cycle::Bool, nohybridladder::Bool, nreject::Integer, nrejectmax::Integer,
    constraints::Vector{TopologyConstraint},
    ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
    γcache::CacheGammaLiNC)

    nmoves = 0
    moveweights = copy(moveweights_LiNC)
    if maxhybrid == 0 # prevent some moves, to avoid proposing non-admissible moves
        moveweights[2] = 0.0 # addhybrid
        moveweights[3] = 0.0 # deletehybrid
        moveweights ./= sum(moveweights)
    end
    while nmoves < maxmoves && nreject < nrejectmax # both should be true to continue
        currLik = obj.loglik
        movechoice = sample(movelist_LiNC, moveweights_LiNC)
        if movechoice == "nni"
            result = nni_LiNC!(obj, no3cycle,  nohybridladder,
                               constraints, ftolRel, ftolAbs, xtolRel, xtolAbs,
                               γcache)
        elseif movechoice == "addhybrid"
            nh = obj.net.numHybrids
            nh == maxhybrid && continue # skip & don't count towards nmoves if enough hybrids in net
            # alternative: switch "add" and "delete" as appropriate (not chosen)
            #              would decrease the weight of NNIs compared to add/delete hybrid
            nh <= maxhybrid || # nh > maxhybrid should never happen
                error("The network has $nh hybrids: more than the allowed max of $maxhybrid")
            result = addhybridedgeLiNC!(obj, currLik, maxhybrid, no3cycle,
                            nohybridladder, constraints, ftolRel, ftolAbs,
                            xtolRel, xtolAbs, γcache)
        elseif movechoice == "deletehybrid"
            obj.net.numHybrids == 0  && continue # skip & don't count towards nmoves if no hybrid in net
            result = deletehybridedgeLiNC!(obj, currLik, maxhybrid,
                        no3cycle, nohybridladder, constraints,
                        ftolRel, ftolAbs, xtolRel, xtolAbs, γcache)
        else # change root (doesn't affect likelihood)
            result = moveroot!(obj.net, constraints)
        end
        # optimize γ's, delete hybrid edges with γ=0
        optimizeallgammas_LiNC!(obj, ftolAbs, γcache, 1) &&
            updateSSM!(obj, true; constraints=constraints)
        # pick a random edge (internal or external), optimize adjancent lengths
        e = Random.rand(obj.net.edge)
        optimizelocalBL_LiNC!(obj, e, ftolRel,ftolAbs,xtolRel,xtolAbs)
        nmoves += 1
        # next: update the number of consecutive rejections ignoring any root move
        #       (which don't change the semi-directed topology or loglik)
        if !isnothing(result) && movechoice != "root"
            if result # move successful and accepted (loglik increase)
                nreject = 0  # reset
            else # move made, rejected, and undone
                nreject += 1
            end
        end
        @debug "$(movechoice) move was " *
          (isnothing(result) ? "not permissible" : (result ? "accepted" : "rejected and undone")) *
          ", $nmoves total moves, $nreject rejected,\nloglik = $(obj.loglik)"
    end
    return nreject
end

"""
    nni_LiNC!(obj::SSM, no3cycle::Bool, nohybridladder::Bool,
              constraints::Vector{TopologyConstraint},
              ftolRel::Float64, ftolAbs::Float64,
              xtolRel::Float64, xtolAbs::Float64,
              γcache::CacheGammaLiNC)

Loop over possible edges for a nearest-neighbor interchange move until one is
found. Performs move and compares the original and modified likelihoods. If the
modified likelihood is greater than the original by `likAbs`, the move
is accepted.

Return true if move accepted, false if move rejected. Return nothing if there
are no nni moves possible in the network.

For arguments, see [`phyLiNC!`](@ref).

Called by [`optimizestructure!`](@ref), which is called by [`phyLiNC!`](@ref).

Note: an RR move does not change the best likelihood. RR means that there's
a hybrid ladder, so it looks like a hard polytomy at the reticulation after
unzipping. Theoretically, we could avoid the re-optimizing the likelihood
accept the move: just change inheritance values to get same likelihood,
and update the tree priors. *Not* done.
"""
function nni_LiNC!(obj::SSM, no3cycle::Bool, nohybridladder::Bool,
                  constraints::Vector{TopologyConstraint},
                  ftolRel::Float64, ftolAbs::Float64,
                  xtolRel::Float64, xtolAbs::Float64,
                  γcache::CacheGammaLiNC)
    currLik = obj.loglik
    edgefound = false
    remainingedges = collect(1:length(obj.net.edge)) # indices in net.edge
    while !edgefound
        isempty(remainingedges) && return nothing
        i_in_remaining = Random.rand(1:length(remainingedges))
        e1 = obj.net.edge[remainingedges[i_in_remaining]]
        undoinfo = nni!(obj.net,e1,nohybridladder,no3cycle,constraints)
        #fixit: for clades, need to update cladeconstraints, then restore below
        if isnothing(undoinfo)
            deleteat!(remainingedges, i_in_remaining)
            continue # edge not found yet, try again
        end
        edgefound = true
        # save displayed trees, priorltw, BLs, and gammas in case we need to undo move
        # todo: why saved here, and not before nni! (for lengths & gammas),
        #       and before the while loop (for displayed trees and priorltw)
        saveddisplayedtree = obj.displayedtree
        savedpriorltw = copy(obj.priorltw)
        savededges = adjacentedges(e1)
        savedlen = [e.length for e in savededges]
        savedgam = [e.gamma for e in savededges]
        updateSSM!(obj)
        optimizelocalgammas_LiNC!(obj, e1, ftolAbs, γcache)
        # don't delete a hybrid edge with γ=0: to be able to undo the NNI
        optimizelocalBL_LiNC!(obj, e1, ftolRel,ftolAbs,xtolRel,xtolAbs)
        if obj.loglik < currLik
            nni!(undoinfo...) # undo move
            obj.displayedtree = saveddisplayedtree # restore displayed trees and weights
            obj.priorltw = savedpriorltw
            obj.loglik = currLik # restore to loglik before move
            for (i,e) in enumerate(savededges) # restore edge lengths and gammas
                e.length = savedlen[i]
                e.gamma  = savedgam[i]
            end
            return false # false means: move was rejected
        else # keep nni move (note: if gammas now = zero, the edge will be removed after move in optimizestructure)
            return true # move was accepted: # rejections will be reset to zero
        end
    end
    return nothing # no NNI were available: none were DAGs or met constraints
end

"""
    addhybridedgeLiNC!(obj::SSM, currLik::Float64, maxhybrid::Int,
        no3cycle::Bool, nohybridladder::Bool,
        constraints::Vector{TopologyConstraint},
        ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
        γcache::CacheGammaLiNC)

Completes checks, adds hybrid in a random location, updates SSM object, and
optimizes branch lengths and gammas locally as part of PhyLiNC optimization.

Return true if accepted add hybrid move. If move not accepted, return false.
If cannot add a hybrid, return nothing.

For arguments, see [`phyLiNC!`](@ref).

Assumptions:
- called by [`optimizestructure!`](@ref)
"""
function addhybridedgeLiNC!(obj::SSM, currLik::Float64, maxhybrid::Int,
    no3cycle::Bool, nohybridladder::Bool,
    constraints::Vector{TopologyConstraint}, ftolRel::Float64,
    ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
    γcache::CacheGammaLiNC)
    #= save displayed trees, priorltw, BLs, and gammas in case we need to undo move
      branch lengths are changed by addhybridedge. We can't use
      adjacentedges because newhybridedge is chosen randomly, so we
      save all branch lengths and gammas here =#
    saveddisplayedtree = obj.displayedtree
    savedpriorltw = copy(obj.priorltw)
    savedlen = [e.length for e in obj.net.edge]
    savedgam = [e.gamma for e in obj.net.edge]
    result = addhybridedge!(obj.net, nohybridladder, no3cycle, constraints;fixroot=true)
    # fixroot=true: to restore edge2 if need be, with deletehybridedge!
    if isnothing(result)
        return nothing
    end
    newhybridnode, newhybridedge = result
    # unzip only at new node and its child edge
    unzipat_canonical!(newhybridnode, getChildEdge(newhybridnode))
    updateSSM!(obj)
    optimizelocalgammas_LiNC!(obj, newhybridedge, ftolRel, γcache)
    if newhybridedge.gamma == 0.0
        deletehybridedge!(obj.net, newhybridedge, false, true) # don't keep nodes; unroot
        obj.displayedtree = saveddisplayedtree # restore displayed trees and tree weights
        obj.priorltw = savedpriorltw
        obj.loglik = currLik # restore to loglik before move
        # todo: why is this needed? optimizelocalgammas_LiNC! leaves obj.loglik up to date, correct?
        for (i,e) in enumerate(obj.net.edge) # restore
            e.gamma = savedgam[i]
        end
        return false
    elseif newhybridedge.gamma == 1.0 # ≃ subtree prune and regraft (SPR) move
        # loglik must be better, yet without any real new reticulation: accept
        deletehybridedge!(obj.net, getMinorParentEdge(newhybridnode), false, true)
        updateSSM!(obj, true; constraints=constraints)
        optimizelocalBL_LiNC!(obj, obj.net.edge[end], ftolRel,ftolAbs,xtolRel,xtolAbs)
        # optimize BL around... which edge? no guarantee that the last edge is the hybrid edge that had gamma=1: it was merged with the edge below.
        #todo: delete call to optimizelocalBL_LiNC! above? would make it symmetric with case gamma=0.
        @debug "addhybrid resulted in SPR move: new hybrid edge had γ=1.0, its partner was deleted"
        return true
    end
    optimizelocalBL_LiNC!(obj, newhybridedge, ftolRel, ftolAbs, xtolRel, xtolAbs)
    if obj.loglik - currLik < likAbsAddHybLiNC # improvement too small or negative: undo
        deletehybridedge!(obj.net, newhybridedge, false, true) # don't keep nodes; unroot
        #= rezip not needed because likelihood unchanged.
        If hybrid edge addition used hybridpartnernew = false, the root was
        changed and direction of edge2 was changed: but that doesn't affect
        the likelihood. We don't restore the same rooted network, but we do
        restore the same semi-directed version. The former displayed trees are fine
        to use, even though they don't use the same root as that in the restored rooted network.
        Further, fixroot=true so hybridpartnernew = true.
        =#
        obj.displayedtree = saveddisplayedtree # restore original displayed trees and weights
        obj.priorltw = savedpriorltw
        obj.loglik = currLik # restore to loglik before move
        for (i,e) in enumerate(obj.net.edge) # restore edges and length
            e.length = savedlen[i]
            e.gamma = savedgam[i]
        end
        return false
    else
        return true
    end
end
"""
    deletehybridedgeLiNC!(obj::SSM, currLik::Float64, maxhybrid::Int,
        no3cycle::Bool, nohybridladder::Bool,
        constraints::Vector{TopologyConstraint},
        ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
        γcache::CacheGammaLiNC)

Deletes a random hybrid edge and updates SSM object as part of
PhyLiNC optimization.

Return true if accepted delete hybrid move. If move not accepted, return false.

For a description of arguments, see [`phyLiNC!`](@ref).

Assumptions:
- called by [`optimizestructure!`](@ref) which does some checks.
"""
function deletehybridedgeLiNC!(obj::SSM, currLik::Float64, maxhybrid::Int,
        no3cycle::Bool, nohybridladder::Bool,
        constraints::Vector{TopologyConstraint}, ftolRel::Float64,
        ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
        γcache::CacheGammaLiNC)

    nh = length(obj.net.hybrid)
    hybridnode = obj.net.hybrid[Random.rand(1:nh)]
    minorhybridedge = getMinorParentEdge(hybridnode)
    #= if type-3 constraints: check that proposed deletion meets constraints
      the constraint's stem edge must be a tree edge -> not disrupted
      species (type 1) or clade type-2 constraints: no problem, because
      the major tree displaying the constraint has major edges only
      type-3 constraints: problem if the tree displaying the constraint
      has the minor hybrid edge to be deleted. Not implemented.
        hybindices = Random.shuffle(1:nh)
        edgenotfound = true
        for hi in hybindices
            hybridnode = obj.net.hybrid[hi]
            minorhybridedge = getMinorParentEdge(hybridnode)
            # edgenotfound = whether any type-3 constraints is violated
            edgenotfound || break
        end
        if edgenotfound # tried all hybrids
            return nothing
        end
    =#
    # set γ to 0 to delete the hybrid edge: makes it easy to undo.
    # update minor edge, and prior log tree weights in obj
    γ0 = minorhybridedge.gamma
    setGamma!(minorhybridedge, 0.0)
    l1mγ = log(1.0-γ0)
    nt, hase = updatecache_hase!(γcache, obj, minorhybridedge.number,
                                getMajorParentEdge(hybridnode).number)
    for it in 1:nt
        @inbounds h = hase[it]
        ismissing(h) && continue # tree has unchanged weight: skip below
        if h  obj.priorltw[it]  = -Inf
        else  obj.priorltw[it] -= l1mγ
        end
    end
    #= warning: after setting γ=0, the 2 other edges connected to the minor
              parent are non-identifiable -> do not optimize them.
     if hybrid ladder (minor parent is hybrid node): not identifiable at all
     otherwise: only their sum is identifiable, but not individual lengths
     so: optimize length of major hybrid edge only =#
    majhyb = getMajorParentEdge(hybridnode)
    len0 = majhyb.length # to restore later if deletion rejected
    optimizeBL_LiNC!(obj, [majhyb], ftolRel,ftolAbs,xtolRel,xtolAbs)
    # don't optimize gammas: because we want to constrain one of them to 0.0
    if obj.loglik - currLik > likAbsDelHybLiNC # -0.1: loglik can decrease for parsimony
        deletehybridedge!(obj.net, minorhybridedge, false, true) # don't keep nodes; unroot
        updateSSM!(obj, true; constraints=constraints) # obj.loglik updated by optimizeBL_LiNC above
        return true
    else # keep hybrid
        majhyb.length = len0
        setGamma!(minorhybridedge, γ0)
        updateSSM_priorltw!(obj) # displayedtree will be correct here
        obj.loglik = currLik # restore to likelihood before move
        return false
    end
end

"""
    updateSSM!(obj::SSM, renumber=false::Bool;
               constraints=TopologyConstraint[]::Vector{TopologyConstraint})

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

Warning:
Does not update the likelihood.

```jldoctest
julia> maxhybrid = 3;

julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple.aln"));

julia> obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69, maxhybrid);

julia> PhyloNetworks.checknetwork_LiNC!(obj.net, maxhybrid, true, true);

julia> using Random; Random.seed!(432);

julia> PhyloNetworks.addhybridedge!(obj.net, obj.net.edge[8], obj.net.edge[1], true, 0.0, 0.4);

julia> writeTopology(obj.net)
"(((B:1.0)#H1:0.1::0.9,(A:1.0)#H2:1.0::0.6):1.5,(C:0.6,#H1:1.0::0.1):1.0,(D:1.25,#H2:0.0::0.4):1.25);"

julia> length(obj.displayedtree) # still as if 1 single reticulation
2

julia> PhyloNetworks.updateSSM!(obj);

julia> length(obj.displayedtree) # now correct 4 displayed treess: 2 reticulations
4
```
"""
function updateSSM!(obj::SSM, renumber=false::Bool;
                   constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    if renumber # traits are in leaf.number order, so leaf nodes not reordered
        resetNodeNumbers!(obj.net; checkPreorder=false, type=:internalonly)
            # preorder not used when type = internalonly
        resetEdgeNumbers!(obj.net, false) # verbose=false
        updateconstraintfields!(constraints, obj.net)
    end
    # extract displayed trees
    obj.displayedtree = displayedTrees(obj.net, 0.0; nofuse=true)
    for tree in obj.displayedtree
        preorder!(tree) # no need to call directEdges! before: already done on net
    end
    # log tree weights: sum log(γ) over edges, for each displayed tree
    updateSSM_priorltw!(obj)
    @debug begin
        all(!ismissing, obj.priorltw) ? "" :
        "one or more inheritance γ's are missing or negative. fix using setGamma!(network, edge)"
    end
    return obj
end

# Warning: requires displayed trees are correct.
         # If they are not, call updatedisplayedtrees!() below
function updateSSM_priorltw!(obj::SSM)
    pltw = obj.priorltw
    resize!(pltw, length(obj.displayedtree))
    for (i,t) in enumerate(obj.displayedtree)
        @inbounds pltw[i] = inheritanceWeight(t)
    end
    return nothing
end

# updates displayed trees when one gamma changes
function updatedisplayedtrees!(trees::Vector, enum::Int, pnum::Int, gammae::Float64, hase::Vector)
    gammap = 1.0 - gammae
    for (i,t) in enumerate(trees)
        h = hase[i]
        ismissing(h) && continue # tree doesn't either e nor partner
        if h
            ei = findfirst(e -> e.number == enum, t.edge)
        else
            ei = findfirst(e -> e.number == pnum, t.edge)
        end
        t.edge[ei].gamma = (h ? gammae : gammap)
    end
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
        forceMinorLength0=true, ultrametric=false)
        # force minor length to 0 to avoid non-identifiability at zippers
        # works well if the true (or "the" best-fit) length of the minor parent
        # edge is less than the true length of the parent edge.
        # (zips all the way up instead of unzipping all the way down, as we do when the child edge = 0)
    for e in net.edge # improve identifiability by not allowing very small lengths
        if e.length < 1.0e-10
            e.length = 0.0001
        end
    end
    if unzip
        unzip_canonical!(net)
    end
    return net
end


"""
    optimizelocalBL_LiNC!(obj::SSM, edge::Edge,
                          ftolRel::Float64, ftolAbs::Float64,
                          xtolRel::Float64, xtolAbs::Float64)

Optimize branch lengths in `net` locally around `edge`. Update all edges that
share a node with `edge` (including itself).
Constrains branch lengths to zero below hybrid nodes.
Return vector of updated `edges`.

Used after `nni!` or `addhybridedge!` moves to update local branch lengths.

For other arguments, see [`phyLiNC!`](@ref).

```jldoctest
julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple.aln"));

julia> obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69);

julia> obj.net.edge[4].length
1.5

julia> using Random; Random.seed!(1234);

julia> PhyloNetworks.optimizelocalBL_LiNC!(obj, obj.net.edge[4], 1e-6,1e-6,1e-2,1e-2);

julia> obj.net.edge[4].length
2.3780835440902453e-20

julia> writeTopology(obj.net; round=true)
"(((A:0.22,(B:1.0)#H1:0.0::0.9):0.0,(C:0.6,#H1:1.0::0.1):0.0):0.0,D:2.0);"

```
"""
function optimizelocalBL_LiNC!(obj::SSM, edge::Edge,
                               ftolRel::Float64, ftolAbs::Float64,
                               xtolRel::Float64, xtolAbs::Float64)
    edges = adjacentedges(edge)
    optimizeBL_LiNC!(obj, edges, ftolRel, ftolAbs, xtolRel, xtolAbs)
    return edges # return all edges from adjacentedges(), even those constrained to 0.0
end

"""
    optimizeBL_LiNC!(obj::SSM, edges::Vector{Edge},
                ftolRel::Float64, ftolAbs::Float64,
                xtolRel::Float64, xtolAbs::Float64,
                maxeval=1000::Int)

Optimize branch lengths for edges in vector `edges`.
Constrains branch lengths to zero below hybrid nodes.

Return vector of updated `edges` Warning: This vector does not include edges
constrained to remain at zero.

For a description of arguments, see [`phyLiNC!`](@ref).

Assumption: None of the branch length are negative.

Warning: always pass a shallow copy of edges using copy().
"""
function optimizeBL_LiNC!(obj::SSM, edges::Vector{Edge},
    ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
    maxeval=1000::Int)
    startlik = discrete_corelikelihood!(obj)
    if !isempty(obj.net.hybrid)
        for i in length(edges):-1:1 # remove constrained edges from edges
            @inbounds e = edges[i]
            if getParent(e).hybrid
                deleteat!(edges, i) # shallow copy!!
            end
        end
    end
    startingvalues = getlengths(edges)
    function loglikfunBL(lengths::Vector{Float64}, grad::Vector{Float64})
        setlengths!(edges, lengths) # set lengths in order of vector `edges`
        res = discrete_corelikelihood!(obj)
        isempty(grad) || error("gradient not implemented")
        return res
    end
    # set-up optimization object for BL parameter
    # no gradient methods:
    # :LN_COBYLA for (non)linear constraits, :LN_BOBYQA for bound constraints
    nparBL = length(edges)
    optBL = NLopt.Opt(:LN_COBYLA, nparBL)
    NLopt.ftol_rel!(optBL,ftolRel) # relative criterion
    NLopt.ftol_abs!(optBL,ftolAbs) # absolute criterion
    NLopt.xtol_rel!(optBL,xtolRel)
    NLopt.xtol_abs!(optBL,xtolAbs)
    defaultinstep = NLopt.initial_step(optBL, getlengths(edges))
    replace!(x -> min(x,0.1), defaultinstep) # 0.1 substitutions / site is kinda large
    NLopt.initial_step!(optBL, defaultinstep)
    NLopt.maxeval!(optBL, maxeval) # max number of iterations
    NLopt.lower_bounds!(optBL, zeros(length(edges)))
    NLopt.upper_bounds!(optBL, fill!(Vector{Float64}(undef, length(edges)), 10.0))
        # max_branch_length = 10 used in IQ-TREE v2.0 (see utils/tools.cpp)
    NLopt.max_objective!(optBL, loglikfunBL)
    fmax, xmax, ret = NLopt.optimize(optBL, getlengths(edges)) # get lengths in order of edges vector
    setlengths!(edges, xmax) # set lengths in order of vector `edges`
    obj.loglik = fmax
    @debug "BL: got $(round(fmax, digits=5)) at BL = $(round.(xmax, digits=5)) after $(optBL.numevals) iterations (return code $(ret))"
    if startlik > obj.loglik
        @debug "The starting likelihood was greater than the post-optimization likelihood.
        Branch lengths will be reassigned to their starting values."
        setlengths!(edges, startingvalues)
        obj.loglik = startlik
    end
    return edges
end

"""
    optimizeallgammas_LiNC!(obj::SSM, ftolAbs::Float64,
                            γcache::CacheGammaLiNC, maxeval=1000::Int)

Optimize all γ's in a network, one by one using [`optimizegamma_LiNC!`]
until `maxeval` γ's have been optimized or until the difference in
log-likelihood falls below `ftolAbs`.

At the end: hybrid edges with γ=0 are deleted (if any).

Output: true if reticulations have been deleted, false otherwise.
If true, `updateSSM!` needs to be called afterwards, with constraints if any.
(Constraints are not known here).
"""
function optimizeallgammas_LiNC!(obj::SSM, ftolAbs::Float64,
                                 γcache::CacheGammaLiNC, maxeval::Int)
    hybnodes = obj.net.hybrid
    nh = length(hybnodes)      # also = obj.net.numHybrids
    if nh==0 return false; end # no gammas to optimize
    hybs = [getMinorParentEdge(h) for h in hybnodes]
    discrete_corelikelihood!(obj) # prerequisite for optimizegamma_LiNC!
    nevals = 0
    ll = obj.loglik
    llnew = +Inf; lldiff = +Inf
    while nevals < maxeval && lldiff > ftolAbs
        for he in hybs
            llnew = optimizegamma_LiNC!(obj, he, ftolAbs, γcache)
        end
        lldiff = llnew - ll
        ll = llnew
        nevals += nh
    end
    reduced = false
    hi = nh
    while hi > 0
        he = getMinorParentEdge(hybnodes[hi])
        if he.gamma == 0.0
            deletehybridedge!(obj.net, he, false, true) # don't keep nodes; unroot
            reduced = true
            nh = length(hybnodes) # normally nh-1, but could be less: deleting
            # one hybrid may delete others indirectly, e.g. if hybrid ladder
            # or if generation of a 2-cycle
        end
        hi -= 1
        if hi>nh hi=nh; end
    end
    # todo: check & delete 3-cycles
    return reduced
end

"""
    optimizelocalgammas_LiNC!(obj::SSM, edge::Edge,
                              ftolAbs::Float64, γcache::CacheGammaLiNC,
                              maxeval=1000::Int)

Optimize γ's in `net` locally around `edge`. Update all edges adjacent to
`edge` (including itself), one by one using [`optimizegamma_LiNC!`]
until `maxeval` γ's have been optimized or until the difference in
log-likelihood falls below `ftolAbs`.
`nothing` is returned.

Used after `nni!` or `addhybridedge!` moves to update local gammas.

Assumptions:
- correct `isChild1` field for `edge` and for hybrid edges
- no in-coming polytomy: a node has 0, 1 or 2 parents, no more

```jldoctest
julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple.aln"));

julia> obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69);

julia> obj.net.edge[3].gamma
0.9

julia> using Random; Random.seed!(1234);

julia> PhyloNetworks.optimizelocalgammas_LiNC!(obj, obj.net.edge[3], 1e-3, PhyloNetworks.CacheGammaLiNC(obj));

julia> obj.net.edge[3].gamma
0.0
````
"""
function optimizelocalgammas_LiNC!(obj::SSM, edge::Edge,
                                   ftolAbs::Float64, γcache::CacheGammaLiNC,
                                   maxeval=1000::Int)
    # get edges that are: hybrid and adjacent to edge
    neighborhybs = Edge[]
    for e in edge.node[1].edge
        e.hybrid || continue # below: e is a hybrid edge
        push!(neighborhybs, e)
    end
    for e in edge.node[2].edge
        e.hybrid || continue
        e === edge && continue
        push!(neighborhybs, e)
    end
    # next: keep minor edges, and replace major hybrid edges
    #       by their minor partner if not already in the list
    for i in length(neighborhybs):-1:1
        e = neighborhybs[i]
        e.isMajor || continue # skip below for minor edges
        p = getPartner(e)     # minor partner
        j = findfirst(x -> x===p, neighborhybs)
        if isnothing(j)
            neighborhybs[i] = p # replace major e by its minor partner
        else
            deleteat!(neighborhybs,i) # delete e
        end
    end
    nh = length(neighborhybs)
    if nh==0
        return nothing
    end
    discrete_corelikelihood!(obj) # prerequisite for optimizegamma_LiNC!
    nevals = 0
    ll = obj.loglik
    llnew = +Inf; lldiff = +Inf
    while nevals < maxeval && lldiff > ftolAbs
        for he in neighborhybs
            llnew = optimizegamma_LiNC!(obj, he, ftolAbs, γcache)
        end
        lldiff = llnew - ll
        ll = llnew
        nevals += nh
    end
    return nothing
end

"""
    optimizegamma_LiNC!(obj::SSM, focusedge::Edge,
                    ftolAbs::Float64, cache::CacheGammaLiNC, maxNR=10::Int)

Optimize γ on a single hybrid `edge` using the Newton-Raphson method
(the log-likelihood is concave).
The new log-likelihood is returned after
updating `obj` and its fields with the new γ.

The search stops if the absolute difference in log-likelihood
between two consecutive iterations is below `ftolAbs`,
or after `maxNR` Newton-Raphson iterations.

Warnings:
- no check that `edge` is hybrid
- `obj._loglikcache` and displayed trees (etc.) are assumed to be up-to-date,
  and is updated alongside the new γ.

Used by [`optimizelocalgammas_LiNC!`](@ref) and
[`optimizeallgammas_LiNC!`](@ref).
"""
function optimizegamma_LiNC!(obj::SSM, focusedge::Edge,
        ftolAbs::Float64, cache::CacheGammaLiNC, maxNR=10::Int)

    ## step 1: prepare vectors constant during the search
    edgenum = focusedge.number
    partner = getPartner(focusedge)
    partnernum = partner.number
    clike = cache.clike # conditional likelihood under focus edge
    clikp = cache.clikp # conditional likelihood under partner edge
    ulik  = obj._sitecache  # unconditional likelihood and more
    fill!(clike, 0.0); fill!(clikp, 0.0)
    nt, hase = updatecache_hase!(cache, obj, edgenum, partnernum)
    γ0 = focusedge.gamma
    if γ0<1e-7 # then prior weight and loglikcachetoo small (-Inf if γ0=0)
        @debug "γ0 too small ($γ0): was changed to 1e-7 prior to optimization"
        γ0 = 1e-7
        setGamma!(focusedge, γ0)
        updatedisplayedtrees!(obj.displayedtree, edgenum, partnernum, γ0, hase)
        updateSSM_priorltw!(obj)
        discrete_corelikelihood!(obj) # to update obj._loglikcache
    elseif γ0>0.9999999
        @debug "γ0 too large ($γ0): was changed to 1 - 1e-7 prior to optimization"
        γ0 = 0.9999999
        setGamma!(focusedge, γ0)
        updatedisplayedtrees!(obj.displayedtree, edgenum, partnernum, γ0, hase)
        updateSSM_priorltw!(obj)
        discrete_corelikelihood!(obj) # to update obj._loglikcache
    end
    # obj._loglikcache[tree,rate,site] = log P(site | tree,rate) + log gamma(tree)
    cadjust = maximum(view(obj._loglikcache, 1:nt,:,:)) # to avoid underflows
    nr = length(obj.ratemodel.ratemultiplier)
    for it in 1:nt # sum likelihood over all displayed trees
        @inbounds h = hase[it]
        ismissing(h) && continue # skip below if tree doesn't have e or partner
        for ir in 1:nr # sum over all rate categories
            if h
                clike .+= exp.(obj._loglikcache[it,ir,:] .- cadjust)
            else
                clikp .+= exp.(obj._loglikcache[it,ir,:] .- cadjust)
            end
        end
    end
    ## step 2: Newton-Raphson
    clike ./= γ0
    clikp ./= 1.0 - γ0
    adjustment = obj.totalsiteweight * (log(nr) - cadjust)
    ll = obj.loglik + adjustment
    noweights = obj.siteweight === nothing
    # evaluate if best γ is at the boundary: 0 or 1
    ulik .= (clike .- clikp) ./ clikp # at γ=0, get derivative of loglik
    llg0 = (noweights ? sum(ulik) : sum(obj.siteweight .* ulik))
    inside01 = true
    if llg0 < 0
        γ = 0.0
        inside01 = false
        ll = (noweights ? sum(log.(clikp)) : sum(obj.siteweight .* log.(clikp)))
        @debug "before NR, saw that gamma = 0 best, gamma assigned to zero"
    else
        ulik .= (clike .- clikp) ./ clike # at γ=1
        llg1 = (noweights ? sum(ulik) : sum(obj.siteweight .* ulik))
        if llg1 > 0
            γ = 1.0
            inside01 = false
            ll = (noweights ? sum(log.(clike)) : sum(obj.siteweight .* log.(clike)))
            @debug "before NR, saw that gamma = 1 best, gamma assigned to 1"
        end
    end
    if inside01
    # use interpolation to get a good starting point? γ = llg0 / (llg0 - llg1) in [0,1]
        γ = γ0
        ulik .= γ .* clike + (1.0-γ) .* clikp
        # ll: rescaled log likelihood. llg = gradient, llh = hessian below
        for istep in 1:maxNR
            # ll from above is also: sum(obj.siteweight .* log.(ulik)))
            ulik .= (clike .- clikp) ./ ulik # gradient of unconditional lik
            llg = (noweights ?  sum(ulik) :  sum(obj.siteweight .* ulik))
            map!(x -> x^2, ulik, ulik) # now ulik = -hessian contributions
            llh = (noweights ? -sum(ulik) : -sum(obj.siteweight .* ulik))
            cγ = γ - llg/llh # candidate γ: will be new γ if inside (0,1)
            if cγ >= 1.0
                γ = γ/2 + 0.5
            elseif cγ <= 0.0
                γ = γ/2
            else
                γ = cγ
            end
            ulik .= γ .* clike + (1.0-γ) .* clikp
            ll_new = (noweights ? sum(log.(ulik)) : sum(obj.siteweight .* log.(ulik)))
            lldiff = ll_new - ll
            ll = ll_new
            lldiff < ftolAbs && break
        end
    end
    ## step 3: update SSM object with new γ
    focusedge.gamma = γ
    partner.gamma = 1.0 - γ
    newmajor = γ > 0.5
    if newmajor != focusedge.isMajor
        focusedge.isMajor = newmajor
        partner.isMajor = !newmajor
        # fixit: check clade constraints. perhaps at the start:
        # if this was problematic for constraints, restrict the search
        # to interval [0,0.5] or [0.5,1] as appropriate
    end
    ll -= adjustment
    obj.loglik = ll
    lγdiff = log(γ) - log(γ0); l1mγdiff = log(1.0-γ) - log(1.0-γ0)
    for it in 1:nt
        @inbounds h = hase[it]
        ismissing(h) && continue # tree has unchanged weight: skip below
        obj.priorltw[it] += (h ? lγdiff : l1mγdiff)
        obj._loglikcache[it,:,:] .+= (h ? lγdiff : l1mγdiff)
    end
    return ll
end

function updatecache_hase!(cache::CacheGammaLiNC, obj::SSM,
                           edgenum::Int, partnernum::Int)
    nt = length(obj.displayedtree) # could change if non-tree-child net
    hase = cache.hase
    resize!(hase, nt)
    for i in 1:nt
        @inbounds t = obj.displayedtree[i]
        @inbounds hase[i] = missing
        for e in t.edge
            if e.number == edgenum
                @inbounds hase[i] = true
                break
            elseif e.number == partnernum
                @inbounds hase[i] = false
                break
            end
        end
    end
    return nt, hase
end
