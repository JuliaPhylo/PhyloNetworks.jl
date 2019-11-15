"""
    StatisticalSubstitutionModel

Subtype of `StatsBase.StatisticalModel`, to fit discrete data to a model
of trait substitution along a network.
See [`fitdiscrete`](@ref) to fit a trait substitution model to discrete data.
It returns an object of type `StatisticalSubstitutionModel`, to which standard
functions can be applied, like `loglikelihood(object)`, `aic(object)` etc.
"""
mutable struct StatisticalSubstitutionModel <: StatsBase.StatisticalModel
    model::SubstitutionModel
    ratemodel::RateVariationAcrossSites #allows rates to vary according to gamma
    net::HybridNetwork
    """ data: trait[i] for leaf with n.number = i
        type Int: for indices of trait labels in getlabels(model)
        allows missing, but not half-ambiguous and not multi-state"""
    trait::Vector{Vector{Union{Missing, Int}}}
    "number of columns in the data: should be equal to `length(trait[i])` for all i"
    nsites::Int
    "vector of weights, one for each column, 1s if missing"
    siteweight::Union{Missings.Missing, Vector{Float64}}
    """
    log of transition probabilities: where e goes from X to Y
    logtrans[i,j,e,r] = P{Y=j|X=i} where i=start_state, j=end_state, e=edge.number, r = rate (if using RateVariationAcrossSites)
    all displayed trees use same edge numbers and same logtrans as in full network
    """
    loglik::Union{Missings.Missing, Float64}
    """
    log of transition probabilities.
    size: k,k, net.numEdges, r where k=nstates(model)
    """
    logtrans::Array{Float64,4}
    "index of active site, for state inference (reconstruction)"
    activesite::Int
    # type based on extracting displayed trees
    displayedtree::Vector{HybridNetwork}
    """
    prior log tree weight: log product of γ's.
    In fit! priorltw = `inheritanceWeight.(trees)`
    This can be missing if any γs are negative.
    """
    priorltw::Vector{Union{Missing,Float64}}
    "posterior log tree weight: P{tree and data}"
    postltw::Vector{Float64}
    """
    partial likelihoods for active trait, at indices [i, n.number or e.number ,t] :
    - forward likelihood: log P{data below node n in tree t given state i at n}
    - direct likelihood:  log P{data below edge e in tree t given i at parent of e}
    - backward likelihood:log P{data at all non-descendants of n and state i at n}
    sizes: k, net.numNodes or net.numEdges, number of displayed trees
    """
    forwardlik::Array{Float64,3} # size: k, net.numNodes, number of displayed trees
    directlik::Array{Float64,3}  # size: k, net.numEdges, number of displayed trees
    backwardlik::Array{Float64,3}# size: k, net.numNodes, number of displayed trees
    "inner (default) constructor: from model, rate model, network, trait and site weights"
    function StatisticalSubstitutionModel(model::SubstitutionModel,
            ratemodel::RateVariationAcrossSites,
            net::HybridNetwork,
            trait::AbstractVector,
            siteweight=missing::Union{Missings.Missing, Vector{Float64}})

        length(trait) > 0 || error("no trait data!")
        nsites = length(trait[1])
        ismissing(siteweight) || length(siteweight) == nsites ||
            error("siteweight must be of same length as the number of traits")
        # T = eltype(getlabels(model))
        # extract displayed trees
        trees = displayedTrees(net, 0.0; keepNodes=true)
        nnodes = length(net.node)
        for tree in trees
            preorder!(tree) # no need to call directEdges! before: already done on net
            length(tree.nodes_changed) == nnodes ||
                error("displayed tree with too few nodes: $(writeTopology(tree))")
            length(tree.edge) == length(net.edge)-net.numHybrids ||
                error("displayed tree with too few edges: $(writeTopology(tree))")
        end
        ntrees = length(trees)
        # log tree weights: sum log(γ) over edges, for each displayed tree
        priorltw = inheritanceWeight.(trees)
        k = nstates(model)
        # fixit: use SharedArray's below to parallelize things
        logtrans   = zeros(Float64, k,k,length(net.edge), length(ratemodel.ratemultiplier))
        forwardlik = zeros(Float64, k, nnodes,           ntrees)
        directlik  = zeros(Float64, k, length(net.edge), ntrees)
        backwardlik= zeros(Float64, k, nnodes,           ntrees)
        postltw    = Vector{Float64}(undef, ntrees)
        new(deepcopy(model), deepcopy(ratemodel),
            net, trait, nsites, siteweight, missing, # missing log likelihood
            logtrans, 1, trees,
            priorltw, postltw, forwardlik, directlik, backwardlik)
    end
end
const SSM = StatisticalSubstitutionModel

StatsBase.loglikelihood(obj::SSM) = obj.loglik
StatsBase.islinear(SSM) = false
StatsBase.dof(obj::SSM) = nparams(obj.model) + nparams(obj.ratemodel)
function Base.show(io::IO, obj::SSM)
    disp =  "$(typeof(obj)):\n"
    disp *= string(obj.model)
    disp *= "$(obj.nsites) traits, $(length(obj.trait)) species\n"
    if obj.ratemodel.ncat != 1
        disp *= "variable rates across sites ~ discretized gamma with\n alpha=$(obj.ratemodel.alpha)"
        disp *= "\n $(obj.ratemodel.ncat) categories"
        disp *= "\n rate multipliers: $(round.(obj.ratemodel.ratemultiplier, digits=5))\n"
    end
    disp *= "on a network with $(obj.net.numHybrids) reticulations"
    if !ismissing(obj.loglik)
        disp *= "\nlog-likelihood: $(round(obj.loglik, digits=5))"
    end
    print(io, disp)
end
# nobs: nsites * nspecies, minus any missing, but ignores correlation between species
# fixit: extend the StatsBase methods
# coef, coefnames, coeftable, confint,
# deviance (not from loglik in StatsBase), nulldeviance (from intercept only. here?)
# nullloglikelihood (from intercept only. here?)
# score (gradient of loglik)
# informationmatrix: Fisher by default, observed info matrix if expected=false
# stderror, vcov, params
# weights

"""
    fitdiscrete(net, model, tipdata)
    fitdiscrete(net, model, RateVariationAcrossSites, tipdata)
    fitdiscrete(net, model, species, traits)
    fitdiscrete(net, model, RateVariationAcrossSites, species, traits)
    fitdiscrete(net, model, dnadata, dnapatternweights)
    fitdiscrete(net, model, RateVariationAcrossSites, dnadata, dnapatternweights)
    fitdiscrete(net, modSymbol, species, traits)
    fitdiscrete(net, modSymbol, dnadata, dnapatternweights)
    
Calculate the maximum likelihood (ML) score of a network or tree given
one or more discrete characters at the tips. Along each edge, transitions
are modelled with a continous time Markov `model`, whose parameters are
estimated (by maximizing the likelihood). At each hybrid node,
the trait is assumed to be inherited from either of the two immediate
parents according to the parents' average genetic contributions
(inheritance γ). The model ignores incomplete lineage sorting.
The algorithm extracts all trees displayed in the network.

Data can given in one of the following:
- `tipdata`: dictionary taxon => state label, for a single trait.
- `tipdata`: data frame for a single trait, in which case the taxon names
  are to appear in column 1 or in a column named "taxon" or "species", and
  trait *labels* are to appear in column 2 or in a column named "trait".
  Here, trait labels should be as they appear in `getlabels(model)`.
- `species`: vector of strings, and `traits`: DataFrame of traits,
  with rows in the order corresponding to the order of species names.
  Again, trait labels should be as they appear in `getlabels(model)`.
  All traits are assumed to follow the same model, with same parameters.
- 'dnadata': the first part of the output of readfastatodna, 
    a dataframe of BioSequence DNA sequences, with taxon in column 1 and a column for each site.
- 'dnapatternweights': the second part of the output of readfastatodna, 
    an array of weights, one weights for each of the site columns. The length of the weight is equal to nsites.
    If using dnapatternweights, must provide dnadata.
- RateVariationAcrossSites: model for rate variation (optional)

Optional arguments (default):
- `optimizeQ` (true): should model rate parameters be fixed, or should they be optimized?
- `optimizeRVAS` (true): should the model optimize the variable rates across sites?
- `optimizeTOPO` (false): should the model optimize network topology and branch lengths?
- `NLoptMethod` (`:LN_COBYLA`, derivative-free) for the optimization algorithm.
  For other options, see the
  [NLopt](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/).
- tolerance values to control when the optimization is stopped:
  `ftolRel` (1e-12), `ftolAbs` (1e-10) on the likelihood, and
  `xtolRel` (1e-10), `xtolAbs` (1e-10) on the model parameters.
- bounds for the alpha parameter of the Gamma distribution of rates across sites:
  `alphamin=0.05`, `alphamax=500`.
- `verbose` (false): if true, more information is output.

# examples:

```jldoctest fitDiscrete_block
julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> m1 = BinaryTraitSubstitutionModel([0.1, 0.1], ["lo", "hi"]);

julia> using DataFrames

julia> dat = DataFrame(species=["C","A","B","D"], trait=["hi","lo","lo","hi"]);

julia> fit1 = fitdiscrete(net, m1, dat)
PhyloNetworks.StatisticalSubstitutionModel:
Binary Trait Substitution Model:
rate lo→hi α=0.27222
rate hi→lo β=0.34981
1 traits, 4 species
on a network with 1 reticulations
log-likelihood: -2.7277

julia> tips = Dict("A" => "lo", "B" => "lo", "C" => "hi", "D" => "hi");

julia> fit2 = fitdiscrete(net, m1, tips; xtolRel=1e-16, xtolAbs=1e-16, ftolRel=1e-16)
PhyloNetworks.StatisticalSubstitutionModel:
Binary Trait Substitution Model:
rate lo→hi α=0.27222
rate hi→lo β=0.34981
1 traits, 4 species
on a network with 1 reticulations
log-likelihood: -2.7277
```

Note that a copy of the network is stored in the fitted object,
but the internal representation of the network may be different in
`fit1.net` and in the original network `net`:

```jldoctest fitDiscrete_block
julia> net = readTopology("(sp1:3.0,(sp2:2.0,(sp3:1.0,sp4:1.0):1.0):1.0);");

julia> using BioSymbols

julia> tips = Dict("sp1" => BioSymbols.DNA_A, "sp2" => BioSymbols.DNA_A, "sp3" => BioSymbols.DNA_G, "sp4" => BioSymbols.DNA_G);

julia> mJC69 = JC69([0.25], false);

julia> fitJC69 = fitdiscrete(net, mJC69, tips)
PhyloNetworks.StatisticalSubstitutionModel:
Jukes and Cantor 69 Substitution Model,
absolute rate version
off-diagonal rates equal to 0.29234/3.
rate matrix Q:
               A       C       G       T
       A       *  0.0974  0.0974  0.0974
       C  0.0974       *  0.0974  0.0974
       G  0.0974  0.0974       *  0.0974
       T  0.0974  0.0974  0.0974       *
1 traits, 4 species
on a network with 0 reticulations
log-likelihood: -4.99274

julia> rv = RateVariationAcrossSites()
Rate Variation Across Sites using Discretized Gamma Model
alpha: 1.0
categories for Gamma discretization: 4
ratemultiplier: [0.14578, 0.51313, 1.07083, 2.27025]

julia> fitdiscrete(net, mJC69, rv, tips; optimizeQ=false, optimizeRVAS=false)
PhyloNetworks.StatisticalSubstitutionModel:
Jukes and Cantor 69 Substitution Model,
absolute rate version
off-diagonal rates equal to 0.25/3.
rate matrix Q:
               A       C       G       T
       A       *  0.0833  0.0833  0.0833
       C  0.0833       *  0.0833  0.0833
       G  0.0833  0.0833       *  0.0833
       T  0.0833  0.0833  0.0833       *
1 traits, 4 species
variable rates across sites ~ discretized gamma with
 alpha=1.0
 4 categories
 rate multipliers: [0.14578, 0.51313, 1.07083, 2.27025]
on a network with 0 reticulations
log-likelihood: -5.2568
```
#TODO add option to allow users to specify root prior (either equal frequencies or stationary frequencies)
for trait models
"""
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel, #tips::Dict no ratemodel version
    tips::Dict; kwargs...)
    ratemodel = RateVariationAcrossSites(1.0, 1)
    fitdiscrete(net, model, ratemodel, tips; kwargs...)
end

#tips::Dict version with ratemodel
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel, ratemodel::RateVariationAcrossSites, 
    tips::Dict; kwargs...)
        species = String[]
    dat = Vector{Int}[] # indices of trait labels
    for (k,v) in tips
        !ismissing(v) || continue
        push!(species, k)
        vi = findfirst(isequal(v), getlabels(model))
        vi !== nothing || error("trait $v not found in model")
        push!(dat, [vi])
    end
    o, net = check_matchtaxonnames!(species, dat, net) # dat[o] would make a shallow copy only
    StatsBase.fit(StatisticalSubstitutionModel, net, model, ratemodel, view(dat, o); kwargs...)
end

#dat::DataFrame, no rate model version
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel, dat::DataFrame; kwargs...)
    ratemodel = RateVariationAcrossSites(1.0, 1)
    fitdiscrete(net, model, ratemodel, dat; kwargs...)
end

#dat::DataFrame with rate model version
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel, ratemodel::RateVariationAcrossSites,
        dat::DataFrame; kwargs...)
    i = findfirst(isequal(:taxon), DataFrames.names(dat))
    
    if i===nothing i = findfirst(isequal(:species), DataFrames.names(dat)); end
    if i===nothing i=1; end # first column if no column named "taxon" or "species"
    j = findfirst(isequal(:trait), DataFrames.names(dat))
    if j===nothing j=2; end
    if i==j
        error("""expecting taxon names in column 'taxon', or 'species' or column 1,
              and trait values in column 'trait' or column 2.""")
    end
    species = dat[:,i]    # modified in place later
    dat = traitlabels2indices(dat[!,j], model)   # vec of vec, indices
    o, net = check_matchtaxonnames!(species, dat, net)
    StatsBase.fit(StatisticalSubstitutionModel, net, model, ratemodel, view(dat, o); kwargs...)
end

#species, dat version, no ratemodel
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel, species::Array{String}, dat::DataFrame; kwargs...)
    ratemodel = RateVariationAcrossSites(1.0, 1)
    fitdiscrete(net, model, ratemodel, species, dat; kwargs...)
end

#species, dat version with ratemodel
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel, ratemodel::RateVariationAcrossSites,
        species::Array{String}, dat::DataFrame; kwargs...)
    dat2 = traitlabels2indices(dat, model) # vec of vec, indices
    o, net = check_matchtaxonnames!(copy(species), dat2, net)
    StatsBase.fit(StatisticalSubstitutionModel, net, model, ratemodel, view(dat2, o); kwargs...)
end

#wrapper: species, dat version with model symbol
function fitdiscrete(net::HybridNetwork, modSymbol::Symbol,  
    species::Array{String}, dat::DataFrame, rvSymbol=:noRV::Symbol; kwargs...)
    rate = startingrate(net)
    labels = learnlabels(modSymbol, dat)
    if modSymbol == :JC69
        model = JC69([1.0], true) # 1.0 instead of rate because relative version
    elseif modSymbol == :HKY85
        model = HKY85([1.0], # transition/transversion rate ratio
                       empiricalDNAfrequencies(dat, repeat([1.], inner=ncol(dat))),
                       true)
    elseif modSymbol == :ERSM
        model = EqualRatesSubstitutionModel(length(labels), rate, labels);
    elseif modSymbol == :BTSM
        model = BinaryTraitSubstitutionModel([rate, rate], labels)
    elseif modSymbol == :TBTSM
        model = TwoBinaryTraitSubstitutionModel([rate, rate, rate, rate, rate, rate, rate, rate], labels)
    else
        error("model $modSymbol is unknown or not implemented yet")
    end

    if rvSymbol == :RV
        rvas = RateVariationAcrossSites(1.0, 4)
    else
        rvas = RateVariationAcrossSites(1.0, 1)
    end
    fitdiscrete(net, model, rvas, species, dat; kwargs...)
end

#dnadata with dnapatternweights version, no ratemodel
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel, dnadata::DataFrame, 
    dnapatternweights::Array{Float64}; kwargs...)
    ratemodel = RateVariationAcrossSites(1.0, 1)
    fitdiscrete(net, model, ratemodel, dnadata, dnapatternweights; kwargs...)
end

#dnadata with dnapatternweights version with ratemodel
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel, ratemodel::RateVariationAcrossSites,
    dnadata::DataFrame, dnapatternweights::Array{Float64}; kwargs...)
    
    dat2 = traitlabels2indices(dnadata[!,2:end], model)
    o, net = check_matchtaxonnames!(dnadata[:,1], dat2, net)
    kwargs = (:siteweights => dnapatternweights, kwargs...)
    StatsBase.fit(StatisticalSubstitutionModel, net, model, ratemodel, view(dat2, o); 
        kwargs...)
end

#wrapper for dna data
function fitdiscrete(net::HybridNetwork, modSymbol::Symbol, dnadata::DataFrame, 
    dnapatternweights::Array{Float64}, rvSymbol=:noRV::Symbol; kwargs...)
    rate = startingrate(net)
    if modSymbol == :JC69
        model = JC69([1.0], true)  # 1.0 instead of rate because relative version (relative = true)
    elseif modSymbol == :HKY85
        model = HKY85([1.0], # transition/transversion rate ratio
                      empiricalDNAfrequencies(view(dnadata, :, 2:ncol(dnadata)), dnapatternweights),
                      true)
    elseif modSymbol == :ERSM
        model = EqualRatesSubstitutionModel(4, rate, [BioSymbols.DNA_A, BioSymbols.DNA_C, BioSymbols.DNA_G, BioSymbols.DNA_T]);
    elseif modSymbol == :BTSM
        error("Binary Trait Substitution Model supports only two trait states, but dna data has four states.")
    elseif modSymbol == :TBTSM
        error("Two Binary Trait Substitution Model does not support dna data: it supports two sets of potentially correlated two trait states.")
    else
        error("model $modSymbol is unknown or not implemented yet")
    end

    if rvSymbol == :RV
        rvas = RateVariationAcrossSites(1.0, 4)
    else
        rvas = RateVariationAcrossSites(1.0, 1)
    end
    
    fitdiscrete(net, model, rvas, dnadata, dnapatternweights; kwargs...)
end

"""
    fit(StatisticalSubstitutionModel, net, model, traits; kwargs...)
    fit!(StatisticalSubstitutionModel; kwargs...)

Internal function called by [`fitdiscrete`](@ref): with same key word arguments `kwargs`.
But dangerous: `traits` should be a vector of vectors as for [`fitdiscrete`](@ref)
**but** here `traits` need to contain the *indices* of trait values corresponding
to the indices in `getlabels(model)`, and species should appear in `traits` in the
order corresponding to the node numbers in `net`.
See [`traitlabels2indices`](@ref) to convert trait labels to trait indices.
**Warning**: does *not* perform checks. [`fitdiscrete`](@ref) calls this function
after doing checks, preordering nodes in the network, making sure nodes have
consecutive numbers, species are matched between data and network etc.
Like, snaq!(), we estimate the network starting from tree (or network) `astraltree`.
"""
function StatsBase.fit(::Type{SSM}, net::HybridNetwork, model::SubstitutionModel, 
    ratemodel::RateVariationAcrossSites, trait::AbstractVector; kwargs...)

    sw = missing
    if haskey(kwargs, :siteweights)
        sw = kwargs[:siteweights]
        kwargs = filter(p -> p.first != :siteweights, kwargs)
    end
    obj = StatisticalSubstitutionModel(model, ratemodel, net, trait, sw)
    fit!(obj; kwargs...)
end

function fit!(obj::SSM; optimizeQ=true::Bool, optimizeRVAS=true::Bool,verbose=false::Bool,
      NLoptMethod=:LD_MMA::Symbol, ftolRel=fRelBL::Float64, ftolAbs=fAbsBL::Float64,
      xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64,
      alphamin=0.05, alphamax=500)
    all(x -> x >= 0.0, [e.length for e in obj.net.edge]) || error("branch lengths should be >= 0")
    all(x -> x >= 0.0, [e.gamma for e in obj.net.edge]) || error("gammas should be >= 0")
    # if optimizeTOPO
    #     obj.net = startingBL!(obj.net, Array(obj.trait), obj.siteweight) #TODO check if Array() conversion needed
    #     # TODO obj.net = updateTopo!(obj.net, liketolAbs, Nfail, data, hmax,
    #     #     ftolRel, ftolAbs, xtolRel, 
    #     #     xtolAbs, verbose, closeN, NmovO,
    #     #     logfile, writelog)
    # end
    if optimizeQ && nparams(obj.model) <1
        @debug "The Q matrix for this model is fixed, so there are no rate parameters to optimize. optimizeQ will be set to false."
        optimizeQ = false
    end
    if optimizeRVAS && (nparams(obj.ratemodel) == 0)
        @debug "Rate model has one rate category, so there are no parameters to optimize. optimizeRVAS will be set to false."
        optimizeRVAS = false
    end
    if !optimizeQ && !optimizeRVAS
        discrete_corelikelihood!(obj)
        verbose && println("loglik = $(loglikelihood(obj)) under fixed parameters, no optimization")
        return obj
    end
    counter = [0]
    if optimizeQ
        function loglikfun(x::Vector{Float64}, grad::Vector{Float64}) # modifies obj
            counter[1] += 1
            setrates!(obj.model, x)
            res = discrete_corelikelihood!(obj)
            verbose && println("loglik: $res, model rates: $x")
            length(grad) == 0 || error("gradient not implemented")
            return res
        end
        # optimize Q under fixed RVAS parameters
        # set-up optimization object for Q
        NLoptMethod=:LN_COBYLA # no gradient
        # :LN_COBYLA for (non)linear constraits, :LN_BOBYQA for bound constraints
        nparQ = nparams(obj.model)
        optQ = NLopt.Opt(NLoptMethod, nparQ)
        NLopt.ftol_rel!(optQ,ftolRel) # relative criterion
        NLopt.ftol_abs!(optQ,ftolAbs) # absolute criterion
        NLopt.xtol_rel!(optQ,xtolRel)
        NLopt.xtol_abs!(optQ,xtolAbs)
        NLopt.maxeval!(optQ,1000) # max number of iterations
        # NLopt.maxtime!(optQ, t::Real)
        NLopt.lower_bounds!(optQ, zeros(Float64, nparQ))
        # fixit: set upper bound depending on branch lengths in network?
        counter[1] = 0
        NLopt.max_objective!(optQ, loglikfun)
        fmax, xmax, ret = NLopt.optimize(optQ, obj.model.rate) # optimization here!
        verbose && println("got $(round(fmax, digits=5)) at $(round.(xmax, digits=5)) after $(counter[1]) iterations (return code $(ret))")
    end
    if optimizeRVAS
        function loglikfunRVAS(alpha::Vector{Float64}, grad::Vector{Float64}) # modifies obj
            counter[1] += 1
            setalpha!(obj.ratemodel, alpha[1]) # obj.ratemodel.alpha is a scalar, alpha is a vector within here
            res = discrete_corelikelihood!(obj)
            verbose && println("loglik: $res, rate variation model shape parameter alpha: $(alpha[1])")
            length(grad) == 0 || error("gradient not implemented")
            return res
        end
        # set-up optimization object for RVAS parameter
        NLoptMethod=:LN_COBYLA # no gradient
        # :LN_COBYLA for (non)linear constraits, :LN_BOBYQA for bound constraints
        nparRVAS = nparams(obj.ratemodel)
        optRVAS = NLopt.Opt(NLoptMethod, nparRVAS)
        NLopt.ftol_rel!(optRVAS,ftolRel) # relative criterion
        NLopt.ftol_abs!(optRVAS,ftolAbs) # absolute criterion
        NLopt.xtol_rel!(optRVAS,xtolRel)
        NLopt.xtol_abs!(optRVAS,xtolAbs)
        NLopt.maxeval!(optRVAS,1000) # max number of iterations
        # NLopt.maxtime!(optRVAS, t::Real)
        NLopt.lower_bounds!(optRVAS, fill(alphamin, (nparRVAS,)) ) # for 0 as lower bound: zeros(Float64, nparRVAS)
        NLopt.upper_bounds!(optRVAS, fill(alphamax, (nparRVAS,)) ) # delete to remove upper bound
        counter[1] = 0
        NLopt.max_objective!(optRVAS, loglikfunRVAS)
        fmax, xmax, ret = NLopt.optimize(optRVAS, [obj.ratemodel.alpha]) # optimization here!
        verbose && println("RVAS: got $(round(fmax, digits=5)) at $(round.(xmax, digits=5)) after $(counter[1]) iterations (return code $(ret))")
    end
    if optimizeQ && optimizeRVAS
        # optimize Q under fixed RVAS parameters: a second time
        counter[1] = 0
        fmax, xmax, ret = NLopt.optimize(optQ, obj.model.rate)
        verbose && println("got $(round(fmax, digits=5)) at $(round.(xmax, digits=5)) after $(counter[1]) iterations (return code $(ret))")
        # optimize RVAS under fixed Q: a second time
        counter[1] = 0
        fmax, xmax, ret = NLopt.optimize(optRVAS, [obj.ratemodel.alpha])
        verbose && println("RVAS: got $(round(fmax, digits=5)) at $(round.(xmax, digits=5)) after $(counter[1]) iterations (return code $(ret))")
    end
    return obj
end


"""
    discrete_corelikelihood!(obj::StatisticalSubstitutionModel; whichtrait=:all)
    discrete_corelikelihood_tree!(obj, t::Integer, traitrange::AbstractArray)

Calculate the likelihood and update `obj.loglik` for discrete characters on a network
(or on a single tree: `t`th tree displayed in the network, for the second form).
Update forward and direct partial likelihoods while doing so.
The algorithm extracts all displayed trees and weights the likelihood under all these trees.
""" 
function discrete_corelikelihood!(obj::SSM; whichtrait=:all::Union{Symbol,Integer})
    if whichtrait == :all
        traitrange = 1:obj.nsites
    elseif isinteger(whichtrait) && whichtrait > 0 && whichtrait <= obj.nsites
        obj.activesite = whichtrait
        traitrange = range(whichtrait, length=1)
    elseif whichtrait == :active
        traitrange = range(obj.activesite, length=1)
    else
        error("'whichtrait' should be :all or :active or an integer in the correct range")
    end
    startingP = P(obj.model, 1.0) #sets t = 1 for starting P for efficency
    for edge in obj.net.edge # update logtrans: same for all displayed trees, all traits
        for i = 1:length(obj.ratemodel.ratemultiplier)
            obj.logtrans[:,:,edge.number, i] = log.(P!(startingP, obj.model, edge.length*obj.ratemodel.ratemultiplier[i])) # element-wise
        end
    end
    for t in 1:length(obj.displayedtree) # calculate P{data | tree t} & store in obj.postltw[t]
        discrete_corelikelihood_tree!(obj, t, traitrange)
    end
    # fixit: paralellize with
    # ll = pmap(t -> discrete_corelikelihood_tree!(obj,t), 1:ntrees)
    obj.postltw .+= obj.priorltw # P{tree t and data} .+= not += to re-use memory
    res = StatsFuns.logsumexp(obj.postltw)
    obj.loglik = res
    obj.postltw .-= res # now P{tree t | data}
    # fixit: write a function to get these posterior probabilities (just take exp.)
    return res
end

@doc (@doc discrete_corelikelihood!) discrete_corelikelihood_tree!
function discrete_corelikelihood_tree!(obj::SSM, t::Integer, traitrange::AbstractArray)
    tree = obj.displayedtree[t]
    # @info "tree: $(writeTopology(tree))"
    forwardlik = view(obj.forwardlik, :,:,t)
    directlik  = view(obj.directlik,  :,:,t)
    k = nstates(obj.model)   # also = size(logtrans,1) if not RateVariationAcrossSites
    nr = length(obj.ratemodel.ratemultiplier)
    fullloglik = 0.0         # full: for all characters
    for ci in traitrange     # ci = character index
        currentloglik = 0.0  # current character
        for iratemultiplier in 1:nr
            fill!(forwardlik, 0.0) # re-initialize for each trait, each iteration
            fill!(directlik,  0.0)
            for ni in reverse(1:length(tree.nodes_changed)) # post-order
                n = tree.nodes_changed[ni]
                nnum = n.number # same n.number across trees for a given node
                if n.leaf # need forwardlik initialized at 0: keep at 0 = log(1) if no data
                    state = obj.trait[nnum][ci] # here: data assumed in a row n.number
                    if !ismissing(state)
                        for i in 1:k
                        forwardlik[i,nnum] = -Inf64 # log(0) = -Inf if i != observed state
                        end
                        forwardlik[state, nnum] = 0.
                    end
                else # forward likelihood = product of direct likelihood over all children edges
                    for e in n.edge
                        n == getParent(e) || continue # to next edge if n is not parent of e
                        forwardlik[:,nnum] .+= view(directlik, :,e.number)
                    end
                end
                if ni==1 # root is first index in nodes changed
                    if typeof(obj.model) == NASM 
                        logprior = log.(stationary(obj.model))
                    else #other trait models
                        logprior = [-log(k) for i in 1:k] # uniform prior at root
                    end
                    loglik = logsumexp(logprior + view(forwardlik, :,nnum)) # log P{data for ci | tree t}
                    if iratemultiplier == 1
                        currentloglik = loglik
                    else # next rate multiplier:
                        currentloglik = logaddexp(currentloglik, loglik) 
                    end
                    break # out of loop over nodes
                end
                # if we keep going, n is not the root
                # calculate direct likelihood on the parent edge of n
                for e in n.edge
                    if n == getChild(e)
                        lt = view(obj.logtrans, :,:,e.number, iratemultiplier)
                        for i in 1:k # state at parent node
                            directlik[i,e.number] = logsumexp(view(lt,i,:) + view(forwardlik,:,nnum))
                        end
                        break # we visited the parent edge: break out of for loop
                    end
                end #loop over edges
            end # of loop over nodes
        end # of loop over rate multipliers
        if !ismissing(obj.siteweight) #if dna data with dna site pattern weights, multiplied here
            currentloglik *= obj.siteweight[ci]
        end
        # add loglik of character ci to full loglik:
        fullloglik += currentloglik #warning: this loglik missing one term - log(length(obj.ratemodel.ratemultiplier)), corrected below
    end #of loop over traits
    obj.postltw[t] = fullloglik - log(nr)*length(traitrange) #logL divided by (#of rates)(# of chars)
    return fullloglik - log(nr)*length(traitrange)
end

"""
    traitlabels2indices(data, model::SubstitutionModel)

Check that the character states in `data` are compatible with (i.e. subset of)
the trait labels in `model`. All columns are used.
`data` can be a DataFrame or a Matrix (multiple traits), or a Vector (one trait).
Return a vector of vectors (one per species) with integer entries,
where each state (label) is replaced by its index in `model`. 
For DNA data, any ambiguous site is treated as missing.
"""
function traitlabels2indices(data::AbstractVector, model::SubstitutionModel)
    A = Vector{Vector{Union{Missings.Missing,Int}}}(undef, 0) # indices of trait labels
    labs = getlabels(model)
    for l in data
        vi = missing
        if !ismissing(l)
            vi = findfirst(isequal(l), getlabels(model)) # value index in model labels
            vi !== nothing || error("trait $l not found in model")
        end
        push!(A, [vi])
    end
    return A
end

function traitlabels2indices(data::Union{AbstractMatrix,AbstractDataFrame},
                             model::SubstitutionModel)
    A = Vector{Vector{Union{Missings.Missing,Int}}}(undef, 0) # indices of trait labels
    labs = getlabels(model)
    isDNA = typeof(labs) == Array{DNA,1}
    ntraits = size(data,2) #starting with spcies column
    for i in 1:size(data,1) # go row by row (which is species by species)
        V = Vector{Union{Missings.Missing,Int}}(undef, ntraits)
        for j in 1:ntraits
            vi = missing # value index
            @inbounds l = data[i,j] # value label
            if !isDNA && ismissing(l)
                vi = missing
            elseif isDNA #else if DNA
                if typeof(l) == String #takes string and converts to a Char so that we can convert to DNA
                    l = Vector{Char}(l)[1]
                end
                l = convert(DNA, l)
            end
            if !ismissing(l)
                vi = findfirst(isequal(l), labs) 
                if vi == nothing
                    #FIXIT ideally, replace isambiguous with isgap and handle ambiguous DNA types
                    #@show BioSymbols.isambiguous(l) #TODO this is false
                    if isDNA #&& BioSymbols.isambiguous(l)
                        vi = missing 
                    else
                        error("trait $l not found in model")
                    end
                end
            end
            V[j] = vi
        end
        push!(A, V)
    end
    return A
end
#TODO for data with multiple traits, add test for 2 traits to test above function
"""
    check_matchtaxonnames!(species, data, net)

Modify `species` and `dat` by removing the species (rows) absent from the network.
Return a new network (`net` is *not* modified) with tips matching those in species:
if some species in `net` have no data, these species are pruned from the network.
The network also has its node names reset, such that leaves have nodes have
consecutive numbers starting at 1, with leaves first.
Used by [`fitdiscrete`](@ref) to build a new [`StatisticalSubstitutionModel`](@ref).
"""
function check_matchtaxonnames!(species::AbstractVector, dat::AbstractVector, net::HybridNetwork)
    # 1. basic checks for dimensions and types
    eltt = eltype(dat)
    @assert eltt <: AbstractVector "traits should be a vector of vectors"
    @assert nonmissingtype(eltype(eltt)) <: Integer "traits should be integers (label indices)"
    @assert !isempty(dat) "empty data vector!"
    ntraits = length(dat[1])
    for d in dat
        @assert length(d)==ntraits "all species should have the same number of traits"
    end
    @assert length(dat) == length(species) "need as many species as rows in trait data"
    # 2. match taxon labels between data and network
    netlab = tipLabels(net)
    ind2notinnet = findall(x -> x ∉ netlab, species) # species not in network
    deleteat!(species, ind2notinnet)
    deleteat!(dat,     ind2notinnet)
    nvalues = [sum(.!ismissing.(d)) for d in dat] # species with completely missing data
    indmissing = findall(x -> x==0, nvalues)
    deleteat!(species, indmissing)
    deleteat!(dat,     indmissing)
    indnotindat = findall(x -> x ∉ species, netlab) # species not in data
    net = deepcopy(net)
    if !isempty(indnotindat)
        @warn "the network contains taxa with no data: those will be pruned"
        for i in indnotindat
            deleteleaf!(net, netlab[i])
        end
    end
    # 3. calculate order of rows to have species with node.number i on ith row
    resetNodeNumbers!(net; checkPreorder=true) # tip species: now with numbers 1:n
    resetEdgeNumbers!(net) # to use edge as indices: 1:numEdges
    netlab = [n.name for n in sort(net.leaf, by = x -> x.number)]
    nspecies = length(netlab)
    o = Vector{Int}(undef, nspecies)
    for i in 1:nspecies
        @inbounds o[i] = something(findfirst(isequal(netlab[i]), species),0)
    end
    count(!iszero, o) == nspecies || # number of non-zeros should be total size
        error("weird: even after pruning, species in network have no data")
    return (o,net)
end

"""
    ancestralStateReconstruction(obj::SSM, trait::Integer)
    ancestralStateReconstruction(obj::SSM)

Estimate the marginal probability of ancestral states for discrete character
number `trait`, or for the active trait if `trait` is unspecified: `obj.activesite`.
The parameters of the [`StatisticalSubstitutionModel`](@ref) object `obj`
must first be fitted using [`fitdiscrete`](@ref), and ancestral state reconstruction
is conditional on the estimated parameters. If these parameters were estimated
using all traits, they are used as is to do ancestral state reconstruction of the
particular `trait` of interest.
number `trait`, or for the active trait if `trait` is unspecified: `obj.activesite`.
**output**: data frame with a first column for the node numbers, a second column for
the node labels, and a column for each possible state: the entries in these columns
give the marginal probability that a given node has a given state.

warnings:
- node numbers and node labels refer to those in `obj.net`, which might
  have a different internal representation of nodes than the original network
  used to build `obj`.
- `obj` is modified: its likelihood fields (forward, directional & backward)
  are updated to make sure that they correspond to the current parameter values
  in `obj.model`, and to the `trait` of interest.

See also [`discrete_backwardlikelihood_tree!`](@ref) to update `obj.backwardlik`.

# examples

```jldoctest
julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> m1 = BinaryTraitSubstitutionModel([0.1, 0.1], ["lo", "hi"]);

julia> using DataFrames

julia> dat = DataFrame(species=["C","A","B","D"], trait=["hi","lo","lo","hi"]);

julia> fit1 = fitdiscrete(net, m1, dat);

julia> asr = ancestralStateReconstruction(fit1)
9×4 DataFrames.DataFrame
│ Row │ nodenumber │ nodelabel │ lo       │ hi       │
│     │ Int64      │ String    │ Float64  │ Float64  │
├─────┼────────────┼───────────┼──────────┼──────────┤
│ 1   │ 1          │ A         │ 1.0      │ 0.0      │
│ 2   │ 2          │ B         │ 1.0      │ 0.0      │
│ 3   │ 3          │ C         │ 0.0      │ 1.0      │
│ 4   │ 4          │ D         │ 0.0      │ 1.0      │
│ 5   │ 5          │ 5         │ 0.286019 │ 0.713981 │
│ 6   │ 6          │ 6         │ 0.319454 │ 0.680546 │
│ 7   │ 7          │ 7         │ 0.168549 │ 0.831451 │
│ 8   │ 8          │ 8         │ 0.76736  │ 0.23264  │
│ 9   │ 9          │ H1        │ 0.782777 │ 0.217223 │

julia> round.(exp.(fit1.postltw), digits=6) # marginal (posterior) probability that the trait evolved on each displayed tree
2-element Array{Float64,1}:
 0.919831
 0.080169

julia> using PhyloPlots

julia> plot(fit1.net, :R, nodeLabel = asr[!,[:nodenumber, :lo]], tipOffset=0.2); # pp for "lo" state
```
"""
ancestralStateReconstruction(obj::SSM) = ancestralStateReconstruction(obj, obj.activesite)
function ancestralStateReconstruction(obj::SSM, trait::Integer)
    # posterior probability of state i at node n: proportional to
    # sum_{trees t} exp( ltw[t] + backwardlik[i,n,t] + forwardlik[i,n,t] )
    bkd = obj.backwardlik
    fill!(bkd, 0.0) # initialize
    discrete_corelikelihood!(obj; whichtrait=trait) # update forward, direct, logtrans, postltw, loglik
    for t in 1:length(obj.displayedtree)
        discrete_backwardlikelihood_tree!(obj, t, trait)
    end
    # fixit: paralellize with
    # ll = pmap(t -> discrete_backwardlikelihood_tree!(obj,t, trait), 1:ntrees)
    k = nstates(obj.model)
    nnodes = length(obj.net.node)
    res = Array{Float64}(undef, k,nnodes)
    frd = obj.forwardlik
    ltw = obj.priorltw
    for i in 1:k
        for n in 1:nnodes
            res[i,n] = logsumexp(ltw + view(bkd, i,n,:) + view(frd, i,n,:))
        end
    end
    ll = obj.loglik
    map!(x -> exp(x - ll), res, res)  # to normalize: p_i / sum(p_j over all states j)
    # alternative syntax: res .= exp.(res .- obj.loglik)
    nodestringlabels = Vector{String}(undef, nnodes)
    for n in obj.net.node
        nodestringlabels[n.number] = (n.name == "" ? string(n.number) : n.name)
    end
    dat = DataFrame(transpose(res), Symbol.(getlabels(obj.model)))
    insertcols!(dat, 1, :nodenumber => collect(1:nnodes), makeunique=true)
    insertcols!(dat, 2, :nodelabel  => nodestringlabels,  makeunique=true)
    return dat
end

"""
    discrete_backwardlikelihood_tree!(obj::SSM, tree::Integer, trait::Integer)

Update `obj.backwardlik`; assume correct forward likelihood, directional likelihood
and transition probabilities.
"""
function discrete_backwardlikelihood_tree!(obj::SSM, t::Integer, trait::Integer)
    tree = obj.displayedtree[t]
    frdlik = view(obj.forwardlik, :,:,t)
    dirlik = view(obj.directlik , :,:,t)
    bkdlik = view(obj.backwardlik,:,:,t)
    k = nstates(obj.model)
    bkwtmp = Vector{Float64}(undef, k) # to hold bkw lik without parent edge transition
    if typeof(obj.model) == NASM 
        logprior = log.(stationary(obj.model))
    else #trait models
        logprior = [-log(k) for i in 1:k] # uniform prior at root
    end
    for ni in 1:length(tree.nodes_changed) # pre-order traversal to calculate backwardlik
        n = tree.nodes_changed[ni]
        nnum = n.number
        if ni == 1 # n is the root
            bkdlik[:,nnum] = logprior
        else
            pe = getMajorParentEdge(n)
            pn = getParent(pe)
            bkwtmp[:] = bkdlik[:,pn.number] # use bktmp's original memory
            for se in pn.edge
                if se != pe && pn == getParent(se) # then se is sister edge to pe
                    bkwtmp .+= view(dirlik, :,se.number)
                end
            end
            lt = view(obj.logtrans, :,:,pe.number)
            for j in 1:k # state at node n
                bkdlik[j,nnum] = logsumexp(bkwtmp + view(lt,:,j))
            end
        end
    end
    return nothing
end

# fixit: new type for two (dependent) binary traits
# need new algorithm for model at hybrid nodes: displayed trees aren't enough
"""
    learnlabels(model::Symbol, dat::DataFrame)

Return unique non-missing values in `dat`, and check that these labels
can be used to construct of substitution model of type `model`.
Examples:

```jldoctest
julia> using DataFrames

julia> dat = DataFrame(trait1 = ["A", "C", "A", missing])
4×1 DataFrame


julia> PhyloNetworks.learnlabels(:BTSM, dat)
2-element Array{String,1}:
 "A"
 "C"

julia> PhyloNetworks.learnlabels(:JC69, dat)
2-element Array{String,1}:
 "A"
 "C"
```
`"""
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
    startingBL!(net::HybridNetwork, trait::Vector{Vector} [, siteweight::Vector])

Calibrate branch lengths in `net` by minimizing the mean squared error
between the JC-adjusted pairwise distance between taxa, and network-predicted
pairwise distances, using [`calibrateFromPairwiseDistances!`](@ref).
`siteweight[k]` gives the weight of site (or site pattern) `k` (default: all 1s).

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
function startingBL!(net::HybridNetwork,
        trait::AbstractVector{Vector{Union{Missings.Missing,Int}}},
        siteweight=ones(length(trait[1]))::AbstractVector{Float64})
    nspecies = net.numTaxa
    M = zeros(Float64, nspecies, nspecies) # pairwise distances, initialized to zeros
    
    # count pairwise differences, then multiply by pattern weight
    ncols = length(trait[1]) # assumption: all species have same # columns
    length(siteweight) == ncols ||
      error("$(length(siteweight)) site weights but $ncols columns in the data")
    for i in 2:nspecies
        species1 = trait[i]
        for j in 1:(i-1)
            species2 = trait[j]
            for col in 1:ncols
                if !(ismissing(species1[col]) || ismissing(species2[col])) && (species1[col] != species2[col])
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
    return net
end

# ## Prep and Wrapper Functions ##

"""
    datatoSSM(net::HybridNetwork, dnadata::DataFrame, modsymbol::Symbol)

Create an SSM object for use in wrapper function. This should include all actions 
that can happen only once. Probably need different versions for different 
kinds of data (snp, amino acids), but works for DNA now.
Call readfastatodna(), startingrate(), startingBL!()
Similar to fitdiscrete()
"""
function datatoSSM(net::HybridNetwork, fastafile::String, modsymbol::Symbol)

    data, siteweights = readfastatodna(fastafile, true)
    model = symboltomodel(net, modsymbol, data, siteweights)
    ratemodel = RateVariationAcrossSites(1.0, 1) #TODO add option for users

    dat2 = traitlabels2indices(view(data, :, 2:ncol(data)), model)
    o, net = check_matchtaxonnames!(data[:,1], dat2, net)
    trait = dat2[o]
    startingBL!(net, trait, siteweights)
    obj = StatisticalSubstitutionModel(model, ratemodel, net, trait, siteweights)
end

"""
    addclades!(net::HybridNetwork, cladesDict::Dict{String, Vector{Int64}}

Add clade designation according to an array of clade dictionaries, where the array
has length equal to the number of species in the non-main clades 
(e.g. outgroups, known clades, etc).
note: all nodes start with clade designation "Main" or 1.
TODO See how this is done in other software (raxML). 
Maybe we force a polytomy at top of outgroup? Let users give us a tree?
"""
function addclades!(net::HybridNetwork, cladesDict::Dict{String,Vector{Int64}})
    #adds clades for leaves not in main clade
    #interate over items in input dictionary 

    #compare keys with net.leaf
    for key in keys(cladesDict)
        for index in 1:length(net.leaf) #net.leaf is an array of the network's leaves
            if net.leaf[index].name == key
                #removes main clade designation
                filter!(x->x≠1,net.leaf[index].clade)
                #adds new clade membership
                for i in 1:length(cladesDict[key])
                    cladetoadd = cladesDict[key][i]
                    push!(net.leaf[index].clade, cladetoadd)
                    #removes duplicates to catche when funct run with same values
                    net.leaf[index].clade = unique(net.leaf[index].clade) 
                end
            end
        end
    end
end

"""
    symboltomodel(network, modsymbol::Symbol, data::DataFrame, siteweights::Vector)

Take a symbol and data to create a relative statistical substitution model using
default site weights.
"""
#TODO is net necessary?
function symboltomodel(net::HybridNetwork, modsymbol::Symbol, data::DataFrame,
        siteweights=repeat([1.], inner=ncol(data))::Vector)
    rate = startingrate(net)
    actualdat = view(data, :, 2:ncol(data))
    labels = learnlabels(modsymbol, actualdat)
    if modsymbol == :JC69
        return JC69([1.0], true) # 1.0 instead of rate because relative version
    elseif modsymbol == :HKY85 # transition/transversion rate ratio 
        return HKY85([1.0], empiricalDNAfrequencies(actualdat, siteweights), true)
    elseif modsymbol == :ERSM
        return EqualRatesSubstitutionModel(length(labels), rate, labels);
    elseif modsymbol == :BTSM
        return BinaryTraitSubstitutionModel([rate, rate], labels)
    elseif modsymbol == :TBTSM
        return TwoBinaryTraitSubstitutionModel([rate, rate, rate, rate, rate, rate, rate, rate], labels)
    else
        error("model $modsymbol is unknown or not implemented yet")
    end
end
