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
    ratemodel::RateVariationAcrossSites
    """stationary for NASM models (log(1/4)), uniform in all other cases"""
    prioratroot::Vector{Float64}
    net::HybridNetwork
    """ data: trait[i] for leaf with n.number = i
        type Int: for indices of trait labels in getlabels(model)
        allows missing, but not half-ambiguous and not multi-state"""
    trait::Vector{Vector{Union{Missing, Int}}}
    "number of columns in the data: should be equal to `length(trait[i])` for all i"
    nsites::Int
    "vector of weights, one for each column, taken as 1s if not provided"
    siteweight::Union{Nothing, Vector{Float64}}
    "total sum of site weights"
    totalsiteweight::Float64
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
    # type based on extracting displayed trees
    displayedtree::Vector{HybridNetwork}
    """
    prior log tree weight: log product of γ's.
    In fit!: priorltw = `inheritanceWeight.(trees)`
    (which returns missing for any negative γ and would cause an error here)
    """
    priorltw::Vector{Float64}
    """
    partial log-likelihoods for active trait(s) given a particular displayed tree
    and a particular rate category, at indices [i, n.number or e.number]:
    - forward likelihood: log P{data below node n (in tree t) given state i at n}
    - direct likelihood:  log P{data below edge e (in tree t) given i at parent of e}
    - backward likelihood:log P{data at all non-descendants of n (in t) and state i at n}
    sizes: k, net.numNodes or net.numEdges.
    """
    # reset for each trait and rate
    forwardlik::Array{Float64,2} # size: k, net.numNodes
    directlik::Array{Float64,2}  # size: k, net.numEdges
    backwardlik::Array{Float64,2}# size: k, net.numNodes
    "log-likelihood of site k"
    _sitecache::Array{Float64,1} # size: nsites
    "log-likelihood of ith displayed tree t & rate category j, of site k"
    _loglikcache::Array{Float64, 3} # size: nsites, nrates, ntrees

    "inner (default) constructor: from model, rate model, network, trait and site weights"
    function StatisticalSubstitutionModel(model::SubstitutionModel,
            ratemodel::RateVariationAcrossSites,
            net::HybridNetwork, trait::AbstractVector,
            siteweight=nothing::Union{Nothing, Vector{Float64}},
            maxhybrid=length(net.hybrid)::Int)
        length(trait) > 0 || error("no trait data!")
        nsites = length(trait[1])
        siteweight === nothing || length(siteweight) == nsites ||
            error("siteweight must be of same length as the number of traits")
        totalsiteweight = (siteweight === nothing ? float(nsites) : sum(siteweight))
        k = nstates(model)
        if typeof(model) <: NASM
            prioratroot = log.(stationary(model)) #refactor to save in obj
        else # other trait models
            prioratroot = [-log(k) for i in 1:k] # uniform prior at root
        end
        # T = eltype(getlabels(model))
        # extract displayed trees
        trees = displayedTrees(net, 0.0; nofuse=true, keeporiginalroot=true)
        for tree in trees
            preorder!(tree) # no need to call directEdges! before: already done on net
        end
        ntrees = 2^maxhybrid
        ntrees >= length(trees) ||
            error("""maxhybrid is too low.
                    Call using maxhybrid >= current number of hybrids""")
        # log tree weights: sum log(γ) over edges, for each displayed tree
        priorltw = inheritanceWeight.(trees)
        all(!ismissing, priorltw) ||
          error("one or more inheritance γ's are missing or negative. fix using setGamma!(network, edge)")
        maxedges = length(net.edge) + 3*(maxhybrid-length(net.hybrid))
        maxnodes = length(net.node) + 2*(maxhybrid-length(net.hybrid))
        logtrans   = zeros(Float64, k,k, maxedges, length(ratemodel.ratemultiplier))
        forwardlik = zeros(Float64, k, maxnodes)
        directlik  = zeros(Float64, k, maxedges)
        backwardlik= zeros(Float64, k, maxnodes)
        _sitecache = Vector{Float64}(undef, nsites)
        _loglikcache = zeros(Float64, nsites, length(ratemodel.ratemultiplier), ntrees)
        new(deepcopy(model), deepcopy(ratemodel), prioratroot,
            net, trait, nsites, siteweight, totalsiteweight, missing, # missing log likelihood
            logtrans, trees,
            priorltw, forwardlik, directlik, backwardlik,_sitecache,_loglikcache)
    end
end
const SSM = StatisticalSubstitutionModel

# fasta constructor: from net, fasta filename, modsymbol, and maxhybrid
# Works for DNA in fasta format. Probably need different versions for
# different kinds of data (snp, amino acids). Similar to fitdiscrete()
"""
    StatisticalSubstitutionModel(model::SubstitutionModel,
            ratemodel::RateVariationAcrossSites,
            net::HybridNetwork, trait::AbstractVector,
            siteweight=nothing::Union{Nothing, Vector{Float64}},
            maxhybrid=length(net.hybrid)::Int)

Inner constructor. Makes a deep copy of the input model, rate model.
Warning: does *not* make a deep copy of the network:
modification of the `object.net` would modify the input `net`.
Assumes that the network has valid gamma values (to extract displayed trees).

    StatisticalSubstitutionModel(net::HybridNetwork, fastafile::String,
            modsymbol::Symbol, rvsymbol=:noRV::Symbol,
            ratecategories=4::Int;
            maxhybrid=length(net.hybrid)::Int)

Constructor from a network and a fasta file.
The model symbol should be one of `:JC69`, `:HKY85`, `:ERSM` or `:BTSM`.
The `rvsymbol` should be as required by [`RateVariationAcrossSites`](@ref).

The network's gamma values are modified if they are missing. After that,
a deep copy of the network is passed to the inner constructor.
"""
function StatisticalSubstitutionModel(net::HybridNetwork, fastafile::String,
        modsymbol::Symbol, rvsymbol=:noRV::Symbol, ratecategories=4::Int;
        maxhybrid=length(net.hybrid)::Int)
    for e in net.edge # check for missing or inappropriate γ values
        if e.hybrid
            e.gamma > 0.0 && continue
            setGamma!(e, (e.isMajor ? 0.6 : 0.4)) # to maintain isMajor as is
        else
            e.gamma == 1.0 || error("tree edge number $(e.number) has γ not 1.0")
        end
    end
    data, siteweights = readfastatodna(fastafile, true)
    model = defaultsubstitutionmodel(net, modsymbol, data, siteweights)
    ratemodel = RateVariationAcrossSites(rvsymbol, ratecategories)
    dat2 = traitlabels2indices(view(data, :, 2:size(data,2)), model)
    # check_matchtaxonnames makes a deep copy of the network
    o, net = check_matchtaxonnames!(data[:,1], dat2, net) # calls resetNodeNumbers, which calls preorder!
    trait = dat2[o]
    obj = StatisticalSubstitutionModel(model, ratemodel, net, trait, siteweights,
                           maxhybrid)
end

StatsBase.loglikelihood(obj::SSM) = obj.loglik
StatsBase.islinear(SSM) = false
StatsBase.dof(obj::SSM) = nparams(obj.model) + nparams(obj.ratemodel)
function Base.show(io::IO, obj::SSM)
    disp =  "$(typeof(obj)):\n"
    disp *= replace(string(obj.model), r"\n" => "\n  ") * "\n"
    if nparams(obj.ratemodel) > 0
        disp *= replace(string(obj.ratemodel), r"\n" => "\n  ") * "\n"
    end
    disp *= "on a network with $(obj.net.numHybrids) reticulations\n"
    print(io, disp)
    showdata(io, obj)
    if !ismissing(obj.loglik)
        print(io, "\nlog-likelihood: $(round(obj.loglik, digits=5))")
    end
end
"""
    showdata(io::IO, obj::SSM, fullsiteinfo=false::Bool)

Return information about the data in an SSM object:
number of species, number or traits or sites, number of distinct patterns,
and more information if `fullsiteinfo` is true:
number sites with missing data only,
number of invariant sites, number of sites with 2 distinct states,
number of parsimony-informative sites (with 2+ states being observed in 2+ tips),
number of sites with some missing data, and
overall proportion of entries with missing data.

Note: Missing is not considered an additional state. For example,
if a site contains some missing data, but all non-missing values take the same
state, then this site is counted in the category "invariant".
"""
function showdata(io::IO, obj::SSM, fullsiteinfo=false::Bool)
    disp =  "data:\n  $(length(obj.trait)) species"
    ns = obj.totalsiteweight
    ns = (isapprox(ns, round(ns), atol=1e-5) ? Int(round(ns)) : ns)
    disp *= (ns ≈ 1 ? "\n  $ns trait" : "\n  $ns sites")
    if !isapprox(obj.nsites, ns, atol=1e-5)
        disp *= "\n  $(obj.nsites) distinct patterns"
    end
    print(io, disp)
    (fullsiteinfo && obj.nsites != 1) || return nothing
    # if more than 1 trait and if the user wants full information:
    nsv = MVector{6,Float64}(undef) # vector to count number of
    # sites with: 0, 1, 2 states, parsimony informative, with 1+ missing value,
    # missing values across all sites.
    fill!(nsv, 0.0)
    text = ["sites with no data", "invariant sites",
            "sites with 2 distinct states", "parsimony-informative sites",
            "sites with 1 or more missing values", "missing values overall"]
    trackstates = zeros(Int, nstates(obj.model)) # states seen at a given site
    ntaxa = length(obj.trait)
    for i in 1:(obj.nsites) # over sites
        sweight = (isnothing(obj.siteweight) ? 1.0 : obj.siteweight[i])
        missone = false
        fill!(trackstates, 0)
        for j in 1:ntaxa # over taxa
            data = obj.trait[j][i]
            if ismissing(data)
                nsv[6] += sweight # total number of missing values
                missone && continue
                nsv[5] += sweight # sites with 1+ missing values
                missone = true
            else # mark state seen
                trackstates[data] += 1 # 1 more taxon
            end
        end
        # add site's weight to appropriate nstates
        nstates = sum(trackstates .> 0)
        if nstates < 3
            nsv[nstates+1] += sweight
        end
        if nstates>1 # are there 2 states observed at 2+ taxa each?
            nstates_2taxa = sum(trackstates .> 1)
            if nstates_2taxa>1
                nsv[4] += sweight
            end
        end
    end
    nsv_r = map(x -> begin y=round(x); (isapprox(y,x,atol=1e-5) ? Int(y) : x); end, nsv)
    for i in 1:5
        print(io, "\n  $(nsv_r[i]) $(text[i]) ($(round(100*nsv[i]/ns, digits=2))%)")
    end
    print(io, "\n  $(round(100*nsv[6]/(ns*ntaxa), digits=2))% $(text[6])")
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
- `dnadata`: the first part of the output of readfastatodna,
  a dataframe of BioSequence DNA sequences, with taxon in column 1 and
  a column for each site.
- `dnapatternweights`: the second part of the output of readfastatodna,
  an array of weights, one weights for each of the site columns.
  The length of the weight is equal to nsites.
  If using dnapatternweights, must provide dnadata.
- RateVariationAcrossSites: model for rate variation (optional)

Optional arguments (default):
- `optimizeQ` (true): should model rate parameters be fixed,
  or should they be optimized?
- `optimizeRVAS` (true): should the model optimize the parameters
  for the variability of rates across sites (α and/or p_invariable)?
- `NLoptMethod` (`:LN_COBYLA`, derivative-free) for the optimization algorithm.
  For other options, see the
  [NLopt](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/).
- tolerance values to control when the optimization is stopped:
  `ftolRel` (1e-12), `ftolAbs` (1e-10) on the likelihood, and
  `xtolRel` (1e-10), `xtolAbs` (1e-10) on the model parameters.
- bounds for the alpha parameter of the Gamma distribution of
  rates across sites: `alphamin=0.05`, `alphamax=50`.
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
on a network with 1 reticulations
data:
  4 species
  1 trait
log-likelihood: -2.7277

julia> tips = Dict("A" => "lo", "B" => "lo", "C" => "hi", "D" => "hi");

julia> fit2 = fitdiscrete(net, m1, tips; xtolRel=1e-16, xtolAbs=1e-16, ftolRel=1e-16)
PhyloNetworks.StatisticalSubstitutionModel:
Binary Trait Substitution Model:
  rate lo→hi α=0.27222
  rate hi→lo β=0.34981
on a network with 1 reticulations
data:
  4 species
  1 trait
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
  off-diagonal rates equal to 0.29233/3.
  rate matrix Q:
                 A       C       G       T
         A       *  0.0974  0.0974  0.0974
         C  0.0974       *  0.0974  0.0974
         G  0.0974  0.0974       *  0.0974
         T  0.0974  0.0974  0.0974       *
on a network with 0 reticulations
data:
  4 species
  1 trait
log-likelihood: -4.99274

julia> rv = RateVariationAcrossSites(alpha=1.0, ncat=4)
Rate variation across sites: discretized Gamma
alpha: 1.0
categories for Gamma discretization: 4
rates: [0.146, 0.513, 1.071, 2.27]

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
Rate variation across sites: discretized Gamma
  alpha: 1.0
  categories for Gamma discretization: 4
  rates: [0.146, 0.513, 1.071, 2.27]
on a network with 0 reticulations
data:
  4 species
  1 trait
log-likelihood: -5.2568

```
fixit: add option to allow users to specify root prior,
using either equal frequencies or stationary frequencies for trait models.
"""
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel,
    tips::Dict; kwargs...) #tips::Dict no ratemodel version
    ratemodel = RateVariationAcrossSites(ncat=1)
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
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel,
    dat::DataFrame; kwargs...)
    ratemodel = RateVariationAcrossSites(ncat=1)
    fitdiscrete(net, model, ratemodel, dat; kwargs...)
end

#dat::DataFrame with rate model version
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel,
    ratemodel::RateVariationAcrossSites, dat::DataFrame; kwargs...)
    i = findfirst(isequal(:taxon), DataFrames.propertynames(dat))
    if i===nothing i = findfirst(isequal(:species), DataFrames.propertynames(dat)); end
    if i===nothing i=1; end # first column if no column "taxon" or "species"
    j = findfirst(isequal(:trait), DataFrames.propertynames(dat))
    if j===nothing j=2; end
    if i==j
        error("""expecting taxon names in column 'taxon', or 'species' or
        column 1, and trait values in column 'trait' or column 2.""")
    end
    species = dat[:,i]    # modified in place later
    dat = traitlabels2indices(dat[!,j], model)   # vec of vec, indices
    o, net = check_matchtaxonnames!(species, dat, net)
    StatsBase.fit(StatisticalSubstitutionModel, net, model, ratemodel, view(dat, o); kwargs...)
end

#species, dat version, no ratemodel
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel,
    species::Array{String}, dat::DataFrame; kwargs...)
    ratemodel = RateVariationAcrossSites(ncat=1)
    fitdiscrete(net, model, ratemodel, species, dat; kwargs...)
end

#species, dat version with ratemodel
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel,
    ratemodel::RateVariationAcrossSites, species::Array{String},
    dat::DataFrame; kwargs...)
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
                       empiricalDNAfrequencies(dat, repeat([1.], inner=size(dat, 2))),
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

    rvas = RateVariationAcrossSites(rvSymbol)
    fitdiscrete(net, model, rvas, species, dat; kwargs...)
end

#dnadata with dnapatternweights version, no ratemodel
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel,
    dnadata::DataFrame, dnapatternweights::Array{Float64}; kwargs...)
    ratemodel = RateVariationAcrossSites(ncat=1)
    fitdiscrete(net, model, ratemodel, dnadata, dnapatternweights; kwargs...)
end

#dnadata with dnapatternweights version with ratemodel
function fitdiscrete(net::HybridNetwork, model::SubstitutionModel,
    ratemodel::RateVariationAcrossSites,dnadata::DataFrame,
    dnapatternweights::Array{Float64}; kwargs...)
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
        model = JC69([1.0], true)  # 1.0 instead of rate because relative version
    elseif modSymbol == :HKY85
        model = HKY85([1.0], # transition/transversion rate ratio
                      empiricalDNAfrequencies(view(dnadata, :, 2:size(dnadata,2)), dnapatternweights),
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

    rvas = RateVariationAcrossSites(rvSymbol)
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
"""
function StatsBase.fit(::Type{SSM}, net::HybridNetwork, model::SubstitutionModel,
    ratemodel::RateVariationAcrossSites, trait::AbstractVector; kwargs...)

    sw = nothing
    if haskey(kwargs, :siteweights)
        sw = kwargs[:siteweights]
        kwargs = filter(p -> p.first != :siteweights, kwargs)
    end
    obj = StatisticalSubstitutionModel(model, ratemodel, net, trait, sw)
    fit!(obj; kwargs...)
end

function fit!(obj::SSM; optimizeQ=true::Bool, optimizeRVAS=true::Bool,
    closeoptim=false::Bool, verbose=false::Bool, maxeval=1000::Int,
    ftolRel=fRelBL::Float64, ftolAbs=fAbsBL::Float64,
    xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64,
    alphamin=alphaRASmin, alphamax=alphaRASmax,
    pinvmin=pinvRASmin, pinvmax=pinvRASmax)

    all(x -> x >= 0.0, [e.length for e in obj.net.edge]) || error("branch lengths should be >= 0")
    all(x -> x >= 0.0, [e.gamma for e in obj.net.edge]) || error("gammas should be >= 0")
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
    if optimizeQ
        function loglikfun(x::Vector{Float64}, grad::Vector{Float64}) # modifies obj
            setrates!(obj.model, x)
            res = discrete_corelikelihood!(obj)
            verbose && println("loglik: $res, model rates: $x")
            length(grad) == 0 || error("gradient not implemented")
            return res
        end
        # optimize Q under fixed RVAS parameters
        # set-up optimization object for Q
        NLoptMethod=:LN_COBYLA # no gradient # :LN_COBYLA for (non)linear constraits, :LN_BOBYQA for bound constraints
        nparQ = nparams(obj.model)
        optQ = NLopt.Opt(NLoptMethod, nparQ)
        NLopt.ftol_rel!(optQ,ftolRel) # relative criterion
        NLopt.ftol_abs!(optQ,ftolAbs) # absolute criterion
        NLopt.xtol_rel!(optQ,xtolRel)
        NLopt.xtol_abs!(optQ,xtolAbs)
        NLopt.maxeval!(optQ, maxeval) # max number of iterations
        # NLopt.maxtime!(optQ, t::Real)
        NLopt.lower_bounds!(optQ, zeros(Float64, nparQ))
        if typeof(obj.model) == HKY85 # set an upper bound on kappa values
            NLopt.upper_bounds!(optQ, fill(kappamax,nparQ))
        end
        NLopt.max_objective!(optQ, loglikfun)
        fmax, xmax, ret = NLopt.optimize(optQ, obj.model.rate)
        setrates!(obj.model, xmax)
        obj.loglik = fmax
        verbose && println("got $(round(fmax, digits=5)) at $(round.(xmax, digits=5)) after $(optQ.numevals) iterations (return code $(ret))")
    end
    if optimizeRVAS
        function loglikfunRVAS(alpha::Vector{Float64}, grad::Vector{Float64})
            setparameters!(obj.ratemodel, alpha)
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
        rvind = getparamindex(obj.ratemodel)
        NLopt.lower_bounds!(optRVAS, [pinvmin,alphamin][rvind] )
        NLopt.upper_bounds!(optRVAS, [pinvmax,alphamax][rvind] )
        NLopt.max_objective!(optRVAS, loglikfunRVAS)
        fmax, xmax, ret = NLopt.optimize(optRVAS, getparameters(obj.ratemodel)) # optimization here!
        setparameters!(obj.ratemodel, xmax)
        obj.loglik = fmax
        verbose && println("RVAS: got $(round(fmax, digits=5)) at $(round.(xmax, digits=5)) after $(optRVAS.numevals) iterations (return code $(ret))")
    end
    if optimizeQ && optimizeRVAS && closeoptim
        # optimize Q under fixed RVAS parameters: a second time
        fmax, xmax, ret = NLopt.optimize(optQ, obj.model.rate)
        setrates!(obj.model, xmax)
        obj.loglik = fmax
        verbose && println("got $(round(fmax, digits=5)) at $(round.(xmax, digits=5)) after $(optQ.numevals) iterations (return code $(ret))")
        # optimize RVAS under fixed Q: a second time
        fmax, xmax, ret = NLopt.optimize(optRVAS, getparameters(obj.ratemodel))
        setparameters!(obj.ratemodel, xmax)
        obj.loglik = fmax
        verbose && println("RVAS: got $(round(fmax, digits=5)) at $(round.(xmax, digits=5)) after $(optRVAS.numevals) iterations (return code $(ret))")
    end
    return obj
end

"""
    update_logtrans(obj::SSM)

Initialize and update `obj.logtrans`, the log transition probabilities
along each edge in the full network.
They are re-used for each displayed tree, which is why edges are not fused
around degree-2 nodes when extracting displayed trees.
"""
function update_logtrans(obj::SSM)
    rates = obj.ratemodel.ratemultiplier
    k = nstates(obj.model)
    Ptmp = MMatrix{k,k,Float64}(undef) # memory to be re-used
    for edge in obj.net.edge # update logtrans: same for all displayed trees, all traits
        enum = edge.number
        len = edge.length
        for i = 1:length(rates)
            obj.logtrans[:,:,enum, i] .= log.(P!(Ptmp, obj.model, len * rates[i]))
        end
    end
end

"""
    update_logtrans(obj::SSM, edge::Edge)

Update the log-transition probabilities associates to one particular `edge`
in the network.
"""
function update_logtrans(obj::SSM, edge::Edge)
    rates = obj.ratemodel.ratemultiplier
    enum = edge.number
    len = edge.length
    for i in 1:length(rates)
        pmat = view(obj.logtrans, :,:,enum,i)
        @inbounds pmat .= log.(P!(pmat, obj.model, len * rates[i]))
    end
end

"""
    discrete_corelikelihood!(obj::StatisticalSubstitutionModel;
                             whichtrait::AbstractVector{Int} = 1:obj.nsites)

Calculate the likelihood and update `obj.loglik` for discrete characters on a network,
calling [`discrete_corelikelihood_trait!`](@ref).
The algorithm extracts all displayed trees and weighs the likelihood under all these trees.
The object's partial likelihoods are updated:
- forward and direct partial likelihoods are re-used, one trait at a time,
- overall likelihoods on each displayed tree, given each rate category and for
  each given site/trait: are cached in `_loglikcache`.
"""
function discrete_corelikelihood!(obj::SSM; whichtrait::AbstractVector{Int} = 1:obj.nsites)
    # fill _loglikcache
    nr = length(obj.ratemodel.ratemultiplier)
    nt = length(obj.displayedtree)
    update_logtrans(obj)
    for t in 1:nt
        for ri in 1:nr
            for ci in whichtrait
                obj._loglikcache[ci,ri,t] = discrete_corelikelihood_trait!(obj,t,ci,ri)
                # conditional: log P(site | rate & tree)
                # note: -Inf expected if 2 tips have different states, but separated by path of total length 0.0
            end
        end
    end
    # aggregate over trees and rates
    lrw = obj.ratemodel.lograteweight
    for ti in 1:nt
        ltprior = obj.priorltw[ti]
        for ri in 1:nr
            obj._loglikcache[:,ri,ti] .+= ltprior + lrw[ri]
            # now unconditional: log P(site & rate & tree)
        end
    end
    for ci in whichtrait
        obj._sitecache[ci] = logsumexp(view(obj._loglikcache, ci,:,1:nt))
    end
    if obj.siteweight !== nothing
        obj._sitecache .*= obj.siteweight
    end
    loglik = sum(obj._sitecache)
    obj.loglik = loglik
    return loglik
end

"""
    discrete_corelikelihood_trait!(obj::SSM, t::Integer, ci::Integer, ri::Integer)

Return the likelihood for tree `t`, trait (character/site) index `ci` and rate category `ri`.
Update & modify the forward & directional log-likelihoods `obj.forwardlik`
and `obj.directlik`, which are indexed by [state, node_number or edge_number].
Used by [`discrete_corelikelihood!`](@ref).

**Preconditions**: `obj.logtrans` updated, edges directed, nodes/edges preordered
"""
function discrete_corelikelihood_trait!(obj::SSM, t::Integer, ci::Integer, ri::Integer)
    forwardlik = obj.forwardlik
    directlik  = obj.directlik
    tree = obj.displayedtree[t]
    k = nstates(obj.model)   # also = size(logtrans,1) if not RateVariationAcrossSites
    fill!(forwardlik, 0.0) # re-initialize for each trait, each iteration
    fill!(directlik,  0.0)
    loglik = 0.
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
            loglik = logsumexp(obj.prioratroot + view(forwardlik, :,nnum)) # log P{data for ci | tree t, rate ri}
            break # out of loop over nodes
        end
        # if we keep going, n is not the root
        # calculate direct likelihood on the parent edge of n
        for e in n.edge
            if n == getChild(e)
                lt = view(obj.logtrans, :,:,e.number, ri)
                for i in 1:k # state at parent node
                    directlik[i,e.number] = logsumexp(view(lt,i,:) + view(forwardlik,:,nnum))
                end
                break # we visited the parent edge: break out of for loop
            end
        end #loop over edges
    end # of loop over nodes
    return loglik
end

"""
    posterior_logtreeweight(obj::SSM, trait = 1)

Return an array A such that A[t] = log of P(tree `t` and trait `trait`)
if a single `trait` is requested, or A[i,t]= log of P(tree `t` and trait `i`)
if `trait` is a vector or range (e.g. `trait = 1:obj.nsites`).
These probabilities are conditional on the model parameters in `obj`.

**Precondition**: `_loglikcache` updated by [`discrete_corelikelihood!`](@ref)

# examples

```jldoctest
julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> m1 = BinaryTraitSubstitutionModel([0.1, 0.1], ["lo", "hi"]); # arbitrary rates

julia> using DataFrames

julia> dat = DataFrame(species=["C","A","B","D"], trait=["hi","lo","lo","hi"]);

julia> fit = fitdiscrete(net, m1, dat); # optimized rates: α=0.27 and β=0.35

julia> pltw = PhyloNetworks.posterior_logtreeweight(fit)
2-element Array{Float64,1}:
 -0.08356519024776699
 -2.523619878044531  

julia> exp.(pltw) # posterior trees probabilities (sum up to 1)
2-element Array{Float64,1}:
 0.9198311206979973 
 0.08016887930200293

julia> round.(exp.(fit.priorltw), digits=4) # the prior tree probabilities are similar here (tiny data set!)
2-element Array{Float64,1}:
 0.9
 0.1
```
"""
function posterior_logtreeweight(obj::SSM, trait = 1)
    # ts[site,tree] = log P(data and tree) at site, integrated over rates
    d = length(size(trait)) # 0 if single trait, 1 if vector of several traits
    ts = dropdims(mapslices(logsumexp, view(obj._loglikcache, trait,:,:),
                            dims=d+1); dims=1)
    if d>0 ts = permutedims(ts); end # now: ts[tree] or ts[tree,site]
    siteliks = mapslices(logsumexp, ts, dims=1) # 1 x ntraits array (or 1-element vector)
    ts .-= siteliks
    return ts
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
                if vi === nothing
                    # ideally, handle ambiguous DNA types optimally
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
    resetNodeNumbers!(net; checkPreorder=true, type=:ape) # tip species: now with numbers 1:n
    resetEdgeNumbers!(net, false) # to use edge as indices: 1:numEdges
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
    ancestralStateReconstruction(obj::SSM, trait::Integer = 1)

Estimate the marginal probability of ancestral states for discrete character
number `trait` (first trait by default).
The parameters of the [`StatisticalSubstitutionModel`](@ref) object `obj`
must first be fitted using [`fitdiscrete`](@ref), and ancestral state reconstruction
is conditional on the estimated parameters. If these parameters were estimated
using all traits, they are used as is, to do ancestral state reconstruction of the
particular `trait` of interest.

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

limitations: the following are not checked.
- Assumes that every node in the large network is also present
  (with descendant leaves) in each displayed tree.
  This is not true if the network is not tree-child...
- Assumes that the root is also in each displayed tree, which
  may not be the case if the root had a hybrid child edge.

See also [`posterior_logtreeweight`](@ref) and
[`discrete_backwardlikelihood_trait!`](@ref) to update `obj.backwardlik`.

# examples

```jldoctest
julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> m1 = BinaryTraitSubstitutionModel([0.1, 0.1], ["lo", "hi"]);

julia> using DataFrames

julia> dat = DataFrame(species=["C","A","B","D"], trait=["hi","lo","lo","hi"]);

julia> fit1 = fitdiscrete(net, m1, dat);

julia> asr = ancestralStateReconstruction(fit1)
9×4 DataFrame
│ Row │ nodenumber │ nodelabel │ lo       │ hi       │
│     │ Int64      │ String    │ Float64  │ Float64  │
├─────┼────────────┼───────────┼──────────┼──────────┤
│ 1   │ 1          │ A         │ 1.0      │ 0.0      │
│ 2   │ 2          │ B         │ 1.0      │ 0.0      │
│ 3   │ 3          │ C         │ 0.0      │ 1.0      │
│ 4   │ 4          │ D         │ 0.0      │ 1.0      │
│ 5   │ 5          │ 5         │ 0.286018 │ 0.713982 │
│ 6   │ 6          │ 6         │ 0.319454 │ 0.680546 │
│ 7   │ 7          │ 7         │ 0.168548 │ 0.831452 │
│ 8   │ 8          │ 8         │ 0.767361 │ 0.232639 │
│ 9   │ 9          │ H1        │ 0.782777 │ 0.217223 │

julia> using PhyloPlots

julia> plot(fit1.net, :R, nodeLabel = asr[!,[:nodenumber, :lo]], tipOffset=0.2); # pp for "lo" state
```
"""
function ancestralStateReconstruction(obj::SSM, trait::Integer = 1)
    # posterior probability of state i at node n: proportional to
    # sum_{tree t, rate r} exp( ltw[t] + backwardll[i,n] given t,r + forwardll[i,n] given t,r ) / nr
    trait <= obj.nsites || error("trait $trait is larger than the number of traits in the data")
    nnodes = length(obj.net.node) # may be smaller than 2nd size of bkd or frd
    update_logtrans(obj)
    bkd = view(obj.backwardlik, :, 1:nnodes)
    frd = view(obj.forwardlik, :, 1:nnodes)
    ltw = obj.priorltw
    res = similar(bkd) # first: hold the cumulative logsumexp of bkd + frd + ltw
    fill!(res, -Inf64)
    nr = length(obj.ratemodel.ratemultiplier)
    lrw = obj.ratemodel.lograteweight
    for t in 1:length(obj.displayedtree)
        ltprior = ltw[t]
        for ri in 1:nr
            # update forward & directional likelihoods
            discrete_corelikelihood_trait!(obj,t,trait,ri)
            # update backward likelihoods
            discrete_backwardlikelihood_trait!(obj,t,ri)
            # P{state i at node n} ∝ bkd[i,n] * frd[i,n] given tree & rate:
            # res = logaddexp(res, ltw[t] + lrw[ri] + bkd + frd)
            broadcast!(logaddexp, res, res, (ltprior + lrw[ri]) .+ bkd .+ frd)
        end
    end
    # normalize the results at each node: p_i / sum(p_j over all states j)
    traitloglik = logsumexp(res[:,1]) # sum p_j at node 1 or at any node = loglikelihood
    res .= exp.(res .- traitloglik)
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
    discrete_backwardlikelihood_trait!(obj::SSM, tree::Integer, ri::Integer)

Update and return the backward likelihood (last argument `backwardlik`)
assuming rate category `ri` and tree index `tree`,
using current forward and backwards likelihoods in `obj`:
these depend on the trait (or site) given to the last call to
`discrete_corelikelihood_trait!`.
Used by `ancestralStateReconstruction`.

**warning**: assume correct transition probabilities.
"""
function discrete_backwardlikelihood_trait!(obj::SSM, t::Integer, ri::Integer)
    backwardlik = obj.backwardlik
    directlik  = obj.directlik
    tree = obj.displayedtree[t]
    k = nstates(obj.model)
    fill!(backwardlik, 0.0) # re-initialize for each trait, each iteration
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
            backwardlik[:,nnum] = logprior
        else
            pe = getMajorParentEdge(n)
            pn = getParent(pe)
            bkwtmp[:] = backwardlik[:,pn.number] # use bktmp's original memory
            for se in pn.edge
                if se != pe && pn == getParent(se) # then se is sister edge to pe
                    bkwtmp .+= view(directlik, :,se.number)
                end
            end
            lt = view(obj.logtrans, :,:,pe.number,ri)
            for j in 1:k # state at node n
                backwardlik[j,nnum] = logsumexp(bkwtmp + view(lt,:,j))
            end
        end
    end
    return backwardlik
end

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
    defaultsubstitutionmodel(network, modsymbol::Symbol, data::DataFrame,
                  siteweights::Vector)

Return a statistical substitution model (SSM) with appropriate state labels
and a rate appropriate for the branch lengths in `net`
(see [`startingrate`](@ref)).
The `data` frame must have the actual trait/site data in columns 2 and up,
as when the species names are in column 1.
For DNA data, the relative rate model is returned, with a
stationary distribution equal to the empirical frequencies.
"""
function defaultsubstitutionmodel(net::HybridNetwork, modsymbol::Symbol, data::DataFrame,
        siteweights=repeat([1.], inner=size(data,2))::AbstractVector)
    rate = startingrate(net)
    actualdat = view(data, :, 2:size(data,2))
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
        return TwoBinaryTraitSubstitutionModel([rate, rate, rate, rate, rate,
                                                rate, rate, rate], labels)
    else
        error("model $modsymbol is unknown or not implemented yet")
    end
end

# fixit: new type for two (dependent) binary traits
# need new algorithm for model at hybrid nodes: displayed trees aren't enough
