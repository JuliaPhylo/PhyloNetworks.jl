"""
    StatisticalSubstitutionModel

Subtype of `StatsBase.StatisticalModel`, to fit discrete data to a model
of trait substitution along a network.
See [`fitDiscrete`](@ref) to fit a trait substitution model to discrete data.
It returns an object of type `StatisticalSubstitutionModel`, to which standard
functions can be applied, like `loglikelihood(object)`, `aic(object)` etc.
"""
mutable struct StatisticalSubstitutionModel{T} <: StatsBase.StatisticalModel
    model::TraitSubstitutionModel{T} # T: type of trait labels
    net::HybridNetwork
    # data: trait[i] for leaf with n.number = i
    #       type Int: for indices of trait labels in model.label
    #       allows missing, but not half-ambiguous and not multi-state
    trait::Vector{Vector{Union{Missings.Missing,Int}}}
    ntraits::Int
    loglik::Union{Missings.Missing, Float64}
    # log of transition probabilities: where e goes from X to Y
    # logtrans[i,j,e] = P{Y=j|X=i} where i=start_state, j=end_state, e=edge.number
    # all displayed trees use same edge numbers and same logtrans as in full network
    logtrans::Array{Float64,3}   # size: k,k, net.numEdges where k=nStates(model)
    activetrait::Int
    # type based on extracting displayed trees
    displayedtree::Vector{HybridNetwork}
    priorltw::Vector{Float64} # prior log tree weight: log product of γ's
    postltw::Vector{Float64}  # posterior log tree weight: P{tree and data}
    # partial likelihoods for active trait, at indices [i, n.number or e.number ,t] :
    # - forward likelihood: log P{data below node n in tree t given state i at n}
    # - direct likelihood:  log P{data below edge e in tree t given i at parent of e}
    # - backward likelihood:log P{data at all non-descendants of n and state i at n}
    forwardlik::Array{Float64,3} # size: k, net.numNodes, number of displayed trees
    directlik::Array{Float64,3}  # size: k, net.numEdges, number of displayed trees
    backwardlik::Array{Float64,3}# size: k, net.numNodes, number of displayed trees
end
const SSM = StatisticalSubstitutionModel{T} where T

StatsBase.loglikelihood(obj::SSM) = obj.loglik
StatsBase.islinear(SSM) = false
StatsBase.dof(obj::SSM) = nparams(obj.model)
function Base.show(io::IO, obj::SSM)
    disp =  "$(typeof(obj)):\n"
    disp *= string(obj.model)
    disp *= "$(obj.ntraits) traits, $(length(obj.trait)) species, "
    disp *= "on a network with $(obj.net.numHybrids) reticulations"
    if !ismissing(obj.loglik)
        disp *= "\nlog-likelihood: $(round(obj.loglik, digits=5))"
    end
    print(io, disp)
end
# nobs: ntraits * nspecies, minus any missing, but ignores correlation between species
# fixit: extend the StatsBase methods
# coef, coefnames, coeftable, confint,
# deviance (not from loglik in StatsBase), nulldeviance (from intercept only. here?)
# nullloglikelihood (from intercept only. here?)
# score (gradient of loglik)
# informationmatrix: Fisher by default, observed info matrix if expected=false
# stderror, vcov, params
# weights

"""
    fitDiscrete(net, model, tipdata)
    fitDiscrete(net, model, species, traits)

Calculate the maximum likelihood (ML) score of a network given
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
  Here, trait labels should be as they appear in `model.label`.
- `species`: vector of strings, and `traits`: DataFrame of traits,
  with rows in the order corresponding to the order of species names.
  Again, trait labels should be as they appear in `model.label`.
  All traits are assumed to follow the same model, with same parameters.

Optional arguments (default):
- `fixedparam` (false): should model rate parameters be fixed, or should they be optimized?
- `NLoptMethod` (`:LN_COBYLA`, derivative-free) for the optimization algorithm.
  For other options, see the
  [NLopt](https://nlopt.readthedocs.io/en/stable/NLopt_Algorithms/).
- tolerance values to control when the optimization is stopped:
  `ftolRel` (1e-12), `ftolAbs` (1e-10) on the likelihood, and
  `xtolRel` (1e-10), `xtolAbs` (1e-10) on the model parameters.
- `verbose` (false): if true, more information is output.

# examples:

```jldoctest fitDiscrete
julia> net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> m1 = BinaryTraitSubstitutionModel([0.1, 0.1], ["lo", "hi"]);

julia> using DataFrames

julia> dat = DataFrame(species=["C","A","B","D"], trait=["hi","lo","lo","hi"]);

julia> fit1 = fitDiscrete(net, m1, dat; fixedparam=true)
PhyloNetworks.StatisticalSubstitutionModel{String}:
Binary Trait Substitution Model:
rate lo→hi α=0.1
rate hi→lo β=0.1
1 traits, 4 species, on a network with 1 reticulations
log-likelihood: -3.10754

julia> PhyloNetworks.fit!(fit1; fixedparam=false)
PhyloNetworks.StatisticalSubstitutionModel{String}:
Binary Trait Substitution Model:
rate lo→hi α=0.27222
rate hi→lo β=0.34981
1 traits, 4 species, on a network with 1 reticulations
log-likelihood: -2.7277

julia> tips = Dict("A" => "lo", "B" => "lo", "C" => "hi", "D" => "hi");

julia> fit2 = fitDiscrete(net, m1, tips; xtolRel=1e-16, xtolAbs=1e-16, ftolRel=1e-16)
PhyloNetworks.StatisticalSubstitutionModel{String}:
Binary Trait Substitution Model:
rate lo→hi α=0.27222
rate hi→lo β=0.34981
1 traits, 4 species, on a network with 1 reticulations
log-likelihood: -2.7277
```

Note that a copy of the network is stored in the fitted object,
but the internal representation of the network may be different in
`fit1.net` and in the original network `net`:

```jldoctest fitDiscrete
julia> [n.number for n in fit2.net.node]
9-element Array{Int64,1}:
 1
 2
 9
 8
 3
 7
 6
 4
 5

julia> [n.number for n in net.node]
9-element Array{Int64,1}:
  1
  2
  3
 -4
  4
 -6
 -3
  5
 -2
```
"""
function fitDiscrete(net, model, tips::Dict; kwargs...)
    species = String[]
    dat = Vector{Int}[] # indices of trait labels
    for (k,v) in tips
        !ismissing(v) || continue
        push!(species, k)
        vi = findfirst(isequal(v), model.label) # value index in model labels
        vi !== nothing || error("trait $v not found in model")
        push!(dat, [vi])
    end
    o, net = check_matchtaxonnames!(species, dat, net)
    # dat[o] would make a shallow copy only
    StatsBase.fit(StatisticalSubstitutionModel, net, model, view(dat, o); kwargs...)
end

function fitDiscrete(net, model, dat::DataFrame; kwargs...)
    i = findfirst(isequal(:taxon), DataFrames.names(dat))
    if i===nothing i = findfirst(isequal(:species), DataFrames.names(dat)); end
    if i===nothing i=1; end # first column if no column named "taxon" or "species"
    j = findfirst(isequal(:trait), DataFrames.names(dat))
    if j===nothing j=2; end
    if i==j
        error("""expecting taxon names in column 'taxon', or 'species' or column 1,
              and trait values in column 'trait' or column 2.""")
    end
    species = copy(dat[i])    # modified in place later
    dat = traitlabels2indices(dat[j], model)   # vec of vec, indices
    o, net = check_matchtaxonnames!(species, dat, net)
    StatsBase.fit(StatisticalSubstitutionModel, net, model, view(dat, o); kwargs...)
end

function fitDiscrete(net, model, species::Array{String}, dat::DataFrame; kwargs...)
    dat2 = traitlabels2indices(dat, model) # vec of vec, indices
    o, net = check_matchtaxonnames!(copy(species), dat2, net)
    StatsBase.fit(StatisticalSubstitutionModel, net, model, view(dat2, o); kwargs...)
end

"""
    fit(StatisticalSubstitutionModel, net, model, traits; kwargs...)
    fit!(StatisticalSubstitutionModel; kwargs...)

Internal function called by [`fitDiscrete`](@ref): with same key word arguments `kwargs`.
But dangerous: `traits` should be a vector of vectors as for [`fitDiscrete`](@ref)
**but** here `traits` need to contain the *indices* of trait values corresponding
to the indices in `model.label`, and species should appear in `traits` in the
order corresponding to the node numbers in `net`.
See [`traitlabels2indices`](@ref) to convert trait labels to trait indices.

**Warning**: does *not* perform checks. [`fitDiscrete`](@ref) calls this function
after doing checks, preordering nodes in the network, making sure nodes have
consecutive numbers, species are matched between data and network etc.
"""
function StatsBase.fit(::Type{SSM}, net::HybridNetwork, model::TraitSubstitutionModel,
    trait::AbstractVector; kwargs...)
    T = eltype(model.label)
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
    k = nStates(model)
    # fixit: use SharedArray's below to parallelize things
    logtrans   = zeros(Float64, k,k,length(net.edge))
    forwardlik = zeros(Float64, k, nnodes,           ntrees)
    directlik  = zeros(Float64, k, length(net.edge), ntrees)
    backwardlik= zeros(Float64, k, nnodes,           ntrees)
    postltw    = Vector{Float64}(undef, ntrees)
    # create new model object then fit:
    fit!(StatisticalSubstitutionModel{T}(deepcopy(model),
            net, trait, length(trait[1]), missing,
            logtrans, 1, trees, priorltw, postltw, forwardlik, directlik, backwardlik);
        kwargs...)
end

function fit!(obj::SSM; fixedparam=false::Bool, verbose=false::Bool,
      NLoptMethod=:LD_MMA::Symbol, ftolRel=fRelBL::Float64, ftolAbs=fAbsBL::Float64,
      xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64)

    if fixedparam
        discrete_corelikelihood!(obj)
        verbose && println("loglik = $(loglikelihood(obj)) under fixed parameters, no optimization")
    else
        # set-up optimization object
        NLoptMethod=:LN_COBYLA # no gradient
        # :LN_COBYLA for (non)linear constraits, :LN_BOBYQA for bound constraints
        npar = nparams(obj.model)
        opt = NLopt.Opt(NLoptMethod, npar)
        NLopt.ftol_rel!(opt,ftolRel) # relative criterion
        NLopt.ftol_abs!(opt,ftolAbs) # absolute criterion
        NLopt.xtol_rel!(opt,xtolRel)
        NLopt.xtol_abs!(opt,xtolAbs)
        NLopt.maxeval!(opt,1000) # max number of iterations
        # NLopt.maxtime!(opt, t::Real)
        NLopt.lower_bounds!(opt, zeros(Float64, npar))
        # fixit: set upper bound depending on branch lengths in network?
        counter = [0]
        function loglikfun(x::Vector{Float64}, grad::Vector{Float64}) # modifies obj
            counter[1] += 1
            obj.model.rate[:] = x
            res = discrete_corelikelihood!(obj)
            verbose && println("loglik: $res, model rates: $x")
            length(grad) == 0 || error("gradient not implemented")
            return res
        end
        NLopt.max_objective!(opt, loglikfun)
        fmax, xmax, ret = NLopt.optimize(opt, obj.model.rate) # optimization here!
        verbose && println("got $(round(fmax, digits=5)) at $(round.(xmax, digits=5)) after $(counter[1]) iterations (return code $(ret))")
    end
    # return fmax,xmax,ret
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
        traitrange = 1:obj.ntraits
    elseif isinteger(whichtrait) && whichtrait > 0 && whichtrait <= obj.ntraits
        obj.activetrait = whichtrait
        traitrange = range(whichtrait, length=1)
    elseif whichtrait == :active
        traitrange = range(obj.activetrait, length=1)
    else
        error("'whichtrait' should be :all or :active or an integer in the correct range")
    end
    for edge in obj.net.edge # update logtrans: same for all displayed trees, all traits
        obj.logtrans[:,:,edge.number] = log.(P(obj.model, edge.length)) # element-wise
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
    # fixit: function to get these posterior probabilities (just take exp.)
    return res
end

@doc (@doc discrete_corelikelihood!) discrete_corelikelihood_tree!
function discrete_corelikelihood_tree!(obj::SSM, t::Integer, traitrange::AbstractArray)
    tree = obj.displayedtree[t]
    # @info "tree: $(writeTopology(tree))"
    forwardlik = view(obj.forwardlik, :,:,t)
    directlik  = view(obj.directlik,  :,:,t)
    k = nStates(obj.model)   # also = size(logtrans,1)
    fullloglik = 0.0
    for ci in traitrange     # ci = character index
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
            logprior = [-log(k) for i in 1:k] # uniform prior at root
            loglik = logsumexp(logprior + view(forwardlik, :,nnum)) # log P{data for ci | tree t}
            fullloglik += loglik # add loglik of character ci
            break # out of loop over nodes
        end
        # if we keep going: n is not the root
        # calculate direct likelihood on the parent edge of n
        for e in n.edge
            if n == getChild(e)
                lt = view(obj.logtrans, :,:,e.number)
                for i in 1:k # state at parent node
                    directlik[i,e.number] = logsumexp(view(lt,i,:) + view(forwardlik,:,nnum))
                end
                break # we visited the parent edge: break out of for loop
            end
        end
      end # of loop over nodes
    end # of loop over traits
    obj.postltw[t] = fullloglik
    return fullloglik
end

"""
    traitlabels2indices(data, model::TraitSubstitutionModel)

Check that the character states in `data` are compatible with (i.e. subset of)
the trait labels in `model`. All columns are used.
`data` can be a DataFrame or a Matrix (multiple traits), or a Vector (one trait).

Return a vector of vectors (one per species) with integer entries,
where each state (label) is replaced by its index in `model`.
"""
function traitlabels2indices(data::AbstractVector, model::TraitSubstitutionModel)
    A = Vector{Vector{Union{Missings.Missing,Int}}}(undef, 0) # indices of trait labels
    labs = model.label
    for l in data
        vi = missing
        if !ismissing(l)
            vi = findfirst(isequal(l), model.label) # value index in model labels
            vi !== nothing || error("trait $l not found in model")
        end
        push!(A, [vi])
    end
    return A
end
function traitlabels2indices(data::Union{AbstractMatrix,DataFrame},
                             model::TraitSubstitutionModel)
    A = Vector{Vector{Union{Missings.Missing,Int}}}(undef, 0) # indices of trait labels
    labs = model.label
    ntraits = size(data,2)
    for i in 1:size(data,1) # go row by row
        V = Vector{Union{Missings.Missing,Int}}(undef, ntraits)
        for j in 1:ntraits
            vi = missing # value index
            @inbounds l = data[i,j] # value label
            if !ismissing(l)
                vi = findfirst(isequal(l), labs)
                vi !== nothing || error("trait $l not found in model")
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
Used by [`fitDiscrete`](@ref) to build a new [`StatisticalSubstitutionModel`](@ref).
"""
function check_matchtaxonnames!(species::AbstractVector, dat::AbstractVector, net::HybridNetwork)
    # 1. basic checks for dimensions and types
    eltt = eltype(dat)
    @assert eltt <: AbstractVector "traits should be a vector of vectors"
    @assert Missings.T(eltype(eltt)) <: Integer "traits should be integers (label indices)"
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
number `trait`, or for the active trait if `trait` is unspecified: `obj.activetrait`.
The parameters of the [`StatisticalSubstitutionModel`](@ref) object `obj`
must first be fitted using [`fitDiscrete`](@ref), and ancestral state reconstruction
is conditional on the estimated parameters. If these parameters were estimated
using all traits, they are used as is to do ancestral state reconstruction of the
particular `trait` of interest.

**output**: data frame with a first column for the node numbers, a second column for
the node labels, and a column for each possible state: the entries in these columns
give the marginal probability that a given node has a given state.

warnings

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

julia> fit1 = fitDiscrete(net, m1, dat);

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
│ 9   │ 9          │ #H1       │ 0.782777 │ 0.217223 │

julia> round.(exp.(fit1.postltw), digits=6) # marginal (posterior) probability that the trait evolved on each displayed tree
2-element Array{Float64,1}:
 0.919831 
 0.080169

julia> using PhyloPlots

julia> plot(fit1.net, :R, nodeLabel = asr[[:nodenumber, :lo]], tipOffset=0.2); # pp for "lo" state
```
"""
ancestralStateReconstruction(obj::SSM) = ancestralStateReconstruction(obj, obj.activetrait)
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
    k = nStates(obj.model)
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
    dat = DataFrame(transpose(res), Symbol.(obj.model.label))
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
    k = nStates(obj.model)
    bkwtmp = Vector{Float64}(undef, k) # to hold bkw lik without parent edge transition
    logprior = [-log(k) for i in 1:k]
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
