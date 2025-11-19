"""
    getnodeages(net)

vector of node ages in pre-order, as in `vec_node`.

*Warnings*: `net` is assumed to
- have been preordered before (to calculate `vec_node`)
- be time-consistent (all paths to the root to a given hybrid have the same length)
- be ultrametric (all leaves have the same age: 0)
"""
function getnodeages(net::HybridNetwork)
    x = Vector{Float64}(undef, length(net.vec_node))
    for i in reverse(1:length(net.vec_node)) # post-order
        n = net.vec_node[i]
        if n.leaf
            x[i] = 0.0
            continue
        end
        for e in n.edge
            if getparent(e) == n # n parent of e
                childnode = getchild(e)
                childIndex = getIndex(childnode, net.vec_node)
                x[i] = x[childIndex] + e.length
                break # use 1st child, ignores all others
            end
        end
    end
    return x
end

"""
    pairwisetaxondistancematrix(net; type=:average, keepInternal=false,
                                checkpreorder=true, nodeAges=[])
    pairwisetaxondistancematrix!(M, net, type, nodeAges)

Return the matrix `M` of pairwise distances between nodes in the network:
- between all nodes (internal and leaves) if `keepInternal=true`,
  in which case the nodes are listed in `M` in the
  order in which they appear in `net.vec_node`
- between taxa only otherwise, in which case the nodes are listed
  in `M` in the order in which they appear in `tiplabels(net)`
  (i.e. same order as in `net.leaf`)

The second form modifies `M` in place, assuming all nodes.

Between any two nodes n1 and n2, we consider all up-down paths
between n1 and n2: paths starting from n1, going up to a common ancestor
then back down to n2. The distance between n1 and n2 associated with an
up-down path is the path length, the sum of its edge lengths. Each path also
has an inheritance weight γ: the product of γs of all edges on the path,
representing the proportion of the genome that took this path.

By default, the distance `type` is `:average` to get a weighted average
distance. This is the average of the up-down paths' lengths,
weighted by the paths' inheritance γs.
It measures the average genetic distance across the genome,
if branch lengths are in substitutions/site.

Another distance type is `:maximum`. This is the maximum length
among all up-down paths between two nodes. For this maximum distance,
the inheritance probabilities γ (or path weights) have no impact.

Other optional arguments:

- `checkpreorder`: if true, `net.vec_node` is updated to get a
  topological ordering of nodes.
- `nodeAges`: if not provided, i.e. empty vector, the network is *not* modified.  
  If provided and non-empty, `nodeAges` should list node ages in the
  pre-order in which nodes are listed in `vec_node` (including leaves),
  and **edge lengths** in `net` **are modified** accordingly.

Providing node ages hence makes the network time consistent: such that
all paths from the root to a given hybrid node have the same length.
If node ages are not provided, the network need not be time consistent.
"""
function pairwisetaxondistancematrix(
    net::HybridNetwork;
    type::Symbol=:average,
    keepInternal::Bool=false,
    checkpreorder::Bool=true,
    nodeAges::Vector{Float64}=Float64[]
)
    net.isrooted || error("net needs to be rooted for preorder recursion")
    if checkpreorder
        preorder!(net)
    end
    nnodes = net.numnodes
    M = zeros(Float64,nnodes,nnodes)
    isempty(nodeAges) || length(nodeAges) == nnodes ||
        error("there should be $nnodes node ages")
    pairwisetaxondistancematrix!(M, net, type, nodeAges)
    if !keepInternal
        M = getTipSubmatrix(M, net)
    end
    return M
end
"""
    getTipSubmatrix(M, net; indexation=:both)

Extract submatrix of M, with rows and/or columns corresponding to
tips in the network, ordered like in `net.leaf`.
In M, rows and/or columns are assumed ordered as in `net.vec_node`.

indexation: one of `:rows`, `:cols` or `:both`: are nodes numbers indexed
in the matrix by rows, by columns, or both? Subsetting is taken accordingly.

"""
function getTipSubmatrix(M::Matrix, net::HybridNetwork; indexation=:both)
    nodenames = [n.name for n in net.vec_node]
    tipind = Int[]
    for l in tiplabels(net)
        push!(tipind, findfirst(isequal(l), nodenames))
    end
    if indexation == :both
        return M[tipind, tipind]
    elseif indexation == :rows
        return M[tipind, :]
    elseif indexation == :cols
        return M[:, tipind]
    else
        error("indexation must be one of :both :rows or :cols")
    end
end

function pairwisetaxondistancematrix!(
    M::Matrix{Float64},
    net::HybridNetwork,
    type::Symbol,
    nodeAges
)
    updatehybrid = get(updatehybrid_distance_funs, type, nothing)
    if isnothing(updatehybrid)
        error("Unknown distance type $type. Supported types are :maximum or :average.")
    end
    traversal_preorder!(net.vec_node, M, # updates M in place
            traversalupdate_default!, # does nothing
            updatetree_pairwisedistancematrix!,
            updatehybrid,
            nodeAges)
end

function updatetree_pairwisedistancematrix!(
    V::Matrix,
    i::Int,
    parindx::Int, # index of parent node
    edge::Edge,
    nodeages, # assumed pre-ordered, as in vec_node
)
    if !isempty(nodeages)
        edge.length = nodeages[parindx] - nodeages[i]
    end
    for j in 1:(i-1)
        V[i,j] = V[parindx,j] + edge.length
        V[j,i] = V[i,j]
    end
    # V[i,i] = zero(eltype(V)) # not needed bc: initialized to 0 and never changes
    return true
end

function average_updatehybrid_pairwisedistancematrix!(
    V::Matrix,
    i::Int,
    parindx::AbstractVector{Int},
    paredge::AbstractVector{Edge},
    nodeages, # node ages should be pre-ordered
)
    if !isempty(nodeages)
        for (pi,pe) in zip(parindx, paredge)
            pe.length = nodeages[pi] - nodeages[i]
        end
    end
    for j in 1:(i-1)
        V[i,j] = zero(eltype(V))
        #= V[i,j] initialized at 0 by pairwisetaxondistancematrix but
        re-initialization needed during optimization in
        calibratefrompairwisedistances! =#
        for (pi,pe) in zip(parindx, paredge)
            V[i,j] += pe.gamma * (V[pi,j] + pe.length)
        end
        V[j,i] = V[i,j]
    end
    return true
end

function max_updatehybrid_pairwisedistancematrix!(
    V::Matrix,
    i::Int,
    parindx::AbstractVector{Int},
    paredge::AbstractVector{Edge},
    nodeages,
)
    if !isempty(nodeages)
        for (pi,pe) in zip(parindx, paredge)
            pe.length = nodeages[pi] - nodeages[i]
        end
    end
    for j in 1:(i-1)
        V[i,j] = maximum(((pi, pe),) -> V[pi, j] + pe.length, zip(parindx, paredge))
        V[j,i] = V[i,j]
    end
    return true
end

const updatehybrid_distance_funs = Dict(
    :maximum => max_updatehybrid_pairwisedistancematrix!,
    :average => average_updatehybrid_pairwisedistancematrix!
)

"""
    pairwisetaxondistance_gradient(net; checkEdgeNumber=true, nodeAges=[])

3-dim array: gradient of pairwise distances between all nodes.
(internal and leaves); gradient with respect to edge lengths
if `nodeAges` is empty; with respect to node ages otherwise.
Assume correct `net.vec_node` (preorder).  
This gradient depends on the network's topology and γ's only,
not on branch lengths or node ages (distances are linear in either).

WARNING: edge numbers need to range between 1 and #edges.
"""
function pairwisetaxondistance_gradient(
    net::HybridNetwork;
    checkEdgeNumber::Bool=true,
    nodeAges::Vector{Float64}=Float64[]
)
    if checkEdgeNumber
      sort([e.number for e in net.edge]) == collect(1:net.numedges) ||
        error("edge numbers must range between 1 and #edges")
    end
    n = (length(nodeAges)==0 ? net.numedges : net.numnodes)
    M = zeros(Float64, net.numnodes, net.numnodes, n)
    traversal_preorder!(net.vec_node, M,
            traversalupdate_default!,  # does nothing
            updatetree_pairwisedistancegrad!,
            updatehybrid_pairwisedistancegrad!,
            nodeAges) # nodeAges assumed pre-ordered, like vec_node
    return M
end
function updatetree_pairwisedistancegrad!(
    V::Array{Float64,3},
    i::Int,
    parindx::Int, # index of parent node
    edge::Edge,
    nodeages,  # assumed pre-ordered
)
    emptyages = isempty(nodeages)
    for j in 1:(i-1)
        if emptyages  # d/d(edge length)
            V[i,j,edge.number] = 1.0
        else          # d/d(node age)
            V[i,j,parindx] = 1.0
            V[i,j,i] = -1.0
        end
        for k in 1:size(V)[3] # brute force...
            V[i,j,k] += V[parindx,j,k]
            V[j,i,k] = V[i,j,k]
        end
    end
    # V[i,i,k] initialized to 0.0 already
    return true
end
function updatehybrid_pairwisedistancegrad!(
    V::Array{Float64,3},
    i::Int,
    parindx::AbstractVector{Int},
    paredge::AbstractVector{Edge},
    nodeages,
)
    emptyages = isempty(nodeages)
    for j in 1:(i-1)
        if emptyages # d/d(edge length)
            for pe in paredge
                V[i,j,pe.number] = pe.gamma
            end
        else         # d/d(node age)
            for (pi,pe) in zip(parindx, paredge)
                V[i,j,pi] = pe.gamma
            end
            V[i,j,i] = - 1.0 # because γ1+γ2 = 1
        end
        for k in 1:size(V)[3] # 1:size(V,3)
            for (pi,pe) in zip(parindx, paredge)
                V[i,j,k] += pe.gamma * V[pi,j,k]
            end
          V[j,i,k] = V[i,j,k]
        end
    end
    return true
end


"""
    calibratefrompairwisedistances!(net, distances::Matrix{Float64},
        taxon_names::Vector{<:AbstractString})

Calibrate the network to match (as best as possible) input
pairwise distances between taxa, such as observed from sequence data.
`taxon_names` should provide the list of taxa, in the same order
in which they they are considered in the `distances` matrix.
The optimization criterion is the sum of squares between the
observed distances, and the distances from the network
(weighted average of tree distances, weighted by γ's).
The network's edge lengths are modified.

Warning: for many networks, mutiple calibrations can fit the pairwise
distance data equally well (lack of identifiability).
This function will output *one* of these equally good calibrations.

optional arguments (default):
- checkpreorder (true)
- forceMinorLength0 (false) to force minor hybrid edges to have a length of 0
- ultrametric (true) to force the network to be
  * time-consistent: all paths from the root to a given node must have the same
    length, so the age of this node is well-defined, and
  * ultrametric: all tips are at the same distance from the root, so have the same age.
- NLoptMethod (`:LD_MMA`) for the optimization algorithm.
  Other options include `:LN_COBYLA` (derivative-free); see NLopt package.
- tolerance values to control when the optimization is stopped:
  ftolRel (1e-12), ftolAbs (1e-10) on the criterion, and
  xtolRel (1e-10), xtolAbs (1e-10) on branch lengths / divergence times.
- verbose (false)
"""
function calibratefrompairwisedistances!(
    net::HybridNetwork,
    D::Array{Float64,2},
    taxNames::Vector{<:AbstractString};
    checkpreorder::Bool=true,
    forceMinorLength0::Bool=false,
    verbose::Bool=false,
    ultrametric::Bool=true,
    NLoptMethod::Symbol=:LD_MMA,
    ftolRel::Float64=fRelBL,
    ftolAbs::Float64=fAbsBL,
    xtolRel::Float64=xRelBL,
    xtolAbs::Float64=xAbsBL
)
    checkpreorder && preorder!(net)
    # fixit: remove root node if of degree 2, and if ultrametric=false
    defaultedgelength = median(D)/(length(net.edge)/2)
    for e in net.edge
        if e.length == -1.0 e.length=defaultedgelength; end
        # get smarter starting values: NJ? fast dating?
    end
    if ultrametric # get all node ages in pre-order
        na = getnodeages(net)
    else na = Float64[]; end
    # get number and indices of edge/nodes to be optimized
    if forceMinorLength0 && !ultrametric
        parind = filter(i -> net.edge[i].ismajor, 1:net.numedges)
        nparams = length(parind) # edges to be optimized: indices
        par = [e.length for e in net.edge][parind] # and lengths
        for i in 1:net.numedges
            if !net.edge[i].ismajor net.edge[i].length=0.0; end
        end
    elseif !forceMinorLength0 && !ultrametric
        nparams = length(net.edge) # number of parameters to optimize
        par = [e.length for e in net.edge]
        parind = 1:nparams
    elseif !forceMinorLength0 && ultrametric
        nparams = net.numnodes - net.numtaxa # internal nodes to be optimized
        parind = filter(i -> !net.vec_node[i].leaf, 1:net.numnodes)
        par = na[parind]
    else # forceMinorLength0 && ultrametric
        nparams = net.numnodes - net.numtaxa - net.numhybrids
        parind = filter(i -> !(net.vec_node[i].leaf || net.vec_node[i].hybrid),
                        1:net.numnodes)
        par = na[parind]
        hybInd = filter(i -> net.vec_node[i].hybrid, 1:net.numnodes)
        hybParentInd = Int[] # index in vec_node of minor parent
        hybGParentI = Int[] # index in 1:nparams of minor (grand-)parent in param list
        for i in hybInd
            n = net.vec_node[i]
            p = getparentminor(n)
            pi = findfirst(n -> n===p, net.vec_node)
            push!(hybParentInd, pi)
            pii = findfirst(isequal(pi), parind)
            while pii===nothing # in case minor parent of n is also hybrid node
                p = getparentminor(p)
                pi = findfirst(n -> n===p, net.vec_node)
                pii = findfirst(isequal(pi), parind)
            end
            push!(hybGParentI, pii)
        end
    end
    # initialize M=dist b/w all nodes, G=gradient (constant)
    M = pairwisetaxondistancematrix(net, keepInternal=true,
            checkpreorder=false, nodeAges=na)
    if !ultrametric && sort([e.number for e in net.edge]) != collect(1:net.numedges)
        for i in 1:net.numedges # renumber edges, needed for G
            net.edge[i].number = i
        end
    end # G assumes edges numbered 1:#edges, if optim edge lengths
    G = pairwisetaxondistance_gradient(net, checkEdgeNumber=false, nodeAges=na) .* 2
    # match order of leaves in input matrix, versus pre-order
    nodenames = [n.name for n in net.vec_node] # pre-ordered
    ntax = length(taxNames)
    tipind = Int[] # pre-order index for leaf #i in dna distances
    for l in taxNames
        i = findfirst(isequal(l), nodenames)
        i !== nothing || error("taxon $l not found in network")
        push!(tipind, i)
    end
    # contraints: to force a parent to be older than its child
    if ultrametric
      numConstraints = length(parind) -1 + net.numhybrids
      # calculate indices in param list of child & (grand-)parent once
      chii = Int[] # non-root internal node, can be repeated: once per constraint
      anii = Int[] # closest ancestor in param list
      ci = 1 # index of constraint
      for i in 2:length(net.vec_node) # 1=root, can skip
        n = net.vec_node[i]
        if n.leaf continue; end # node ages already bounded by 0
        if n.hybrid && forceMinorLength0          # get index in param list of
          nii = hybGParentI[findfirst(isequal(i), hybInd)] # minor grand-parent (same age)
        else
          nii = findfirst(isequal(i), parind)
        end
        for e in n.edge
          if getchild(e) == n # n child of e
            p = getparent(e)  # parent of n
            if forceMinorLength0 && n.hybrid && !e.ismajor
                continue; end # p and n at same age already
            pi = findfirst(isequal(p.number), [no.number for no in net.vec_node])
            if forceMinorLength0 && p.hybrid
              pii = hybGParentI[findfirst(isequal(pi), hybInd)]
            else
              pii = findfirst(isequal(pi), parind)
            end
            push!(chii, nii)
            push!(anii, pii)
            verbose && println("node $(net.vec_node[parind[nii]].number) constrained by age of parent $(net.vec_node[parind[pii]].number)")
            ci += 1
          end
        end
      end
      length(chii) == numConstraints ||
        error("incorrect number of node age constraints: $numConstraints")
      function ageConstraints(result, nodeage, grad)
        if length(grad) > 0 # grad: nparams x nConstraints: ∂cj/∂xi = grad[i,j]
            fill!(grad, 0.0)
        end
        for j in 1:numConstraints
          nii = chii[j]; pii = anii[j]
          result[j] = nodeage[nii] - nodeage[pii] # jth constraint: cj ≤ 0
          if length(grad) > 0 # nparams x nConstraints
            grad[nii, j] =  1.
            grad[pii, j] = -1.
          end
        end
      end
    end
    opt = NLopt.Opt(NLoptMethod,nparams) # :LD_MMA to use gradient
    # :LN_COBYLA for (non)linear constraits, :LN_BOBYQA for bound constraints
    NLopt.maxeval!(opt,1000) # max iterations
    NLopt.ftol_rel!(opt,ftolRel)
    NLopt.ftol_abs!(opt,ftolAbs)
    NLopt.xtol_rel!(opt,xtolRel)
    NLopt.xtol_abs!(opt,xtolAbs)
    # NLopt.maxtime!(opt, t::Real)
    NLopt.lower_bounds!(opt, zeros(nparams))
    if ultrametric
      NLopt.inequality_constraint!(opt,ageConstraints,fill(0.0,numConstraints))
    end
    counter = [0]
    function obj(x::Vector{Float64}, grad::Vector{Float64})
        verbose && println("mismatch objective, BL = $(x)")
        counter[1] += 1
        # update edge lengths or node ages
        if ultrametric # update na using x, in place
            for i in 1:nparams # na=0 at leaves already
                na[parind[i]] = x[i]
            end
            if forceMinorLength0 # set hybrid age to minor parent age
                for i in 1:net.numhybrids
                    na[hybInd[i]] = na[hybParentInd[i]] # pre-order important
                end
            end
        else # not ultrametric: optimize branch lengths
            for i in 1:nparams # update network
                net.edge[parind[i]].length = x[i]
            end
        end
        # update distances in M, in place
        pairwisetaxondistancematrix!(M,net,:average,na)
        ss = 0.0 # sum of squares between M and observed distances
        for i in 2:ntax; for j in 1:(i-1)
            ss += (M[tipind[i],tipind[j]]-D[i,j])^2
        end; end
        if length(grad) > 0 # sum_ij 2 * dM_ij/dx_t * (M_ij-D_ij)
          for t in 1:nparams grad[t] = 0.0; end;
          for i in 2:ntax; for j in 1:(i-1);
            for t in 1:nparams
              G[tipind[i],tipind[j],parind[t]] != 0.0 || continue
              grad[t] += G[tipind[i],tipind[j],parind[t]] *
                        (M[tipind[i],tipind[j]]-D[i,j])
            end
            if ultrametric && forceMinorLength0
              for i in 1:net.numhybrids # na[hybrid] set to na[minor parent]
                grad[hybGParentI[i]] += G[tipind[i],tipind[j],hybInd[i]] *
                        (M[tipind[i],tipind[j]]-D[i,j])
              end
            end
          end; end
        end
        return ss
    end
    NLopt.min_objective!(opt,obj)
    fmin, xmin, ret = NLopt.optimize(opt,par) # optimization here!
    verbose && println("got $(round(fmin, digits=5)) at $(round.(xmin, digits=5)) after $(counter[1]) iterations (return code $(ret))")
    return fmin,xmin,ret
end

"""
    hammingdistancematrix(
        trait::AbstractVector{Vector{Union{Missings.Missing,Int}}},
        traitweight::AbstractVector{<:Real}=ones(length(trait[1]));
        scaled::Bool=true
    )

Matrix of pairwise Hamming distances between taxa:
proportion (or number) of traits at which 2 taxa differ, possibly weighted.

Arguments:
- `trait`: vector of vectors of integers.
  `trait[i]` has the data for taxon number `i`. All vectors `trait[i]` should
  be of the same length: same number of traits across taxa.
  `trait[i][j] = k` means that taxon `i` has state index `k` for trait `j`.
- `traitweight[k]` gives the weight of trait `k` (default: all 1s)
- `scaled`: if `false`, the output distance is the number (or total weight) of
  traits at which 2 taxa differ. If `scaled` is `true`, the output distance is
  a proportion: the un-scaled distance divided by the number of traits
  (or sum of trait weights).

Missing data:
For each pair, the Hamming distances is calculated ignoring any trait in
which one of 2 values is missing. The total number of differences is then
divided by the total weight of all sites, ignoring that some of them may
have been missing for the pair. This is an inexact rescaling of the
Hamming distance, assuming a small proportion of missing values for each pair.
"""
function hammingdistancematrix(
    trait::AbstractVector{Vector{T}},
    traitweight::AbstractVector{<:Real}=ones(length(trait[1]));
    scaled::Bool=true
) where T
    ntax = length(trait)
    M = zeros(Float64, ntax, ntax)
    ncols = length(trait[1])
    all(length.(trait) .== ncols) ||
        error("all taxa should have the same number of traits")
    length(traitweight) == ncols ||
        error("$(length(traitweight)) weights but $ncols columns (traits) in the data")
    # count pairwise differences, then multiply by pattern weight
    for i in 2:ntax
        species1 = trait[i]
        for j in 1:(i-1)
            species2 = trait[j]
            for col in 1:ncols
                if !(ismissing(species1[col]) || ismissing(species2[col])) &&
                    (species1[col] != species2[col])
                    M[i, j] += traitweight[col]
                end
            end
            M[j,i] = M[i,j]
        end
    end
    if scaled # normalize to get proportion
        M ./= sum(traitweight)
    end
    return M
end

"""
    distancecorrection_JC!(M::AbstractMatrix, nstates::Integer; scalar=1.01)

Jukes-Cantor (JC) corrected distance matrix, to estimate pairwise evolutionary
distances from scaled Hamming distances (proportion of sites at which the pair
differs). `M` should contain non-negative values, such as Hamming distances
from [`hammingdistancematrix`](@ref).
The JC correction is defined as:
`- 0.75 log(1 - M/0.75)` if there are 4 states, and more generally
`- dmax log(1 - M/dmax)` where `dmax = (n-1/n)` on `n` states.

Theoretically, scaled Hamming distances from the Jukes-Cantor model are < 0.75
on four states, or < `(n-1)/n` on `n` states.

If `M` has some values ≥ 0.75 (or `dmax` more generally), then `M` does not fit
the Jukes-Cantor model, and `log(1 - M/0.75)` would give negative values.
To avoid negative corrected distances, and numerical issues if some M values
are close to 0.75, the JC correction here applies a following stronger
normalization. `M` is modified in place to contain:

`d_JC = - 0.75 log(1 - M/max{0.75, m*1.01})`

for 4 states (replace 0.75 by (n-1)/n for n states),
where `m` is the maximum observed distance in `M`.

The `scalar=1.01` value above is a keyword argument, to adjust the strength
of normalization.
"""
function distancecorrection_JC!(
    M::AbstractMatrix,
    nstates::Integer;
    scalar::Real=1.01
)
    maxdist = (nstates-1)/nstates
    M ./= max(maxdist, maximum(M)*scalar) # values ≤ 0.9901: log(1-M) well defined
    M .= - maxdist .* log.( 1.0 .- M)     # 0 -> -0, not nice
    for I in eachindex(M)
        if M[I] == -0
            M[I] = 0
        end
    end
    return M
end

"""
    startingBL!(net::HybridNetwork,
                trait::AbstractVector{Vector{Union{Missings.Missing,Int}}},
                siteweight::AbstractVector{<:Real}=ones(length(trait[1])))

Calibrate branch lengths in `net` by minimizing the mean squared error
between the JC-adjusted pairwise distance between taxa, and network-predicted
pairwise distances, using [`calibratefrompairwisedistances!`](@ref).
The network is *not* forced to be time-consistent nor ultrametric.
To avoid one source of non-identifiability, the network is "zipped" by
forcing minor hybrid edges to have length 0.
Finally, any edge length smaller than 1.0e-10 is reset to 0.0001,
to avoid the 0 boundary.

Traits and weights:
- `trait[i]` is for leaf with `node.number = i` in `net`, and
  `trait[i][j] = k` means that leaf number `i` has state index `k` for trait `j`.
  These indices are those used in a substitution model
  by [PhyloTraits](https://juliaphylo.github.io/PhyloTraits.jl/stable/man/simulate_discrete/#Discrete-trait-simulation):
  kth value of `PhyloTraits.getlabels(model)`.
- `siteweight[k]` gives the weight of site (or site pattern) `k` (default: all 1s).

See [`hammingdistancematrix`](@ref) for how missing data are treated to
calculate the Hamming distance.  
See [`distancecorrection_JC!`](@ref) for the Jukes-Cantor correction.
For this correction, the number of traits is taken as the maximum trait
state index (denoted `k` above).
"""
function startingBL!(
    net::HybridNetwork,
    trait::AbstractVector{Vector{Union{Missings.Missing,Int}}},
    siteweight::AbstractVector{<:Real}=ones(Float64,length(trait[1]))
)
    dhat = hammingdistancematrix(trait, siteweight)
    net.numtaxa == size(dhat,1) ||
        error("$(net.numtaxa) in net, yet $(join(size(dhat),"×")) distances")
    nstates = mapreduce(x -> maximum(skipmissing(x)), max, trait)
    distancecorrection_JC!(dhat, nstates) # default scalar=1.01
    taxonnames = [net.leaf[i].name for i in sortperm([n.number for n in net.leaf])]
    # taxon names: to tell the calibration that row i of dhat is for taxonnames[i]
    # trait[i][j] = trait j for taxon at node number i: 'node.number' = i
    calibratefrompairwisedistances!(net, dhat, taxonnames,
        forceMinorLength0=true, ultrametric=false)
        #= force minor length to 0 to avoid non-identifiability at zippers.
        works well if the true (or "the" best-fit) minor parent edge is shorter
        than that of the major parent edge: to zip-up =#
    for e in net.edge # avoid starting at the boundary
        if e.length < 1.0e-10
            e.length = 0.0001
        end
    end
    return net
end
