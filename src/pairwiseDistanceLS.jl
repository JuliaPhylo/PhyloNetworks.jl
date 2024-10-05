"""
    getNodeAges(net)

vector of node ages in pre-order, as in `nodes_changed`.

*Warnings*: `net` is assumed to
- have been preordered before (to calculate `nodes_changed`)
- be time-consistent (all paths to the root to a given hybrid have the same length)
- be ultrametric (all leaves have the same age: 0)
"""
function getNodeAges(net::HybridNetwork)
    x = Vector{Float64}(undef, length(net.nodes_changed))
    for i in reverse(1:length(net.nodes_changed)) # post-order
        n = net.nodes_changed[i]
        if n.leaf
            x[i] = 0.0
            continue
        end
        for e in n.edge
            if getparent(e) == n # n parent of e
                childnode = getchild(e)
                childIndex = getIndex(childnode, net.nodes_changed)
                x[i] = x[childIndex] + e.length
                break # use 1st child, ignores all others
            end
        end
    end
    return x
end

"""
    pairwiseTaxonDistanceMatrix(net; keepInternal=false,
                                checkPreorder=true, nodeAges=[])
    pairwiseTaxonDistanceMatrix!(M, net, nodeAges)

Return the matrix `M` of pairwise distances between nodes in the network:
- between all nodes (internal and leaves) if `keepInternal=true`,
  in which case the nodes are listed in `M` in the
  order in which they appear in `net.nodes_changed`
- between taxa only otherwise, in which case the nodes are listed
  in `M` in the order in which they appear in `tipLabels(net)`
  (i.e. same order as in `net.leaf`)

The second form modifies `M` in place, assuming all nodes.

The distance between the root and a given hybrid node (to take an example)
is the weighted average of path lengths from the root to that node,
where each path is weighted by the product of γs of all edges on that path.
This distance measures the average genetic distance across the genome,
if branch lengths are in substitutions/site.

optional arguments:

- `checkPreorder`: if true, `net.nodes_changed` is updated to get a
  topological ordering of nodes.
- `nodeAges`: if not provided, i.e. empty vector, the network is *not* modified.  
  If provided and non-empty, `nodeAges` should list node ages in the
  pre-order in which nodes are listed in `nodes_changed` (including leaves),
  and **edge lengths** in `net` **are modified** accordingly.

Providing node ages hence makes the network time consistent: such that
all paths from the root to a given hybrid node have the same length.
If node ages are not provided, the network need not be time consistent.
"""
function pairwiseTaxonDistanceMatrix(net::HybridNetwork;
            keepInternal=false::Bool, checkPreorder=true::Bool,
            nodeAges=Float64[]::Vector{Float64})
    net.isRooted || error("net needs to be rooted for preorder recursion")
    if(checkPreorder)
        preorder!(net)
    end
    nnodes = net.numNodes
    M = zeros(Float64,nnodes,nnodes)
    isempty(nodeAges) || length(nodeAges) == nnodes ||
        error("there should be $nnodes node ages")
    pairwiseTaxonDistanceMatrix!(M,net,nodeAges)
    if !keepInternal
        M = getTipSubmatrix(M, net)
    end
    return M
end
"""
    getTipSubmatrix(M, net; indexation=:both)

Extract submatrix of M, with rows and/or columns corresponding to
tips in the network, ordered like in `net.leaf`.
In M, rows and/or columns are assumed ordered as in `net.nodes_changed`.

indexation: one of `:rows`, `:cols` or `:both`: are nodes numbers indexed
in the matrix by rows, by columns, or both? Subsetting is taken accordingly.

"""
function getTipSubmatrix(M::Matrix, net::HybridNetwork; indexation=:both)
    nodenames = [n.name for n in net.nodes_changed]
    tipind = Int[]
    for l in tipLabels(net)
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

function pairwiseTaxonDistanceMatrix!(M::Matrix{Float64},net::HybridNetwork,nodeAges)
    traversal_preorder!(net.nodes_changed, M, # updates M in place
            traversalupdate_default!, # does nothing
            updateTreePairwiseTaxonDistanceMatrix!,
            updateHybridPairwiseTaxonDistanceMatrix!,
            nodeAges)
end

function updateTreePairwiseTaxonDistanceMatrix!(
    V::Matrix,
    i::Int,
    parindx::Int, # index of parent node
    edge::Edge,
    nodeages, # assumed pre-ordered, as in nodes_changed
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

function updateHybridPairwiseTaxonDistanceMatrix!(
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
        for (pi,pe) in zip(parindx, paredge)
            V[i,j] += pe.gamma * (V[pi,j] + pe.length)
        end
        V[j,i] = V[i,j]
    end
    return true
end

"""
    pairwiseTaxonDistanceGrad(net; checkEdgeNumber=true, nodeAges=[])

3-dim array: gradient of pairwise distances between all nodes.
(internal and leaves); gradient with respect to edge lengths
if `nodeAges` is empty; with respect to node ages otherwise.
Assume correct `net.nodes_changed` (preorder).  
This gradient depends on the network's topology and γ's only,
not on branch lengths or node ages (distances are linear in either).

WARNING: edge numbers need to range between 1 and #edges.
"""
function pairwiseTaxonDistanceGrad(net::HybridNetwork;
        checkEdgeNumber=true::Bool, nodeAges=Float64[]::Vector{Float64})
    if checkEdgeNumber
      sort([e.number for e in net.edge]) == collect(1:net.numEdges) ||
        error("edge numbers must range between 1 and #edges")
    end
    n = (length(nodeAges)==0 ? net.numEdges : net.numNodes)
    M = zeros(Float64, net.numNodes, net.numNodes, n)
    traversal_preorder!(net.nodes_changed, M,
            traversalupdate_default!,  # does nothing
            updateTreePairwiseTaxonDistanceGrad!,
            updateHybridPairwiseTaxonDistanceGrad!,
            nodeAges) # nodeAges assumed pre-ordered, like nodes_changed
    return M
end
function updateTreePairwiseTaxonDistanceGrad!(
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
function updateHybridPairwiseTaxonDistanceGrad!(
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
    calibrateFromPairwiseDistances!(net, distances::Matrix{Float64},
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
- checkPreorder (true)
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
function calibrateFromPairwiseDistances!(net::HybridNetwork,
      D::Array{Float64,2}, taxNames::Vector{<:AbstractString};
      checkPreorder=true::Bool, forceMinorLength0=false::Bool, verbose=false::Bool,
      ultrametric=true::Bool, NLoptMethod=:LD_MMA::Symbol,
      ftolRel=fRelBL::Float64, ftolAbs=fAbsBL::Float64,
      xtolRel=xRelBL::Float64, xtolAbs=xAbsBL::Float64)

    checkPreorder && preorder!(net)
    # fixit: remove root node if of degree 2, and if ultrametric=false
    defaultedgelength = median(D)/(length(net.edge)/2)
    for e in net.edge
        if e.length == -1.0 e.length=defaultedgelength; end
        # get smarter starting values: NJ? fast dating?
    end
    if ultrametric # get all node ages in pre-order
        na = getNodeAges(net)
    else na = Float64[]; end
    # get number and indices of edge/nodes to be optimized
    if forceMinorLength0 && !ultrametric
        parind = filter(i -> net.edge[i].isMajor, 1:net.numEdges)
        nparams = length(parind) # edges to be optimized: indices
        par = [e.length for e in net.edge][parind] # and lengths
        for i in 1:net.numEdges
            if !net.edge[i].isMajor net.edge[i].length=0.0; end
        end
    elseif !forceMinorLength0 && !ultrametric
        nparams = length(net.edge) # number of parameters to optimize
        par = [e.length for e in net.edge]
        parind = 1:nparams
    elseif !forceMinorLength0 && ultrametric
        nparams = net.numNodes - net.numTaxa # internal nodes to be optimized
        parind = filter(i -> !net.nodes_changed[i].leaf, 1:net.numNodes)
        par = na[parind]
    else # forceMinorLength0 && ultrametric
        nparams = net.numNodes - net.numTaxa - net.numHybrids
        parind = filter(i -> !(net.nodes_changed[i].leaf || net.nodes_changed[i].hybrid),
                        1:net.numNodes)
        par = na[parind]
        hybInd = filter(i -> net.nodes_changed[i].hybrid, 1:net.numNodes)
        hybParentInd = Int[] # index in nodes_changed of minor parent
        hybGParentI = Int[] # index in 1:nparams of minor (grand-)parent in param list
        for i in hybInd
            n = net.nodes_changed[i]
            p = getparentminor(n)
            pi = findfirst(n -> n===p, net.nodes_changed)
            push!(hybParentInd, pi)
            pii = findfirst(isequal(pi), parind)
            while pii===nothing # in case minor parent of n is also hybrid node
                p = getparentminor(p)
                pi = findfirst(n -> n===p, net.nodes_changed)
                pii = findfirst(isequal(pi), parind)
            end
            push!(hybGParentI, pii)
        end
    end
    # initialize M=dist b/w all nodes, G=gradient (constant)
    M = pairwiseTaxonDistanceMatrix(net, keepInternal=true,
            checkPreorder=false, nodeAges=na)
    if !ultrametric && sort([e.number for e in net.edge]) != collect(1:net.numEdges)
        for i in 1:net.numEdges # renumber edges, needed for G
            net.edge[i].number = i
        end
    end # G assumes edges numbered 1:#edges, if optim edge lengths
    G = pairwiseTaxonDistanceGrad(net, checkEdgeNumber=false, nodeAges=na) .* 2
    # match order of leaves in input matrix, versus pre-order
    nodenames = [n.name for n in net.nodes_changed] # pre-ordered
    ntax = length(taxNames)
    tipind = Int[] # pre-order index for leaf #i in dna distances
    for l in taxNames
        i = findfirst(isequal(l), nodenames)
        i !== nothing || error("taxon $l not found in network")
        push!(tipind, i)
    end
    # contraints: to force a parent to be older than its child
    if ultrametric
      numConstraints = length(parind) -1 + net.numHybrids
      # calculate indices in param list of child & (grand-)parent once
      chii = Int[] # non-root internal node, can be repeated: once per constraint
      anii = Int[] # closest ancestor in param list
      ci = 1 # index of constraint
      for i in 2:length(net.nodes_changed) # 1=root, can skip
        n = net.nodes_changed[i]
        if n.leaf continue; end # node ages already bounded by 0
        if n.hybrid && forceMinorLength0          # get index in param list of
          nii = hybGParentI[findfirst(isequal(i), hybInd)] # minor grand-parent (same age)
        else
          nii = findfirst(isequal(i), parind)
        end
        for e in n.edge
          if getchild(e) == n # n child of e
            p = getparent(e)  # parent of n
            if forceMinorLength0 && n.hybrid && !e.isMajor
                continue; end # p and n at same age already
            pi = findfirst(isequal(p.number), [no.number for no in net.nodes_changed])
            if forceMinorLength0 && p.hybrid
              pii = hybGParentI[findfirst(isequal(pi), hybInd)]
            else
              pii = findfirst(isequal(pi), parind)
            end
            push!(chii, nii)
            push!(anii, pii)
            verbose && println("node $(net.nodes_changed[parind[nii]].number) constrained by age of parent $(net.nodes_changed[parind[pii]].number)")
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
                for i in 1:net.numHybrids
                    na[hybInd[i]] = na[hybParentInd[i]] # pre-order important
                end
            end
        else # not ultrametric: optimize branch lengths
            for i in 1:nparams # update network
                net.edge[parind[i]].length = x[i]
            end
        end
        # update distances in M, in place
        pairwiseTaxonDistanceMatrix!(M,net,na)
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
              for i in 1:net.numHybrids # na[hybrid] set to na[minor parent]
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
    startingBL!(net::HybridNetwork,
                trait::AbstractVector{Vector{Union{Missings.Missing,Int}}},
                siteweight=ones(length(trait[1]))::AbstractVector{Float64})

Calibrate branch lengths in `net` by minimizing the mean squared error
between the JC-adjusted pairwise distance between taxa, and network-predicted
pairwise distances, using [`calibrateFromPairwiseDistances!`](@ref).
The network is *not* forced to be time-consistent nor ultrametric.
To avoid one source of non-identifiability, the network is "zipped" by
forcing minor hybrid edges to have length 0.

`siteweight[k]` gives the weight of site (or site pattern) `k` (default: all 1s).

Assumptions:

- all species have the same number of traits (sites): `length(trait[i])` constant
- `trait[i]` is for leaf with `node.number = i` in `net`, and
  `trait[i][j] = k` means that leaf number `i` has state index `k` for trait `j`.
  These indices are those used in a substitution model (see `PhyloTraits.jl`):
  kth value of `getlabels(model)`.

Other:
- Hamming distances are calculated for each pair, ignoring any site in which
  one of 2 values is missing. The total number of differences is then divided
  by the total weight of all sites, ignoring that some of them may have been
  missing for the pair. This is an inexact rescaling of the hamming distance,
  assuming a small proportion of missing values for each pair.
- Theoretically, Hamming distances are < 0.75 with four states,
  or < (n-1)/n for n states.
  If not, all pairwise hamming distances are scaled by `.75/(m*1.01)` where `m`
  is the maximum observed hamming distance, to make them all < 0.75.
  The JC correction calculates: d_JC = - 0.75 log(1-d_hamming/0.75) for n=4 states
- At the end, any edge length smaller than 1.0e-10 is reset to 0.0001,
  to avoid the 0 boundary.
"""
function startingBL!(
    net::HybridNetwork,
    trait::AbstractVector{Vector{Union{Missings.Missing,Int}}},
    siteweight=ones(length(trait[1]))::AbstractVector{Float64}
)
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
    # taxon names: to tell the calibration that row i of dhat is for taxonnames[i]
    # trait[i][j] = trait j for taxon at node number i: 'node.number' = i
    calibrateFromPairwiseDistances!(net, dhat, taxonnames,
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
