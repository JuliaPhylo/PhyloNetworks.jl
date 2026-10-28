#####################################################
## Tree NJ functions
#####################################################

function check_distance_matrix(D::Matrix{<:Real})
    size(D, 1) > 1 || throw(DomainError(D, "Too few nodes"))
    if any(D .< 0)
        throw(DomainError(D, "Negative entries in distance matrix"))
    end
    if any(diag(D) .!= 0.0)
        throw(DomainError(D, "Diagonal entry not 0"))
    end
    return nothing
end

# Find the pair of indices that minimize the Q criterion:
# Q_ij = (n-2) * D_ij - sum_{a} (D_ia + D_ja)
# and compute new edge lengths
function Q_findpair(D::Matrix{<:Real}, active_size::Int)
    # compute Q matrix and find pair with minimum Q
    sums = sum(D[1:active_size, 1:active_size], dims = 1)
    cur = Inf
    min_index = (0, 0)
    for i = 1:active_size
        for j = i:active_size
            if j != i
                qij = (active_size-2) * D[i,j] - sums[i] - sums[j]
                if qij < cur
                    cur = qij
                    min_index = (i, j)
                end
            end
        end
    end

    (i, j) = min_index
    return i, j, sums
end

function nj_findbranchlengths(D::Matrix{<:Real}, active_size::Int, sums::Matrix{<:Real}, i, j, force_nonnegative_edges::Bool, neglenp::Int)

    dik = D[i,j] / 2 + (sums[i] - sums[j]) / (2 * (active_size - 2))
    djk = D[i,j] - dik

    if dik < 0.0
        neglenp += 1
        if force_nonnegative_edges dik = 0.0; end
    end
    if djk < 0.0
        neglenp += 1
        if force_nonnegative_edges djk = 0.0; end
    end

    return dik, djk, neglenp
end

# Modify in place D and net to join tips i and j
function join_ij!(D::Matrix{<:Real}, net::HybridNetwork, active_nodes::Array{Node,1}, active_size::Int, i::Int, j::Int, dik::Real, djk::Real)
    # create new edges and node, update tree
    edgenum = net.numedges
    eik = Edge(edgenum + 1, dik)
    ejk = Edge(edgenum + 2, djk)
    node_k = Node(net.numnodes+1, false, false, [eik, ejk])
    node_i = active_nodes[i]
    node_j = active_nodes[j]
    setNode!(eik, Node[node_i, node_k]) # correct direction, bc ischild1 true by default
    setNode!(ejk, Node[node_j, node_k])
    setEdge!(node_i, eik)
    setEdge!(node_j, ejk)
    pushEdge!(net, eik)
    pushEdge!(net, ejk)
    pushNode!(net, node_k)

    # Update distances to new node k, placed in position i
    if i == active_size # swap with j instead, if i is last active position
        (i,j) = (j,i)
    end
    for l in 1:active_size
        if !(l in [i, j])
            D[l,i] = (D[l,i] + D[l,j] - D[i,j]) / 2
            D[i,l] = D[l,i]
        end
    end
    # swap j with last active position using whole column/row operations
    if j != active_size
        D[:, j], D[:, active_size] = D[:, active_size], D[:, j]
        D[j, :], D[active_size, :] = D[active_size, :], D[j, :]
        active_nodes[j] = active_nodes[active_size]
    end
    # inactivate last active position: okay bc k replaces i != active_size
    active_size -= 1
    active_nodes[i] = node_k

    return active_size
end

function join_12!(net::HybridNetwork, active_nodes::Array{Node,1}, d12::Real)
    node1 = active_nodes[1]
    node2 = active_nodes[2]
    newedge = Edge(net.numedges+1, d12)
    setNode!(newedge, [node1, node2]) # ischild1 = true by default: nodes are [child, parent]
    if node1.number == net.numnodes
        newedge.ischild1 = false # direct the edge: last created node -> other node
    end
    setEdge!(node1, newedge)
    setEdge!(node2, newedge)
    pushEdge!(net, newedge)
end

"""
    nj!(D::Matrix{Float64}, names::AbstractVector{<:AbstractString}=String[];
        force_nonnegative_edges::Bool=false)

Construct a phylogenetic tree from the input distance matrix and
vector of names (as strings), using the Neighbour-Joinging algorithm
([Satou & Nei 1987](https://doi.org/10.1093/oxfordjournals.molbev.a040454)).
The order of the `names` argument should match that of the row (and
column) of the matrix `D`.
With `force_nonnegative_edges` being `true`, any negative edge length
is changed to 0.0 (with a message).

Warning: `D` is modified.
"""
function nj!(
    D::Matrix{Float64},
    names::AbstractVector{<:AbstractString}=String[];
    force_nonnegative_edges::Bool=false
)
    check_distance_matrix(D)
    n = size(D, 1)    # number of species
    active_size = n   # tracks the current active size of matrix
    if isempty(names)
        names = string.(1:n)
    end
    # create empty network with n unconnected leaf nodes
    nodes = map(function(i)
                node = Node(i, true)
                node.name = names[i]
                return node
                end,
                1:n)
    net = HybridNetwork(nodes, Edge[])
    active_nodes = copy(nodes)
    neglenp = 0  # number of negative edge lengths

    while active_size > 2
        # compute Q matrix and find pair with minimum Q
        i, j, sums = Q_findpair(D, active_size)
        # compute the length of the edges
        dik, djk, neglenp = nj_findbranchlengths(D, active_size, sums, i, j, force_nonnegative_edges, neglenp)
        # update tree
        active_size = join_ij!(D, net, active_nodes, active_size, i, j, dik, djk)
    end

    # base case: 2 leaves
    if D[1,2] < 0.0
        neglenp += 1
        if force_nonnegative_edges D[1,2] = 0.0; end
    end
    join_12!(net, active_nodes, D[1,2])
    net.rooti = net.numnodes # root = last created node, which is internal
    net.isrooted = true # because rooti is valid, and edges are directed away from it
    if neglenp > 0
        infostr = (force_nonnegative_edges ?
                   "$neglenp branch(es) of negative length, reset to 0" :
                   "$neglenp branch(es) of negative length" )
        @info infostr
    end
    return net
end

"""
    nj(D::DataFrame; force_nonnegative_edges::Bool=false)

Construct a tree from a distance matrix by neighbor joining, where
`D` is a `DataFrame` of the distance matrix, with taxon names taken
from the header of the data frame.
The rows are assumed to correspond to tips in the tree in the same order
as they do in columns.
With `force_nonnegative_edges` being `true`, any negative edge length
is changed to 0.0 (with a message).

For the algorithm, see
[Satou & Nei 1987](https://doi.org/10.1093/oxfordjournals.molbev.a040454).

See [`nj!`](@ref) for using a matrix as input.
"""
function nj(D::DataFrame; force_nonnegative_edges::Bool=false)
    nj!(Matrix{Float64}(D), names(D); force_nonnegative_edges=force_nonnegative_edges)
end

#####################################################
## F2 and f4 fonctions - diagonal case
#####################################################

# Un-normalized f4 statistic for tips (a,b,c,d) in D
function getf4(D::Matrix{<:Real}, a::Int, b::Int, c::Int, d::Int)
    return D[a,d] + D[b,c] - D[a,c] - D[b,d]
end

function getf4(D::Matrix{<:Real}, t::Vector{Int})
    return D[t[1],t[4]] + D[t[2],t[3]] - D[t[1],t[3]] - D[t[2],t[4]]
end

# pair of f4
# If (ab|cd) verify the 4 point condition, we should have:
# f4(D,a,b,c,d) = 0, f4(D,a,c,b,d) > 0
function getpairf4(D::Matrix{<:Real}, a::Int, b::Int, c::Int, d::Int)
    return getf4(D,a,b,c,d), getf4(D,a,c,b,d)
end

function getpairf4(D::Matrix{<:Real}, t::Vector{Int})
    return getf4(D,t), getf4(D,t[[1,3,2,4]])
end

# lrt statistic assuming V = 1/nloci * I
# If (ab|cd) verify the 4 point condition, we should have:
# f4(D,a,b,c,d) = 0
# f4(D,a,c,b,d) > 0
function getlrtidentity(D::Matrix{<:Real}, ind_A::Vector{Int}, ind_B::Vector{Int}, nloci::Int)
    n_A = size(ind_A, 1)
    n_B = size(ind_B, 1)
    lrt = 0.0
    i = 0
    for a in 1:n_A
        for b in a:n_A
            if (b != a)
                for c in 1:n_B
                    for d in c:n_B
                        if (d != c)
                            f1, f2 = getpairf4(D, ind_A[a], ind_A[b], ind_B[c], ind_B[d])
                            lrt += (f1 - 0.0)^2 + (f2 - max(0, f2))^2
                            i += 1
                        end
                    end
                end
            end
        end
    end
    return lrt * nloci, i
end

# Mixture of chi2 distribution
# Assuming a diagonal V matrix
# weight[i] = Prob[Z has exactly i positive components] with Z ~ N(0, V)
function mixtureChi2Diagonal(r)
    weights = pdf(Binomial(r, 0.5), r:-1:0)
    return MixtureModel(Chisq, r:2*r, weights)
end

function cutblobpvalue(D::Matrix{<:Real}, clus_ind::BitVector, nloci::Int)
    cutblobpvalue(D, findall(clus_ind), findall(.!clus_ind), nloci::Int)
end

function cutblobpvalue(D::Matrix{<:Real}, ind_A::Vector{Int}, ind_B::Vector{Int}, nloci::Int)
    lrt, r = getlrtidentity(D, ind_A, ind_B, nloci)
    return ccdf(mixtureChi2Diagonal(r), lrt)
end

#####################################################
## F2 and f4 fonctions - non-diagonal case
#####################################################

function check_variance_matrix(Qf2::Matrix{<:Real}, n::Int)
    if size(Qf2, 1) > 0
        expdim = n * (n-1) ÷ 2
        if (size(Qf2, 1) != expdim) || (size(Qf2, 2) != expdim)
            throw(DomainError(Qf2, "Qf2 should be a square matrix of dimension n(n+1)/2 with n the number of nodes"))
        end
    end
    return nothing
end

# get the index of the f2 value in the Q matrix
# Q matrix is sorted in lexicographic order
function get_index_f2(i::Int, j::Int, n::Int)
    if (i < j)
        return (i-1) * n - (i-1)*i ÷ 2 + j - i
    else
        return get_index_f2(j,i,n)
    end
end

# get all the indexes in Q
# indf4 is nf4 times 4 the indices of the f4 values
function get_f4_ind(indf4::Matrix{Int}, n::Int)
    f(x::Vector) = [get_index_f2(x[1],x[4], n), get_index_f2(x[2],x[3], n), get_index_f2(x[1],x[3], n), get_index_f2(x[2],x[4], n)]
    return mapslices(f, indf4, dims = 2)
end
function get_covar_f4(Qmatrix::Matrix{Float64}, indf4::Matrix{Int}, n::Int)
    dimv = size(indf4, 1)
    v = zeros(dimv, dimv)
    # indexes
    indf2 = get_f4_ind(indf4, n)
    indsign = [1, 1, -1, -1]
    # var
    for iv in 1:dimv
        for jv in iv:dimv
            for (i,ia) in enumerate(indf2[iv,:])
                for (j,jb) in enumerate(indf2[jv,:])
                    v[iv,jv] += Qmatrix[ia,jb] * indsign[i]*indsign[j]
                end
            end
            v[jv,iv] = v[iv,jv]
        end
    end
    return v
end

# lrt statistic assuming V = Q
# If (ab|cd) verify the 4 point condition, we should have:
# f4(D,a,b,c,d) = 0
# f4(D,a,c,b,d) > 0
function getlrtQ(D::Matrix{<:Real}, Q::Matrix{<:Real}, ind_A::Vector{Int}, ind_B::Vector{Int}, active_size::Int, lambda::Real=1e-10)
    n_A = size(ind_A, 1)
    n_B = size(ind_B, 1)
    dimf4 = (n_A * (n_A - 1) ÷ 2) * (n_B * (n_B - 1) ÷ 2)
    f4vec = zeros(2 * dimf4)
    indf4 = zeros(Int, 2 * dimf4, 4)
    i = 0
    for a in 1:n_A
        for b in a:n_A
            if (b != a)
                for c in 1:n_B
                    for d in c:n_B
                        if (d != c)
                            i += 1
                            indf4[i,:] = [ind_A[a], ind_A[b], ind_B[c], ind_B[d]]
                            indf4[dimf4 + i,:] = indf4[i,[1,3,2,4]]
                            f1, f2 = getpairf4(D, indf4[i,:])
                            f4vec[i] = f1
                            f4vec[dimf4 + i] = f2
                        end
                    end
                end
            end
        end
    end
    V = get_covar_f4(Q, indf4, active_size)

    # TODO: find a better formula for regularization ?

    # Simple regularization
    Vprime = V + lambda .* tr(V) .* I(size(V,1)) 
    # in admixtools, they add diag_f3 * trace(Q) * Identity matrix, with diag_f3 = 1e-5 by default
    # See https://github.com/uqrmaie1/admixtools/blob/472d04e98ac8aa70e5c9cbe65722c7e5dfbe8200/R/qpgraph.R#L524
    # and https://github.com/uqrmaie1/admixtools/blob/472d04e98ac8aa70e5c9cbe65722c7e5dfbe8200/R/qpgraph.R#L204
    # in poolfstat, they just add Q.lambda, default to zero
    # see https://github.com/cran/poolfstat/blob/47d0fe5e1bc867d27362a5bbd6edddd632d3ea7d/R/fit.graph.R#L19
    Vinv = inv(Symmetric(Vprime))

    # # EB simple estimator from Montoya 2026
    # Vprime = (nloci * V + I(size(V,1))) / (nloci - 2 + size(V,1)+1)
    # Vinv = inv(Symmetric(Vprime))

    # # EB from Colucci 2015
    # p = size(V,1)
    # nu = p + 1
    # Vprime = (nloci * V + I(p)) / (nloci + nu - p)
    # Vinv = inv(Symmetric(Vprime))

    # # Haff 1980 EB estimator
    # # S ~ W_p(Sigma, k), k = nloci - 1, p = dim(V)
    # k = nloci - 1
    # p = size(V, 1)
    # S = k * V 
    # # Sigma^-1 ~ W(1/gamma I, k'), k' = p+1
    # kprime = p + 1
    # # 3.4:  a = 1/(k+k'-p-1)
    # a = 1 / (k + kprime - p - 1)
    # # 3.7:  pk / (k+k') * 1/tr(S^-1 C), C = I
    # gam = (p * k) / (k + kprime) * 1 / sum(1/eigen(S).values)
    # Vprime = a * (S + gam * I(p))
    # Vinv = inv(Symmetric(Vprime))

    # # Efron Moris 1976
    # k = nloci - 1
    # p = size(V,1)
    # S = k * V 
    # Sinv = inv(Symmetric(S))
    # Vinv = (k - p - 1) * Sinv + (p^2 - p - 2)/tr(S) * I(p)



    # d = DAQP.Model()
    # f = - Vinv * f4vec
    # A = Matrix(1.0I, 2*dimf4, 2*dimf4)
    # # blower = zeros(Float64, 2 * dimf4)
    # # bupper = [zeros(dimf4); 1e16 .* ones(Float64, dimf4)]
    # blower = -1e16 .* ones(Float64, 2 * dimf4)
    # bupper = 1e16 .* ones(Float64, 2 * dimf4)
    # sense = zeros(Cint, 2 * dimf4)
    # # sense[1:dimf4] .= 5 # equality constraint
    # DAQP.setup(d,Vinv,f,A,bupper,blower,sense)
    # DAQP.settings(d, Dict(:iter_limit =>2000))
    # x,fval,exitflag,info = DAQP.solve(d)
    # 2 * fval + transpose(f4vec) * Vinv * f4vec
    # ff = copy(f4vec)
    # ff[1:dimf4] .=0
    # transpose(f4vec - ff) * Vinv * (f4vec - ff)

    d = COSMO.Model()
    Imat = Matrix(1.0I, dimf4, dimf4)
    # first dimf4 entries equal to zero
    constraint1 = COSMO.Constraint(hcat(Imat, 0.0*Imat), zeros(dimf4), COSMO.ZeroSet)
    # last dimf4 entries greater than zero
    constraint2 = COSMO.Constraint(hcat(0.0Imat, Imat), zeros(dimf4), COSMO.Nonnegatives)
    # Model
    P = Vinv
    q = - Vinv * f4vec
    COSMO.assemble!(d, P, q, [constraint1; constraint2])
    res = COSMO.optimize!(d)
    return 2 * res.obj_val - transpose(f4vec) * q, Vprime

    # transpose(res.x) * Vinv * res.x - 2 * transpose(f4vec) * Vinv * res.x
    # transpose(res.x) * P * res.x + 2 * transpose(q)  * res.x
    # transpose(f4vec - res.x) * Vinv * (f4vec - res.x)
end

function cutblobpvalue(D::Matrix{<:Real}, Q::Matrix{<:Real}, clus_ind::BitVector, active_size::Int, lambda::Real)
    cutblobpvalue(D, Q, findall(clus_ind), findall(.!clus_ind), active_size, lambda)
end

function cutblobpvalue(D::Matrix{<:Real}, Q::Matrix{<:Real}, ind_A::Vector{Int}, ind_B::Vector{Int}, active_size::Int, lambda::Real)
    lrt, V = getlrtQ(D, Q, ind_A, ind_B, active_size, lambda)
    # return ccdf(mixtureF(r, nloci), lrt)
    return ccdf(mixtureChi2(V), lrt)
end

# Mixture of f distribution
# n is the number of blocks in jacknife estimation
function mixtureF(r::Int, n::Int)
    if r > 1
        return MixtureModel([(r - 1) / (n - r + 1) * FDist(r - 1, n - r + 1), r / (n - r) * FDist(r, n - r)], [0.5, 0.5])
    else # Dirac
        return MixtureModel([Normal(0.0, 0.0), r / (n - r) * FDist(r, n - r)], [0.5, 0.5])
    end
end

# Mixture of chi2 distribution
# Assuming a diagonal V matrix
# weight[i] = Prob[Z has exactly i positive components] with Z ~ N(0, V)
function mixtureChi2(V::Matrix{<:Real}, nmc::Int=10000*size(V,1))
    r = size(V, 1) ÷ 2
    ## Compute A = V22 - V12' * inv(V11) * V12
    ## Naive approach
    # V11 = view(V, 1:r, 1:r)
    # V22 = view(V, (r+1):(2*r), (r+1):(2*r))
    # V12 = view(V, 1:r, (r+1):(2*r))
    # A = Symmetric(V22 - transpose(V12) * (V11 \ V12))
    ## Using block Cholesky, cf https://math.stackexchange.com/questions/1384769/proving-that-l-22l-22t-s-is-the-schur-complement-of-a-cholesky-factorizat
    L_V = cholesky(V).L
    L_A = view(L_V, (r+1):(2*r), (r+1):(2*r))
    ## Other approach : cholesky of smaller matrices
    # L_1 = cholesky(view(V, 1:r, 1:r)).L
    # L_12t = L_1 \ view(V, 1:r, (r+1):(2*r))
    # A = view(V, (r+1):(2*r), (r+1):(2*r)) - L_12t' * L_12t
    # L_A = cholesky(A).L

    ## Draw Z and count positives
    # Z = L_A * rand(MvNormal(zeros(r), I), nmc)
    Z = L_A * reshape(rand(Normal(0,1), nmc * r), r, nmc)
    npos = mapslices(x -> sum(x .> 0.0), Z, dims = 1)
    weights = zeros(Int, r + 1)
    for val in npos
        weights[val + 1] += 1
    end
    weights = weights / nmc
    return MixtureModel(Chisq, r:2*r, reverse(weights))
end

#####################################################
## Tree of blob NJ functions
#####################################################
# Find distance between clusters of tips:
# average distances over all the tips in the clusters.
# Note: diagonal elements are left to zero.
function dist_clusters!(Dclust::Matrix{<:Real}, D::Matrix{<:Real}, active_size::Int, clusters::Vector{Int})
    uniclust = unique(clusters[1:active_size])
    nclust = size(uniclust, 1)
    for (ic1,c1) in enumerate(uniclust)
        in_c1 = findall(clusters[1:active_size] .== c1)
        for ic2 in ic1:nclust
            if ic2 != ic1
                c2 = uniclust[ic2]
                in_c2 = findall(clusters[1:active_size] .== c2)
                Dclust[ic1,ic2] = sum(D[in_c1,in_c2]) / size(in_c1, 1) / size(in_c2, 1)
                Dclust[ic2,ic1] = Dclust[ic1,ic2]
            end
        end
    end
end

# Merge two clusters in one
# New cluster inherits from the smallest cluster number
# clusters is a vector of n entries, with entry i giving the cluster number of node i
function clusterize!(clusters::Vector{Int}, c1::Int, c2::Int)
    newclust = min(c1, c2)
    clusters[clusters .== c1] .= newclust
    clusters[clusters .== c2] .= newclust
end

# Join all the tips of a blob into a single tip
# Tips of the blob form a star tree from the new tip
# Branch lengths are arbitrarily set to 1.0
# Distance of the new node to all other nodes is the average distance to the nodes of the blob
function join_blob!(D::Matrix{<:Real}, net::HybridNetwork, active_nodes::Array{Node,1}, active_size::Int, clusters::Vector{Int}, clus_ind::BitVector)
    # create new edges and node, update tree
    # create a star tree from new node node_k
    clus_ind = findall(clus_ind)
    n_in_clus = size(clus_ind, 1)
    edgenum = net.numedges
    newEdges = Vector{Edge}(undef, n_in_clus)
    for (i,a) in enumerate(clus_ind)
        eak = Edge(edgenum + i, 1.0)
        setEdge!(active_nodes[a], eak)
        newEdges[i] = eak
    end
    node_k = Node(net.numnodes+1, false, false, newEdges)
    for (i,a) in enumerate(clus_ind)
        setNode!(newEdges[i], Node[active_nodes[a], node_k]) # correct direction, bc ischild1 true by default
        pushEdge!(net, newEdges[i])
    end
    pushNode!(net, node_k)

    # Update distances to new node k, placed in position i
    i = clus_ind[1]
    for l in 1:active_size
        if !(l in clus_ind)
            D[l,i] = sum(D[l,clus_ind]) / n_in_clus
            D[i,l] = D[l,i]
        end
    end
    # remove all nodes of the blob from the scope
    for j in reverse(clus_ind[2:n_in_clus])
        # swap j with last active position using whole column/row operations
        if j != active_size
            D[:, j], D[:, active_size] = D[:, active_size], D[:, j]
            D[j, :], D[active_size, :] = D[active_size, :], D[j, :]
            active_nodes[j] = active_nodes[active_size]
            clusters[j] = clusters[active_size]
        end
        active_size -= 1
    end
    # inactivate last active position: okay bc k replaces i != active_size
    active_nodes[i] = node_k

    return active_size
end

function update_covar_f2(Qmatrix::Matrix{Float64}, clus_ind::Vector{Int}, i::Int, l::Int, active_size::Int)
    n_in_clus = size(clus_ind, 1)
    # index to update
    il = get_index_f2(i, l, active_size)
    # update
    for a in 1:active_size
        for b in a:active_size
            if (b != a)
                ab = get_index_f2(a, b, active_size)
                Qmatrix[il, ab] = 0.0
                for k in clus_ind
                    lk = get_index_f2(l, k, active_size)
                    Qmatrix[il, ab] += Qmatrix[lk, ab] ./ n_in_clus
                end
                Qmatrix[ab, il] = Qmatrix[il, ab]  
            end
        end
    end
end

function reorder_covar_f2(Qmatrix::Matrix{Float64}, j::Int, active_size::Int)
    for a in 1:active_size
        if (a != j && a != active_size)
            aj = get_index_f2(a, j, active_size)
            aac = get_index_f2(a, active_size, active_size)
            Qmatrix[:, aj], Qmatrix[:, aac] = Qmatrix[:, aac], Qmatrix[:, aj]
            Qmatrix[aj, :], Qmatrix[aac, :] = Qmatrix[aac, :], Qmatrix[aj, :]  
        end
    end
end

function join_blob!(D::Matrix{<:Real}, Qf2::Matrix{<:Real}, net::HybridNetwork, active_nodes::Array{Node,1}, active_size::Int, clusters::Vector{Int}, clus_ind::BitVector)
    # create new edges and node, update tree
    # create a star tree from new node node_k
    clus_ind = findall(clus_ind)
    n_in_clus = size(clus_ind, 1)
    edgenum = net.numedges
    newEdges = Vector{Edge}(undef, n_in_clus)
    for (i,a) in enumerate(clus_ind)
        eak = Edge(edgenum + i, 1.0)
        setEdge!(active_nodes[a], eak)
        newEdges[i] = eak
    end
    node_k = Node(net.numnodes+1, false, false, newEdges)
    for (i,a) in enumerate(clus_ind)
        setNode!(newEdges[i], Node[active_nodes[a], node_k]) # correct direction, bc ischild1 true by default
        pushEdge!(net, newEdges[i])
    end
    pushNode!(net, node_k)

    # Update distances to new node k, placed in position i
    i = clus_ind[1]
    for l in 1:active_size
        if !(l in clus_ind)
            D[l,i] = sum(D[l,clus_ind]) / n_in_clus
            D[i,l] = D[l,i]
            update_covar_f2(Qf2, clus_ind, i, l, active_size)
        end
    end
    # remove all nodes of the blob from the scope
    for j in reverse(clus_ind[2:n_in_clus])
        # swap j with last active position using whole column/row operations
        if j != active_size
            D[:, j], D[:, active_size] = D[:, active_size], D[:, j]
            D[j, :], D[active_size, :] = D[active_size, :], D[j, :]
            reorder_covar_f2(Qf2, j, active_size)
            active_nodes[j] = active_nodes[active_size]
            clusters[j] = clusters[active_size]
        end
        active_size -= 1
    end
    # inactivate last active position: okay bc k replaces i != active_size
    active_nodes[i] = node_k

    return active_size
end

function nj_treeofblobs!(
    D::Matrix{Float64},
    names::AbstractVector{<:AbstractString}=String[];
    force_nonnegative_edges::Bool=false,
    threshold::Float64=0.01,
    nloci::Int=1000
)
    nj_treeofblobs!(D, zeros(0,0), names; force_nonnegative_edges = force_nonnegative_edges, threshold = threshold, nloci = nloci)
end

"""
    nj_treeofblobs!(D::Matrix{Float64}, names::AbstractVector{<:AbstractString}=String[];
        force_nonnegative_edges::Bool=false, epsilon::Float64=sqrt(eps(Float64))

Construct a tree of blob from the input distance matrix and
vector of names (as strings), using an adaptation of the the Neighbour-Joinging algorithm.
The order of the `names` argument should match that of the row (and
column) of the matrix `D`.
All branch lengths are arbitrarilly set to 1.0.
epsilon is the tolerance value for calling blobs.

Warning: `D` is modified.
See [`nj`](@ref) for the classical neighbor joining algorithm.
"""
function nj_treeofblobs!(
    D::Matrix{Float64},
    Qf2::Matrix{Float64},
    names::AbstractVector{<:AbstractString}=String[];
    force_nonnegative_edges::Bool=false,
    threshold::Float64=0.01,
    nloci::Int=1000,
    lambda::Real=1e-10
)
    check_distance_matrix(D)
    n = size(D, 1)    # number of species
    active_size = n   # tracks the current active size of matrix
    if isempty(names)
        names = string.(1:n)
    end
    check_variance_matrix(Qf2, n)
    clusters = collect(1:n)
    nclust = n
    Dclust = copy(D)
    # create empty network with n unconnected leaf nodes
    nodes = map(function(i)
                node = Node(i, true)
                node.name = names[i]
                return node
                end,
                1:n)
    net = HybridNetwork(nodes, Edge[])
    active_nodes = copy(nodes)
    neglenp = 0  # number of negative edge lengths

    while active_size > 2
        # compute Q matrix and find pair with minimum Q
        dist_clusters!(Dclust, D, active_size, clusters)
        ic1, ic2, _ = Q_findpair(Dclust, nclust)
        uniclust = unique(clusters[1:active_size])
        c1 = uniclust[ic1]
        c2 = uniclust[ic2]
        clus_ind = (clusters[1:active_size] .== c1) .|| (clusters[1:active_size] .== c2)
        if ((active_size - sum(clus_ind)) >= 2)
            if size(Qf2, 1) == 0
                blobpval = cutblobpvalue(D, clus_ind, nloci)
            else
                blobpval = cutblobpvalue(D, Qf2, clus_ind, active_size, lambda)
            end
        else # check that there are more than 2 tips left : otherwise, error when finding pairs
            blobpval = 1
            clusters[1:active_size] .= c1
            clus_ind = clusters[1:active_size] .== c1
        end
        clusterize!(clusters, c1, c2)
        nclust = nclust - 1
        # compute the length of the edges
        if (blobpval >= threshold)
            if size(Qf2, 1) == 0
                active_size = join_blob!(D, net, active_nodes, active_size, clusters, clus_ind)
            else
                active_size = join_blob!(D, Qf2, net, active_nodes, active_size, clusters, clus_ind)
            end
        end
    end

    # base case: 2 leaves
    if active_size >= 2
        join_12!(net, active_nodes, 1.0)
    end
    net.rooti = net.numnodes # root = last created node, which is internal
    net.isrooted = true # because rooti is valid, and edges are directed away from it
    if neglenp > 0
        infostr = (force_nonnegative_edges ?
                   "$neglenp branch(es) of negative length, reset to 0" :
                   "$neglenp branch(es) of negative length" )
        @info infostr
    end
    return net
end

"""
    nj_treeofblobs(D::DataFrame; force_nonnegative_edges::Bool=false, epsilon::Float64=sqrt(eps(Float64))

Construct a tree of blob from a distance matrix by neighbor joining, where
`D` is a `DataFrame` of the distance matrix, with taxon names taken
from the header of the data frame.
The rows are assumed to correspond to tips in the tree in the same order
as they do in columns.
All branch lengths are arbitrarilly set to 1.0.
epsilon is the tolerance value for calling blobs.

See [`nj_treeofblobs!`](@ref) for using a matrix as input,
and [`nj`](@ref) for the classical neighbor joining algorithm.
"""
function nj_treeofblobs(D::DataFrame; force_nonnegative_edges::Bool=false, threshold::Float64=0.01, nloci::Int=1000, lambda::Real=1e-5)
    nj_treeofblobs!(Matrix{Float64}(D), names(D); force_nonnegative_edges=force_nonnegative_edges, threshold = threshold, nloci = nloci, lambda = lambda)
end