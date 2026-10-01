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

# Un-normalized f4 statitic for tips (a,b,c,d) in D
function getf4(D::Matrix{<:Real}, a::Int, b::Int, c::Int, d::Int)
    return D[a,d] + D[b,c] - D[a,c] - D[b,d]
end

# If (ab|cd) verify the 4 point condition, we should have:
# getf4(D,a,b,c,d) = 0, getf4(D,a,c,b,d) > 0 and getf4(D,a,d,b,c) > 0
function getminf4(D::Matrix{<:Real}, a::Int, b::Int, c::Int, d::Int)
    return min(getf4(D,a,b,c,d), getf4(D,a,c,b,d), getf4(D,a,d,b,c))
end

function getminf4(D::Matrix{<:Real}, t::Tuple)
    return getminf4(D, t[1], t[2], t[3], t[4])
end

# Get the min of getminf4 over all possible combinations (a,b,c,d),
# with a,b in clus_ind, and c,d not in clus_ind
function minblobf4(D::Matrix{<:Real}, clus_ind::BitVector)
    minblobf4(D, findall(clus_ind), findall(.!clus_ind))
end

function minblobf4(D::Matrix{<:Real}, ind_A::Vector{Int}, ind_B::Vector{Int})
    # comb_AB is a 4 dim vector, with upper diagonal having the
    # all the (a,b,c,d) with a<b in ind_A and c<d in ind_B
    # (There might be better / more efficient ways ?)
    comb_AB = collect(Iterators.product(ind_A, ind_A, ind_B, ind_B))
    minblob = Inf
    for a in 1:(size(ind_A, 1))
        for b in a:(size(ind_A, 1))
            if (b != a)
                for c in 1:(size(ind_B, 1))
                    for d in c:(size(ind_B, 1))
                        if (d != c)
                            tmp = getminf4(D, comb_AB[a,b,c,d])
                            if tmp < minblob
                                minblob = tmp
                            end
                        end
                    end
                end
            end
        end
    end
    return minblob
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
    names::AbstractVector{<:AbstractString}=String[];
    force_nonnegative_edges::Bool=false,
    epsilon::Float64=sqrt(eps(Float64))
)
    check_distance_matrix(D)
    n = size(D, 1)    # number of species
    active_size = n   # tracks the current active size of matrix
    if isempty(names)
        names = string.(1:n)
    end
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
            minblob = minblobf4(D, clus_ind)
        else # check that there are more than 2 tips left : otherwise, error when finding pairs
            minblob = 1
            clusters[1:active_size] .= c1
            clus_ind = clusters[1:active_size] .== c1
        end
        clusterize!(clusters, c1, c2)
        nclust = nclust - 1
        # compute the length of the edges
        if (minblob >= -epsilon)
            active_size = join_blob!(D, net, active_nodes, active_size, clusters, clus_ind)
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
function nj_treeofblobs(D::DataFrame; force_nonnegative_edges::Bool=false, epsilon::Float64=sqrt(eps(Float64)))
    nj_treeofblobs!(Matrix{Float64}(D), names(D); force_nonnegative_edges=force_nonnegative_edges, epsilon = epsilon)
end