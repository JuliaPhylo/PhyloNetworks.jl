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
        # connect the nodes, compute the length of the edge
        (i, j) = min_index
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
        # create new edges and node, update tree
        edgenum = net.numEdges
        eik = Edge(edgenum + 1, dik)
        ejk = Edge(edgenum + 2, djk)
        node_k = Node(net.numNodes+1, false, false, [eik, ejk])
        node_i = active_nodes[i]
        node_j = active_nodes[j]
        setNode!(eik, Node[node_i, node_k])
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
    end

    # base case: 2 leaves
    if D[1,2] < 0.0
        neglenp += 1
        if force_nonnegative_edges D[1,2] = 0.0; end
    end
    node1 = active_nodes[1]
    node2 = active_nodes[2]
    newedge = Edge(net.numEdges+1, D[1,2])
    setNode!(newedge, [node1, node2]) # isChild1 = true by default: nodes are [child, parent]
    if node1.number == net.numNodes
        newedge.isChild1 = false # direct the edge: last created node -> other node
    end
    setEdge!(node1, newedge)
    setEdge!(node2, newedge)
    pushEdge!(net, newedge)

    net.root = net.numNodes # root = last created node, which is internal
    if neglenp > 0
        infostr = (force_nonnegative_edges ?
                   "$neglenp branches had negative lengths, reset to 0" :
                   "$neglenp branches have negative lengths" )
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
