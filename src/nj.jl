function check_distance_matrix(D::Matrix{<:Real})
    size(D, 1) > 1 || throw(DomainError(D, "Too few nodes"))
    if any(D .< 0)
        throw(DomainError(D, "Negative entries in distance matrix"))
    end
    if any(diag(D) .!= 0.0)
        throw(DomainError(D, "Diagonal entry not 0"))
    end
end

@doc (@doc nj) nj!
"""
    nj(df::DataFrame)
    nj!(D::Matrix{<:Real}, names::AbstractVector{String}=String[])

The `nj!` function takes the distance matrix and the vector of names
(as strings) as input, and constructs the corresponding tree.  The
order of the `names` argument should match that of the row (and
column) of the matrix `D`.

"""

function nj!(D::Matrix{<:Real}, names::AbstractVector{String}=String[])

    check_distance_matrix(D)
    n = size(D, 1)              # number of species

    # when no names arg is supplied
    if isempty(names)
        names = fill("", n)
    end

    # create empty network with n unconnected leaf nodes
    nodes = map(function(i)
                node = Node(i, true)
                node.name = names[i]
                return node
                end,
                1:n)

    net = HybridNetwork(nodes, Edge[])
    # an array of Node s.t. active_nodes[i] would correspond to the
    # ith entry in distance matrix D at each iteration
    active_nodes = nodes

    while n > 2
        # compute Q matrix and find min
        # warning: here i and j are indeces for D, not actual id's
        sums = sum(D, dims = 1)
        cur = Inf
        min_index = (0, 0)
        for i = 1:n
            for j = i:n
                if j != i
                    qij = (n-2) * D[i,j] - sums[i] - sums[j]
                    if qij < cur
                        cur = qij
                        min_index = (i, j)
                    end
                end
            end
        end

        # connect the nodes, compute the length of the edge
        (i, j) = min_index
        dik = D[i,j] / 2 + (sums[i] - sums[j]) / (2 * (n - 2))
        djk = D[i,j] - dik

        # create new edges and node, update tree
        edgenum = net.numEdges
        eik = Edge(edgenum + 1, dik)
        ejk = Edge(edgenum + 2, djk)
        node_k = Node(net.numNodes+1, false, false, [eik, ejk]) # new node
        node_i = active_nodes[i]
        node_j = active_nodes[j]
        setNode!(eik, Node[node_i, node_k])
        setNode!(ejk, Node[node_j, node_k])
        setEdge!(node_i, eik)
        setEdge!(node_j, ejk)
        pushEdge!(net, eik)
        pushEdge!(net, ejk)
        pushNode!(net, node_k)

        # update map and D
        # replace D[l, i] with D[l, k], delete D[ , j]
        for l in 1:n
            if !(l in [i j])
                D[l, i] = (D[l,i] + D[j,l] - D[i,j]) / 2
                D[i, l] = D[l, i]
            end
        end
        newidx = filter!(u->u!=j, 1:n) # index 1:n\{j}
        D = view(D, newidx, newidx)

        # update active_nodes
        active_nodes[i] = node_k
        active_nodes = view(active_nodes, newidx)

        n = n - 1
    end

    # base case
    node1 = active_nodes[1]
    node2 = active_nodes[2]
    newedge = Edge(net.numEdges+1, D[1,2])
    setNode!(newedge, [node1, node2])
    setEdge!(node1, newedge)
    setEdge!(node2, newedge)
    pushEdge!(net, newedge)
    return net
end

"""
    `nj(D::DataFrame)`

Construct a tree from a distance matrix by neighbor joining, where
`D` is a `DataFrame` of the distance matrix, with taxon names taken
from the header of the data frame.

"""
function nj(D::DataFrame)
    nj!(convert(Matrix, D), string.(names(D)))
end
