function check_distance_matrix(D::Matrix{Float64})
    size(D, 1) > 1 || throw(DomainError(D, "Too few nodes"))
    if any(D .< 0)
        throw(DomainError(D, "Negative entries in distance matrix"))
    end
    if any(diag(D) .!= 0.0)
        throw(DomainError(D, "Diagonal entry not 0"))
    end
end

function nj!(D::Matrix{Float64})
    check_distance_matrix(D)
    n = size(D, 1)              # number of species
    # create empty network with n unconnected leaf nodes
    nodes = [ Node(i, true) for i = 1:n ]
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
        newidx = filter!(u->u!=j, collect(1:n))
        D = view(D, newidx, newidx)

        # update active_nodes
        active_nodes[i] = node_k

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

D = float([0 5 9 9 8; 5 0 10 10 9; 9 10 0 8 7; 9 10 8 0 3; 8 9 7 3 0])
