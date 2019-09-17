# first prototype a recursive version
# update the network once
function nj_rec!(D::Matrix{Float64}, # distance matrix
                 id::Vector{Int64}, # mapping for distance matrix
                 net::HybridNetwork    # network from last iteration
                 )
    n = length(id)
    # base case: if n = 2, connect nodes
    if n == 2
        node1 = net.node[getIndexNode(id[1], net)]
        node2 = net.node[getIndexNode(id[2], net)]
        newedge = Edge(net.numEdges+1, D[1,2])
        setNode!(newedge, [node1, node2])
        setEdge!(node1, newedge)
        setEdge!(node2, newedge)
        pushEdge!(net, newedge)
        return net
    end

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

    # connect the nodes, compute the length of the edge, update tree
    (i, j) = min_index
    new_node_id = maximum(id) + 1 # call this node k
    dik = D[i,j] / 2 + (sums[i] - sums[j]) / (2 * (n - 2))
    djk = D[i,j] - dik

    # create new edges and node, update network
    edgenum = net.numEdges
    eik = Edge(edgenum + 1, dik)
    ejk = Edge(edgenum + 2, djk)
    node_k = Node(new_node_id, false, false, [eik, ejk])
    node_i = net.node[getIndexNode(id[i], net)]
    node_j = net.node[getIndexNode(id[j], net)]
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
    D = D[newidx, newidx]
    id = filter!(u->u!=id[j], id) # delete node_j from id
    id[i] = new_node_id

    # continue recursively
    return nj_rec!(D, id, net)
end

function nj_proto!(D::Matrix{Float64})
    n = size(D, 1)              # number of species
    # create empty network with n unconnected leaf nodes
    nodes = Node[]
    for i in 1:n
        node = Node(i, true)
        push!(nodes, node)
    end
    net = HybridNetwork(nodes, Edge[])
    return nj_rec!(D, collect(1:n), net)
end

function nj!(D::Matrix{Float64})
    n = size(D, 1)              # number of species
    # create empty network with n unconnected leaf nodes
    nodes = Node[]
    for i in 1:n
        node = Node(i, true)
        push!(nodes, node)
    end
    net = HybridNetwork(nodes, Edge[])
    id = collect(1:n)

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

        # connect the nodes, compute the length of the edge, update tree
        (i, j) = min_index
        new_node_id = maximum(id) + 1 # call this node k
        dik = D[i,j] / 2 + (sums[i] - sums[j]) / (2 * (n - 2))
        djk = D[i,j] - dik

        # create new edges and node, update network
        edgenum = net.numEdges
        eik = Edge(edgenum + 1, dik)
        ejk = Edge(edgenum + 2, djk)
        node_k = Node(new_node_id, false, false, [eik, ejk])
        node_i = net.node[getIndexNode(id[i], net)]
        node_j = net.node[getIndexNode(id[j], net)]
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
        D = D[newidx, newidx]
        id = filter!(u->u!=id[j], id) # delete node_j from id
        id[i] = new_node_id

        n = n - 1
    end

    # base case
    node1 = net.node[getIndexNode(id[1], net)]
    node2 = net.node[getIndexNode(id[2], net)]
    newedge = Edge(net.numEdges+1, D[1,2])
    setNode!(newedge, [node1, node2])
    setEdge!(node1, newedge)
    setEdge!(node2, newedge)
    pushEdge!(net, newedge)
    return net
end

D = float([0 5 9 9 8; 5 0 10 10 9; 9 10 0 8 7; 9 10 8 0 3; 8 9 7 3 0])
