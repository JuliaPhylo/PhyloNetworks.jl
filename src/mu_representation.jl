"""
    EdgeMuRepresentation

A representation of the μ vector for edges in a network.
The μ vector is a vector of integers representing the number of paths from the edge to each leaf.
The order of the leaves in the μ vector is determined by the provided tiplabels vector.

# Arguments

- `mu_map`: A dictionary mapping each edge to its μ vector
- `tip_labels`: A vector of strings representing the labels of the leaves
"""
struct EdgeMuRepresentation
    "Dictionary: maps edge number => μ vector for that edge"
    mu_map::Dict{Int, Tuple{Vararg{Tuple{Vararg{Int}}}}}
    "vector of tip (leaf) labels: listed in the same order as in the μ vector"
    tip_labels::Vector{String}
end

function Base.show(io::IO, obj::EdgeMuRepresentation)
    disp = "$(typeof(obj))\n$(length(obj.tip_labels)) taxa, $(length(obj.mu_map)) edges,"
    disp *= "\ntaxon order in μ vectors:\n" * string(obj.tip_labels)
    disp *= "\nedge number => μ vector map:\n" * string(obj.mu_map)
    println(io, disp)
end

"""
    @Override ==(a::EdgeMuRepresentation, b::EdgeMuRepresentation)

Check if two EdgeMuRepresentation objects are equal using the distance function.
"""
function ==(a::EdgeMuRepresentation, b::EdgeMuRepresentation)
    return edge_vec_distance(a,b)==0
end

"""
    NodeMuRepresentation

A representation of the μ vector for nodes in a network.
The μ vector is a vector of integers representing the number of paths from the node to each leaf.
The order of the leaves in the μ vector is determined by the provided tiplabels vector.

# Arguments
- `mu_map`: A dictionary mapping each node to its μ vector
- `tip_labels`: A vector of strings representing the labels of the leaves
"""
struct NodeMuRepresentation
    "Dictionary: maps node number => μ vector for that node"
    mu_map::Dict{Int, Tuple{Vararg{Int}}}
    "vector of tip (leaf) labels: listed in the same order as in the μ vector"
    tip_labels::Vector{String}
end
function Base.show(io::IO, obj::NodeMuRepresentation)
    disp = "$(typeof(obj))\n$(length(obj.tip_labels)) taxa, $(length(obj.mu_map)) nodes,"
    disp *= "\ntaxon order in μ vectors:\n" * string(obj.tip_labels)
    disp *= "\nnode number => μ vector map:\n" * string(obj.mu_map)
    println(io, disp)
end

"""
    @Override ==(a::EdgeMuRepresentation, b::EdgeMuRepresentation)
    Check if two NodeMuRepresentation objects are equal using the distance function.
"""
function ==(a::NodeMuRepresentation, b::NodeMuRepresentation)
    return node_vec_distance(a,b)==0
end


"""
    node_mu_vectors(Network, tiplabels, preprocess))

Warning: Network cannot have multiple roots

Get the μ vector for each node in the network.
The μ vector is a vector of integers representing the number of paths from the node to each leaf.
The order of the leaves in the μ vector is determined by the provided tiplabels vector.

# Arguments
- `Network`: A HybridNetwork object
- `tiplabels`: A vector of strings representing the labels of the leaves
- `preprocess`: A boolean indicating whether to preprocess the network (default: true)
# Returns
- `nodemurepresentation`:  An NodeMuRepresentation object containing the μ vectors for each node in the network
"""
function node_mu_vectors(Net::HybridNetwork, labels::Vector{String}, preprocess=true )

    Ntiplabels = tiplabels(Net)

    # warn if labels in the network are not in the provided labels vector
    for l in Ntiplabels
        if !(l in labels)
            @warn("Leaf label $(l) not in the provided labels. will be missing from μ vectors")
        end
    end

    #check if labels are repeated in the provided labels vector or in the network
    if length(labels) != length(Set(labels))
        error("Tip Labels are repeated in the provided labels vector.")
    end
    if length(Ntiplabels) != length(Set(Ntiplabels))
        error("Tip Labels are repeated in the network.")
    end


    # preprocess the network if needed
    if preprocess
        directedges!(Net)
        preorder!(Net)
    end

    # dictionary to store the mu values for each node
    mu_map = Dict{Int,Vector{Int}}([node.number => zeros(length(Net.leaf)) for node in Net.node])

    # dictionary to store the mapping of leaves to their indices
    label_map = Dict{String, Int}(leaf => i for (i, leaf) in enumerate(labels))

    # use vec_node indices

    for currnode in reverse(Net.vec_node)

        # if leaf, only one way to reach a leaf
        if currnode.leaf
            if haskey(label_map, currnode.name)
                # get the index of the leaf in the labels vector
                index = label_map[currnode.name]
                mu_map[currnode.number][index]= 1
            end
        # add the number of paths to the leaves of the child to parent
        else
            for curredge in currnode.edge
                c = getchild(curredge)
                if c ==  currnode
                    continue
                end
                mu_map[currnode.number] .+= mu_map[c.number]
            end
        end
    end 
    node_mu_dict = Dict{Int, Tuple{Vararg{Int}}}([node.number => Tuple(mu_map[node.number]) for node in Net.node])
    

    nodemurepresentation = NodeMuRepresentation(node_mu_dict, labels)
    return nodemurepresentation       
end

"""
    edge_mu_vectors(Network, tiplabels, preprocess)

Warning: Network cannot have multiple roots

Get the μ vector for each edge in the network.
The μ vector is a vector of integers representing the number of paths from the edge to each leaf.
The order of the leaves in the μ vector is determined by the provided tiplabels vector.

# Arguments
- `Network`: A HybridNetwork object
- `tiplabels`: A vector of strings representing the labels of the leaves
- `preprocess`: A boolean indicating whether to preprocess the network (default: true)
# Returns
- `edgemurepresentation`: An EdgeMuRepresentation object containing the μ vectors for each edge in the network
"""
function edge_mu_vectors(Net::HybridNetwork, labels::Vector{String}, preprocess=true)

    # Get the mu values for the nodes
    nodemu = node_mu_vectors(Net,labels, preprocess).mu_map
    edge_mu = Dict{Int, Tuple{Vararg{Tuple{Vararg{Int64}}}}}()
    rho = collect(nodemu[getroot(Net).number])
    for curredge in Net.edge
        if !curredge.hybrid && curredge.containroot && !getchild(curredge).leaf
            edge_mu[curredge.number] = (nodemu[getchild(curredge).number],Tuple(rho - collect(nodemu[getparent(curredge).number])))
        else
            edge_mu[curredge.number] = (nodemu[getchild(curredge).number],)
        end
    end

    # Create a mu representation object
    edgemurepresentation = EdgeMuRepresentation(edge_mu, labels)

    return edgemurepresentation
end

"""
    node_vec_distance(Net1::NodeMuRepresentation, Net2::NodeMuRepresentation)

Calculate the distance between two node based μ vectors.

# Arguments
- `Net1`: A NodeMuRepresentation object
- `Net2`: A NodeMuRepresentation object

# Returns
- `dist`: The distance between the two node based μ vectors.
"""
function node_vec_distance(Net1::NodeMuRepresentation, Net2::NodeMuRepresentation)
    #convert to not use multiet
    m1 = Multiset(collect(values(Net1.mu_map)))
    m2 = Multiset(collect(values(Net2.mu_map)))
    d1 = setdiff(m1, m2)
    d2 = setdiff(m2, m1)
    dist = length(d1) + length(d2)
    return dist
end

"""
    distance_node(Net1::HybridNetwork, Net2::HybridNetwork)

Distance between two networks based on their node based μ vectors.

Warning: Networks cannot have multiple roots.

# Arguments
- `Net1`: A HybridNetwork object
- `Net2`: A HybridNetwork object

# Returns
- `dist`: The distance between the two networks based on their node based μ vectors.
"""
function network_node_mu_distance(Net1::HybridNetwork, Net2::HybridNetwork, preprocess=true)
    # make sure mu vectors are consistent
    if length(tiplabels(Net1)) > length(tiplabels(Net2))
        labels = tiplabels(Net1)
    else
        labels = tiplabels(Net2)
    end

    # Get the mu values for the nodes
    nodemu1 = node_mu_vectors(Net1,labels, preprocess)
    nodemu2 = node_mu_vectors(Net2,labels, preprocess)
    return node_vec_distance(nodemu1, nodemu2)
    
end

"""
    edge_vec_distance(Net1::EdgeMuRepresentation, Net2::EdgeMuRepresentation)

Calculate the distance between two edge based μ vectors.

# Arguments
- `Net1`: A EdgeMuRepresentation object
- `Net2`: A EdgeMuRepresentation object

# Returns
- `dist`: The distance between the two edge based μ vectors.
"""
function edge_vec_distance(Net1::EdgeMuRepresentation, Net2::EdgeMuRepresentation)

    m1 = Multiset(collect(values(Net1.mu_map)))
    m2 = Multiset(collect(values(Net2.mu_map)))
    d1 = setdiff(m1, m2)
    d2 = setdiff(m2, m1)
    dist = length(d1) + length(d2)
    return dist
end
"""
    distance_edge(Net1::HybridNetwork, Net2::HybridNetwork)

Distance between two networks based on their edge based μ vectors.

Warning: Networks cannot have multiple roots.

# Arguments
- `Net1`: A HybridNetwork object
- `Net2`: A HybridNetwork object

# Returns
- `dist`: The distance between the two networks based on their edge based μ vectors.
"""
function network_edge_mu_distance(Net1::HybridNetwork, Net2::HybridNetwork, preprocess=true)
    # make sure mu vectors are consistent
    if length(tiplabels(Net1)) > length(tiplabels(Net2))
        labels = tiplabels(Net1)
    else
        labels = tiplabels(Net2)
    end

    # Get the mu values for the nodes
    edgemu1 = edge_mu_vectors(Net1,labels, preprocess)
    edgemu2 = edge_mu_vectors(Net2,labels, preprocess)
    return edge_vec_distance(edgemu1, edgemu2)
    
end