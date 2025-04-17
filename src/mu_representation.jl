using PhyloNetworks, DataStructures
import Base: ==
using Multisets

"""
    edge_mu_rep

    A representation of the μ vector for edges in a network.
    The μ vector is a vector of integers representing the number of paths from the edge to each leaf.
    The order of the leaves in the μ vector is determined by the provided tiplabels vector.

    # Arguments
    - `mu_map`: A dictionary mapping each edge to its μ vector
    - `tip_labels`: A vector of strings representing the labels of the leaves
"""
struct edge_mu_rep
    mu_map::Dict{Int, Vector{Int}}   # Node => μ vector
    tip_labels::Vector{String}        # Tip labels (leaf names)
end
"""
    @Override ==(a::edge_mu_rep, b::edge_mu_rep)
    Check if two edge_mu_rep objects are equal using the distance function.
"""
# function ==(a::edge_mu_rep, b::edge_mu_rep)
#     return a.mu_map == b.mu_map && a.tip_labels == b.tip_labels
# end

"""
    node_mu_rep

    A representation of the μ vector for nodes in a network.
    The μ vector is a vector of integers representing the number of paths from the node to each leaf.
    The order of the leaves in the μ vector is determined by the provided tiplabels vector.

    # Arguments
    - `mu_map`: A dictionary mapping each node to its μ vector
    - `tip_labels`: A vector of strings representing the labels of the leaves
"""
struct node_mu_rep
    mu_map::Dict{Int, Vector{Int}}   # Node => μ vector
    tip_labels::Vector{String}        # Tip labels (leaf names)
end
"""
    @Override ==(a::edge_mu_rep, b::edge_mu_rep)
    Check if two node_mu_rep objects are equal using the distance function.
"""
function ==(a::node_mu_rep, b::node_mu_rep)
    return node_mu_distance(a,b)==0
end


"""
    node_mu(Network, tiplabels, preprocess))

Warning: Network cannot have multiple roots

Get the μ vector for each node in the network.
The μ vector is a vector of integers representing the number of paths from the node to each leaf.
The order of the leaves in the μ vector is determined by the provided tiplabels vector.

# Arguments
- `Network`: A HybridNetwork object
- `tiplabels`: A vector of strings representing the labels of the leaves
- `preprocess`: A boolean indicating whether to preprocess the network (default: true)
# Returns
- `mu_map`: A dictionary mapping each node to its μ vector
"""
function node_mu(N::HybridNetwork, labels::Vector{String}, preprocess=true )

    #TODO: map node number instead of edge? : cant access thennode from node number unless you build a map
    #TODO: instead of a map use of of the unused node variables and just return the network : maybe later confirm algo works first
    #TODO: get index for leaves :  done

    Ntiplabels = tipLabels(N)

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
        directEdges!(N)
        preorder!(N)
    end

    # dictionary to store the mu values for each node
    mu_map = Dict{Int,Vector{Int}}([node.number => zeros(length(N.leaf)) for node in N.node])

    # dictionary to store the mapping of leaves to their indices
    label_map = Dict{String, Int}(leaf => i for (i, leaf) in enumerate(labels))

    # use vec_node indices

    for currnode in N.vec_node
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

    mu = node_mu_rep(mu_map, labels)
    return mu       
end

"""
    edge_mu(Network, tiplabels, preprocess)

    Warning: Network cannot have multiple roots

    Get the μ vector for each edge in the network.
    The μ vector is a vector of integers representing the number of paths from the edge to each leaf.
    The order of the leaves in the μ vector is determined by the provided tiplabels vector.

    # Arguments
    - `Network`: A HybridNetwork object
    - `tiplabels`: A vector of strings representing the labels of the leaves
    - `preprocess`: A boolean indicating whether to preprocess the network (default: true)
    # Returns
    - `mu`: A MuRep object containing the μ vectors for each edge in the network
"""
function edge_mu(N::HybridNetwork, labels::Vector{String}, preprocess=true)

    # Get the mu values for the nodes
    nodemu = node_mu(N,labels, preprocess).mu_map
    # edge_mu = Dict{Int, Vector{Vector{Int}}}(edge.number => Vector{Vector{Int}}() for edge in N.edge)
    edge_mu = Dict{Int, Vector{Tuple{Int}}}
    rho = nodemu[getroot(net).number]
    for curredge in N.edge
        if !curredge.hybrid && curredge.containroot && !getchild(curredge).leaf
            edge_mu[curredge.number] = [Tuple(nodemu[getchild(curredge).number]),Tuple(rho - nodemu[getparent(curredge).number])]
        else
            edge_mu[curredge.number] = [Tuple(nodemu[getchild(curredge).number])]
        end
    end

    # Create a mu representation object
    mu = edge_mu_rep(edge_mu, labels)

    return mu
end

"""
    node_mu_distance(N1::node_mu_rep, N2::node_mu_rep)

    Calculate the distance between two μ vectors.

    # Arguments
    - `N1`: A MuRep object
    - `N2`: A MuRep object

    # Returns
    - `dist`: The distance between the two μ vectors.
"""
function node_mu_distance(N1::node_mu_rep, N2::node_mu_rep)
    #convert to not use multiet
    m1 = Multiset(collect(values(N1.mu_map)))
    m2 = Multiset(collect(values(N2.mu_map)))
    d1 = setdiff(m1, m2)
    d2 = setdiff(m2, m1)
    dist = length(d1) + length(d2)
    return dist
end

"""
    distance_node(N1::HybridNetwork, N2::HybridNetwork)

    Calculate the distance between two networks based on their μ vectors.

    Warning: Networks cannot have multiple roots.
    
    # Arguments
    - `N1`: A HybridNetwork object
    - `N2`: A HybridNetwork object

    # Returns
    - `dist`: The distance between the two networks based on their μ vectors.
"""
function distance_node(N1::HybridNetwork, N2::HybridNetwork, preprocess=true)
    # make sure mu vectors are consistent
    if length(tipLabels(N1)) > length(tipLabels(N2))
        labels = tipLabels(N1) 
    else
        labels = tipLabels(N2)
    end

    # Get the mu values for the nodes
    nodemu1 = node_mu(N1,labels, preprocess)
    nodemu2 = node_mu(N2,labels, preprocess)
    return node_mu_distance(nodemu1, nodemu2)
    
end

#TODO: Distance using edge mu