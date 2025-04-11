using PhyloNetworks, DataStructures
import Base: ==
using Multisets


#TODO: override the == operator for the MuRep struct
struct MuRep
    mu_map::Dict{Node, Vector{Int}}   # Node => μ vector
    tip_labels::Vector{String}        # Tip labels (leaf names)
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
    mu_map = Dict{Node,Vector{Int}}([node => zeros(length(N.leaf)) for node in N.node])

    # dictionary to store the mapping of leaves to their indices
    label_map = Dict{String, Int}(leaf => i for (i, leaf) in enumerate(labels))



    for currnode in N.node
        # if leaf, only one way to reach a leaf
        if currnode.leaf
            if currnode.name in label_map
                # get the index of the leaf in the labels vector
                index = label_map[currnode.name]
                mu_map[currnode][index]= 1
            end
        # add the number of paths to the leaves of the child to parent
        else
            for curredge in currnode.edge
                c = getchild(curredge)
                if c ==  currnode
                    continue
                end
                mu_map[currnode] .+= mu_map[c]
            end
        end
    end

    

    return mu_map       
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

# TODO:make sure code can be adapted easily to multiple root components
function edge_mu(N::HybridNetwork, labels::Vector{String}, preprocess=true)

    # Get the mu values for the nodes
    nodemu = node_mu(N,labels, preprocess)
    edge_mu = Dict{Edge, Vector{Vector{Int}}}(edge => Vector{Vector{Int}}() for edge in N.edge)
    rho = nodemu[N.node[N.rooti]] # root node
    for curredge in N.edge
        if !curredge.hybrid && curredge.containroot && !getchild(curredge).leaf
            # switch to get child
            edge_mu[curredge] = [nodemu[getchild(curredge)],rho - nodemu[getparent(curredge)]]
        else
            edge_mu[curredge] = [nodemu[getchild(curredge)]]
        end
    end

    # Create a mu representation object
    mu = MuRep(edge_mu, labels)

    return mu
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
function distance_node(N1::HybridNetwork, N2::HybridNetwork)
    # make sure mu vectors are consistent
    if length(tipLabels(N1)) > length(tipLabels(N2))
        labels = tipLabels(N1) 
    else
        labels = tipLabels(N2)
    end

    # Get the mu values for the nodes
    nodemu1 = node_mu(N1,labels, preprocess)
    nodemu2 = node_mu(N2,labels, preprocess)

    m1 = Multiset(Tuple(v) for  v in nodemu1.mu_map.values())
    m2 = Multiset(Tuple(v) for  v in nodemu2.mu_map.values())
    d1 = setdiff(m1, m2)
    d2 = setdiff(m2, m1)
    dist = len(d1) + len(d2)


    return dist
end

#TODO: Distance using edge mu