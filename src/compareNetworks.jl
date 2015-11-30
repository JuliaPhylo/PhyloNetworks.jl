# functions to compare unrooted networks
# Claudia & Cecile November 2015

"""
`deleteHybridThreshold!(net::HybridNetwork,gamma::Float64)`

Deletes from a network all hybrid edges with heritability below a threshold gamma.

- if gamma<0.5: deletes     minor hybrid edges with gamma value <  threshold
- if gamma=0.5: deletes all minor hybrid edges (i.e gamma value <= threshold)

Warning: assumes that the network has branch lengths (for now) and correct 'isMajor' attributes.
"""
function deleteHybridThreshold!(net::HybridNetwork,gamma::Float64)
    gamma <= 0.5 || error("deleteHybridThreshold! called with gamma = $(gamma)>0.5")
    for(i = net.numHybrids:-1:1)
    # starting from last because net.hybrid changes as hybrids are removed. Empty range if 0 hybrids.
        hybedges = hybridEdges(net.hybrid[i]) # vector of 3 edges: major, minor, tree edge
        # remove minor hybrid edge if: gamma < threshold OR threshold=0.5=gamma
        # warning: no check that the minor edge actually has gamma <= 0.5.
        if(hybedges[2].gamma < gamma || gamma == 0.5)
            deleteHybrid!(net.hybrid[i],net,true,false)
            # warning: deleteHybrid assumes non-missing branch lengths.
            # note: deleteHybridizationUpdate! requires a level-1 network with 
            # all attributes updated: inCycle, containRoot, etc., as after being 
            # read by readTopologyLevel1.
        end
    end
end

# extracts the two networks that simplify a given network at a given hybrid node:
#          deleting either one or the other parent hybrid edge.
# the original network is modified: minor edge removed.
# returns one HybridNetwork object (the network with major edge removed)
function displayedNetworks!(net::HybridNetwork, node::Node)
    node.hybrid || error("will not extract networks from tree node $(node.number)")
    ind = getIndex(node,net)
    netmin = deepcopy(net)
    hybedges = hybridEdges(node)
    majorgamma = hybedges[1].gamma
    deleteHybrid!(   net.node[ind],net   ,true,false)
    hybedges = hybridEdges(netmin.node[ind]) # hybrid edges in netmin
    setGamma!(hybedges[2],majorgamma) # set major gamma to minor edge (to delete old major = new minor)
    deleteHybrid!(netmin.node[ind],netmin,true,false)
    return netmin
end

"""
`displayedTrees(net::HybridNetwork, gamma::Float64)`

Warning: for now, requires the network to have non-missing positive branch lengths and no missing heritabilities.

Extracts all trees displayed in a network, following hybrid edges
with heritability >= gamma threshold (or >0.5 if threshold=0.5)
and ignoring any hybrid edge with heritability lower than gamma.
Returns an array of trees, as HybridNetwork objects.
"""
function displayedTrees(net0::HybridNetwork, gamma::Float64)
    trees = HybridNetwork[]
    net = deepcopy(net0)
    deleteHybridThreshold!(net,gamma)
    displayedTrees!(trees,net)
    return trees # should have length 2^net.numHybrids
end

"""
`majorTree(net::HybridNetwork)`

Warning: for now, requires the network to have non-missing positive branch lengths and no missing heritabilities.

Extracts the major tree displayed in a network, keeping the major edge and dropping the minor edge at each hybrid node.
Returns a HybridNetwork object.
"""
majorTree(net::HybridNetwork) = displayedTrees(net,0.5)[1]


# expands current list of trees, with trees displayed in a given network
function displayedTrees!(trees::Array{HybridNetwork,1},net::HybridNetwork)
    if (isTree(net))
        push!(trees, net)
    else
        netmin = displayedNetworks!(net,net.hybrid[1])
        displayedTrees!(trees,net)
        displayedTrees!(trees,netmin)
    end
end

"""
`minorTreeAt(net::HybridNetwork, hybindex::Int64)`

Warning: for now, requires the network to have non-missing positive branch lengths and no missing heritabilities.

Extracts the tree displayed in the network, following the major hybrid edge at each hybrid node, except at the ith hybrid node (i=hybindex), where the minor hybrid edge is kept instead of the major hybrid edge.
"""
function minorTreeAt(net::HybridNetwork, hybindex::Int64)
    hybindex <= length(net.hybrid) || error("network has fewer hybrid nodes than index $(hybindex).")
    tree = deepcopy(net)
    hybedges = hybridEdges(tree.hybrid[hybindex])
    majorgamma = hybedges[1].gamma
    setGamma!(hybedges[2],majorgamma) # set major gamma to minor edge (to delete old major = new minor)
    deleteHybrid!(tree.hybrid[hybindex],tree,true,false) # major edge at hybrid removed.
    return majorTree(tree) # all remaining minor edges removed: now it's a tree.
end


"""
`hardwiredClusterDistance(net1::HybridNetwork, net2::HybridNetwork, rooted::Bool)`

Warning: the current function is only implemented for trees, i.e. networks with 0 hybridizations.

Takes 2 networks and returns their hardwired cluster distance, that is, the number of hardwired clusters found in one network and not in the other. Note that this is not a distance per se on the full space of hybrid networks: there are pairs of different networks for which this measure is 0. But it is a distance on some network subspaces.

If the 2 networks are trees, this is the Robinson-Foulds distance.
If rooted=false, the trees are considered unrooted.
"""
function hardwiredClusterDistance(net1::HybridNetwork, net2::HybridNetwork, rooted::Bool)
    (net1.numHybrids == 0 && net2.numHybrids == 0) || error("hardwiredClusterDistance not implemented yet for non-tree networks.")
    taxa = sort!(ASCIIString[net1.leaf[i].name for i in 1:net1.numTaxa])
    length(setdiff(taxa, ASCIIString[net2.leaf[i].name for i in 1:net2.numTaxa])) == 0 ||
        error("net1 and net2 do not share the same taxon set. Please prune networks first.")
    nTax = length(taxa)
    M1 = tree2Matrix(net1, taxa, rooted=rooted)
    M2 = tree2Matrix(net2, taxa, rooted=rooted)
    dis = 0

    for (i1=1:size(M1)[1])
        found = false
        m1 = 1-M1[i1,2:(nTax+1)]
        for (i2=1:size(M2)[1])
            if (M1[i1,2:(nTax+1)] == M2[i2,2:(nTax+1)] ||
                  ( !rooted && m1 == M2[i2,2:(nTax+1)])     )
                found = true
                break
            end
        end
        if (!found)
            dis += 1
        end
    end # (size(M1)[1] - dis) edges have been found in net2, dis edges have not.
    # so size(M2)[1] - (size(M1)[1] - dis) edges in net2 are not in net1.
    dis + dis + size(M2)[1] - size(M1)[1]
end
