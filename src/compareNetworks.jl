# functions to compare unrooted networks
# Claudia & Cecile November 2015

"""
`deleteHybridThreshold!(net::HybridNetwork,gamma::Float64)`

Deletes from a network all hybrid edges with heritability below a threshold gamma.

- if gamma<0.5: deletes     minor hybrid edges with gamma value <  threshold
- if gamma=0.5: deletes all minor hybrid edges (i.e gamma value <= threshold)
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
`displayedTrees(net::HybridNetwork; gamma=0.0::Float64)`

Extracts all trees displayed in a network, following hybrid edges
with heritability >= gamma threshold (or >0.5 if threshold=0.5)
and ignoring any hybrid edge with heritability lower than gamma.
Returns an array of trees, as HybridNetwork objects.
gamma threshold = 0.0 by default.
"""
function displayedTrees(net0::HybridNetwork;gamma=0.0::Float64)
    trees = HybridNetwork[]
    net = deepcopy(net0)
    deleteHybridThreshold!(net,gamma)
    displayedTrees!(trees,net)
    return trees # should have length 2^net.numHybrids
end

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


