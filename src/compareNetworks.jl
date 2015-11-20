# functions to compare unrooted networks
# Claudia November 2015

# function to delete hybrid edges below a gamma threshold in a given network
function deleteHybridThreshold!(net::HybridNetwork,gamma::Float64)
    for(n in net.hybrid)
        hybedges = hybridEdges(n)
        if(hybedges[2].gamma < gamma) #lower gamma than threshold in minor edge
            deleteHybridizationUpdate!(net,n,false,false)
        end
    end
end

# function to extract the two trees corresponding to a
# given hybridization from a network
# returns vector of HybridNetwork object (with two trees)
function extractTrees(net::HybridNetwork, node::Node)
    node.hybrid || error("cannot extract trees from a hybridization on a tree node $(node.number)")
    ind = getIndex(node,net)
    treemaj = deepcopy(net)
    treemin = deepcopy(net)
    hybedges = hybridEdges(node)
    majorgamma = hybedges[1].gamma
    deleteHybridizationUpdate!(treemaj,treemaj.node[ind],false,false)
    hybedges = hybridEdges(treemin.node[ind]) #hybrid edges in treemin
    setGamma!(hybedges[2],majorgamma) #set major gamma to minor edge (to delete minor)
    deleteHybridizationUpdate!(treemin,treemin.node[ind],false,false)
    return treemaj,treemin
end

# function to extract all tree from a network that have hybrid edge
# with gamma greater than threshold
function extractTrees(net0::HybridNetwork,gamma::Float64)
    if(!isTree(net0))
        net = deepcopy(net0)
        deleteHybridThreshold!(net,gamma)
        trees = HybridNetwork[]
        if(!isTree(net)) #net is still network
            for(n in net.hybrid)
                treemaj,treemin = extractTrees(net,n)
                push!(trees, treemaj)
                push!(trees, treemin)
            end
        end
        return trees
    end
end



