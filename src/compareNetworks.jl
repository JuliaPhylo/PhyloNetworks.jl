# functions to compare unrooted networks
# Claudia & Cecile November 2015

# Traverses a tree postorder and modifies matrix M
# with edges as rows and species as columns (see tree2Matrix)
# S should be sorted --fixit: are you sure?
function traverseTree2Matrix!(node::Node, edge::Edge, ie::Vector{Int64}, M::Matrix{Int}, S::Union{Vector{ASCIIString},Vector{Int64}})
    child = getOtherNode(edge,node) # must not be a leaf
    indedge = ie[1]
    M[indedge,1] = edge.number
    ie[1] += 1 # mutable array: to modify edge index 'ie' outside function scope
    for(e in child.edge) #postorder traversal
        if(!isEqual(e,edge)) # assumes a tree here
            grandchild = getOtherNode(e,child)
            if (grandchild.leaf)
                indsp= 0
                try
                    indsp = getIndex(grandchild.name,S)
                catch
                    error("leaf $(grandchild.name) not in species list $(S)")
                end
                M[indedge,indsp+1] = 1 #indsp+1 bc first column is edge numbers
            else
                inde = ie[1];
                # inde = getIndex(e.number,M[:,1])
                traverseTree2Matrix!(child,e,ie,M,S)
                M[indedge,2:size(M,2)] |= M[inde,2:size(M,2)]
            end
        end
    end
end

# takes a tree and a list of species as input,
# and produces a matrix M with edges as rows and species as columns:
# Mij=1 if species j is descendant of edge i, 0 ow.
# allows for missing taxa:
# Mij=0 if species not present in tree. This is handled in calculateObsCFAll with sameTaxa function
function tree2Matrix(T::HybridNetwork, S::Union{Vector{ASCIIString},Vector{Int64}}; rooted=true::Bool)
    length(T.hybrid)==0 || error("tree2Matrix only works on trees. Network has $(T.numHybrids) hybrid nodes.")
    # sort!(S) # why sort 'taxa', again and again for each tree? Benefits?
    ne = length(T.edge)-T.numTaxa # number of internal branch lengths
    if (T.node[T.root].leaf)      # the root is a leaf: the 1 edge stemming from the root is an external edge
        ne += 1                   # and will need to get a row in the matrix, to be deleted later.
    end
    M = zeros(Int,ne,length(S)+1)
    # M[:,1] = sort!([e.number for e in T.edge])
    ie = [1] # index of next edge to be documented: row index in M
    for(e in T.node[T.root].edge)
        child = getOtherNode(e,T.node[T.root])
        if (!child.leaf)
            traverseTree2Matrix!(T.node[T.root],e,ie,M,S)
        end
    end
    if (!rooted && length(T.node[T.root].edge)<3)
        # remove first row of M: 1st edge at the root, duplicated edge
        # if the tree is to be considered as unrooted, or just leaf.
        M = M[2:size(M,1),:] # makes a copy, too bad.
    end
    return M
end

"""
`hardwiredClusters(net::HybridNetwork, S::Union{Vector{ASCIIString},Vector{Int64}})`


Returns a matrix describing all the hardwired clusters in a network.
Warning: clusters are rooted, so the root must be correct.
Warning: allows for missing taxa, with entries all 0.

Each row corresponds to one internal edge, that is, external edges are excluded.
If the root is a leaf node, the external edge to that leaf is included (first row).
Both parent hybrid edges to a given hybrid node only contribute a single row (they share the same hardwired cluster). 

- first column: edge number
- next columns: 0/1 values. 1=descendant of edge, 0=not a descendant, or missing taxon.
- last column:  10/11 values. 10=tree edge, 11=hybrid edge
"""
function hardwiredClusters(net::HybridNetwork, S::Union{Vector{ASCIIString},Vector{Int64}})
    ne = length(net.edge)-net.numTaxa # number of internal branch lengths
    ne -= length(net.hybrid)          # to remove duplicate rows for the 2 parent edges of each hybrid
    if (net.node[net.root].leaf)      # root is leaf: the 1 edge stemming from the root is an external edge
        ne += 1                       #               but needs to be included still (rooted clusters).
    end
    M = zeros(Int,ne,length(S)+2)
    ie = [1] # index of next edge to be documented: row index in M
    for(e in net.node[net.root].edge)
        hardwiredClusters!(net.node[net.root],e,ie,M,S)
    end
    return M
end

function hardwiredClusters!(node::Node, edge::Edge, ie::Vector{Int64}, M::Matrix{Int},
                            S::Union{Vector{ASCIIString},Vector{Int64}})
    child = getOtherNode(edge,node)

    !child.leaf || return 0 # do nothing if child is a leaf.

    if (edge.hybrid) # edge is a hybrid. Need to find its partner.
        (edge.isChild1 ? edge.node[1] == child : edge.node[2] == child) || error(
        "inconsistency during network traversal: node $(child.number) supposed to be child of hybrid edge $(edge.number), inconsistent with isChild1.")
        partner = nothing
        partnerVisited = true
        indpartner = 0
        for (e in child.edge)
            if (e.hybrid && e != edge && (e.isChild1 ? e.node[1] == child : e.node[2] == child))
                partner = e
                try
                    indpartner = getIndex(partner.number,M[:,1])
                catch
                    partnerVisited = false # will need to continue traversal
                end
                break # partner hybrid edge was found
            end
        end
        partner != nothing || error("partner hybrid edge not found for edge $(edge.number), child $(child.number)")
        !partnerVisited || return indpartner
    end

    indedge = ie[1]
    M[indedge,1] = edge.number
    M[indedge,end] = edge.hybrid ? 11 : 10
    ie[1] += 1 # mutable array

    for (e in child.edge) # postorder traversal
        if (e != edge && (!edge.hybrid || e!=partner)) # do not go back to (either) parent edge.
            grandchild = getOtherNode(e,child)
            if (grandchild.leaf)
                indsp = 0
                try
                    indsp = getIndex(grandchild.name,S)
                catch
                    error("leaf $(grandchild.name) not in species list $(S)")
                end
                M[indedge,indsp+1] = 1 #indsp+1 because first column is edge numbers
            else
                inde = hardwiredClusters!(child,e,ie,M,S)
                M[indedge,2:end-1] |= M[inde,2:end-1]
            end
        end
    end
    return indedge
end



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
`displayedNetworkAt!(net::HybridNetwork, node::Node)`

Warning: for now, requires the network to have non-missing positive branch lengths and no missing heritabilities.

Deletes all the minor hybrid edges, except at input node. The network is left with a single hybridization, and otherwise displays the same major tree as before.
"""
function displayedNetworkAt!(net::HybridNetwork, node::Node)
    node.hybrid || error("will not extract network from tree node $(node.number)")
    for(i = net.numHybrids:-1:1)
    # starting from last because net.hybrid changes as hybrids are removed. Empty range if 0 hybrids.
        net.hybrid[i] != node || continue
        deleteHybrid!(net.hybrid[i],net,true,false)
    end
end


"""
`hardwiredClusterDistance(net1::HybridNetwork, net2::HybridNetwork, rooted::Bool)`

Warning: the current function is only implemented for trees, i.e. networks with 0 hybridizations.

Takes 2 networks and returns their hardwired cluster distance, that is, the number of hardwired clusters found in one network and not in the other. Note that this is not a distance per se on the full space of hybrid networks: there are pairs of different networks for which this measure is 0. But it is a distance on some network subspaces.

If the 2 networks are trees, this is the Robinson-Foulds distance.
If rooted=false, the trees are considered unrooted.
"""
function hardwiredClusterDistance(net1::HybridNetwork, net2::HybridNetwork, rooted::Bool)
    bothtrees = (net1.numHybrids == 0 && net2.numHybrids == 0)
    rooted || bothtrees || error("unrooted hardwired cluster distance not defined for non-tree networks.")
    taxa = sort!(ASCIIString[net1.leaf[i].name for i in 1:net1.numTaxa])
    length(setdiff(taxa, ASCIIString[net2.leaf[i].name for i in 1:net2.numTaxa])) == 0 ||
        error("net1 and net2 do not share the same taxon set. Please prune networks first.")
    nTax = length(taxa)
    if (bothtrees) # even if rooted, different treatment at the root if root=leaf
        M1 = tree2Matrix(net1, taxa, rooted=rooted)
        M2 = tree2Matrix(net2, taxa, rooted=rooted)
    else
        M1 = hardwiredClusters(net1, taxa) # last row: 10/11 if tree/hybrid edge.
        M2 = hardwiredClusters(net2, taxa)
        #println("M1="); print(M1); println("\nM2="); print(M2); println("\n");
    end
    dis = 0

    for (i1=1:size(M1)[1])
        found = false
        m1 = 1-M1[i1,2:end] # going to the end: i.e. we want to match a tree edge with a tree edge
        for (i2=1:size(M2)[1])                                  # and hybrid edge with hybrid edge
            if (M1[i1,2:end] == M2[i2,2:end] ||
                  ( !rooted && m1 == M2[i2,2:end])     )
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
