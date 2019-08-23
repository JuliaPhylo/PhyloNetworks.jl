# functions to compare unrooted networks
# Claudia & Cecile November 2015

# Traverses a tree postorder and modifies matrix M
# with edges as rows and species as columns (see tree2Matrix)
# S should be sorted --fixit: are you sure?
function traverseTree2Matrix!(node::Node, edge::Edge, ie::Vector{Int}, M::Matrix{Int}, S::Union{Vector{String},Vector{Int}})
    child = getOtherNode(edge,node) # must not be a leaf
    indedge = ie[1]
    M[indedge,1] = edge.number
    ie[1] += 1 # mutable array: to modify edge index 'ie' outside function scope
    for e in child.edge #postorder traversal
        if(!isEqual(e,edge)) # assumes a tree here
            grandchild = getOtherNode(e,child)
            if grandchild.leaf
                indsp = findfirst(isequal(grandchild.name), S)
                indsp != nothing || error("leaf $(grandchild.name) not in species list $(S)")
                M[indedge,indsp+1] = 1 #indsp+1 bc first column is edge numbers
            else
                inde = ie[1];
                traverseTree2Matrix!(child,e,ie,M,S)
                M[indedge,2:size(M,2)] .|= M[inde,2:size(M,2)]
            end
        end
    end
end

# takes a tree and a list of species as input,
# and produces a matrix M with edges as rows and species as columns:
# Mij=1 if species j is descendant of edge i, 0 ow.
# allows for missing taxa:
# Mij=0 if species not present in tree. This is handled in calculateObsCFAll with sameTaxa function
function tree2Matrix(T::HybridNetwork, S::Union{Vector{String},Vector{Int}}; rooted=true::Bool)
    length(T.hybrid)==0 || error("tree2Matrix only works on trees. Network has $(T.numHybrids) hybrid nodes.")
    # sort!(S) # why sort 'taxa', again and again for each tree? Benefits?
    ne = length(T.edge)-T.numTaxa # number of internal branch lengths
    if (T.node[T.root].leaf)      # the root is a leaf: the 1 edge stemming from the root is an external edge
        ne += 1                   # and will need to get a row in the matrix, to be deleted later.
    end
    M = zeros(Int,ne,length(S)+1)
    # M[:,1] = sort!([e.number for e in T.edge])
    ie = [1] # index of next edge to be documented: row index in M
    for e in T.node[T.root].edge
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
`hardwiredClusters(net::HybridNetwork, S::Union{Vector{String},Vector{Int}})`


Returns a matrix describing all the hardwired clusters in a network.
Warnings: Clusters are rooted, so the root must be correct.
          Allows for missing taxa, with entries all 0.

Each row corresponds to one internal edge, that is, external edges are excluded.
If the root is a leaf node, the external edge to that leaf is included (first row).
Both parent hybrid edges to a given hybrid node only contribute a single row (they share the same hardwired cluster).

- first column: edge number
- next columns: 0/1 values. 1=descendant of edge, 0=not a descendant, or missing taxon.
- last column:  10/11 values. 10=tree edge, 11=hybrid edge
"""
function hardwiredClusters(net::HybridNetwork, S::Union{Vector{String},Vector{Int}})
    ne = length(net.edge)-net.numTaxa # number of internal branch lengths
    ne -= length(net.hybrid)          # to remove duplicate rows for the 2 parent edges of each hybrid
    if (net.node[net.root].leaf)      # root is leaf: the 1 edge stemming from the root is an external edge
        ne += 1                       #               but needs to be included still (rooted clusters).
    end
    M = zeros(Int,ne,length(S)+2)
    ie = [1] # index of next edge to be documented: row index in M
    for e in net.node[net.root].edge
        hardwiredClusters!(net.node[net.root],e,ie,M,S)
    end
    return M
end

function hardwiredClusters!(node::Node, edge::Edge, ie::Vector{Int}, M::Matrix{Int},
                            S::Union{Vector{String},Vector{Int}})
    child = getOtherNode(edge,node)

    !child.leaf || return 0 # do nothing if child is a leaf.

    if (edge.hybrid) # edge is a hybrid. Need to find its partner.
        (edge.isChild1 ? edge.node[1] == child : edge.node[2] == child) || error(
        "inconsistency during network traversal: node $(child.number) supposed to be child of hybrid edge $(edge.number), inconsistent with isChild1.")
        partner = nothing
        partnerVisited = true
        indpartner = 0
        for e in child.edge
            if (e.hybrid && e != edge && (e.isChild1 ? e.node[1] == child : e.node[2] == child))
                partner = e
                indpartner = findfirst(isequal(partner.number), M[:,1])
                if isnothing(indpartner)
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

    for e in child.edge # postorder traversal
        if (e != edge && (!edge.hybrid || e!=partner)) # do not go back to (either) parent edge.
            grandchild = getOtherNode(e,child)
            if grandchild.leaf
                indsp = findfirst(isequal(grandchild.name), S)
                indsp != nothing || error("leaf $(grandchild.name) not in species list $(S)")
                M[indedge,indsp+1] = 1 #indsp+1 because first column is edge numbers
            else
                inde = hardwiredClusters!(child,e,ie,M,S)
                M[indedge,2:end-1] .|= M[inde,2:end-1]
            end
        end
    end
    return indedge
end


"""
    hardwiredCluster(edge::Edge,taxa::Union{Vector{String},Vector{Int}})
    hardwiredCluster!(v::Vector{Bool},edge::Edge,taxa::Union{Vector{String},Vector{Int}})
    hardwiredCluster!(v::Vector{Bool},edge::Edge,taxa::Union{Vector{String},Vector{Int}},
                      visited::Vector{Int})

Calculate the hardwired cluster of `node`, coded a vector of booleans:
true for taxa that are descendent of nodes, false for other taxa (including missing taxa).

The node should belong in a rooted network for which isChild1 is up-to-date.
Run directEdges! beforehand. This is very important, otherwise one might enter an infinite loop,
and the function does not test for this.

visited: vector of node numbers, of all visited nodes.

# Examples:
```jldoctest
julia> net5 = "(A,((B,#H1),(((C,(E)#H2),(#H2,F)),(D)#H1)));" |> readTopology |> directEdges! ;

julia> taxa = net5 |> tipLabels # ABC EF D
6-element Array{String,1}:
 "A"
 "B"
 "C"
 "E"
 "F"
 "D"

julia> hardwiredCluster(net5.edge[12], taxa) # descendants of 12th edge = CEF
6-element Array{Bool,1}:
 0
 0
 1
 1
 1
 0
```
"""
function hardwiredCluster(edge::Edge,taxa::Union{Vector{String},Vector{Int}})
    v = zeros(Bool,length(taxa))
    hardwiredCluster!(v,edge,taxa)
    return v
end

hardwiredCluster!(v::Vector{Bool},edge::Edge,taxa::Union{Vector{String},Vector{Int}}) =
    hardwiredCluster!(v,edge,taxa,Int[])

function hardwiredCluster!(v::Vector{Bool},edge::Edge,taxa::Union{Vector{String},Vector{Int}},
                           visited::Vector{Int})
    n = getChild(edge)
    if n.leaf
        j = findall(isequal(n.name), taxa)
        length(j)==1 || error("taxon $(n.name) was not found in taxon list, or more than once")
        v[j[1]]=true
        return nothing
    end
    if n.number in visited
        return nothing  # n was already visited: exit. avoid infinite loop is isChild1 was bad.
    end
    push!(visited, n.number)
    for ce in n.edge
        if n == getParent(ce)
            hardwiredCluster!(v,ce,taxa,visited)
        end
    end
    return nothing
end


"""
    descendants(edge::Edge)

Return the node numbers of all the descendants of a given edge.

The node should belong in a rooted network for which isChild1 is up-to-date.
Run directEdges! beforehand. This is very important, otherwise one might enter an infinite loop,
and the function does not test for this.

## Examples
```jldoctest
julia> net5 = "(A,((B,#H1),(((C,(E)#H2),(#H2,F)),(D)#H1)));" |> readTopology |> directEdges! ;

julia> PhyloNetworks.descendants(net5.edge[12]) # descendants of 12th
7-element Array{Int64,1}:
 -6
 -7
  4
  6
  5
 -9
  7
```
"""
function descendants(edge::Edge)
    visited = Int[]
    descendants!(edge, visited)
    return visited
end

function descendants!(edge::Edge, visited::Vector{Int})
    n = getChild(edge)
    if n.number in visited
        return nothing  # n was already visited: exit. avoid infinite loop is isChild1 was bad.
    end
    push!(visited, n.number)
    for ce in n.edge
        if n == getParent(ce)
            descendants!(ce, visited)
        end
    end
    return nothing
end



"""
    deleteHybridThreshold!(net::HybridNetwork, threshold::Float64, keepNodes=false)

Deletes from a network all hybrid edges with heritability below a threshold gamma.
Returns the network.

- if threshold<0.5: delete minor hybrid edges with γ < threshold
  (or with a missing γ, for any threshold > -1.0)
- if threshold=0.5: delete all minor hybrid edges (i.e normally with γ < 0.5, if γ non-missing)
- `keepNodes`: if true, keep all original nodes; delete edges only.

Warnings:

- by default, `keepNodes` is false, and partner hybrid edges have their γ changed to 1.0.
  If `keepNodes` is true: the γ's of partner hybrid edges are unchanged.
- assumes correct isMajor attributes.
"""
function deleteHybridThreshold!(net::HybridNetwork, gamma::Float64, keepNodes=false::Bool)
    gamma <= 0.5 || error("deleteHybridThreshold! called with gamma = $(gamma)>0.5")
    for i = net.numHybrids:-1:1
    # starting from last because net.hybrid changes as hybrids are removed. Empty range if 0 hybrids.
        e = getMinorParentEdge(net.hybrid[i])
        # remove minor edge e if γ < threshold OR threshold=0.5
        # warning: no check if γ and isMajor are in conflict
        if e.gamma < gamma || gamma == 0.5 # note: γ=-1 if missing, so < gamma threshold
            # deleteHybrid!(net.hybrid[i],net,true,false) # requires non-missing edge lengths
            # deleteHybridizationUpdate! requires level-1 network with corresponding attributes
            deleteHybridEdge!(net, e, keepNodes) # does not update inCycle, containRoot, etc.
        end
    end
    return net
end

"""
    displayedNetworks!(net::HybridNetwork, node::Node, keepNode=false)

Extracts the two networks that simplify a given network at a given hybrid node:
deleting either one or the other parent hybrid edge.
If `keepNodes` is true, all original nodes are kept in both networks.

- the original network is modified: the minor edge removed.
- returns one HybridNetwork object: the network with the major edge removed
"""
function displayedNetworks!(net::HybridNetwork, node::Node, keepNodes=false::Bool)
    node.hybrid || error("will not extract networks from tree node $(node.number)")
    ind = findfirst(x -> x===node, net.node)
    ind !== nothing || error("node $(node.number) was not found in net")
    netmin = deepcopy(net)
    emin = getMinorParentEdge(node)
    deleteHybridEdge!(net   , emin, keepNodes)  # *no* update of inCycle, containRoot, etc.
    emaj = getMajorParentEdge(netmin.node[ind]) # hybrid node & edge in netmin
    deleteHybridEdge!(netmin, emaj, keepNodes)
    return netmin
end

"""
    displayedTrees(net::HybridNetwork, gamma::Float64; keepNodes=false::Bool)

Extracts all trees displayed in a network, following hybrid edges
with heritability >= γ threshold (or >0.5 if threshold=0.5)
and ignoring any hybrid edge with heritability lower than γ.
Returns an array of trees, as HybridNetwork objects.

`keepNodes`: if true, keep all nodes during hybrid edge removal.

Warnings:

- if `keepNodes` is true: the retained partner hybrid edges have
  their γ values unchanged, but their `isMajor` is changed to true
- assume correct `isMajor` attributes.
"""
function displayedTrees(net0::HybridNetwork, gamma::Float64; keepNodes=false::Bool)
    trees = HybridNetwork[]
    net = deepcopy(net0)
    deleteHybridThreshold!(net,gamma,keepNodes)
    displayedTrees!(trees,net,keepNodes)
    return trees # should have length 2^net.numHybrids
end

"""
    inheritanceWeight(tree::HybridNetwork)

Return the *log* inheritance weight of a network or tree
(as provided by [`displayedTrees`](@ref) with `keepNodes` = true for instance).
For a tree displayed in a network, its inheritance weight is the log of the product
of γ's of all edges retained in the tree. To avoid underflow, the log is calculated:
i.e. sum of log(γ) across retained edges.

If any edge has a negative γ, it is assumed to mean that its γ is missing,
and the function returns `missing`.

# Example

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);");

julia> trees = displayedTrees(net,0.0; keepNodes=true);

julia> PhyloNetworks.inheritanceWeight.(trees)
2-element Array{Float64,1}:
 -0.105361
 -2.30259 
```
"""
function inheritanceWeight(tree::HybridNetwork)
    ltw = 0.0
    for e in tree.edge
        e.gamma != 1.0 || continue
        if e.gamma < 0.0 return missing; end # any negative γ assumed -1.0: missing
        ltw += log(e.gamma)
    end
    return ltw
end

"""
`majorTree(net::HybridNetwork)`

Warning: assumes correct isMajor attributes.

Extracts the major tree displayed in a network, keeping the major edge and dropping the minor edge at each hybrid node.
Returns a HybridNetwork object.
"""
majorTree(net::HybridNetwork) = displayedTrees(net,0.5)[1]


# expands current list of trees, with trees displayed in a given network
function displayedTrees!(trees::Array{HybridNetwork,1},net::HybridNetwork, keepNodes=false::Bool)
    if isTree(net)
        # warning: no update of edges' containRoot (true) or edges' and nodes' inCycle (-1)
        push!(trees, net)
    else
        netmin = displayedNetworks!(net, net.hybrid[1], keepNodes)
        displayedTrees!(trees, net, keepNodes)
        displayedTrees!(trees, netmin, keepNodes)
    end
end

"""
    minorTreeAt(net::HybridNetwork, hybindex::Integer, keepNodes=false)

Extract the tree displayed in the network, following the major hybrid edge
at each hybrid node, except at the ith hybrid node (i=`hybindex`),
where the minor hybrid edge is kept instead of the major hybrid edge.
If `keepNodes` is true, all nodes are kept during edge removal.

Warning: assume correct `isMajor` fields.
"""
function minorTreeAt(net::HybridNetwork, hybindex::Integer, keepNodes=false::Bool)
    hybindex <= length(net.hybrid) || error("network has fewer hybrid nodes than index $(hybindex).")
    tree = deepcopy(net)
    hybedges = hybridEdges(tree.hybrid[hybindex])
    deleteHybridEdge!(tree, hybedges[2], keepNodes)
    # majorgamma = hybedges[1].gamma
    # setGamma!(hybedges[2],majorgamma) # set major gamma to minor edge (to delete old major = new minor)
    # deleteHybrid!(tree.hybrid[hybindex],tree,true,false) # major edge at hybrid removed.
    return majorTree(tree) # all remaining minor edges removed: now it's a tree.
end

"""
    displayedNetworkAt!(net::HybridNetwork, node::Node, keepNodes=false)

Delete all the minor hybrid edges, except at input node. The network is left
with a single hybridization, and otherwise displays the same major tree as before.
If `keepNodes` is true, all nodes are kept during edge removal.

Warning: assume correct `isMajor` fields.
"""
function displayedNetworkAt!(net::HybridNetwork, node::Node, keepNodes=false::Bool)
    node.hybrid || error("will not extract network from tree node $(node.number)")
    for i = net.numHybrids:-1:1
    # starting from last because net.hybrid changes as hybrids are removed. Empty range if 0 hybrids.
        net.hybrid[i] != node || continue
        emin = getMinorParentEdge(net.hybrid[i])
        deleteHybridEdge!(net, emin, keepNodes)
    end
end


"""
`hardwiredClusterDistance(net1::HybridNetwork, net2::HybridNetwork, rooted::Bool)`

Takes 2 networks and returns their hardwired cluster distance, that is,
the number of hardwired clusters found in one network and not in the other.
Note that this is not a distance per se on the full space of hybrid networks:
there are pairs of different networks for which this measure is 0.
But it is a distance on some network subspaces.

If the 2 networks are trees, this is the Robinson-Foulds distance.
If rooted=false, the trees are considered unrooted.
"""
function hardwiredClusterDistance(net1::HybridNetwork, net2::HybridNetwork, rooted::Bool)
    bothtrees = (net1.numHybrids == 0 && net2.numHybrids == 0)
    rooted || bothtrees || error("unrooted hardwired cluster distance not defined for non-tree networks.")
    taxa = sort!(String[net1.leaf[i].name for i in 1:net1.numTaxa])
    length(setdiff(taxa, String[net2.leaf[i].name for i in 1:net2.numTaxa])) == 0 ||
        error("net1 and net2 do not share the same taxon set. Please prune networks first.")
    nTax = length(taxa)
    if bothtrees # even if rooted, different treatment at the root if root=leaf
        M1 = tree2Matrix(net1, taxa, rooted=rooted)
        M2 = tree2Matrix(net2, taxa, rooted=rooted)
    else
        M1 = hardwiredClusters(net1, taxa) # last row: 10/11 if tree/hybrid edge.
        M2 = hardwiredClusters(net2, taxa)
        #println("M1="); print(M1); println("\nM2="); print(M2); println("\n");
    end
    dis = 0

    for i1=1:size(M1)[1]
        found = false
        m1 = 1 .- M1[i1,2:end] # going to the end: i.e. we want to match a tree edge with a tree edge
        for i2=1:size(M2)[1]                                  # and hybrid edge with hybrid edge
            if (M1[i1,2:end] == M2[i2,2:end] ||
                  ( !rooted && m1 == M2[i2,2:end])     )
                found = true
                break
            end
        end
        if !found
            dis += 1
        end
    end # (size(M1)[1] - dis) edges have been found in net2, dis edges have not.
    # so size(M2)[1] - (size(M1)[1] - dis) edges in net2 are not in net1.
    dis + dis + size(M2)[1] - size(M1)[1]
end
