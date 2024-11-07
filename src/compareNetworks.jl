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
function tree2Matrix(
    T::HybridNetwork,
    S::Union{Vector{String},Vector{Int}};
    rooted::Bool=true
)
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
        if !isleaf(child)
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
    hardwiredClusters(net::HybridNetwork, taxon_labels)


Returns a matrix describing all the hardwired clusters in a network, with
taxa listed in same order as in `taxon_labels` to describe their membership
in each cluster. Allows for missing taxa, with entries all 0.

Warnings:
- clusters are rooted, so the root must be correct.
- each hybrid node is assumed to have exactly 2 parents (no more).

Each row corresponds to one internal edge, that is, external edges are excluded.
If the root is a leaf node, the external edge to that leaf is included (first row).
Both parent hybrid edges to a given hybrid node only contribute a single row (they share the same hardwired cluster).

- first column: edge number
- next columns: 0/1. 1=descendant of edge, 0=not a descendant, or missing taxon.
- last column:  10/11 values. 10=tree edge, 11=hybrid edge

See also [`hardwiredClusterDistance`](@ref) and [`hardwiredCluster`](@ref).
"""
function hardwiredClusters(net::HybridNetwork, S::Union{AbstractVector{String},AbstractVector{Int}})
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

function hardwiredClusters!(node::Node, edge::Edge, ie::AbstractVector{Int}, M::Matrix{Int},
                            S::Union{AbstractVector{String},AbstractVector{Int}})
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
        partner !== nothing || error("partner hybrid edge not found for edge $(edge.number), child $(child.number)")
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
    hardwiredCluster(edge::Edge, taxa::Union{AbstractVector{String},AbstractVector{Int}})
    hardwiredCluster!(v::Vector{Bool}, edge::Edge, taxa)
    hardwiredCluster!(v::Vector{Bool}, edge::Edge, taxa, visited::Vector{Int})

Calculate the hardwired cluster of `edge`, coded as a vector of booleans:
true for taxa that are descendent of the edge, false for other taxa (including missing taxa).

The edge should belong in a rooted network for which `isChild1` is up-to-date.
Run `directEdges!` beforehand. This is very important, otherwise one might enter an infinite loop,
and the function does not test for this.

visited: vector of node numbers, of all visited nodes.

# Examples:
```jldoctest
julia> net5 = "(A,((B,#H1),(((C,(E)#H2),(#H2,F)),(D)#H1)));" |> readnewick |> directEdges! ;

julia> taxa = net5 |> tipLabels # ABC EF D
6-element Vector{String}:
 "A"
 "B"
 "C"
 "E"
 "F"
 "D"

julia> hardwiredCluster(net5.edge[12], taxa) # descendants of 12th edge = CEF
6-element Vector{Bool}:
 0
 0
 1
 1
 1
 0
```

See also [`hardwiredClusterDistance`](@ref) and [`hardwiredClusters`](@ref)
"""
function hardwiredCluster(edge::Edge,taxa::Union{AbstractVector{String},AbstractVector{Int}})
    v = zeros(Bool,length(taxa))
    hardwiredCluster!(v,edge,taxa)
    return v
end

hardwiredCluster!(v::Vector{Bool},edge::Edge,taxa::Union{AbstractVector{String},AbstractVector{Int}}) =
    hardwiredCluster!(v,edge,taxa,Int[])

function hardwiredCluster!(v::Vector{Bool},edge::Edge,taxa::Union{AbstractVector{String},AbstractVector{Int}},
                           visited::Vector{Int})
    n = getchild(edge)
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
        if n === getparent(ce)
            hardwiredCluster!(v,ce,taxa,visited)
        end
    end
    return nothing
end


"""
    descendants(edge::Edge, internal::Bool=false)

Return the node numbers of the descendants of a given edge:
all descendant nodes if `internal` is true (internal nodes and tips),
or descendant tips only otherwise (defaults).

`edge` should belong in a rooted network for which `isChild1` is up-to-date.
Run `directEdges!` beforehand. This is very important, otherwise one might enter an infinite loop,
and the function does not test for this.

## Examples
```jldoctest
julia> net5 = "(A,((B,#H1),(((C,(E)#H2),(#H2,F)),(D)#H1)));" |> readnewick |> directEdges! ;

julia> PhyloNetworks.descendants(net5.edge[12], true) # descendants of 12th edge: all of them
7-element Vector{Int64}:
 -6
 -7
  4
  6
  5
 -9
  7

julia> PhyloNetworks.descendants(net5.edge[12]) # descendant leaves only
3-element Vector{Int64}:
 4
 5
 7
```
"""
function descendants(edge::Edge, internal::Bool=false)
    visited = Int[]
    des = Int[]
    descendants!(des, visited, edge, internal)
    return des
end

function descendants!(des::Vector{Int}, visited::Vector{Int}, edge::Edge, internal::Bool=false)
    n = getchild(edge)
    if n.hybrid # only need to check previous visits for hybrid nodes
        n.number in visited && return nothing
        push!(visited, n.number)
    end
    if internal || n.leaf
        push!(des, n.number)
    end
    for ce in n.edge
        if isparentof(n, ce)
            descendants!(des, visited, ce, internal)
        end
    end
    return nothing
end

"""
    isdescendant(des:Node, anc::Node)

Return true if `des` is a strict descendant of `anc`, using `isChild1` fields
to determine the direction of edges. See [`isdescendant_undirected`](@ref)
for a version that does not use `isChild1`.
"""
function isdescendant(des::Node, anc::Node)
    visited = Int[]
    for e in anc.edge
        anc !== getchild(e) || continue # skip parents of anc
        if isdescendant!(visited, des, e)
            return true
        end
    end
    return false
end
function isdescendant!(visited::Vector{Int}, des::Node, e::Edge)
    n = getchild(e)
    if n == des
        return true
    end
    if n.hybrid # only need to check previous visits for hybrid nodes
        if n.number in visited # n & its descendants were already visited: exit
            return false
        end
        push!(visited, n.number)
    end
    for ce in n.edge
        if isparentof(n, ce)
            if isdescendant!(visited, des, ce) return true; end
        end
    end
    return false
end

"""
    isdescendant_undirected(des:Node, ancestor::Node, parentedge)

Return `true` if `des` is a strict descendant of `ancestor` when starting
from edge `parentedge` and going towards `ancestor` onward, regardless
of the field `isChild1` of tree edges; `false` otherwise.

This is useful to know how descendant relationships would change as a result
of reverting the direction of a tree edge, without actually
modifying the direction (`isChild1`) of any edge.

`parentedge` should be connected to `ancestor` (not checked).
The direction of hybrid edges is respected (via `isChild1`), that is,
the traversal does not go from the child to the parent of a hybrid edge.
"""
function isdescendant_undirected(des::Node, anc::Node, parentedge::Edge)
    visited = Int[]
    isdescendant_undirected!(visited, des, anc, parentedge)
end
function isdescendant_undirected!(visited::Vector{Int}, des::Node, anc::Node, parentedge::Edge)
    for e in anc.edge
        e !== parentedge || continue # do not go back up where we came from
        !e.hybrid || getparent(e) === anc || continue # do not go back up a parent hybrid edge of anc
        n = getOtherNode(e, anc)
        if n === des
            return true
        end
        if n.hybrid
            !(n.number in visited) || continue # skip to next edge is n already visited
            push!(visited, n.number)
        end
        if isdescendant_undirected!(visited, des, n, e)
            return true
        end
    end
    return false
end


"""
    ladderpartition(tree::HybridNetwork)

For each node in `tree`, calculate the clade below each child edge of the node,
and each clade moving up the "ladder" from the node to the root. The output is a
tuple of 2 vectors (node) of vector (clade) of vectors (taxon in clade):
`below,above`. More specifically, for node number `n`,
`below[n]` is generally of vector of 2 clades: one for the left child and one
for the right child of the node (unless the node is of degree 2 or is a polytomy).
`above[n]` contains the grade of clades above node number `n`.

WARNING: assumes that
1. node numbers and edge numbers can be used as indices, that is,
   be all distinct, positive, covering exactly 1:#nodes and 1:#edges.
2. edges are corrected directed (`isChild1` is up-to-date) and
   nodes have been pre-ordered already (field `nodes_changed` up-to-date).

# examples
```jldoctest
julia> tree = readnewick("(O,A,((B1,B2),(E,(C,D))));");

julia> PhyloNetworks.resetNodeNumbers!(tree; checkPreorder=true, type=:postorder)

julia> printNodes(tree)
node leaf  hybrid hasHybEdge name inCycle edges'numbers
1    true  false  false      O    -1      1   
2    true  false  false      A    -1      2   
3    true  false  false      B1   -1      3   
4    true  false  false      B2   -1      4   
8    false false  false           -1      3    4    5   
5    true  false  false      E    -1      6   
6    true  false  false      C    -1      7   
7    true  false  false      D    -1      8   
9    false false  false           -1      7    8    9   
10   false false  false           -1      6    9    10  
11   false false  false           -1      5    10   11  
12   false false  false           -1      1    2    11  

julia> below, above = PhyloNetworks.ladderpartition(tree);

julia> below
12-element Vector{Vector{Vector{Int64}}}:
 [[1]]                      
 [[2]]                      
 [[3]]                      
 [[4]]                      
 [[5]]                      
 [[6]]                      
 [[7]]                      
 [[3], [4]]                 
 [[6], [7]]                 
 [[5], [6, 7]]              
 [[3, 4], [5, 6, 7]]        
 [[1], [2], [3, 4, 5, 6, 7]]

julia> for n in 8:12
         println("clades below node ", n, ": ", join(below[n], " "))
       end
clades below node 8: [3] [4]
clades below node 9: [6] [7]
clades below node 10: [5] [6, 7]
clades below node 11: [3, 4] [5, 6, 7]
clades below node 12: [1] [2] [3, 4, 5, 6, 7]

julia> above[8:12] # clades sister to and above nodes 8 through 12:
5-element Vector{Vector{Vector{Int64}}}:
 [[5, 6, 7], [1], [2]]
 [[5], [3, 4], [1], [2]]
 [[3, 4], [1], [2]]     
 [[1], [2]]             
 []                     
```
"""
function ladderpartition(net::HybridNetwork)
    nnodes = length(net.node)
    nleaf  = length(net.leaf)
    for n in net.node
        n.number > 0 && n.number <= nnodes || error("node numbers must be in 1 - #nodes")
    end
    sort!([n.number for n in net.leaf]) == collect(1:nleaf) || error("leaves must come first")
    below = Vector{Vector{Vector{Int}}}(undef, nnodes)
    above = Vector{Vector{Vector{Int}}}(undef, nnodes)
    for nni in nnodes:-1:1 # post-order (not pre-order) traversal
        nn = net.nodes_changed[nni]
        above[nn.number] = Vector{Vector{Int}}(undef,0) # initialize
        if nn.leaf
            below[nn.number] = [[nn.number]]
            continue
        end
        below[nn.number]  = Vector{Vector{Int}}(undef,0) # initialize
        !nn.hybrid || error("ladder partitions not implemented for non-tree networks")
        children = [getchild(e) for e in nn.edge]
        filter!(n -> n!=nn, children)
        for cc in children
            allbelowc = union(below[cc.number]...)
            push!(below[nn.number], allbelowc)
            for cc2 in children
                cc2 !=cc || continue
                push!(above[cc2.number], allbelowc)
            end
        end
    end
    # so far: above[n] contains the clades *sister* to node number n, only.
    #         add those above = any above n's parent, using pre-order this time.
    for nni in 2:nnodes # avoid 1: it's the root, no parent, nothing to update
        nn = net.nodes_changed[nni]
        pn = getparent(nn).number # major parent number
        for clade in above[pn]  push!(above[nn.number], clade); end
    end
    return below,above
end

"""
    deleteHybridThreshold!(net::HybridNetwork, threshold::Float64,
                           nofuse=false, unroot=false, multgammas=false,
                           keeporiginalroot=false)

Deletes from a network all hybrid edges with heritability below a threshold gamma.
Returns the network.

- if threshold<0.5: delete minor hybrid edges with γ < threshold
  (or with a missing γ, for any threshold > -1.0)
- if threshold=0.5: delete all minor hybrid edges (i.e normally with γ < 0.5, if γ non-missing)
- `nofuse`: if true, do not fuse edges and keep original nodes.
- `unroot`: if false, the root will not be deleted if it becomes of degree 2.
- `multgammas`: if true, the modified edges have γ values equal to the
  proportion of genes that the extracted subnetwork represents. For an edge `e`
  in the modified network, the inheritance γ for `e` is the product of γs
  of all edges in the original network that have been merged into `e`.
-`keeporiginalroot`: if true, the root will be retained even if of degree 1.

Warnings:

- by default, `nofuse` is false, partner hybrid edges are fused with their child edge
  and have their γ changed to 1.0.
  If `nofuse` is true: the γ's of partner hybrid edges are unchanged.
- assumes correct `isMajor` fields, and correct `isChild1` fields to update `containRoot`.
"""
function deleteHybridThreshold!(
    net::HybridNetwork,
    gamma::Float64,
    nofuse::Bool=false,
    unroot::Bool=false,
    multgammas::Bool=false,
    keeporiginalroot::Bool=false
)
    gamma <= 0.5 || error("deleteHybridThreshold! called with gamma = $(gamma)>0.5")
    for i = net.numHybrids:-1:1
    # starting from last because net.hybrid changes as hybrids are removed. Empty range if 0 hybrids.
        i > lastindex(net.hybrid) && continue # removing 1 hybrid could remove several, if non-tree child net
        e = getparentedgeminor(net.hybrid[i])
        # remove minor edge e if γ < threshold OR threshold=0.5
        # warning: no check if γ and isMajor are in conflict
        if e.gamma < gamma || gamma == 0.5 # note: γ=-1 if missing, so < gamma threshold
            # deleteHybrid!(net.hybrid[i],net,true,false) # requires non-missing edge lengths
            # deleteHybridizationUpdate! requires level-1 network with corresponding attributes
            deletehybridedge!(net, e, nofuse, unroot, multgammas, true, keeporiginalroot)
        end
    end
    return net
end

"""
    displayedNetworks!(net::HybridNetwork, node::Node, keepNode=false,
                       unroot=false, multgammas=false, keeporiginalroot=false)

Extracts the two networks that simplify a given network at a given hybrid node:
deleting either one or the other parent hybrid edge.
If `nofuse` is true, the original edges (and nodes) are kept in both networks,
provided that they have one or more descendant leaves.
If `unroot` is true, the root will be deleted if it becomes of degree 2.
If `keeporiginalroot` is true, the root is retained even if it is of degree 1.

- the original network is modified: the minor edge removed.
- returns one HybridNetwork object: the network with the major edge removed
"""
function displayedNetworks!(
    net::HybridNetwork,
    node::Node,
    nofuse::Bool=false,
    unroot::Bool=false,
    multgammas::Bool=false,
    keeporiginalroot::Bool=false
)
    node.hybrid || error("will not extract networks from tree node $(node.number)")
    ind = findfirst(x -> x===node, net.node)
    ind !== nothing || error("node $(node.number) was not found in net")
    netmin = deepcopy(net)
    emin = getparentedgeminor(node)
    deletehybridedge!(net   , emin, nofuse, unroot, multgammas, true, keeporiginalroot)
    emaj = getparentedge(netmin.node[ind]) # hybrid node & edge in netmin
    deletehybridedge!(netmin, emaj, nofuse, unroot, multgammas, true, keeporiginalroot)
    return netmin
end

"""
    displayedTrees(net::HybridNetwork, gamma::Float64; nofuse::Bool=false,
                   unroot::Bool=false, multgammas::Bool=false,
                   keeporiginalroot::Bool=false)

Extracts all trees displayed in a network, following hybrid edges
with heritability >= γ threshold (or >0.5 if threshold=0.5)
and ignoring any hybrid edge with heritability lower than γ.
Returns an array of trees, as HybridNetwork objects.

`nofuse`: if true, do not fuse edges (keep degree-2 nodes) during hybrid edge removal.  
`unroot`: if false, the root will not be deleted if it becomes of degree 2 unless
  keeporiginalroot is true.  
`multgammas`: if true, the edges in the displayed trees have γ values
  equal to the proportion of genes that the edge represents, even though all
  these edges are tree edges. The product of all the γ values across all edges
  is the proportion of genes that the tree represents. More specifically,
  edge `e` in a given displayed tree has γ equal to the product of γs
  of all edges in the original network that have been merged into `e`.  
`keeporiginalroot`: if true, keep root even if of degree 1.

Warnings:

- if `nofuse` is true: the retained partner hybrid edges have
  their γ values unchanged, but their `isMajor` is changed to true
- assume correct `isMajor` attributes.
"""
function displayedTrees(
    net0::HybridNetwork,
    gamma::Float64;
    nofuse::Bool=false,
    unroot::Bool=false,
    multgammas::Bool=false,
    keeporiginalroot::Bool=false
)
    trees = HybridNetwork[]
    net = deepcopy(net0)
    deleteHybridThreshold!(net,gamma,nofuse,unroot, multgammas, keeporiginalroot)
    displayedTrees!(trees,net,nofuse,unroot, multgammas, keeporiginalroot)
    return trees # should have length 2^net.numHybrids
end

"""
    inheritanceWeight(tree::HybridNetwork)

Return the *log* inheritance weight of a network or tree
(as provided by [`displayedTrees`](@ref) with `nofuse` = true for instance).
For a tree displayed in a network, its inheritance weight is the log of the product
of γ's of all edges retained in the tree. To avoid underflow, the log is calculated:
i.e. sum of log(γ) across retained edges.

If any edge has a negative γ, it is assumed to mean that its γ is missing,
and the function returns `missing`.

# Example

```julia-repl
julia> net = readnewick("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);");

julia> trees = displayedTrees(net,0.0; nofuse=true);

julia> PhyloNetworks.inheritanceWeight.(trees)
2-element Vector{Float64}:
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
    majorTree(net::HybridNetwork; nofuse::Bool=false, unroot::Bool=false,
              keeporiginalroot::Bool=false)

Extract the major tree displayed in a network, keeping the major edge
and dropping the minor edge at each hybrid node.

`nofuse`: if true, edges and degree-2 nodes are retained during edge removal.
Otherwise, at each reticulation the child edge (below the hybrid node)
is retained: the major hybrid edge is fused with it.

`unroot`: is true, the root will be deleted if it becomes of degree 2.

`keeporiginalroot`: the network's root is kept even if it becomes of degree 1.

Warnings:

- if `nofuse` is true: the hybrid edges that are retained (without fusing)
  have their γ values unchanged, but their `isMajor` is changed to true
- assume correct `isMajor` attributes.
"""
majorTree(
    net::HybridNetwork;
    nofuse::Bool=false,
    unroot::Bool=false,
    keeporiginalroot::Bool=false
) = displayedTrees(net,0.5; nofuse=nofuse, unroot=unroot,
                   keeporiginalroot=keeporiginalroot)[1]


# expands current list of trees, with trees displayed in a given network
function displayedTrees!(
    trees::Array{HybridNetwork,1},
    net::HybridNetwork,
    nofuse::Bool=false,
    unroot::Bool=false,
    multgammas::Bool=false,
    keeporiginalroot::Bool=false
)
    if isTree(net)
        # warning: no update of edges' containRoot (true) or edges' and nodes' inCycle (-1)
        push!(trees, net)
    else
        netmin = displayedNetworks!(net, net.hybrid[1], nofuse, unroot, multgammas, keeporiginalroot)
        displayedTrees!(trees, net, nofuse, unroot, multgammas, keeporiginalroot)
        displayedTrees!(trees, netmin, nofuse, unroot, multgammas, keeporiginalroot)
    end
end

"""
    minorTreeAt(net::HybridNetwork, hybindex::Integer, nofuse=false, unroot::Bool=false)

Extract the tree displayed in the network, following the major hybrid edge
at each hybrid node, except at the ith hybrid node (i=`hybindex`),
where the minor hybrid edge is kept instead of the major hybrid edge.
If `nofuse` is true, edges are not fused (degree-2 nodes are kept).
If `unroot` is true, the root will be deleted if it becomes of degree 2.

Warning: assume correct `isMajor` fields.
"""
function minorTreeAt(
    net::HybridNetwork,
    hybindex::Integer,
    nofuse::Bool=false,
    unroot::Bool=false
)
    hybindex <= length(net.hybrid) || error("network has fewer hybrid nodes than index $(hybindex).")
    tree = deepcopy(net)
    hybedge = getparentedge(tree.hybrid[hybindex])   # major parent
    deletehybridedge!(tree, hybedge, nofuse, unroot) # delete major hybrid edge at reticulation of interest
    return majorTree(tree; nofuse=nofuse, unroot=unroot) # all remaining minor edges removed: now it's a tree.
end

"""
    displayedNetworkAt!(net::HybridNetwork, node::Node, nofuse=false,
                        unroot=false, multgammas=false)

Delete all the minor hybrid edges, except at input node. The network is left
with a single hybridization, and otherwise displays the same major tree as before.
If `nofuse` is true, edges are not fused (degree-2 nodes are kept).

Warning: assume correct `isMajor` fields.
"""
function displayedNetworkAt!(
    net::HybridNetwork,
    node::Node,
    nofuse::Bool=false,
    unroot::Bool=false,
    multgammas::Bool=false
)
    node.hybrid || error("will not extract network from tree node $(node.number)")
    for i = net.numHybrids:-1:1
    # starting from last because net.hybrid changes as hybrids are removed. Empty range if 0 hybrids.
        net.hybrid[i] != node || continue
        emin = getparentedgeminor(net.hybrid[i])
        deletehybridedge!(net, emin, nofuse, unroot, multgammas)
    end
end


"""
    hardwiredClusterDistance(net1::HybridNetwork, net2::HybridNetwork, rooted::Bool)

Hardwired cluster distance between the topologies of `net1` and `net2`, that is,
the number of hardwired clusters found in one network and not in the other
(with multiplicity, see below).

If the 2 networks are trees, this is the Robinson-Foulds distance.
If rooted=false, then both networks are considered as semi-directed.

Networks are assumed bicombining (each hybrid has exactly 2 parents, no more).

## Dissimilarity vs distance

This is *not* a distance per se on the full space of phylogenetic networks:
there are pairs of distinct networks for which this dissimilarity is 0.
But it is a distance on some classes of networks, such as the class of
tree-child networks that are "normal" (without shortcuts), or the class of
tree-child networks that can be assigned node ages such that hybrid edges
have length 0 and tree edges have non-negative lengths. See
[Cardona, Rossello & Valiente (2008)](https://doi.org/10.1016/j.mbs.2007.11.003),
[Cardona, Llabres, Rossello & Valiente (2008)](https://doi.org/10.1109/TCBB.2008.70),
and [Huson, Rupp, Scornavacca (2010)](https://doi.org/10.1017/CBO9780511974076).

## Example

```jldoctest
julia> net1 = readnewick("(t6,(t5,((t4,(t3,((t2,t1))#H1)),#H1)));");

julia> taxa = sort(tipLabels(net1)); # t1 through t6, sorted alphabetically

julia> # using PhyloPlots; plot(net1, showedgenumber=true);

julia> # in matrix below: column 1: edge number. last column: tree (10) vs hybrid (11) edge
       # middle columns: for 'taxa': t1,...t6. 1=descendant, 0=not descendant
       hardwiredClusters(net1, taxa)
6×8 Matrix{Int64}:
 13  1  1  1  1  1  0  10
 12  1  1  1  1  0  0  10
 10  1  1  1  1  0  0  10
  9  1  1  1  0  0  0  10
  8  1  1  0  0  0  0  11
  7  1  1  0  0  0  0  10

julia> net2 = readnewick("(t6,(t5,((t4,(t3)#H1),(#H1,(t1,t2)))));");

julia> hardwiredClusters(net2, taxa)
6×8 Matrix{Int64}:
 13  1  1  1  1  1  0  10
 12  1  1  1  1  0  0  10
  6  0  0  1  1  0  0  10
  5  0  0  1  0  0  0  11
 11  1  1  1  0  0  0  10
 10  1  1  0  0  0  0  10

julia> hardwiredClusterDistance(net1, net2, true) # true: as rooted networks
4
```

## What is a hardwired cluster?

Each edge in a network is associated with its *hardwired cluster*, that is,
the set of all its descendant taxa (leaves). The set of hardwired cluster
of a network is the set of its edges' hardwired clusters. The dissimilarity
`d_hard` defined in [Huson, Rupp, Scornavacca (2010)](https://doi.org/10.1017/CBO9780511974076)
is the number of hardwired clusters that are in one network but not in the other.

This implementation is a slightly more discriminative version of `d_hard`, where
each cluster is counted with multiplicity and annotated with its edge's hybrid
status, as follows:
- External edges are not counted (they are tree edges to a leaf, shared by all
  phylogenetic networks).
- A cluster is counted for each edge for which it's the hardwired cluster.
- At a given hybrid node, both hybrid partner edges have the same cluster,
  so this cluster is only counted once for both partners.
- A given cluster is matched between the two networks only if it's the cluster
  from a tree edge in both networks, or from a hybrid edge in both networks.

In the example above, `net1` has a shortcut (hybrid edge 11) resulting in 2 tree
edges (12 and 10) with the same cluster {t1,t2,t3,t4}. So cluster {t1,t2,t3,t4}
has multiplicity 2 in `net1`. `net2` also has this cluster, but only associated
with 1 tree edge, so this cluster contributes (2-1)=1 towards the hardwired cluster
distance between the two networks. The distance of 4 corresponds to these 4 clusters:
- {t1,t2,t3,t4}: twice in net1, once in net2
- {t3,t4}: absent in net1, once in net2
- {t1,t2}: twice in net1 (from a hybrid edge & a tree edge), once in net2
- {t3}: absent in net1 (because external edges are not counted),
  once in net2 (from a hybrid edge).

Degree-2 nodes cause multiple edges to have the same cluster, so counting
clusters with multiplicity distinguishes a network with extra degree-2 nodes
from the "same" network after these nodes have been suppressed
(e.g. with [`PhyloNetworks.fuseedgesat!`](@ref) or [`PhyloNetworks.shrinkedge!`](@ref)).

## Networks as semi-directed

If `rooted` is false and one of the phylogenies is not a tree (1+ reticulations),
then all degree-2 nodes are removed before comparing the hardwired clusters,
and the minimum distance is returned over all possible ways to root the
networks at internal nodes.

See also: [`hardwiredClusters`](@ref), [`hardwiredCluster`](@ref)
"""
function hardwiredClusterDistance(net1::HybridNetwork, net2::HybridNetwork, rooted::Bool)
    bothtrees = (net1.numHybrids == 0 && net2.numHybrids == 0)
    rooted || bothtrees ||
        return hardwiredClusterDistance_unrooted(net1, net2) # tries all roots, but removes degree-2 nodes
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
    n2ci = collect(1:size(M2, 1)) # cluster indices
    for i1 in 1:size(M1,1)
        found = false
        m1 = 1 .- M1[i1,2:end] # going to the end: i.e. we want to match a tree edge with a tree edge
                                # and hybrid edge with hybrid edge
        for j in length(n2ci):-1:1 # check only unmatched cluster indices, in reverse
            i2 = n2ci[j]
            if (M1[i1,2:end] == M2[i2,2:end] ||
                  ( !rooted && m1 == M2[i2,2:end])     )
                found = true
                deleteat!(n2ci, j) # a cluster can be repeated
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


"""
    hardwiredClusterDistance_unrooted(net1::HybridNetwork, net2::HybridNetwork)

Miminum hardwired cluster dissimilarity between the two networks, considered as
unrooted (or semi-directed). This dissimilarity is defined as the minimum
rooted distance, over all root positions that are compatible with the direction
of hybrid edges.
Called by [`hardwiredClusterDistance`](@ref).

To avoid repeating identical clusters, all degree-2 nodes
are deleted before starting the comparison.
Since rooting the network at a leaf creates a root node of degree 2 and
an extra cluster, leaves are excluded from possible rooting positions.
"""
function hardwiredClusterDistance_unrooted(net1::HybridNetwork, net2::HybridNetwork)
    return hardwiredClusterDistance_unrooted!(deepcopy(net1), deepcopy(net2))
end
function hardwiredClusterDistance_unrooted!(net1::HybridNetwork, net2::HybridNetwork)
    #= fixit: inefficient function, because r1 * r2 "M" matrices of
      hardwiredClusters() are calculated, where ri = # root positions in neti.
      Rewrite to calculate only r1 + r2 M's.
    =#
    removedegree2nodes!(net1) # because re-rooting would remove them in an
    removedegree2nodes!(net2) # unpredictable order
    # find all permissible positions for the root
    net1roots = [n.number for n in net1.node if !n.leaf]
    #= disallow the root at a leaf: adding a degree-2 node adds a cluster
       that could be artificially matched to a cluster from a degree-3 node
       sister to a hybrid edge, when a the leaf edge is the donor. =#
    for i in length(net1roots):-1:1 # reverse order, to delete some of them
        try
            rootatnode!(net1, net1roots[i])
            # tricky: rootatnode adds a degree-2 node if i is a leaf,
            #         and delete former root node if it's of degree 2.
        catch e
            isa(e, RootMismatch) || rethrow(e)
            deleteat!(net1roots, i)
        end
    end
    net2roots = [n.number for n in net2.node if !n.leaf]
    for i in length(net2roots):-1:1
        try
            rootatnode!(net2, net2roots[i])
        catch e
            isa(e, RootMismatch) || rethrow(e)
            deleteat!(net2roots, i)
        end
    end
    bestdissimilarity = typemax(Int)
    bestns = missing
    for n1 in net1roots
        rootatnode!(net1, n1)
        for n2 in net2roots
            rootatnode!(net2, n2)
            diss = hardwiredClusterDistance(net1, net2, true) # rooted = true now
            if diss < bestdissimilarity
                bestns = (n1, n2)
                bestdissimilarity = diss
            end
        end
    end
    # @info "best root nodes: $bestns"
    # warning: original roots (and edge directions) NOT restored
    return bestdissimilarity
end
