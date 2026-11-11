"""
    MuTag

μ-tag to tag μ-vectors, one of `:r` for root components, `:t` for tree edges,
`:h` for hybrid edges, or `:i` for edges incident to (but not in)
a root component, for when we place the root along the edge and follow it
in its opposite direction.
"""
struct MuTag
    tag::Int8
end
# integer that happens to match in-degree for non-incident tag
const mutag_r = MuTag( 0) # root node. integer that matches in-degree
const mutag_t = MuTag( 1) # tree edges / nodes
const mutag_h = MuTag( 2) # hybrid edges / nodes
const mutag_i = MuTag( 3) # incident to root component. ordered last by `isless`
function Base.show(io::IO, obj::MuTag)
    tag = obj.tag
    s = (tag==0 ? ":r" : (tag==1 ? ":t" : (tag==2 ? ":h" : ":i")))
    print(io, s)
end
Base.isless(t1::MuTag, t2::MuTag) = isless(t1.tag, t2.tag)

"""
    MuVector

μ-vector for a single node
(see [Cardona et al. 2024](https://doi.org/10.1109/TCBB.2024.3361390))
separated into
- μ0: number of paths starting at the node or edge and ending at any hybrid node
- μ-vector counting the number of paths from that node to each specific leaf.

This type does not store the leaf labels nor the order in which leaves are
positioned in μ-vectors.
"""
struct MuVector
    "μ0: number of paths to any hybrid node"
    mu_hybs::Int
    "μ-vector: number of paths to tips"
    mu_tips::Vector{Int}
end
function Base.show(io::IO, obj::MuVector)
    s  = "μ0=$(obj.mu_hybs) μ=$(obj.mu_tips)"
    print(io, s)
end
function ==(m1::MuVector, m2::MuVector)
    return m1.mu_hybs == m2.mu_hybs && m1.mu_tips == m2.mu_tips
end
#= comparing μ-vectors: lexicographic order on (μ-vector to leaves..., μ0)
   such that ignoring μ0 does not affect the order (see option useμ0)
!!! info "version info"
    In PhyloNetworks v1.4.0, the order was reversed: μ0 was compared first,
    aka using the lexicographic order on (μ0, μ-vectors to leaves...)
=#
function Base.isless(m1::MuVector, m2::MuVector)
    sgn = cmp(m1.mu_tips, m2.mu_tips) # -1, 0 if equal, 1
    if sgn == 0
        return isless(m1.mu_hybs, m2.mu_hybs)
    end
    return sgn == -1
end
# isless_μtips: to ignore μ0, used for option `useμ0=false`
isless_μtips(m1::MuVector, m2::MuVector) = isless(m1.mu_tips, m2.mu_tips)
has_0μentries_at(m::MuVector, i) = all(m.mu_tips[i] .== 0)
has_0μentries_at(ms::Tuple{MuVector,MuVector}, i) =
    has_0μentries_at(ms[1],i) && has_0μentries_at(ms[2],i)

"""
    TagMuVector

Tagged μ-vector, which wraps a tag [`MuTag`](@ref)
and a "simple" μ-vector [`MuVector`](@ref).
See [Maxfield, Xu & Ané 2025](https://doi.org/10.1109/TCBBIO.2025.3534780) and
[Ribas, Pons & Ané (2026)](https://arxiv.org/abs/2605.06243).
"""
struct TagMuVector
    "tag"
    mu_tag::MuTag
    "μ0 to hybrids & μ-vector to tips"
    mu_vec::MuVector
end
function Base.show(io::IO, obj::TagMuVector)
    print(io, "tag=$(obj.mu_tag) ")
    show(io, obj.mu_vec)
end
function ==(m1::TagMuVector, m2::TagMuVector)
    return m1.mu_tag == m2.mu_tag && m1.mu_vec == m2.mu_vec
end
# comparing tagged μ-vectors: compare the tag first, then the μ-vector
function Base.isless(m1::TagMuVector, m2::TagMuVector)
    if m1.mu_tag == m2.mu_tag
        return isless(m1.mu_vec, m2.mu_vec)
    end
    return isless(m1.mu_tag, m2.mu_tag)
end
isless_μtips(m1::TagMuVector, m2::TagMuVector) = isless_μtips(m1.mu_vec, m2.mu_vec)

# To compare μ-entries with 1 and 2 tagged μ-vectors: use the 1st μ-vector
# to ignore the :i-tagged μ-vector (2nd from incident edges): see option
# `useitag=false`. Not meant to be used for edges in the root component.
Base.isless(m1::TagMuVector, m2::Tuple{TagMuVector,TagMuVector}) =
    isless(m1, m2[1])
Base.isless(m1::Tuple{TagMuVector,TagMuVector}, m2::TagMuVector) =
    isless(m1[1], m2)
isless_μtips(m1::TagMuVector, m2::Tuple{TagMuVector,TagMuVector}) =
    isless_μtips(m1, m2[1])
isless_μtips(m1::Tuple{TagMuVector,TagMuVector}, m2::TagMuVector) =
    isless_μtips(m1[1], m2)

has_0μentries_at(m::TagMuVector, i) = has_0μentries_at(m.mu_vec, i)
has_0μentries_at(ms::Tuple{TagMuVector,TagMuVector}, i) =
    has_0μentries_at(ms[1],i) && has_0μentries_at(ms[2],i)

abstract type MuRepresentation end

tiplabels(m::MuRepresentation) = m.tiplabels

"""
    NodeMuRepresentation

Representation of a rooted network with one [`MuVector`](@ref) per node.
This entry contains the node's μ-vector: a vector of integers counting the
number of paths from that node to each leaf.
`tiplabels` gives the order of leaves in μ-vectors.

[`node_murepresentation`](@ref) constructs such an object, with a μ-vector
for each node. Leaf nodes may be included, but are excluded by default because
leaves have trivial μ-vectors (all 0s but one 1) that match across any network
on the same leaf set. Including leaves is useful if networks don't have the
same taxa, or for building a edge-based μ-representation.

See [Cardona et al. (2024)](https://doi.org/10.1109/TCBB.2024.3361390).
"""
struct NodeMuRepresentation <: MuRepresentation
    "vector of tip (leaf) labels: listed in the same order as in μ-vectors"
    tiplabels::Vector{String}
    "Dictionary: maps node number => [`MuVector`](@ref) for that node"
    mu_map::Dict{Int, MuVector}
    "sorted vector of μ-vectors"
    mu_vec::Vector{MuVector}
end
function Base.show(io::IO, obj::NodeMuRepresentation)
    labs = tiplabels(obj)
    disp = "$(typeof(obj))\n$(length(labs)) taxa in μ-vectors: $labs\n"
    disp *= "$(length(obj.mu_map)) nodes, map node number => μ0 to hybrids and μ-vector to taxa:"
    for (k,v) in obj.mu_map disp *= "\n  $k => $v"; end
    print(io, disp)
end

# constructor: builds sorted vector from map, re-using memory for each μ-vector
function NodeMuRepresentation(labels, mu_map)
    mu_vec = collect(values(mu_map))
    sort!(mu_vec)
    return NodeMuRepresentation(labels, mu_map, mu_vec)
end

function has_0μentries_at(m::NodeMuRepresentation, i)
    for md in values(m.mu_map)
        has_0μentries_at(md, i) || return false
    end
    return true
end

"""
    indexin_check0μentries(m1::MuRepresentation, m2::MuRepresentation)

Check that for any label in `m1` that is missing in `m2`, the corresponding
number of paths in `m1`'s entries are all 0 -- and vice versa, then return
- `nothing` if the check failed,
- otherwise if the check passed: two vectors of indices `(o1,o2)` such that
  `tiplabels(m1)[o1] == tiplabels(m2)[o2]` is the intersection of the label sets.
"""
function indexin_check0μentries(m1::MuRepresentation, m2::MuRepresentation)
    l1 = tiplabels(m1); l2 = tiplabels(m2)
    o1 = collect(1:length(l1))
    o2 = indexin(l1, l2) # order to use for μ-vectors in m2: l2[o2] == l1
    # if a label is in l1 but not in m2: all its entries should be 0 in m1
    in1not2_i = findall(isnothing, o2) # indices in o2 and in m1
    if !isempty(in1not2_i)
        has_0μentries_at(m1, in1not2_i) || return nothing
        deleteat!(o2, in1not2_i)
        deleteat!(o1, in1not2_i)
    end
    if length(l2) + length(in1not2_i) > length(l1)
        # then some labels in are l2 but not in l1: find their indices in l2
        in2not1_i = findall(!in(l1), l2)
        has_0μentries_at(m2, in2not1_i) || return nothing
    end
    return (o1, o2)
end

"""
    ==(m1::NodeMuRepresentation, m2::NodeMuRepresentation)

Equality of two [`NodeMuRepresentation`](@ref) objects, requiring that:
- If `m1` has labels that `m2` doesn't have, then all corresponding values
  for these labels in `m1` should be 0 -- as if a label is tracked in `m1`
  but absent from the original network, and not tracked in `m2`.
  Vice versa: labels in `m2` but not in `m1` should have 0 values in `m2`.
- `m1` and `m2` should have the same number of μ-vectors
- The μ-vectors should have the same subvectors corresponding to the labels
  shared by `m1` and `m2`.
"""
function ==(m1::NodeMuRepresentation, m2::NodeMuRepresentation)
    # order for μ-vectors: tiplabels(m1)[o1] == tiplabels(m2)[o2], others are 0
    o12 = indexin_check0μentries(m1, m2)
    isnothing(o12) && return false
    o1, o2 = o12
    length(m1.mu_vec) == length(m2.mu_vec) || return false
    # below: assumes that labels are sorted in both m1 and m2, such that
    # lexicographic sorting of mu-vectors gives the same order in m1 and m2
    for (μ1, μ2) in zip(m1.mu_vec, m2.mu_vec)
        μ1.mu_hybs == μ2.mu_hybs || return false
        μ1.mu_tips[o1] == μ2.mu_tips[o2] || return false
    end
    return true
end


"""
    EdgeMuRepresentation

Representation of a semidirected network with one μ-entry for each edge,
as described in [Ribas, Pons & Ané (2026)](https://arxiv.org/abs/2605.06243).
See also [Maxfield, Xu & Ané 2025](https://doi.org/10.1109/TCBBIO.2025.3534780),
which does not have the extension to count the number of paths to hybrid nodes,
or to distinguish directed edges where the root might be placed.
[`edge_murepresentation`](@ref) constructs such an object.
Each μ-entry is a set of 1 or 2 tagged μ-vectors.

1. If an edge is in the root component, that is, both of its incident nodes can
   root the network, then its μ-entry is a set of 2 tree-tagged μ-vectors,
   one corresponding to each direction of the edge.
2. If only the edge's parent may be used to root the network but not its child,
   then this child is either a hybrid or a leaf. If, in addition, the edge's
   parent is *not* the only possible root (rooting along the edge would give a
   different rooted network), then the edge's μ-entry has 2 tagged μ-vectors:
   * one corresponding to the edge's direction, tagged by the type of its child
     node (hybrid or tree), and
   * one tagged 'incident', corresponding to the edge in its opposite direction,
     as if it was subdivided into 2 segments by placing the root along the edge.
3. Otherwise, we may not root the network along the edge: its child cannot root
   the root and its parent cannot either (or can but is the only possible root).
   Then the edge's μ-entry is 1 μ-vector tagged by
   the edge type (tree or hybrid) and corresponds to the edge's direction,
   which is independent of where the network is rooted.

All μ-vectors include the number of paths to hybrids, defined for rooted networks
by [Cardona et al. (2024)](https://doi.org/10.1109/TCBB.2024.3361390).
See [`MuVector`](@ref) and [`TagMuVector`](@ref).

!!! info "version info"
    In PhyloNetworks v1.4.0:
    - external edges to a leaf, which are tree edges with a trivial μ-vector,
      were not included in the μ-representation;
    - the 'tree' tag of edges in the root component was not stored;
    - the μ-entry for a hybrid edge adjacent to the root component had a single
      μ-vector, tagged as 'incident', which departs from Maxfield, Xu & Ané (2025)
      (which tag it as 'hybrid') and from Ribas, Pons & Ané (2026), which build
      it with 2 tagged μ-vectors (one with the 'hybrid' tag and one with the
      'incident' tag).

"""
struct EdgeMuRepresentation <: MuRepresentation
    "vector of tip (leaf) labels: listed in the same order as in μ-vectors"
    tiplabels::Vector{String}
    "root μ-vector, assuming the network has a single root component"
    mu_root::MuVector
    """dictionary for edges in a root component:
       edge number => ((:t, μ-vector1), (:t, μ-vector2))"""
    mumap_rootcomp::Dict{Int, Tuple{TagMuVector,TagMuVector}}
    """dictionary for directed edges incident to a root component:
       edge number => ((:t or :h, μ-vector1), (:i, μ-vector2))"""
    mumap_incident::Dict{Int, Tuple{TagMuVector,TagMuVector}}
    """dictionary for directed edges that may not contain the root:
       edge number => tagged μ-vector"""
    mumap_directed::Dict{Int, TagMuVector}
    "sorted vector of sorted μ-entries in a root component"
    muvec_rootcomp::Vector{Tuple{TagMuVector,TagMuVector}}
    "sorted vector of sorted μ-entries incident to a root component"
    muvec_incident::Vector{Tuple{TagMuVector,TagMuVector}}
    "sorted vector of μ-entries in the directed part"
    muvec_directed::Vector{TagMuVector}
end
function Base.show(io::IO, obj::EdgeMuRepresentation)
    labs = tiplabels(obj)
    s = "$(typeof(obj))\n$(length(labs)) taxa in μ-vectors: $labs\n"
    s *= "root μ-vector: $(obj.mu_root)\n"
    s *= "map edge number => μ-entry:"
    s *= "\n$(length(obj.mumap_rootcomp)) edges in the root component"
    for (k,v) in obj.mumap_rootcomp
        s *= "\n  $k => $(v[1]); $(v[2])"
    end
    s *= "\n$(length(obj.mumap_incident)) edges incident to the root component"
    for (k,v) in obj.mumap_incident
        s *= "\n  $k => $(v[1]); $(v[2])"
    end
    s *= "\n$(length(obj.mumap_directed)) directed edges that may not contain the root"
    for (k,v) in obj.mumap_directed
        s *= "\n  $k => $v"
    end
    print(io, s)
end

# constructor: builds sorted vectors from dictionaries
#              re-using memory for each μ-vector
# root component edges: assumes all have 2 'tree' tags
# incident edges: assumes that the μ-vector tagged 'incident' comes second
function EdgeMuRepresentation(labels, mu_root, map_root, map_inc, map_dir)
    vec_root = Vector{Tuple{TagMuVector,TagMuVector}}()
    for (v1,v2) in values(map_root)
        if v1.mu_vec < v2.mu_vec # assumes both have the same tag (:t)
            push!(vec_root, (v1,v2))
        else
            push!(vec_root, (v2,v1))
        end
    end
    sort!(vec_root)
    vec_inc = collect(values(map_inc)) # assumes (:?, :i) with :? < :i, so swap not needed
    sort!(vec_inc)
    vec_dir = collect(values(map_dir))
    sort!(vec_dir)
    return EdgeMuRepresentation(labels, mu_root, map_root, map_inc, map_dir,
                                vec_root, vec_inc, vec_dir)
end

function has_0μentries_at(m::EdgeMuRepresentation, i)
    has_0μentries_at(m.mu_root, i) || return false
    for mr in m.muvec_rootcomp
        has_0μentries_at(mr, i) || return false
    end
    for mr in m.muvec_incident
        has_0μentries_at(mr, i) || return false
    end
    for md in m.muvec_directed
        has_0μentries_at(md, i) || return false
    end
    return true
end

"""
    ==(m1::EdgeMuRepresentation, m2::EdgeMuRepresentation)

Equality of two [`EdgeMuRepresentation`](@ref) objects, requiring that:
- If `m1` has labels that `m2` doesn't have, then all corresponding entries
  for these labels in `m1` should be 0 -- as if a label is tracked in `m1`
  but absent from the original network, and not tracked in `m2`.
  Vice versa: labels in `m2` but not in `m1` should have 0 entries in `m2`.
- `m1` and `m2` should have the same number of μ-entries of the same type
  (either in, incident to the root component, or directed with the same tags).
- All tagged μ-entries should have the same μ-vectors corresponding to the labels
  shared by `m1` and `m2`.

Assumes that when present, the :i-tagged μ-vector is second.
"""
function ==(m1::EdgeMuRepresentation, m2::EdgeMuRepresentation)
    # order for μ-vectors: tiplabels(m1)[o1] == tiplabels(m2)[o2], others are 0
    o12 = indexin_check0μentries(m1, m2)
    isnothing(o12) && return false
    o1, o2 = o12
    m1.mu_root.mu_hybs == m2.mu_root.mu_hybs || return false
    m1.mu_root.mu_tips[o1] == m2.mu_root.mu_tips[o2] || return false
    length(m1.muvec_rootcomp) == length(m2.muvec_rootcomp) || return false
    length(m1.muvec_incident) == length(m2.muvec_incident) || return false
    length(m1.muvec_directed) == length(m2.muvec_directed) || return false
    # below: assumes that labels are sorted in both m1 and m2, such that
    # lexicographic sorting of mu-vectors gives the same order in m1 and m2
    for (ms1, ms2) in zip(m1.muvec_rootcomp, m2.muvec_rootcomp)
        ms1[1].mu_vec.mu_hybs == ms2[1].mu_vec.mu_hybs || return false
        ms1[1].mu_vec.mu_tips[o1] == ms2[1].mu_vec.mu_tips[o2] || return false
        # no need to check second vector: bc 1 root and same root μ-vector
    end
    for (ms1, ms2) in zip(m1.muvec_incident, m2.muvec_incident)
        ms1[1].mu_tag == ms2[1].mu_tag || return false
        ms1[1].mu_vec.mu_hybs == ms2[1].mu_vec.mu_hybs || return false
        ms1[1].mu_vec.mu_tips[o1] == ms2[1].mu_vec.mu_tips[o2] || return false
        # no need to check second vector (:i tag), because 1 root μ-vector
    end
    for (ms1, ms2) in zip(m1.muvec_directed, m2.muvec_directed)
        ms1.mu_tag == ms2.mu_tag || return false # same tags
        ms1.mu_vec.mu_hybs == ms2.mu_vec.mu_hybs || return false # same μ-vectors
        ms1.mu_vec.mu_tips[o1] == ms2.mu_vec.mu_tips[o2] || return false
    end
    return true
end


"""
    node_murepresentation(net::HybridNetwork,
        labels::AbstractVector{<:AbstractString}=sort!(tiplabels(net));
        preorder::Bool=true,
        includeleaves::Bool=false)

[`NodeMuRepresentation`](@ref) object for network `net` considered as rooted,
including a mapping from each non-leaf node in `net` to its μ-vector: vector
integers representing the number of paths from the node to each leaf,
with leaves ordered as in `labels`.

`net` is assumed to have a single root.  
`preorder`: whether to preprocess the network with [`directedges!`](@ref) and
[`preorder!`](@ref).  
`includeleaves`: if `true` then the trivial μ-vectors of leaf nodes are included.

Assumptions about tip labels:
- technical: `labels` should be sorted alphabetically
- `labels` should have no repeats, otherwise an error is thrown
- `net` should have unique tip labels, otherwise an error is thrown
- tip labels in `net` should all appear in `labels`, otherwise a warning is sent
  and μ-vectors will be missing the entry for the tips not listed in `labels`
- if `label[i]` is not a tip in `net`, then all μ-vector will have a count of 0
  at index `i`, with no warning.

Used by [`mudistance_rooted`](@ref), and [`edge_murepresentation`](@ref).

## example

```jldoctest
julia> net = readnewick("(((C,(A)#H1),(#H1,E)),O);");

julia> # using PhyloPlots; plot(net, shownodenumber=true);

julia> PhyloNetworks.node_murepresentation(net) # internal nodes only
PhyloNetworks.NodeMuRepresentation
4 taxa in μ-vectors: ["A", "C", "E", "O"]
5 nodes, map node number => μ0 to hybrids and μ-vector to taxa:
  -3 => μ0=2 μ=[2, 1, 1, 0]
  -6 => μ0=1 μ=[1, 0, 1, 0]
  -4 => μ0=1 μ=[1, 1, 0, 0]
  -2 => μ0=2 μ=[2, 1, 1, 1]
  3 => μ0=1 μ=[1, 0, 0, 0]
```
"""
function node_murepresentation(
    net::HybridNetwork,
    labels::AbstractVector{<:AbstractString}=sort!(tiplabels(net));
    preorder::Bool=true,
    includeleaves::Bool=false,
)
    rn = getroot(net)
    rn.leaf && error("please re-root the network: it's leaf $(rn.name)")
    issorted(labels) || error("labels should be sorted. consider using sort() or sort!().")
    allunique_sorted(labels) || error("some input tip labels are repeated.")
    # map: label => index in μ-vectors, for quick access later
    label_map = Dict{String, Int}(l => i for (i,l) in enumerate(labels))
    net_labels = tiplabels(net)
    allunique(net_labels) ||
        error("the network is a MUL-net: has multiple tips with the same label.")
    # warn if a tip in the network is missing from `labels`
    for l in net_labels
        haskey(label_map, l) || @warn("""
            leaf $l in the network is not in the input list of labels:
            its coordinate will be missing from μ-vectors (count of paths to that leaf).""")
    end
    if preorder
        directedges!(net)
        preorder!(net)
    end
    nlabs = length(labels)
    mu_tips = Dict(node.number => zeros(Int,nlabs) for node in net.node)
    # number of paths ending at any hybrid node: initialize by simply counting
    # trivial paths (h) that start & end at h, where h is a hybrid node.
    mu_hybs = Dict(node.number => (node.hybrid ? 1 : 0) for node in net.node)

    # postorder traversal from the tips to the root: reverse of preorder
    for ni in length(net.vec_node):-1:1
        currnode = net.vec_node[ni]
        currnum = currnode.number
        if currnode.leaf # current node can only reach itself
            if haskey(label_map, currnode.name)
                index = label_map[currnode.name]
                mu_tips[currnum][index] = 1
            end
            # other entries remain as initialized at 0
            # mu_hybs[currnum]=1 already if this leaf is also a hybrid (!)
            continue
        end
        for curredge in currnode.edge
            c = getchild(curredge)
            c === currnode && continue # skip parent edge(s)
            chnum = c.number
            mu_tips[currnum] .+= mu_tips[chnum]
            mu_hybs[currnum]  += mu_hybs[chnum]
        end
    end 
    mu_dict = Dict{Int, MuVector}()
    for no in net.node
        no.leaf && !includeleaves && continue # exclude trivial μ-vectors of leaves
        num = no.number
        push!(mu_dict, num => MuVector(mu_hybs[num], mu_tips[num]))
    end
    return NodeMuRepresentation(labels, mu_dict)
end
# check for unique values, assuming the vector is sorted
function allunique_sorted(v::AbstractVector)
    length(v) < 2 && return true
    x1 = v[1]
    for i in 2:length(v)
        x2 = v[i]
        x2 == x1 && return false
        x1 = x2
    end
    return true
end

"""
    edge_murepresentation(net::HybridNetwork,
        labels::AbstractVector{<:AbstractString}=sort!(tiplabels(net));
        preorder::Bool=true,
        suppressroot::Bool=true)

[`EdgeMuRepresentation`](@ref) object for network `net`, considered as a
semidirected network, as described in
[Ribas, Pons & Ané (2026)](https://arxiv.org/abs/2605.06243).
See also the "basic" edge-based μ-representation first described in
[Maxfield, Xu & Ané 2025](https://doi.org/10.1109/TCBBIO.2025.3534780).
This representation maps each edge in `net` to its μ-entry: a set of
1 or 2 tagged μ-vectors).
A tagged μ-vector has a tag and a simple μ-vector.
The simple μ-vector contains the number of paths from the edge to each leaf,
with leaves ordered as in `labels`. It also contains
the number of paths from the edge to a (any) hybrid node.

- `preorder`: whether to preprocess the network (direct edges away from the root
  and calculate a pre-ordering of nodes).
- `suppressroot`: whether to exclude one internal tree edge incident to the root,
  in case the root is incident to exactly 2 edges which are not both hybrid edges.
  In that case, re-rooting the network would suppress that internal tree edge,
  hence its associated μ-entry. This option does *not* modify the network.

`net` is assumed to have a single root component.

See [`node_murepresentation`](@ref) for assumptions about `labels` and tip labels
in `net`.

Used by [`mudistance_semidirected`](@ref).

!!! info "version info"
    In PhyloNetworks v1.4.0 external edges to a leaf were excluded, because
    - they have trivial μ-vectors under the common assumption that leaves are
      incident to a single edge (that is, leaves must not be hybrid nodes), and
    - those incident to the root were given 1 tagged μ-vector, not 2.

## example

```jldoctest
julia> net = readnewick("(((C,(A)#H1),(#H1,E)),O);");

julia> # using PhyloPlots; plot(net, showedgenumber=true);

julia> PhyloNetworks.edge_murepresentation(net)
PhyloNetworks.EdgeMuRepresentation
4 taxa in μ-vectors: ["A", "C", "E", "O"]
root μ-vector: μ0=2 μ=[2, 1, 1, 1]
map edge number => μ-entry:
2 edges in the root component
  4 => tag=:t μ0=1 μ=[1, 1, 0, 0]; tag=:t μ0=1 μ=[1, 0, 1, 1]
  7 => tag=:t μ0=1 μ=[1, 0, 1, 0]; tag=:t μ0=1 μ=[1, 1, 0, 1]
5 edges incident to the root component
  5 => tag=:h μ0=1 μ=[1, 0, 0, 0]; tag=:i μ0=1 μ=[1, 1, 1, 1]
  6 => tag=:t μ0=0 μ=[0, 0, 1, 0]; tag=:i μ0=2 μ=[2, 1, 0, 1]
  9 => tag=:t μ0=0 μ=[0, 0, 0, 1]; tag=:i μ0=2 μ=[2, 1, 1, 0]
  3 => tag=:h μ0=1 μ=[1, 0, 0, 0]; tag=:i μ0=1 μ=[1, 1, 1, 1]
  1 => tag=:t μ0=0 μ=[0, 1, 0, 0]; tag=:i μ0=2 μ=[2, 0, 1, 1]
1 directed edges that may not contain the root
  2 => tag=:t μ0=0 μ=[1, 0, 0, 0]

```
"""
function edge_murepresentation(
    net::HybridNetwork,
    labels::AbstractVector{<:AbstractString}=sort!(tiplabels(net));
    preorder::Bool=true,
    suppressroot::Bool=true,
)
    nodemu = node_murepresentation(net, labels; preorder=preorder,
        includeleaves=true).mu_map
    rn = getroot(net) # node_murepresentation already checked it's not a leaf
    re = rn.edge
    rdeg = length(re)
    if (rdeg==0 && !rn.leaf) || (rdeg==1 && !isexternal(re[1]))
        error("""root with at most 1 child edge:
        please treat the network as a rooted network, or suppress the root,
        e.g. with deleteleaf!(yournetworkname, $(net.rooti); index=true) .""")
    end
    deg2root = rdeg==2
    resolvedroot = deg2root && ( # the only alternative placement for the root
        (re[1].hybrid && re[2].hybrid) || # is along either of its 2 child edges
        (isexternal(re[1]) && isexternal(re[2]))) # but they give the same net
    rho = nodemu[rn.number]
    rμ0 = rho.mu_hybs
    rμv = rho.mu_tips
    mumap_root = Dict{Int, Tuple{TagMuVector,TagMuVector}}()
    mumap_inc  = Dict{Int, Tuple{TagMuVector,TagMuVector}}()
    mumap_dir  = Dict{Int, TagMuVector}()
    for ee in net.edge
        enum = ee.number
        cn = getchild(ee)
        cmu = nodemu[cn.number]
        tmu = TagMuVector((ee.hybrid ? mutag_h : mutag_t), cmu)
        if ee.containroot
            if resolvedroot && getparent(ee) === rn
                push!(mumap_dir, enum => tmu)
                continue
            end
            mu_o = MuVector(rμ0 - cmu.mu_hybs, rμv - cmu.mu_tips)
            if !ee.hybrid && !cn.leaf # inside root component
                push!(mumap_root, enum => (tmu, TagMuVector(mutag_t, mu_o)))
            else # insident, place the :i-tagged μ-vector second
                push!(mumap_inc,  enum => (tmu, TagMuVector(mutag_i, mu_o)))
            end
        else
            push!(mumap_dir, enum => tmu)
        end
    end
    if suppressroot && deg2root && !resolvedroot
        # then exclude an edge in the root component: non-external tree edge
        i1,i2 = (re[1].number < re[2].number ? (1,2) : (2,1))
        e = re[i2] # exclude edge with larger number, if both are eligible
        if e.hybrid || isexternal(e)
            e = re[i1]
            e.hybrid && error("root adjacent to 1 leaf and 1 hybrid")
        end
        delete!(mumap_root, e.number)
    end
    return EdgeMuRepresentation(labels, rho, mumap_root, mumap_inc, mumap_dir)
end

"""
Symmetric distance between two vectors, assumed sorted in ascending order:
number of elements in one vector but not in the order, each element considered
as many times as it appears (with multiplicity). The vectors are *not* mutated.
"""
function symmetricdistance_sorted(
    v1::Vector{T},
    v2::Vector{T},
    lt=isless,
) where T
    d = 0
    l1 = length(v1); l2 = length(v2)
    l1==0 && return l2
    l2==0 && return l1
    x1 = v1[l1]; x2 = v2[l2] # last: largest
    while true
        sgn = cmp(lt, x1, x2)
        if sgn == 1 # x1>x2
            l1 -= 1; d += 1
            if l1>0
                x1 = v1[l1]
            else break; end
        elseif sgn == -1 # x2>x1
            l2 -= 1; d += 1
            if l2>0
                x2 = v2[l2]
            else break; end
        else
            l1 -= 1; l2 -= 1
            if l1>0 && l2>0
                x1 = v1[l1]; x2 = v2[l2]
            else break; end
        end
    end
    return d + l1 + l2 # at least one of them is 0
end
# to ignore i-tagged μ-vectors on the fly, without allocating memory
function symmetricdistance_sorted_ignoreitag(
    v1i::Vector{Tuple{T,T}},
    v1d::Vector{T},
    v2i::Vector{Tuple{T,T}},
    v2d::Vector{T},
    lt=isless,
) where T
    d = 0
    l1i = length(v1i); l1d = length(v1d); l1 = l1i + l1d
    l2i = length(v2i); l2d = length(v2d); l2 = l2i + l2d
    l1 == 0 && return l2
    l2 == 0 && return l1
    cat1, x1 = nextmaxcat(v1i, v1d, l1i, l1d)
    cat2, x2 = nextmaxcat(v2i, v2d, l2i, l2d)
    while true
        sgn = cmp(lt, x1, x2)
        if sgn == 1 # x1>x2
            l1 -= 1; d += 1
            if cat1 == :i; l1i -= 1; else l1d -= 1; end
            if l1>0
                cat1, x1 = nextmaxcat(v1i, v1d, l1i, l1d)
            else break; end
        elseif sgn == -1 # x2>x1
            l2 -= 1; d += 1
            if cat2 == :i; l2i -= 1; else l2d -= 1; end
            if l2>0
                cat2, x2 = nextmaxcat(v2i, v2d, l2i, l2d)
            else break; end
        else
            l1 -= 1; l2 -= 1
            if cat1 == :i; l1i -= 1; else l1d -= 1; end
            if cat2 == :i; l2i -= 1; else l2d -= 1; end
            if l1>0 && l2>0
                cat1, x1 = nextmaxcat(v1i, v1d, l1i, l1d)
                cat2, x2 = nextmaxcat(v2i, v2d, l2i, l2d)
            else break; end
        end
    end
    return d + l1 + l2 # at least one of them is 0
end
function nextmaxcat(vi, vd, li, ld)
    if li == 0
        return (:d, vd[ld])
    elseif ld == 0 || lt(vd[ld], vi[li])
        return (:i, vi[li])
    end
    return (:d, vd[ld])
end

"""
    mudistance(m1::NodeMuRepresentation, m2::NodeMuRepresentation, useμ0::Bool)

Distance between two node-based μ-representations: number of μ-vectors,
considered with their multiplicity, in one but not in the other.
If `useμ0` is false, then only the number of paths to leaves are compared,
otherwise the number of paths to hybrid nodes is also compared.

Warning: this method assumes, without checking, that `m1` and `m2` have the same
labels (and in the same order), for efficiency.
"""
function mudistance(
    m1::NodeMuRepresentation,
    m2::NodeMuRepresentation,
    useμ0::Bool,
)
    # algorithm based on `mu_vec` being already sorted, in ascending order
    lt = useμ0 ? isless : isless_μtips
    return symmetricdistance_sorted(m1.mu_vec, m2.mu_vec, lt)
end

"""
    mudistance_rooted(net1::HybridNetwork, net2::HybridNetwork;
        labels=sort!(union(tiplabels(net1), tiplabels(net2))),
        preorder::Bool=true,
        useμ0::Bool=true)
    mudistance_rooted(nets::AbstractVector{HybridNetwork}; kwargs...)

Distance or matrix of distances between the networks, considered as rooted,
based on their node-based μ-representation: number of nodes in one network
whose μ-vector does not match with that of a node in the other network.
The μ-vector of a node has n+1 coordinates: one for each tip taxon, and one
for hybrids. Each value counts the number of paths starting at the node and
ending at a specific leaf, or ending at a hybrid node (any hybrid node).
With `useμ0=false`, the number of paths to hybrids is ignored in the comparison.

See [Cardona, Rossello & Valiente (2009)](https://doi.org/10.1109/TCBB.2007.70270)
for the original μ-distance (correspondint to `useμ0=false`) and
[Cardona et al. (2024)](https://doi.org/10.1109/TCBB.2024.3361390) for the
extension using paths ending at hybrids.

The distance implemented here considers one μ-vector per non-leaf node:
without ignoring any non-leaf node (including degree-2 nodes if any), and
ignoring leaves whose μ-vectors are trivial.

Assumption: networks have a single root.

See also [`node_murepresentation`](@ref).

## keyword arguments

- `preorder` indicates whether to preprocess each network with
  [`directedges!`](@ref) and [`preorder!`](@ref).
- `labels`: tip labels that each μ-vector will track, by default the union of
  taxon labels across both networks. If a taxon is absent from `net1` (say) but
  tracked because present in `net2`, then there are 0 paths from any node in
  `net1` to that taxon.

If the networks have different tip labels, consider pruning leaves that are
not shared between them beforehand, using [`deleteleaf!`](@ref),
to compare the subnetworks on their shared leaf set.

If a non-default `label` set is used, such as the sorted shared taxon set
`sort!(intersect(tiplabels(net1), tiplabels(net2)))`, then the distance may drop,
because the number of paths to non-shared taxa will be ignored.
But *all* nodes still count towards the distance, including degree-2 nodes
that would be suppressed in subnetworks.
"""
function mudistance_rooted(
    net1::HybridNetwork,
    net2::HybridNetwork;
    labels::AbstractVector{<:AbstractString}=sort!(union(tiplabels(net1), tiplabels(net2))),
    preorder::Bool=true,
    useμ0::Bool=true,
)
    nodemu1 = node_murepresentation(net1, labels; preorder=preorder)
    nodemu2 = node_murepresentation(net2, labels; preorder=preorder)
    return mudistance(nodemu1, nodemu2, useμ0)
end
function mudistance_rooted(
    nets::AbstractVector{HybridNetwork};
    preorder::Bool=true,
    useμ0::Bool=true,
)
    labs = sort!(tiplabels(nets)) # union
    nn = length(nets)
    M = zeros(Int, nn, nn)
    mr = [node_murepresentation(n, labs; preorder=preorder) for n in nets]
    for (i1, m1) in enumerate(mr)
        for i2 in 1:(i1-1)
            d = mudistance(m1, mr[i2], useμ0)
            M[i1,i2] = d; M[i2,i1] = d
        end
    end
    return M
end

"""
    mudistance(m1::EdgeMuRepresentation, m2::EdgeMuRepresentation,
               useμ0::Bool, useitag::Bool, userootμ::Bool)

Distance between two edge-based μ-representations: number of μ-entries,
considered with their multiplicities, in one but not in the other.
- If `useμ0` is false, then only the number of paths to leaves are compared,
  otherwise the number of paths to hybrid nodes is also compared.
- If `useitag` is false, then the μ-vectors with an `:i` tag in μ-entries that
  have such an `:i`-tagged μ-vector are ignored.
  Only their other μ-vector are compared, treating such edges like those
  with a fixed direction and on which the root cannot be placed.
- `userootμ` controls whether the root μ-vector in included in the
  set of μ-entries being compared. It is ignored by default, because different
  root μ-vectors for networks with a single root would cause all μ-entries in
  their root components to be different.

Warning: this method assumes, without checking, that `m1` and `m2` have the same
labels (and in the same order), for efficiency.
"""
function mudistance(
    m1::EdgeMuRepresentation,
    m2::EdgeMuRepresentation,
    useμ0::Bool,
    useitag::Bool,
    userootμ::Bool,
)
    lt = useμ0 ? isless : isless_μtips
    d  = symmetricdistance_sorted(m1.muvec_rootcomp, m2.muvec_rootcomp, lt)
    if useitag
        d += symmetricdistance_sorted(m1.muvec_directed, m2.muvec_directed, lt)
        d += symmetricdistance_sorted(m1.muvec_incident, m2.muvec_incident, lt)
    else
        d += symmetricdistance_sorted_ignoreitag(
            m1.muvec_incident, m1.muvec_directed,
            m2.muvec_incident, m2.muvec_directed, lt)
    end
    if userootμ
        d += (useμ0 ? m1.mu_root != m2.mu_root :
              m1.mu_root.mu_tips != m2.mu_root.mu_tips)
    end
    return d
end

"""
    mudistance_semidirected(net1::HybridNetwork, net2::HybridNetwork; ...)
    mudistance_semidirected(nets::AbstractVector{HybridNetwork}; ...)

Distance or matrix of distances between the networks, considered as semidirected,
based on their edge-based μ-representation: number of edges in one network
whose μ-entry does not match with that of an edge in the other network,
as defined by [Ribas, Pons & Ané (2026)](https://arxiv.org/abs/2605.06243).

The "simple" μ-vector of an edge, assigned some specific direction,
has n+1 coordinates: one for hybrids and one for each tip taxon.
Each value counts the number of paths starting with the edge (directed as
specified) and ending at a (any) hybrid node, or ending at a specific leaf.

Each μ-entry includes 1 or 2 μ-vector(s) and tag(s) indicating the edge type
(e.g. hybrid). The number of μ-vector(s) depends on whether the root can be
placed along the edge or not.

If the networks have different tip labels, the union of all their labels
is considered to build μ-vectors, and the resulting distance must be positive.
Consider pruning leaves that are not shared between the two networks beforehand,
using [`deleteleaf!`](@ref), to compare the subnetworks on their shared leaf set.

Assumption: networks have a single root.

See also [`edge_murepresentation`](@ref).

# keyword arguments

- `preorder=true` indicates whether to preprocess each network with
  [`directedges!`](@ref) and [`preorder!`](@ref).
- `suppressroot=true` to exclude one internal tree edge incident to the root,
  in case the root is incident to exactly 2 edges which are not both hybrid edges.
  In that case, re-rooting the network would suppress that internal tree edge,
  hence its associated μ-entry. This option does *not* modify the network.
- `labels`: tip labels that each μ-vector will track, by default the union of
  taxon labels across both networks, sorted alphabetically.
  If a taxon is absent from `net1` (say) but tracked because present in `net2`,
  then there are 0 paths from any edge in `net1` to that taxon.

If the networks have different tip labels, consider pruning leaves that are
not shared between them beforehand, using [`deleteleaf!`](@ref),
to compare the subnetworks on their shared leaf set.

If a non-default `label` set is used, such as the sorted shared taxon set
`sort!(intersect(tiplabels(net1), tiplabels(net2)))`, then the distance may drop,
as the number of paths to non-shared taxa will be ignored. But *all* edges still
count towards the distance, whereas edges below degree-2 nodes would be
suppressed in subnetworks.

Other keyword arguments are to ignore some properties of the network's
μ-representation, which cause a shorter dissimilarity, and also may lead
to some distinct (non-isomorphic) networks to have dissimilarity 0.
- `useμ0=true`: set to `false` to ignore the the number of paths to some (any)
  hybrid in each μ-vector.
- `useitag=true`: set to `false` to ignore the μ-vectors with an `:i` tag,
  in μ-entries that have such an `:i`-tagged μ-vector.
- `userootμ=false` to ignore the root μ-entry.
  This default is made because we assume here that networks have a single root.
  If two single-rooted networks have a different root μ-entry, they also have
  different μ-entries for all edges that can contain their root.

The μ-dissimilarity in
[Maxfield, Xu & Ané 2025](https://doi.org/10.1109/TCBBIO.2025.3534780)
corresponds to `userootμ=true, useitag=false, useμ0=false`. The μ-dissimilarity
in [Ribas, Pons & Ané 2026](https://arxiv.org/abs/2605.06243) corresponds to
`userootμ=true, useitag=true, useμ0=true`.

!!! info "version info"
    In PhyloNetworks v1.4.0, this distance was closer to the original μ-distance
    of [Maxfield, Xu & Ané 2025](https://doi.org/10.1109/TCBBIO.2025.3534780),
    although it was not identical, because the μ-entry for a hybrid edge adjacent
    to the root component had a single μ-vector, tagged as 'incident',
    and by considering "extended" μ-vectors with the counts of paths ending at
    hybrids, as defined on rooted networks by
    [Cardona et al. (2024)](https://doi.org/10.1109/TCBB.2024.3361390).

"""
function mudistance_semidirected(
    net1::HybridNetwork,
    net2::HybridNetwork;
    labels::AbstractVector{<:AbstractString}=sort!(union(tiplabels(net1), tiplabels(net2))),
    useμ0::Bool=true,
    useitag::Bool=true,
    userootμ::Bool=false,
    kwargs...
)
    m1 = edge_murepresentation(net1, labels; kwargs...)
    m2 = edge_murepresentation(net2, labels; kwargs...)
    return mudistance(m1, m2, useμ0, useitag, userootμ)
end
function mudistance_semidirected(
    nets::AbstractVector{HybridNetwork};
    labels::AbstractVector{<:AbstractString}=sort!(tiplabels(nets)),
    useμ0::Bool=true,
    useitag::Bool=true,
    userootμ::Bool=false,
    kwargs...
)
    nn = length(nets)
    M = zeros(Int, nn, nn)
    mr = [edge_murepresentation(n, labels; kwargs...) for n in nets]
    for (i1, m1) in enumerate(mr)
        for i2 in 1:(i1-1)
            d = mudistance(m1, mr[i2], useμ0, useitag, userootμ)
            M[i1,i2] = d; M[i2,i1] = d
        end
    end
    return M
end
