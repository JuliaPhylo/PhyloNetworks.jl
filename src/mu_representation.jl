struct MuTag
    tag::Int8
end
# integer that happens to match in-degree
const mutag_r = MuTag( 0) # root node. integer that matches in-degree
const mutag_t = MuTag( 1) # tree edges / nodes
const mutag_h = MuTag( 2) # hybrid edges / nodes
const mutag_i = MuTag(-1) # incident to root component
function Base.show(io::IO, obj::MuTag)
    tag = obj.tag
    s = (tag==0 ? "root" : (tag==1 ? "tree" : (tag==2 ? "hybrid" : "incident")))
    print(io, s)
end
isless(t1::MuTag, t2::MuTag) = isless(t1.tag, t2.tag)

"""
    MuVector

μ-vector for a single node (see [Cardona et al. 2024](https://10.1109/TCBB.2024.3361390))
separated into
- μ0: number of paths starting at the node or edge and ending at any hybrid node
- μ-vector counting the number of paths from that node to each leaf.

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
    s  = "    μ0 to hybrids: " * string(obj.mu_hybs)
    s *= "\n    μ-vector to tips:" * string(obj.mu_tips)
    s *= "\n"
    println(io, disp)
end
# comparing μ-vectors: lexicographic order using
# 1. μ0, number of paths to hybrids
# 2. μ-vector itself, lexicographically
function isless(m1::MuVector, m2::MuVector)
    if m1.mu_hybs == m2.mu_hybs
        return isless(m1.mu_tips, m2.mu_tips)
    end
    return isless(m1.mu_hybs, m2.mu_hybs)
end
has_0μentries_at(m::MuVector, i) = all(m.mu_tips[i] .== 0)
function has_0μentries_at(ms::Tuple{MuVector,MuVector}, i)
    b = has_0μentries_at(ms[1],i) && has_0μentries_at(ms[2],i)
    return b
end

abstract type MuRepresentation end

tiplabels(m::MuRepresentation) = m.tiplabels

"""
    NodeMuRepresentation

Representation of a rooted network with one [`MuVector`](@ref) per non-leaf node.
This entry contains the node's μ-vector: a vector of integers counting the
number of paths from that node to each leaf.
`tiplabels` gives the order of leaves in μ-vectors.

Leaf nodes have trivial μ-vector, hence are not included in this node-based
μ-representation.

See [Cardona et al. 2024](https://doi.org/10.1109/TCBB.2024.3361390).
"""
struct NodeMuRepresentation <: MuRepresentation
    "vector of tip (leaf) labels: listed in the same order as in μ-vectors"
    tiplabels::Vector{String}
    "Dictionary: maps node number => [`MuVector`](@ref) object for that node"
    mu_map::Dict{Int, MuVector}
    "sorted vector of μ-vectors"
    mu_vec::Vector{MuVector}
end
function Base.show(io::IO, obj::NodeMuRepresentation)
    labs = tiplabels(obj)
    disp = "$(typeof(obj))\n$(length(labs)) taxa, $(length(obj.mu_map)) nodes,"
    disp *= "\ntaxon order in μ-vectors:\n" * string(labs)
    disp *= "\nmap node number => μ-vector:\n" * string(obj.mu_map)
    println(io, disp)
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

function indexin_check0μentries(m1::MuRepresentation, m2::MuRepresentation)
    l1 = tiplabels(m1); l2 = tiplabels(m2)
    o2 = indexin(l1, l2) # order to use for μ-vectors in m2: l2[o2] == l1
    # if a label is in l1 but not in m2: all its entries should be 0 in m1
    in1not2_i = findall(isnothing, o2) # indices in o2 and in m1
    if !isempty(in1not2_i)
        has_0μentries_at(m1, in1not1_i) || return nothing
        deleteat!(o2, in1not2_i)
    end
    o1 = collect(1:length(l1))
    if length(l2) + length(in1not2_i) > length(l1)
        # then some labels in are l2 but not in l1: find their indices in l2
        in2not1_i = findall(!in(l1), l2)
        has_0μentries_at(m2, in2not1_i) || return nothing
        deleteat!(o1, in2not1_i)
    end
    return (o1, o2)
end

"""
    ==(m1::NodeMuRepresentation, m2::NodeMuRepresentation)

Equality of two [`NodeMuRepresentation`](@ref) objects, requiring that:
- If `m1` has labels that `m2` doesn't have, then all corresponding entries
  for these labels in `m1` should be 0 -- as if a label is tracked in `m1`
  but absent from the original network, and not tracked in `m2`.
  Vice versa: labels in `m2` but not in `m1` should have 0 entries in `m2`.
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
    for (μ1, μ2) in zip(m1.mu_vec, m2.mu_vec)
        μ1[o1] == μ2[o2] || return false
    end
    return true
end


"""
    EdgeMuRepresentation

Representation of a semidirected network with one μ-entry for each edge, see
[Maxfield, Xu & Ané 2025](https://doi.org/10.1109/TCBBIO.2025.3534780).

- If an edge has a fixed direction regardless of where the network is rooted,
  its entry is 1 tagged μ-vector.
  The tag says whether the edge is a tree e}dge or hybrid edge.
- If an edge is in the root component, that is, both of its incident nodes can
  root the network, then its entry is a set of 2 tagged μ-vectors (each with a
  tree tag because the edge must be a tree edge).
- If an edge is not in the root component, but is incident to it (one of its
  incident nodes can be the root), then the root may be placed along the edge,
  and its entry is a set of 2 tagged μ-vectors. The tag is the edge type for
  the μ-vector corresponding to the edge's direction.

Each μ-vector is a vector of integers representing the number of paths from the
edge to each leaf. `tiplabels` gives the order of leaves in μ-vectors.
"""
struct EdgeMuRepresentation <: MuRepresentation
    "vector of tip (leaf) labels: listed in the same order as in μ-vectors"
    tiplabels::Vector{String}
    "root μ-vector, assuming the network has a single root component"
    mu_root::MuVector
    "dictionary for root components: edge number => (μ-vector1, μ-vector2)"
    mumap_rootcomp::Dict{Int, Tuple{MuVector,MuVector}}
    "dictionary for fixed direction: edge number => (μ-vector, tag)"
    mumap_directed::Dict{Int, Tuple{MuTag,MuVector}}
    "sorted vector of sorted μ-entries in the root component"
    muvec_rootcomp::Vector{Tuple{MuVector,MuVector}}
    "sorted vector of μ-entries in the directed part"
    muvec_directed::Vector{Tuple{MuTag,MuVector}}
end
function Base.show(io::IO, obj::EdgeMuRepresentation)
    labs = tiplabels(obj)
    s = "$(typeof(obj))\n$(length(labs)) taxa, $(length(obj.mu_map)) edges,"
    s *= "\ntaxon order in μ-vectors:\n" * string(labs)
    s *= "root μ-vector:\n" * string(obj.mu_root)
    s *= "\nmap edge number => μ-entry, in root component:\n"
    for (k,v) in mumap_rootcomp
        s *= "tree edge $k:\n" * string(v[1]) * string(v[2])
    end
    s *= "\nmap edge number => μ-entry, in directed part:\n"
    for (k,v) in mumap_rootcomp
        s *= "edge $k, tagged $(v[1]):\n" * string(v[2])
    end
    println(io, s)
end

# constructor: builds sorted vectors from dictionaries
#              re-using memory for each μ-vector
function EdgeMuRepresentation(labels, mu_root, map_rootcomp, map_directed)
    vec_root = Vector{Tuple{MuVector,MuVector}}()
    for (v1,v2) in values(map_rootcomp)
        if v1 < v2
            push!(vec_root, (v1,v2))
        else
            push!(vec_root, (v2,v1))
        end
    end
    sort!(vec_root)
    vec_dir = collect(values(map_directed))
    sort!(vec_dir)
    return EdgeMuRepresentation(labels, mu_root, map_rootcomp, map_directed,
                                vec_root, vec_dir)
end

function has_0μentries_at(m::EdgeMuRepresentation, i)
    has_0μentries_at(m.mu_root, i) || return false
    for mr in m.muvec_rootcomp
        has_0μentries_at(mr, i) || return false
    end
    for md in m.muvec_directed
        has_0μentries_at(md[2], i) || return false
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
  (either in the root component, or directed with the same tags).
- All μ-entries should have the same subvectors corresponding to the labels
  shared by `m1` and `m2`.
"""
function ==(m1::EdgeMuRepresentation, m2::EdgeMuRepresentation)
    # order for μ-vectors: tiplabels(m1)[o1] == tiplabels(m2)[o2], others are 0
    o12 = indexin_check0μentries(m1, m2)
    isnothing(o12) && return false
    o1, o2 = o12
    m1.mu_root[o1] == m2.mu_root[o2] || return false
    length(m1.muvec_rootcomp) == length(m2.muvec_rootcomp) || return false
    length(m1.muvec_directed) == length(m2.muvec_directed) || return false
    for (ms1, ms2) in zip(m1.muvec_rootcomp, m2.muvec_rootcomp)
        ms1[1][o1] == ms2[1][o2] || return false
        # no need to check second vector: bc 1 root and same root μ-vector
        # ms1[2][o1] == ms2[2][o2] || return false
    end
    for (ms1, ms2) in zip(m1.muvec_directed, m2.muvec_directed)
        ms1[1] == ms2[1] || return false         # same tags
        ms1[2][o1] == ms2[2][o2] || return false # same μ-vectors
    end
    return true
end


"""
    node_murepresentation(net::HybridNetwork,
        labels::AbstractVector{<:AbstractString},
        preorder::Bool=true)

[`NodeMuRepresentation`](@ref) object for network `net` considered as rooted,
including a mapping from each non-leaf node in `net` to its μ-vector: vector
integers representing the number of paths from the node to each leaf,
with leaves ordered as in `labels`.

`net` is assumed to have a single root.

`preorder`: boolean indicating whether to preprocess the network with
[`directedges!`](@ref) and [`preorder!`](@ref).

Assumptions about tip labels:
- `labels` should have no repeats, otherwise an error is thrown
- `net` should have unique tip labels, otherwise an error is thrown
- tip labels in `net` should all appear in `labels`, otherwise a warning is sent
  and μ-vectors will be missing the entry for the tips not listed in `labels`
- if `label[i]` is not a tip in `net`, then all μ-vector will have a count of 0
  at index `i`, with no warning.
"""
function node_murepresentation(
    net::HybridNetwork,
    labels::AbstractVector{<:AbstractString},
    preorder::Bool=true
)
    allunique(labels) || error("some input tip labels are repeated.")
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
    for currnode in reverse(net.vec_node)
        currnum = currnode.number
        if currnode.leaf # current node can only reach itself
            if haskey(label_map, currnode.name)
                index = label_map[currnode.name]
                mu_tips[currnum][index] = 1
            end
            # other entries remain as initialized at 0
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
        no.leaf && continue # do *not* store the trivial μ-vectors of leaves
        num = no.number
        push!(mu_dict, num => MuVector(mu_tips[num], mu_hybs[num]))
    end
    return NodeMuRepresentation(labels, mu_dict)
end

"""
    edge_murepresentation(net::HybridNetwork,
        labels::AbstractVector{<:AbstractString},
        preorder::Bool=true)

[`EdgeMuRepresentation`](@ref) object for network `net`, considered as a
semidirected network. This representation maps each internal edge in `net` to
its μ-entry (1 or 2 tagged μ-vectors). A μ-vector contains the number of paths
from the edge to each leaf, with leaves ordered as in `labels`. It also contains
the number of paths from the edge to a (any) hybrid node.

External edges are edges whose child is a leaf, and have trivial μ-vectors under
the common assumption that leaves are incident to a single edge (that is, leaves
must not be hybrid nodes).

`net` is assumed to have a single root component.

`preorder`: boolean indicating whether to preprocess the network (direct edges
away from the root and calculate a pre-ordering of nodes).

See [`node_murepresentation`](@ref) for assumptions about `labels` and tip labels
in `net`.
"""
function edge_murepresentation(
    net::HybridNetwork,
    labels::AbstractVector{<:AbstractString},
    preorder::Bool=true
)
    nodemu = node_murepresentation(net, labels, preorder).mu_map
    rho = nodemu[getroot(net).number]
    rμ0 = rho.mu_hybs
    rμv = rho.mu_tips
    mumap_rootcomp::Dict{Int, Tuple{MuVector, MuVector}}()
    mumap_directed::Dict{Int, Tuple{MuVector, MuTag}}()
    for ee in net.edge
        enum = ee.number
        cn = getchild(ee)
        if cn.leaf
            ee.hybrid && @warn("edge $enum is hybrid yet has a leaf child: its μ-entry won't be stored")
            continue # skip this edge
        end
        cmu = nodemu[cn.number]
        if !ee.hybrid && ee.containroot
            mu_otherdirection = MuVector(rμ0 - cmu.mu_hybs, rμv - cmu.mu_tips)
            push!(mumap_rootcomp, enum => (cmu, mu_otherdirection))
        else
            tag = (ee.containroot ? mutag_i : (ee.hybrid ? mutag_h : mutag_t))
            push!(mumap_directed, enum => (cmu, tag))
        end
    end
    return EdgeMuRepresentation(rho, mumap_rootcomp, mumap_directed, labels)
end

"""
    mudistance(m1::NodeMuRepresentation, m2::NodeMuRepresentation)

Calculate the distance between two node-based μ-representations: number of
μ-vectors in one but not in the other.

Warning: this method assumes, without checking, that `m1` and `m2` have the same
labels (and in the same order), for efficiency.
"""
function mudistance(m1::NodeMuRepresentation, m2::NodeMuRepresentation)
    #convert to not use multiet
    m1 = Multiset(collect(values(net1.mu_map)))
    m2 = Multiset(collect(values(net2.mu_map)))
    d1 = setdiff(m1, m2)
    d2 = setdiff(m2, m1)
    dist = length(d1) + length(d2)
    return dist
end

"""
    mudistance_rooted(net1::HybridNetwork, net2::HybridNetwork, preorder::Bool=true)

Distance between two networks, considered as rooted networks, based on their
node-based μ-representation: number of nodes in one network whose μ-vector does
not match with that of a node in the other network.
See [Cardona et al. 2024](https://10.1109/TCBB.2024.3361390).

If `net1` and `net2` have different tip labels, the union of all their labels
is considered to build μ-vectors, and the resulting distance must be positive.

Consider pruning leaves that are not shared between the two networks beforehand,
using [`deleteleaf!`](@ref), to compare the subnetworks on their shared leaf set.

Assumption: networks have a single root.
"""
function mudistance_rooted(
    net1::HybridNetwork,
    net2::HybridNetwork,
    preorder::Bool=true
)
    labels = union(tiplabels(net1), tiplabels(net2))
    nodemu1 = node_murepresentation(net1, labels, preorder)
    nodemu2 = node_murepresentation(net2, labels, preorder)
    return mudistance(nodemu1, nodemu2)
end

"""
    mudistance(m1::EdgeMuRepresentation, m2::EdgeMuRepresentation)

Calculate the distance between two edge-based μ-representations: number of
μ-entries in one but not in the other.

Warning: this method assumes, without checking, that `m1` and `m2` have the same
labels (and in the same order), for efficiency.
"""
function mudistance(m1::EdgeMuRepresentation, m2::EdgeMuRepresentation)
    m1 = Multiset(collect(values(net1.mu_map)))
    m2 = Multiset(collect(values(net2.mu_map)))
    d1 = setdiff(m1, m2)
    d2 = setdiff(m2, m1)
    dist = length(d1) + length(d2)
    return dist
end
"""
    mudistance_semidirected(net1::HybridNetwork, net2::HybridNetwork, preorder::Bool=true)

Distance between two networks, considered as semidirected, based on their
edge-based μ-representation: number of edges in one network whose μ-entry does
not match with that of an edge in the other network.
See [Maxfield, Xu & Ané 2025](https://doi.org/10.1109/TCBBIO.2025.3534780).

If `net1` and `net2` have different tip labels, the union of all their labels
is considered to build μ-vectors, and the resulting distance must be positive.

Consider pruning leaves that are not shared between the two networks beforehand,
using [`deleteleaf!`](@ref), to compare the subnetworks on their shared leaf set.

Assumption: networks have a single root.
"""
function mudistance_semidirected(
    net1::HybridNetwork,
    net2::HybridNetwork,
    preorder::Bool=true,
)
    labels = union(tiplabels(net1), tiplabels(net2))
    m1 = edge_murepresentation(net1, labels, preorder)
    m2 = edge_murepresentation(net2, labels, preorder)
    return mudistance(m1, m2)
end
