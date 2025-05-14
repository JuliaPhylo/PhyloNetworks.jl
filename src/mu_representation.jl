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
end
function Base.show(io::IO, obj::NodeMuRepresentation)
    labs = tiplabels(obj)
    disp = "$(typeof(obj))\n$(length(labs)) taxa, $(length(obj.mu_map)) nodes,"
    disp *= "\ntaxon order in μ-vectors:\n" * string(labs)
    disp *= "\nmap node number => μ-vector:\n" * string(obj.mu_map)
    println(io, disp)
end


"""
    @Override ==(a::NodeMuRepresentation, b::NodeMuRepresentation)
    Check if two NodeMuRepresentation objects are equal using the distance function.
"""
function ==(a::NodeMuRepresentation, b::NodeMuRepresentation)
    return node_vec_distance(a,b)==0
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

"""
    @Override ==(a::EdgeMuRepresentation, b::EdgeMuRepresentation)

Check if two EdgeMuRepresentation objects are equal using the distance function.
"""
function ==(a::EdgeMuRepresentation, b::EdgeMuRepresentation)
    return edge_vec_distance(a,b)==0
end

"""
    node_murepresentation(net::HybridNetwork,
        labels::AbstractVector{<:AbstractString},
        preprocess::Bool=true)

[`NodeMuRepresentation`](@ref) object for network `net` considered as rooted,
including a mapping from each non-leaf node in `net` to its μ-vector: vector
integers representing the number of paths from the node to each leaf,
with leaves ordered as in `labels`.

`net` is assumed to have a single root.

`preprocess`: boolean indicating whether to preprocess the network with
`directedges!` and `preorder!`.

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
    preprocess::Bool=true
)
    allunique(labels) || error("some input tip labels are repeated.")
    # map: label => index in μ-vectors, for quick access later
    label_map = Dict{String, Int}(l => i for (i,l) in enumerate(labels))
    Ntiplabels = tiplabels(net)
    allunique(Ntiplabels) ||
        error("the network is a MUL-net: has multiple tips with the same label.")
    # warn if a tip in the network is missing from `labels`
    for l in Ntiplabels
        haskey(label_map, l) || @warn("""
            leaf $l in the network is not in the input list of labels:
            its coordinate will be missing from μ-vectors (count of paths to that leaf).""")
    end
    if preprocess
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
    return NodeMuRepresentation(mu_dict, labels)
end

"""
    edge_murepresentation(net::HybridNetwork,
        labels::AbstractVector{<:AbstractString},
        preprocess::Bool=true)

[`EdgeMuRepresentation`](@ref) object for network `net`, considered as a
semidirected network. This representation maps each internal edge in `net` to
its μ-entry (1 or 2 tagged μ-vectors). A μ-vector contains the number of paths
from the edge to each leaf, with leaves ordered as in `labels`. It also contains
the number of paths from the edge to a (any) hybrid node.

External edges are edges whose child is a leaf, and have trivial μ-vectors under
the common assumption that leaves are incident to a single edge (that is, leaves
must not be hybrid nodes).

`net` is assumed to have a single root component.

`preprocess`: boolean indicating whether to preprocess the network (direct edges
away from the root and calculate a pre-ordering of nodes).

See [`node_murepresentation`](@ref) for assumptions about `labels` and tip labels
in `net`.
"""
function edge_murepresentation(
    net::HybridNetwork,
    labels::AbstractVector{<:AbstractString},
    preprocess::Bool=true
)
    nodemu = node_murepresentation(net, labels, preprocess).mu_map
    rho = nodemu[getroot(net).number]
    mumap_rootcomp::Dict{Int, Tuple{MuVector, MuVector}}()
    mumap_directed::Dict{Int, Tuple{MuVector, MuTag}}()
    for ee in net.edge
        enum = ee.number
        cn = getchild(ee)
        if cn.leaf
            ee.hybrid && @warn("edge $enum is hybrid yet has a leaf child: its μ-entry won't be stored")
            continue # skip this edge
        end
        if !ee.hybrid && ee.containroot
            pnum = getparent(ee).number
            push!(mumap_rootcomp, enum => (nodemu[cn.number], rho - nodemu[pnum]))
        else
            tag = (ee.containroot ? mutag_i : (ee.hybrid ? mutag_h : mutag_t))
            push!(mumap_directed, enum => (nodemu[cn.number], tag))
        end
    end
    return EdgeMuRepresentation(rho, mumap_rootcomp, mumap_directed, labels)
end

"""
    node_vec_distance(net1::NodeMuRepresentation, net2::NodeMuRepresentation)

Calculate the distance between two node-based μ-vectors.

# Arguments
- `net1`: A NodeMuRepresentation object
- `net2`: A NodeMuRepresentation object

# Returns
- `dist`: The distance between the two node-based μ-vectors.
"""
function node_vec_distance(net1::NodeMuRepresentation, net2::NodeMuRepresentation)
    #convert to not use multiet
    m1 = Multiset(collect(values(net1.mu_map)))
    m2 = Multiset(collect(values(net2.mu_map)))
    d1 = setdiff(m1, m2)
    d2 = setdiff(m2, m1)
    dist = length(d1) + length(d2)
    return dist
end

"""
    distance_node(net1::HybridNetwork, net2::HybridNetwork)

Distance between two networks based on their node-based μ-vectors.

Warning: Networks cannot have multiple roots.

# Arguments
- `net1`: A HybridNetwork object
- `net2`: A HybridNetwork object

# Returns
- `dist`: The distance between the two networks based on their node-based μ-vectors.
"""
function network_node_mu_distance(net1::HybridNetwork, net2::HybridNetwork, preprocess=true)
    # make sure mu vectors are consistent
    if length(tiplabels(net1)) > length(tiplabels(net2))
        labels = tiplabels(net1)
    else
        labels = tiplabels(net2)
    end

    # Get the mu values for the nodes
    nodemu1 = node_murepresentation(net1,labels, preprocess)
    nodemu2 = node_murepresentation(net2,labels, preprocess)
    return node_vec_distance(nodemu1, nodemu2)
    
end

"""
    edge_vec_distance(net1::EdgeMuRepresentation, net2::EdgeMuRepresentation)

Calculate the distance between two edge-based μ-vectors.

# Arguments
- `net1`: A EdgeMuRepresentation object
- `net2`: A EdgeMuRepresentation object

# Returns
- `dist`: The distance between the two edge-based μ-vectors.
"""
function edge_vec_distance(net1::EdgeMuRepresentation, net2::EdgeMuRepresentation)

    m1 = Multiset(collect(values(net1.mu_map)))
    m2 = Multiset(collect(values(net2.mu_map)))
    d1 = setdiff(m1, m2)
    d2 = setdiff(m2, m1)
    dist = length(d1) + length(d2)
    return dist
end
"""
    distance_edge(net1::HybridNetwork, net2::HybridNetwork)

Distance between two networks based on their edge-based μ-vectors.

Warning: networks cannot have multiple roots.

# Arguments
- `net1`: A HybridNetwork object
- `net2`: A HybridNetwork object

# Returns
- `dist`: The distance between the two networks based on their edge-based μ-vectors.
"""
function network_edge_mu_distance(net1::HybridNetwork, net2::HybridNetwork, preprocess=true)
    # make sure mu vectors are consistent
    if length(tiplabels(net1)) > length(tiplabels(net2))
        labels = tiplabels(net1)
    else
        labels = tiplabels(net2)
    end

    # Get the mu values for the nodes
    edgemu1 = edge_murepresentation(net1,labels, preprocess)
    edgemu2 = edge_murepresentation(net2,labels, preprocess)
    return edge_vec_distance(edgemu1, edgemu2)
    
end