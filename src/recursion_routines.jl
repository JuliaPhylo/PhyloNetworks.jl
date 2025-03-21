"""
    traversal_preorder(nodes, init_function, root_function, tree_node_function,
                      hybrid_node_function, parameters...)
    traversal_preorder!(nodes, output_array, root_function, tree_node_function,
                       hybrid_node_function, parameters...)
    preorder_nodeupdate!(node_index, nodes, output_array, root_function, tree_node_function,
                   hybrid_node_function, parameters...)

Generic tool to apply a pre-order (or topological ordering) algorithm.
Used by `sharedpathmatrix` and by `pairwisetaxondistancematrix`, for example.
`preorder_nodeupdate!` is a helper that calls the root / tree node / hybrid node
function as appropriate.

output: array object `output_array`.

arguments:
- `nodes`: array of nodes in the network, pre-ordered, typically the internal
  `net.vec_node` after applying `preorder!(net)`
- `output_array`: output object, of type `AbstractArray`, named `V` below
- `init_function`: to initialize the output array, taking `(nodes, parameters...)`
  as arguments
- `root_function`: to do whatever needs to be done to `V` at the root, using
  `(V, rootnode_index, parameters...)` as arguments
- `tree_node_function`: to do whatever needs to be done to `V` at a tree node,
  using `(V, treenode_index, parentnode_index, parentedge, parameters...)`
  as arguments
- `hybrid_node_function`: to do whatever needs to be done to `V` at a hybrid
  node, using arguments
  `(V, hybnode_index, parentnodes_index_vector, parentedges_vector, parameters...)`

The last 3 functions should return a boolean. If true: traversal continues.
If false: traversal is stopped, that is, the next node is not processed.

See also [`traversal_postorder`](@ref PhyloNetworks.traversal_postorder), and
[`traversalupdate_default!`](@ref PhyloNetworks.traversalupdate_default!)
for a default function that does nothing to `V` and returns true,
with an adequate signature to be used here.
"""
function traversal_preorder(
    nodes::Vector{Node},
    init::Function,
    updateRoot::Function,
    updateTree::Function,
    updateHybrid::Function,
    params...
)
    M = init(nodes, params...)
    traversal_preorder!(nodes, M, updateRoot, updateTree, updateHybrid, params...)
end

@doc (@doc traversal_preorder) traversal_preorder!
function traversal_preorder!(
    nodes::Vector{Node},
    M::AbstractArray,
    updateRoot::Function,
    updateTree::Function,
    updateHybrid::Function,
    params...
)
    for i in 1:length(nodes) # nodes should be listed in topological order
        preorder_nodeupdate!(i, nodes, M, updateRoot, updateTree, updateHybrid, params...) || break
    end
    return M
end

@doc (@doc traversal_preorder) preorder_nodeupdate!
function preorder_nodeupdate!(
    i::Int,
    nodes::Vector{Node},
    V::AbstractArray,
    updateRoot::Function,
    updateTree::Function,
    updateHybrid::Function,
    params...
)
    parnode = Node[] # will be the vector of parent nodes (empty, size 1 or 2)
    paredge = Edge[]
    parindx = Int[]  # index(es) of parent node(s) in `nodes`
    keepgoing=true # true/false returned from each recursion function to determine if we can stop
    for e in nodes[i].edge
        ischildof(nodes[i], e) || continue
        push!(paredge, e)
        pn = getparent(e)
        push!(parnode, pn)
        push!(parindx, getIndex(pn, nodes))
    end
    if isempty(parnode) # nodes[i] is root
        keepgoing=updateRoot(V, i, params...)
    elseif length(parnode) == 1 # nodes[i] is a tree node
        keepgoing=updateTree(V, i, parindx[1], paredge[1], params...)
    elseif length(parnode) > 1 # nodes[i] is a hybrid node
        for e in paredge
            e.hybrid || error("edge $(e.number), parent of node $(nodes[i].number), should be hybrid")
        end
        keepgoing = updateHybrid(V, i, parindx, paredge, params...)
    end
    return keepgoing
end

"""
    traversalupdate_default!(::AbstractArray, ::Int, args...)

Returns `true`. With its signature, this function can be used as a default
update function at any node (root/tree/hybrid) in [`traversal_preorder`](@ref PhyloNetworks.traversal_preorder)
or [`traversal_postorder`](@ref PhyloNetworks.traversal_postorder),
if we need a function that does nothing and keeps the traversal going.
"""
function traversalupdate_default!(::AbstractArray, ::Int, params...)
    return true
end

"""
    traversal_postorder(nodes, init_function,
        tip_function, internalnode_function, parameters...)
    postorder_nodeupdate!(node_index, nodes, output_array,
        tip_function, internalnode_function, parameters...)

Generic tool to apply a post-order (or reverse topological ordering) algorithm,
acting on a matrix where rows & columns correspond to nodes.
Used by [`descendencematrix`](@ref).

output: matrix `output_array`.

arguments:
- `nodes`: array of nodes in the network, pre-ordered, typically the internal
  `net.vec_node` after applying `preorder!(net)`.
  This array is traversed in reverse order.
- `output_array`: output object, of type `Matrix`, named `V` below
- `init_function`: to initialize the output array, taking `(nodes, parameters...)`
  as arguments
- `tip_function`: to do whatever needs to be done to `V` at a leaf node,
  using `(V, tip_index, parameters...)` as arguments
- `internalnode_function`: to do whatever needs to be done to `V` at an internal
  node, using arguments
  `(V, node_index, childnodes_index_vector, childedges_vector, parameters...)`

The last 3 functions should return a boolean. If true: traversal continues.
If false: traversal is stopped, that is, the next node is not processed.

See also [`traversal_preorder`](@ref PhyloNetworks.traversal_preorder), and
[`traversalupdate_default!`](@ref PhyloNetworks.traversalupdate_default!)
for a default function that does nothing to `V` and returns true,
with an adequate signature to be used here.
"""
function traversal_postorder(
    nodes::Vector{Node},
    init::Function,
    updatetip::Function,
    updateinternalnode::Function,
    params...
)
    n = length(nodes)
    M = init(nodes, params...)
    for i in n:-1:1 #sorted list of nodes
        postorder_nodeupdate!(i, nodes, M, updatetip, updateinternalnode, params...) || break
    end
    return M
end
@doc (@doc traversal_postorder) postorder_nodeupdate!
function postorder_nodeupdate!(
    i::Int,
    nodes::Vector{Node},
    V::Matrix,
    updatetip::Function,
    updateinternalnode::Function,
    params...
)
    chnode = Node[] # will be the vector of child nodes (empty, size 1 or 2)
    chedge = Edge[]
    chindx = Int[]  # index(es) of child node(s) in `nodes`
    keepgoing=true # true/false returned from each recursion function to determine if we can stop
    for e in nodes[i].edge
        isparentof(nodes[i], e) || continue
        push!(chedge, e)
        cn = getchild(e)
        push!(chnode, cn)
        push!(chindx, getIndex(cn, nodes))
    end
    if isempty(chnode) # nodes[i] is a tip
        keepgoing = updatetip(V, i, params...)
    else
        keepgoing = updateinternalnode(V, i, chindx, chedge, params...)
    end
    return keepgoing
end
