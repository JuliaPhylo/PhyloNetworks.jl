"""
    recursion_preorder(nodes, init_function, root_function, tree_node_function,
                      hybrid_node_function, parameters...)
    recursion_preorder!(nodes, AbstractArray, root_function, tree_node_function,
                       hybrid_node_function, parameters...)
    updatePreOrder(index, nodes, updated_matrix, root_function, tree_node_function,
                   hybrid_node_function, parameters...)

Generic tool to apply a pre-order (or topological ordering) algorithm.
Used by `sharedPathMatrix` and by `pairwiseTaxonDistanceMatrix`.
"""
function recursion_preorder(
    nodes::Vector{Node},
    init::Function,
    updateRoot::Function,
    updateTree::Function,
    updateHybrid::Function,
    params...
)
    M = init(nodes, params...)
    recursion_preorder!(nodes, M, updateRoot, updateTree, updateHybrid, params...)
end
@doc (@doc recursion_preorder) recursion_preorder!
function recursion_preorder!(
    nodes::Vector{Node},
    M::AbstractArray,
    updateRoot::Function,
    updateTree::Function,
    updateHybrid::Function,
    params...
)
    for i in 1:length(nodes) # nodes should be listed in topological order
        updatePreOrder!(i, nodes, M, updateRoot, updateTree, updateHybrid, params...)
    end
    return M
end

@doc (@doc recursion_preorder) updatePreOrder!
function updatePreOrder!(
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
    for e in nodes[i].edge
        ischildof(nodes[i], e) || continue
        push!(paredge, e)
        pn = getparent(e)
        push!(parnode, pn)
        push!(parindx, getIndex(pn, nodes))
    end
    if isempty(parnode) # nodes[i] is root
        updateRoot(V, i, params...)
    elseif length(parnode) == 1 # nodes[i] is a tree node
        updateTree(V, i, parindx[1], paredge[1], params...)
    elseif length(parnode) > 1 # nodes[i] is a hybrid node
        for e in paredge
            e.hybrid || error("edge $(e.number), parent of node $(nodes[i].number), should be hybrid")
        end
        updateHybrid(V, i, parindx, paredge, params...)
    end
end

function updateRecursion_default!(::AbstractArray, ::Int, params...)
    return nothing
end

"""
    recursion_postorder(net::HybridNetwork, checkPreorder::Bool,
                       init_function, tip_function, node_function,
                       indexation="b", parameters...)
    recursion_postorder(nodes, init_function, tip_function, node_function,
                       parameters...)
    updatePostOrder!(index, nodes, updated_matrix, tip_function, node_function,
                    parameters...)

Generic tool to apply a post-order (or topological ordering) algorithm,
acting on a matrix where rows & columns correspond to nodes.
Used by `descendenceMatrix`.
"""
function recursion_postorder(
    nodes::Vector{Node},
    init::Function,
    updateTip::Function,
    updateNode::Function,
    params...
)
    n = length(nodes)
    M = init(nodes, params...)
    for i in n:-1:1 #sorted list of nodes
        updatePostOrder!(i, nodes, M, updateTip, updateNode, params...)
    end
    return M
end
@doc (@doc recursion_postorder) updatePostOrder!
function updatePostOrder!(
    i::Int,
    nodes::Vector{Node},
    V::Matrix,
    updateTip::Function,
    updateNode::Function,
    params...
)
    chnode = Node[] # will be the vector of child nodes (empty, size 1 or 2)
    chedge = Edge[]
    chindx = Int[]  # index(es) of child node(s) in `nodes`
    for e in nodes[i].edge
        isparentof(nodes[i], e) || continue
        push!(chedge, e)
        cn = getchild(e)
        push!(chnode, cn)
        push!(chindx, getIndex(cn, nodes))
    end
    if isempty(chnode) # nodes[i] is a tip
        updateTip(V, i, params...)
    else
        updateNode(V, i, chindx, chedge, params...)
    end
end
