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
    parent = getparents(nodes[i]) # vector of nodes (empty, size 1 or 2)
    if isempty(parent) #nodes[i] is root
        updateRoot(V, i, params...)
    elseif length(parent) == 1 #nodes[i] is tree
        parentIndex = getIndex(parent[1],nodes)
        edge = getConnectingEdge(nodes[i],parent[1])
        updateTree(V, i, parentIndex, edge, params...)
    elseif length(parent) == 2 #nodes[i] is hybrid
        parentIndex1 = getIndex(parent[1],nodes)
        parentIndex2 = getIndex(parent[2],nodes)
        edge1 = getConnectingEdge(nodes[i],parent[1])
        edge2 = getConnectingEdge(nodes[i],parent[2])
        edge1.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[1].number) should be a hybrid egde")
        edge2.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[2].number) should be a hybrid egde")
        updateHybrid(V, i, parentIndex1, parentIndex2, edge1, edge2, params...)
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
    children = getchildren(nodes[i]) # vector of nodes (empty, size 1 or 2)
    if(isempty(children)) #nodes[i] is a tip
        updateTip(V, i, params)
    else
        childrenIndex = [getIndex(n, nodes) for n in children]
        edges = [getConnectingEdge(nodes[i], c) for c in children]
        updateNode(V, i, childrenIndex, edges, params...)
    end
end
