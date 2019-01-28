# functions for trait evolution on network
# Claudia & Paul Bastide: November 2015

###############################################################################
###############################################################################
## Function to traverse the network in the pre-order, updating a matrix
###############################################################################
###############################################################################

# Matrix with rows and/or columns in topological order of the net.
"""
    MatrixTopologicalOrder

Matrix associated to an [`HybridNetwork`](@ref) sorted in topological order.

The following functions and extractors can be applied to it: [`tipLabels`](@ref), `obj[:Tips]`, `obj[:InternalNodes]`, `obj[:TipsNodes]` (see documentation for function [`getindex(::MatrixTopologicalOrder, ::Symbol)`](@ref)).

Functions [`sharedPathMatrix`](@ref) and [`simulate`](@ref) return objects of this type.

The `MatrixTopologicalOrder` object has fields: `V`, `nodeNumbersTopOrder`, `internalNodeNumbers`, `tipNumbers`, `tipNames`, `indexation`.
Type in "?MatrixTopologicalOrder.field" to get documentation on a specific field.
"""
struct MatrixTopologicalOrder
    "V: the matrix per se"
    V::Matrix # Matrix in itself
    "nodeNumbersTopOrder: vector of nodes numbers in the topological order, used for the matrix"
    nodeNumbersTopOrder::Vector{Int} # Vector of nodes numbers for ordering of the matrix
    "internalNodeNumbers: vector of internal nodes number, in the original net order"
    internalNodeNumbers::Vector{Int} # Internal nodes numbers (original net order)
    "tipNumbers: vector of tips numbers, in the origial net order"
    tipNumbers::Vector{Int} # Tips numbers (original net order)
    "tipNames: vector of tips names, in the original net order"
    tipNames::Vector # Tips Names (original net order)
    """
    indexation: a string giving the type of matrix `V`:
    -"r": rows only are indexed by the nodes of the network
    -"c": columns only are indexed by the nodes of the network
    -"b": both rows and columns are indexed by the nodes of the network
    """
    indexation::AbstractString # Are rows ("r"), columns ("c") or both ("b") indexed by nodes numbers in the matrix ?
end

function Base.show(io::IO, obj::MatrixTopologicalOrder)
    println(io, "$(typeof(obj)):\n$(obj.V)")
end

# docstring already in descriptive.jl
function tipLabels(obj::MatrixTopologicalOrder)
    return obj.tipNames
end

# This function takes an init and update funtions as arguments
# It does the recursion using these functions on a preordered network.
function recursionPreOrder(net::HybridNetwork,
                           checkPreorder=true::Bool,
                           init=identity::Function,
                           updateRoot=identity::Function,
                           updateTree=identity::Function,
                           updateHybrid=identity::Function,
                           indexation="b"::AbstractString,
                           params...)
    net.isRooted || error("net needs to be rooted to get matrix of shared path lengths")
    if(checkPreorder)
        preorder!(net)
    end
    M = recursionPreOrder(net.nodes_changed, init, updateRoot, updateTree, updateHybrid, params)
    # Find numbers of internal nodes
    nNodes = [n.number for n in net.node]
    nleaf = [n.number for n in net.leaf]
    deleteat!(nNodes, indexin(nleaf, nNodes))
    MatrixTopologicalOrder(M, [n.number for n in net.nodes_changed], nNodes, nleaf, [n.name for n in net.leaf], indexation)
end

"""
    recursionPreOrder(nodes, init_function, root_function, tree_node_function,
                      hybrid_node_function, parameters)
    recursionPreOrder!(nodes, AbstractArray, root_function, tree_node_function,
                       hybrid_node_function, parameters)
    updatePreOrder(index, nodes, updated_matrix, root_function, tree_node_function,
                   hybrid_node_function, parameters)

Generic tool to apply a pre-order (or topological ordering) algorithm.
Used by `sharedPathMatrix` and by `pairwiseTaxonDistanceMatrix`.
"""
function recursionPreOrder(nodes::Vector{Node},
                           init::Function,
                           updateRoot::Function,
                           updateTree::Function,
                           updateHybrid::Function,
                           params)
    M = init(nodes, params)
    recursionPreOrder!(nodes, M, updateRoot, updateTree, updateHybrid, params)
end
@doc (@doc recursionPreOrder) recursionPreOrder!
function recursionPreOrder!(nodes::Vector{Node},
                           M::AbstractArray,
                           updateRoot::Function,
                           updateTree::Function,
                           updateHybrid::Function,
                           params)
    for i in 1:length(nodes) #sorted list of nodes
        updatePreOrder!(i, nodes, M, updateRoot, updateTree, updateHybrid, params)
    end
    return M
end

# Update on the network
# Takes three function as arguments : updateRoot, updateTree, updateHybrid
@doc (@doc recursionPreOrder) updatePreOrder!
function updatePreOrder!(i::Int,
                         nodes::Vector{Node},
                         V::AbstractArray, updateRoot::Function,
                         updateTree::Function,
                         updateHybrid::Function,
                         params)
    parent = getParents(nodes[i]) #array of nodes (empty, size 1 or 2)
    if(isempty(parent)) #nodes[i] is root
        updateRoot(V, i, params)
    elseif(length(parent) == 1) #nodes[i] is tree
        parentIndex = getIndex(parent[1],nodes)
        edge = getConnectingEdge(nodes[i],parent[1])
        updateTree(V, i, parentIndex, edge, params)
    elseif(length(parent) == 2) #nodes[i] is hybrid
        parentIndex1 = getIndex(parent[1],nodes)
        parentIndex2 = getIndex(parent[2],nodes)
        edge1 = getConnectingEdge(nodes[i],parent[1])
        edge2 = getConnectingEdge(nodes[i],parent[2])
        edge1.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[1].number) should be a hybrid egde")
        edge2.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[2].number) should be a hybrid egde")
        updateHybrid(V, i, parentIndex1, parentIndex2, edge1, edge2, params)
    end
end

## Same, but in post order (tips to root)
function recursionPostOrder(net::HybridNetwork,
                            checkPreorder=true::Bool,
                            init=identity::Function,
                            updateTip=identity::Function,
                            updateNode=identity::Function,
                            indexation="b"::AbstractString,
                            params...)
    net.isRooted || error("net needs to be rooted to get matrix of shared path lengths")
    if(checkPreorder)
        preorder!(net)
    end
    M = recursionPostOrder(net.nodes_changed, init, updateTip, updateNode, params)
    # Find numbers of internal nodes
    nNodes = [n.number for n in net.node]
    nleaf = [n.number for n in net.leaf]
    deleteat!(nNodes, indexin(nleaf, nNodes))
    MatrixTopologicalOrder(M, [n.number for n in net.nodes_changed], nNodes, nleaf, [n.name for n in net.leaf], indexation)
end

"""
    recursionPostOrder(nodes, init_function, tip_function, node_function,
                       parameters)
    updatePostOrder(index, nodes, updated_matrix, tip_function, node_function,
                    parameters)

Generic tool to apply a post-order (or topological ordering) algorithm.
Used by `descendenceMatrix`.
"""
function recursionPostOrder(nodes::Vector{Node},
                            init::Function,
                            updateTip::Function,
                            updateNode::Function,
                            params)
    n = length(nodes)
    M = init(nodes, params)
    for i in n:-1:1 #sorted list of nodes
        updatePostOrder!(i, nodes, M, updateTip, updateNode, params)
    end
    return M
end
@doc (@doc recursionPostOrder) updatePostOrder!
function updatePostOrder!(i::Int,
                          nodes::Vector{Node},
                          V::Matrix,
                          updateTip::Function,
                          updateNode::Function,
                          params)
    children = getChildren(nodes[i]) #array of nodes (empty, size 1 or 2)
    if(isempty(children)) #nodes[i] is a tip
        updateTip(V, i, params)
    else
        childrenIndex = [getIndex(n, nodes) for n in children]
        edges = [getConnectingEdge(nodes[i], c) for c in children]
        updateNode(V, i, childrenIndex, edges, params)
    end
end

# Extract the right part of a matrix in topological order
# Tips : submatrix corresponding to tips
# InternalNodes : submatrix corresponding to internal nodes
# TipsNodes : submatrix nTips x nNodes of interactions
# !! Extract sub-matrices in the original net nodes numbers !!
# function Base.getindex(obj::MatrixTopologicalOrder, d::Symbol)
#   if d == :Tips # Extract rows and/or columns corresponding to the tips
#       mask = indexin(obj.tipNumbers, obj.nodeNumbersTopOrder)
#       obj.indexation == "b" && return obj.V[mask, mask] # both columns and rows are indexed by nodes
#       obj.indexation == "c" && return obj.V[:, mask] # Only the columns
#       obj.indexation == "r" && return obj.V[mask, :] # Only the rows
#   end
#   if d == :InternalNodes # Idem, for internal nodes
#       mask = indexin(obj.internalNodeNumbers, obj.nodeNumbersTopOrder)
#       obj.indexation == "b" && return obj.V[mask, mask]
#       obj.indexation == "c" && return obj.V[:, mask]
#       obj.indexation == "r" && return obj.V[mask, :]
#   end
#   if d == :TipsNodes
#       maskNodes = indexin(obj.internalNodeNumbers, obj.nodeNumbersTopOrder)
#       maskTips = indexin(obj.tipNumbers, obj.nodeNumbersTopOrder,)
#       obj.indexation == "b" && return obj.V[maskTips, maskNodes]
#       obj.indexation == "c" && error("Both rows and columns must be net
#       ordered to take the submatrix tips vs internal nodes.")
#       obj.indexation == "r" && error("Both rows and columns must be net
#       ordered to take the submatrix tips vs internal nodes.")
#   end
#   d == :All && return obj.V
# end

# If some tips are missing, treat them as "internal nodes"
"""
    getindex(obj, d,[ indTips, nonmissing])

Getting submatrices of an object of type [`MatrixTopologicalOrder`](@ref).

# Arguments
* `obj::MatrixTopologicalOrder`: the matrix from which to extract.
* `d::Symbol`: a symbol precising which sub-matrix to extract. Can be:
  * `:Tips` columns and/or rows corresponding to the tips
  * `:InternalNodes` columns and/or rows corresponding to the internal nodes
  * `:TipsNodes` columns corresponding to internal nodes, and row to tips (works only is indexation="b")
* `indTips::Vector{Int}`: optional argument precising a specific order for the tips (internal use).
* `nonmissing::BitArray{1}`: optional argument saying which tips have data (internal use).

"""
function Base.getindex(obj::MatrixTopologicalOrder,
                       d::Symbol,
                       indTips=collect(1:length(obj.tipNumbers))::Vector{Int},
                       nonmissing=trues(length(obj.tipNumbers))::BitArray{1})
    if d == :Tips # Extract rows and/or columns corresponding to the tips with data
        maskTips = indexin(obj.tipNumbers, obj.nodeNumbersTopOrder)
        maskTips = maskTips[indTips]
        maskTips = maskTips[nonmissing]
        obj.indexation == "b" && return obj.V[maskTips, maskTips] # both columns and rows are indexed by nodes
        obj.indexation == "c" && return obj.V[:, maskTips] # Only the columns
        obj.indexation == "r" && return obj.V[maskTips, :] # Only the rows
    end
    if d == :InternalNodes # Idem, for internal nodes
        maskNodes = indexin(obj.internalNodeNumbers, obj.nodeNumbersTopOrder)
        maskTips = indexin(obj.tipNumbers, obj.nodeNumbersTopOrder)
        maskTips = maskTips[indTips]
        maskNodes = [maskNodes; maskTips[.!nonmissing]]
        obj.indexation == "b" && return obj.V[maskNodes, maskNodes]
        obj.indexation == "c" && return obj.V[:, maskNodes]
        obj.indexation == "r" && return obj.V[maskNodes, :]
    end
    if d == :TipsNodes
        maskNodes = indexin(obj.internalNodeNumbers, obj.nodeNumbersTopOrder)
        maskTips = indexin(obj.tipNumbers, obj.nodeNumbersTopOrder)
        maskTips = maskTips[indTips]
        maskNodes = [maskNodes; maskTips[.!nonmissing]]
        maskTips = maskTips[nonmissing]
        obj.indexation == "b" && return obj.V[maskTips, maskNodes]
        obj.indexation == "c" && error("""Both rows and columns must be net
                                       ordered to take the submatrix tips vs internal nodes.""")
        obj.indexation == "r" && error("""Both rows and columns must be net
                                       ordered to take the submatrix tips vs internal nodes.""")
    end
    d == :All && return obj.V
end

###############################################################################
###############################################################################
## Functions to compute the variance-covariance between Node and its parents
###############################################################################
###############################################################################
"""
    vcv(net::HybridNetwork; model="BM"::AbstractString, 
                            corr=false::Bool,
                            checkPreorder=true::Bool)

This function computes the variance covariance matrix between the tips of the
network, assuming a Brownian model of trait evolution (with unit variance).
If optional argument `corr` is set to `true`, then the correlation matrix is returned instead.

The function returns a `DataFrame` object, with columns named by the tips of the network.

The calculation of the covariance matrix requires a pre-ordering of nodes to be fast.
If `checkPreorder` is true (default), then [`preorder!`](@ref) is run on the network beforehand.
Otherwise, the network is assumed to be already in pre-order.

This function internally calls [`sharedPathMatrix`](@ref), that computes the variance
matrix between all the nodes of the network.

# Examples
```jldoctest
julia> tree_str = "(((t2:0.14,t4:0.33):0.59,t3:0.96):0.14,(t5:0.70,t1:0.18):0.90);";

julia> tree = readTopology(tree_str);

julia> C = vcv(tree)
5×5 DataFrames.DataFrame
│ Row │ t2      │ t4      │ t3      │ t5      │ t1      │
│     │ Float64 │ Float64 │ Float64 │ Float64 │ Float64 │
├─────┼─────────┼─────────┼─────────┼─────────┼─────────┤
│ 1   │ 0.87    │ 0.73    │ 0.14    │ 0.0     │ 0.0     │
│ 2   │ 0.73    │ 1.06    │ 0.14    │ 0.0     │ 0.0     │
│ 3   │ 0.14    │ 0.14    │ 1.1     │ 0.0     │ 0.0     │
│ 4   │ 0.0     │ 0.0     │ 0.0     │ 1.6     │ 0.9     │
│ 5   │ 0.0     │ 0.0     │ 0.0     │ 0.9     │ 1.08    │

```
The following block needs `ape` to be installed (not run):
```julia
julia> using RCall # Comparison with ape vcv function

julia> R"ape::vcv(ape::read.tree(text = \$tree_str))"
RCall.RObject{RCall.RealSxp}
     t2   t4   t3  t5   t1
t2 0.87 0.73 0.14 0.0 0.00
t4 0.73 1.06 0.14 0.0 0.00
t3 0.14 0.14 1.10 0.0 0.00
t5 0.00 0.00 0.00 1.6 0.90
t1 0.00 0.00 0.00 0.9 1.08

```

The covariance can also be calculated on a network
(for the model, see for Bastide et al. 2018)
```jldoctest
julia> net = readTopology("((t1:1.0,#H1:0.1::0.30):0.5,((t2:0.9)#H1:0.2::0.70,t3:1.1):0.4);");

julia> C = vcv(net)
3×3 DataFrames.DataFrame
│ Row │ t1      │ t2      │ t3      │
│     │ Float64 │ Float64 │ Float64 │
├─────┼─────────┼─────────┼─────────┤
│ 1   │ 1.5     │ 0.15    │ 0.0     │
│ 2   │ 0.15    │ 1.248   │ 0.28    │
│ 3   │ 0.0     │ 0.28    │ 1.5     │
```
"""
function vcv(net::HybridNetwork;
             model="BM"::AbstractString,
             corr=false::Bool,
             checkPreorder=true::Bool)
    @assert (model == "BM") "The 'vcv' function only works for a BM process (for now)."
    V = sharedPathMatrix(net; checkPreorder=checkPreorder)
    C = V[:Tips]
    corr && StatsBase.cov2cor!(C, sqrt.(LinearAlgebra.diag(C)))
    Cd = convert(DataFrame, C)
    names!(Cd, map(Symbol, V.tipNames))
    return(Cd)
end


"""
    sharedPathMatrix(net::HybridNetwork; checkPreorder=true::Bool)

This function computes the shared path matrix between all the nodes of a
network. It assumes that the network is in the pre-order. If checkPreorder is
true (default), then it runs function `preoder` on the network beforehand.

Returns an object of type [`MatrixTopologicalOrder`](@ref).

"""
function sharedPathMatrix(net::HybridNetwork;
                          checkPreorder=true::Bool)
    recursionPreOrder(net,
                      checkPreorder,
                      initsharedPathMatrix,
                      updateRootSharedPathMatrix!,
                      updateTreeSharedPathMatrix!,
                      updateHybridSharedPathMatrix!,
                      "b")
end

function updateRootSharedPathMatrix!(V::AbstractArray, i::Int, params)
    return
end


function updateTreeSharedPathMatrix!(V::Matrix,
                                     i::Int,
                                     parentIndex::Int,
                                     edge::Edge,
                                     params)
    for j in 1:(i-1)
        V[i,j] = V[j,parentIndex]
        V[j,i] = V[j,parentIndex]
    end
    V[i,i] = V[parentIndex,parentIndex] + edge.length
end

function updateHybridSharedPathMatrix!(V::Matrix,
                                       i::Int,
                                       parentIndex1::Int,
                                       parentIndex2::Int,
                                       edge1::Edge,
                                       edge2::Edge,
                                       params)
    for j in 1:(i-1)
        V[i,j] = V[j,parentIndex1]*edge1.gamma + V[j,parentIndex2]*edge2.gamma
        V[j,i] = V[i,j]
    end
    V[i,i] = edge1.gamma*edge1.gamma*(V[parentIndex1,parentIndex1] + edge1.length) + edge2.gamma*edge2.gamma*(V[parentIndex2,parentIndex2] + edge2.length) + 2*edge1.gamma*edge2.gamma*V[parentIndex1,parentIndex2]
end

function initsharedPathMatrix(nodes::Vector{Node}, params)
    n = length(nodes)
    return(zeros(Float64,n,n))
end

###############################################################################
###############################################################################
## Functions to compute the network matrix
###############################################################################
###############################################################################
"""
    descendenceMatrix(net::HybridNetwork; checkPreorder=true::Bool)

This function computes the inciednce matrix between all the nodes of a
network. It assumes that the network is in the pre-order. If checkPreorder is
true (default), then it runs function `preoder` on the network beforehand.

Returns an object of type [`MatrixTopologicalOrder`](@ref).

"""
function descendenceMatrix(net::HybridNetwork;
                         checkPreorder=true::Bool)
    recursionPostOrder(net,
                       checkPreorder,
                       initDescendenceMatrix,
                       updateTipDescendenceMatrix!,
                       updateNodeDescendenceMatrix!,
                       "r")
end

function updateTipDescendenceMatrix!(V::Matrix,
                                   i::Int,
                                   params)
    return
end

function updateNodeDescendenceMatrix!(V::Matrix,
                                    i::Int,
                                    childrenIndex::Vector{Int},
                                    edges::Vector{Edge},
                                    params)
    for j in 1:length(edges)
        V[:,i] .+= edges[j].gamma .* V[:,childrenIndex[j]]
    end
end

function initDescendenceMatrix(nodes::Vector{Node}, params)
    n = length(nodes)
    return(Matrix{Float64}(LinearAlgebra.I, n, n)) # identity matrix
end

###############################################################################
###############################################################################
## Function to get the regressor out of a shift
###############################################################################
###############################################################################
"""
    regressorShift(node::Vector{Node}, net::HybridNetwork; checkPreorder=true::Bool)

`regressorShift(edge::Vector{Edge}, net::HybridNetwork; checkPreorder=true::Bool)`

Compute the regressor vectors associated with shifts on edges that are above nodes
`node`, or on edges `edge`, on a network `net`. It uses function [`descendenceMatrix`](@ref), so
`net` might be modified to sort it in a pre-order.
Return a `DataFrame` with as many rows as there are tips in net, and a column for
each shift, each labelled according to the pattern shift_{number_of_edge}. It has
an aditional column labelled `tipNames` to allow easy fitting afterward (see example).

# Examples
```jldoctest
julia> net = readTopology("(A:2.5,((B:1,#H1:0.5::0.4):1,(C:1,(D:0.5)#H1:0.5::0.6):1):0.5);");

julia> preorder!(net)

julia> using PhyloPlots

julia> plot(net, :RCall, showNodeNumber=true); # to locate nodes

julia> nodes_shifts = indexin([1,-5], [n.number for n in net.node]) # Put a shift on edges ending at nodes 1 and -5
2-element Array{Union{Nothing, Int64},1}:
 1
 7

julia> params = ParamsBM(10, 0.1, ShiftNet(net.node[nodes_shifts], [3.0, -3.0],  net))
ParamsBM:
Parameters of a BM with fixed root:
mu: 10
Sigma2: 0.1

There are 2 shifts on the network:
──────────────────────────
  Edge Number  Shift Value
──────────────────────────
          8.0         -3.0
          1.0          3.0
──────────────────────────

julia> using Random; Random.seed!(2468); # sets the seed for reproducibility

julia> sim = simulate(net, params); # simulate a dataset with shifts

julia> using DataFrames # to handle data frames

julia> dat = DataFrame(trait = sim[:Tips], tipNames = sim.M.tipNames)
4×2 DataFrames.DataFrame
│ Row │ trait   │ tipNames │
│     │ Float64 │ String   │
├─────┼─────────┼──────────┤
│ 1   │ 13.392  │ A        │
│ 2   │ 9.55741 │ B        │
│ 3   │ 7.17704 │ C        │
│ 4   │ 7.88906 │ D        │

julia> dfr_shift = regressorShift(net.node[nodes_shifts], net) # the regressors matching the shifts.
4×3 DataFrames.DataFrame
│ Row │ shift_1 │ shift_8 │ tipNames │
│     │ Float64 │ Float64 │ String   │
├─────┼─────────┼─────────┼──────────┤
│ 1   │ 1.0     │ 0.0     │ A        │
│ 2   │ 0.0     │ 0.0     │ B        │
│ 3   │ 0.0     │ 1.0     │ C        │
│ 4   │ 0.0     │ 0.6     │ D        │

julia> dfr = join(dat, dfr_shift, on=:tipNames); # join data and regressors in a single dataframe

julia> using StatsModels # for statistical model formulas

julia> fitBM = phyloNetworklm(@formula(trait ~ shift_1 + shift_8), dfr, net) # actual fit
StatsModels.DataFrameRegressionModel{PhyloNetworkLinearModel,Array{Float64,2}}

Formula: trait ~ 1 + shift_1 + shift_8

Model: BM

Parameter(s) Estimates:
Sigma2: 0.0112618

Coefficients:
────────────────────────────────────────────────────
             Estimate  Std.Error   t value  Pr(>|t|)
────────────────────────────────────────────────────
(Intercept)   9.48238   0.327089  28.9902     0.0220
shift_1       3.9096    0.46862    8.34279    0.0759
shift_8      -2.4179    0.422825  -5.71843    0.1102
────────────────────────────────────────────────────
Log Likelihood: 1.8937302027
AIC: 4.2125395947

```

# See also
[`phyloNetworklm`](@ref), [`descendenceMatrix`](@ref), [`regressorHybrid`](@ref).
"""
function regressorShift(node::Vector{Node},
                        net::HybridNetwork; checkPreorder=true::Bool)
    T = descendenceMatrix(net; checkPreorder=checkPreorder)
    regressorShift(node, net, T)
end

function regressorShift(node::Vector{Node},
                        net::HybridNetwork,
                        T::MatrixTopologicalOrder)
    ## Get the descendence matrix for tips
    T_t = T[:Tips]
    ## Get the indices of the columns to keep
    ind = zeros(Int, length(node))
    for i in 1:length(node)
        !node[i].hybrid || error("Shifts on hybrid edges are not allowed")
        ind[i] = getIndex(node[i], net.nodes_changed)
    end
    df = DataFrame(T_t[:, ind])
    ## Get the names of the columns
    eNum = [getMajorParentEdgeNumber(n) for n in net.nodes_changed[ind]]
    # function tmp_fun(x::Int)
    #     if x<0
    #         return(Symbol("shift_m$(-x)"))
    #     else
    #         return(Symbol("shift_$(x)"))
    #     end
    # end
    function tmp_fun(x::Int)
        return(Symbol("shift_$(x)"))
    end
    names!(df, [tmp_fun(num) for num in eNum])
    df[:tipNames]=T.tipNames
    return(df)
end

function regressorShift(edge::Vector{Edge},
                        net::HybridNetwork; checkPreorder=true::Bool)
    childs = [getChild(ee) for ee in edge]
    return(regressorShift(childs, net; checkPreorder=checkPreorder))
end

regressorShift(edge::Edge, net::HybridNetwork; checkPreorder=true::Bool) = regressorShift([edge], net; checkPreorder=checkPreorder)
regressorShift(node::Node, net::HybridNetwork; checkPreorder=true::Bool) = regressorShift([node], net; checkPreorder=checkPreorder)

"""
    regressorHybrid(net::HybridNetwork; checkPreorder=true::Bool)

Compute the regressor vectors associated with shifts on edges that imediatly below
all hybrid nodes of `net`. It uses function [`descendenceMatrix`](@ref) through
a call to [`regressorShift`](@ref), so `net` might be modified to sort it in a pre-order.
Return a `DataFrame` with as many rows as there are tips in net, and a column for
each hybrid, each labelled according to the pattern shift_{number_of_edge}. It has
an aditional column labelled `tipNames` to allow easy fitting afterward (see example).

This function can be used to test for heterosis.

# Examples
```jldoctest
julia> using DataFrames # Needed to handle data frames.

julia> net = readTopology("(A:2.5,((B:1,#H1:0.5::0.4):1,(C:1,(D:0.5)#H1:0.5::0.6):1):0.5);");

julia> preorder!(net)

julia> using PhyloPlots

julia> plot(net, :RCall, showNodeNumber=true); # to locate nodes: node 5 is child of hybrid node

julia> nodes_hybrids = indexin([5], [n.number for n in net.node]) # Put a shift on edges below hybrids
1-element Array{Union{Nothing, Int64},1}:
 5

julia> params = ParamsBM(10, 0.1, ShiftNet(net.node[nodes_hybrids], [3.0],  net))
ParamsBM:
Parameters of a BM with fixed root:
mu: 10
Sigma2: 0.1

There are 1 shifts on the network:
──────────────────────────
  Edge Number  Shift Value
──────────────────────────
          6.0          3.0
──────────────────────────


julia> using Random; Random.seed!(2468); # sets the seed for reproducibility

julia> sim = simulate(net, params); # simulate a dataset with shifts

julia> dat = DataFrame(trait = sim[:Tips], tipNames = sim.M.tipNames)
4×2 DataFrames.DataFrame
│ Row │ trait   │ tipNames │
│     │ Float64 │ String   │
├─────┼─────────┼──────────┤
│ 1   │ 10.392  │ A        │
│ 2   │ 9.55741 │ B        │
│ 3   │ 10.177  │ C        │
│ 4   │ 12.6891 │ D        │

julia> dfr_hybrid = regressorHybrid(net) # the reressors matching the hybrids.
4×3 DataFrames.DataFrame
│ Row │ shift_6 │ tipNames │ sum     │
│     │ Float64 │ String   │ Float64 │
├─────┼─────────┼──────────┼─────────┤
│ 1   │ 0.0     │ A        │ 0.0     │
│ 2   │ 0.0     │ B        │ 0.0     │
│ 3   │ 0.0     │ C        │ 0.0     │
│ 4   │ 1.0     │ D        │ 1.0     │

julia> dfr = join(dat, dfr_hybrid, on=:tipNames); # join data and regressors in a single dataframe

julia> using StatsModels

julia> fitBM = phyloNetworklm(@formula(trait ~ shift_6), dfr, net) # actual fit
StatsModels.DataFrameRegressionModel{PhyloNetworkLinearModel,Array{Float64,2}}

Formula: trait ~ 1 + shift_6

Model: BM

Parameter(s) Estimates:
Sigma2: 0.041206

Coefficients:
────────────────────────────────────────────────────
             Estimate  Std.Error   t value  Pr(>|t|)
────────────────────────────────────────────────────
(Intercept)  10.064     0.277959  36.2068     0.0008
shift_6       2.72526   0.315456   8.63912    0.0131
────────────────────────────────────────────────────
Log Likelihood: -0.7006021946
AIC: 7.4012043891

```

# See also
[`phyloNetworklm`](@ref), [`descendenceMatrix`](@ref), [`regressorShift`](@ref).
"""
function regressorHybrid(net::HybridNetwork; checkPreorder=true::Bool)
    childs = [getChildren(nn)[1] for nn in net.hybrid]
    dfr = regressorShift(childs, net; checkPreorder=checkPreorder)
    dfr[:sum] = vec(sum(Matrix(dfr[findall(names(dfr) .!= :tipNames)]), dims=2))
    return(dfr)
end

###############################################################################
###############################################################################
## Types for params process
###############################################################################
###############################################################################

# Abstract type of all the (future) types (BM, OU, ...)
abstract type ParamsProcess end

# Type for shifts
"""
    ShiftNet

Shifts associated to a [`HybridNetwork`](@ref) sorted in topological order.
Its `shift` field is a vector of shift values, one for each node,
corresponding to the shift on the parent edge of the node
(which makes sense for tree nodes only: they have a single parent edge).

Two `ShiftNet` objects on the same network can be concatened with `*`.

`ShiftNet(node::Vector{Node}, value::AbstractVector, net::HybridNetwork; checkPreorder=true::Bool)`

Constructor from a vector of nodes and associated values. The shifts are located
on the edges above the nodes provided. Warning, shifts on hybrid edges are not
allowed.

`ShiftNet(edge::Vector{Edge}, value::AbstractVector, net::HybridNetwork; checkPreorder=true::Bool)`

Constructor from a vector of edges and associated values.
Warning, shifts on hybrid edges are not allowed.

Extractors: [`getShiftEdgeNumber`](@ref), [`getShiftValue`](@ref)
"""
struct ShiftNet
    shift::Vector{Float64}
    net::HybridNetwork
end

# Default
ShiftNet(net::HybridNetwork) = ShiftNet(zeros(length(net.node)), net)

function ShiftNet(node::Vector{Node}, value::AbstractVector,
                  net::HybridNetwork; checkPreorder=true::Bool)
    if length(node) != length(value)
        error("The vector of nodes/edges and of values must be of the same length.")
    end
    if checkPreorder
        preorder!(net)
    end
    obj = ShiftNet(net)
    for i in 1:length(node)
        !node[i].hybrid || error("Shifts on hybrid edges are not allowed")
        ind = findfirst(x -> x===node[i], net.nodes_changed)
        obj.shift[ind] = value[i]
    end
    return(obj)
end

# Construct from edges and values
function ShiftNet(edge::Vector{Edge}, value::AbstractVector,
                  net::HybridNetwork; checkPreorder=true::Bool)
    childs = [getChild(ee) for ee in edge]
    return(ShiftNet(childs, value, net; checkPreorder=checkPreorder))
end

ShiftNet(edge::Edge, value::Float64, net::HybridNetwork; checkPreorder=true::Bool) = ShiftNet([edge], [value], net; checkPreorder=checkPreorder)
ShiftNet(node::Node, value::Float64, net::HybridNetwork; checkPreorder=true::Bool) = ShiftNet([node], [value], net; checkPreorder=checkPreorder)

"""
    shiftHybrid(value::Vector{T} where T<:Real, net::HybridNetwork; checkPreorder=true::Bool)

Construct an object [`ShiftNet`](@ref) with shifts on all the edges below
hybrid nodes, with values provided. The vector of values must have the
same length as the number of hybrids in the network.

"""
function shiftHybrid(value::Vector{T} where T<:Real,
                     net::HybridNetwork; checkPreorder=true::Bool)
    if length(net.hybrid) != length(value)
        error("You must provide as many values as the number of hybrid nodes.")
    end
    childs = [getChildren(nn)[1] for nn in net.hybrid]
    return(ShiftNet(childs, value, net; checkPreorder=checkPreorder))
end
shiftHybrid(value::Real, net::HybridNetwork; checkPreorder=true::Bool) = shiftHybrid([value], net; checkPreorder=checkPreorder)

"""
    getShiftEdgeNumber(shift::ShiftNet)

Get the edge numbers where the shifts are located, for an object [`ShiftNet`](@ref).
"""
function getShiftEdgeNumber(shift::ShiftNet)
    nodInd = findall(!iszero, shift.shift)
    [getMajorParentEdgeNumber(n) for n in shift.net.nodes_changed[nodInd]]
end
function getMajorParentEdgeNumber(n::Node)
    try
        getMajorParentEdge(n).number
    catch
        -1
    end
end
"""
    getShiftValue(shift::ShiftNet)

Get the values of the shifts, for an object [`ShiftNet`](@ref).
"""
function getShiftValue(shift::ShiftNet)
    shift.shift[shift.shift .!= 0]
end

function shiftTable(shift::ShiftNet)
    sv = getShiftValue(shift)
    CoefTable(hcat(getShiftEdgeNumber(shift), sv),
              ["Edge Number", "Shift Value"],
              fill("", length(sv)))
end

function Base.show(io::IO, obj::ShiftNet)
    println(io, "$(typeof(obj)):\n",
            shiftTable(obj))
end

function Base.:*(sh1::ShiftNet, sh2::ShiftNet)
    isEqual(sh1.net, sh2.net) || error("Shifts to be concatenated must be defined on the same network.")
    length(sh1.shift) == length(sh2.shift) || error("Shifts to be concatenated must have the same length.")
    shiftNew = zeros(length(sh1.shift))
    for i in 1:length(sh1.shift)
        if iszero(sh1.shift[i])
            shiftNew[i] = sh2.shift[i]
        elseif iszero(sh2.shift[i])
            shiftNew[i] = sh1.shift[i]
        elseif sh1.shift[i] == sh2.shift[i]
            shiftNew[i] = sh1.shift[i]
        else
            error("The two shifts vectors you provided affect the same edges, so I cannot choose which one you want.")
        end
    end
    return(ShiftNet(shiftNew, sh1.net))
end

# function Base.:(==)(sh1::ShiftNet, sh2::ShiftNet)
#     isEqual(sh1.net, sh2.net) || return(false)
#     sh1.shift == sh2.shift || return(false)
#     return(true)
# end

"""
    ParamsBM <: ParamsProcess

Type for a BM process on a network. Fields are `mu` (expectation),
`sigma2` (variance), `randomRoot` (whether the root is random, default to `false`),
and `varRoot` (if the root is random, the variance of the root, defalut to `NaN`).

"""
mutable struct ParamsBM <: ParamsProcess
    mu::Real # Ancestral value or mean
    sigma2::Real # variance
    randomRoot::Bool # Root is random ? default false
    varRoot::Real # root variance. Default NaN
    shift::Union{ShiftNet, Missing} # shifts
end
# Constructor
ParamsBM(mu::Real, sigma2::Real) = ParamsBM(mu, sigma2, false, NaN, missing) # default values
ParamsBM(mu::Real, sigma2::Real, net::HybridNetwork) = ParamsBM(mu, sigma2, false, NaN, ShiftNet(net)) # default values
ParamsBM(mu::Real, sigma2::Real, shift::ShiftNet) = ParamsBM(mu, sigma2, false, NaN, shift) # default values

function anyShift(params::ParamsBM)
    if ismissing(params.shift) return(false) end
    for v in params.shift.shift
        if v != 0 return(true) end
    end
    return(false)
end

function Base.show(io::IO, obj::ParamsBM)
    disp =  "$(typeof(obj)):\n"
    pt = paramstable(obj)
    if obj.randomRoot
        disp = disp * "Parameters of a BM with random root:\n" * pt
    else
        disp = disp * "Parameters of a BM with fixed root:\n" * pt
    end
    println(io, disp)
end

function paramstable(obj::ParamsBM)
    disp = "mu: $(obj.mu)\nSigma2: $(obj.sigma2)"
    if obj.randomRoot
        disp = disp * "\nvarRoot: $(obj.varRoot)"
    end
    if anyShift(obj)
        disp = disp * "\n\nThere are $(length(getShiftValue(obj.shift))) shifts on the network:\n"
        disp = disp * "$(shiftTable(obj.shift))"
    end
    return(disp)
end


###############################################################################
###############################################################################
## Simulation Function
###############################################################################
###############################################################################

"""
    TraitSimulation

Result of a trait simulation on an [`HybridNetwork`](@ref) with function [`simulate`](@ref).

The following functions and extractors can be applied to it: [`tipLabels`](@ref), `obj[:Tips]`, `obj[:InternalNodes]` (see documentation for function [`getindex(::TraitSimulation, ::Symbol)`](@ref)).

The `TraitSimulation` object has fields: `M`, `params`, `model`.
"""
struct TraitSimulation
    M::MatrixTopologicalOrder
    params::ParamsProcess
    model::AbstractString
end

function Base.show(io::IO, obj::TraitSimulation)
    disp = "$(typeof(obj)):\n"
    disp = disp * "Trait simulation results on a network with $(length(obj.M.tipNames)) tips, using a $(obj.model) model, with parameters:\n"
    disp = disp * paramstable(obj.params)
    println(io, disp)
end

# docstring already in descriptive.jl
function tipLabels(obj::TraitSimulation)
    return tipLabels(obj.M)
end


"""
    simulate(net::HybridNetwork, params::ParamsProcess, checkPreorder=true::Bool)

Simulate traits on `net` using the parameters `params`. For now, only
parameters of type [`ParamsBM`](@ref) (Brownian Motion) are accepted.

The simulation using a recursion from the root to the tips of the network,
therefore, a pre-ordering of nodes is needed. If `checkPreorder=true` (default),
[`preorder!`](@ref) is called on the network beforehand. Otherwise, it is assumed
that the preordering has already been calculated.

Returns an object of type [`TraitSimulation`](@ref),
which has a matrix with two rows:
row 1 for the trait expectations at all the nodes, and
row 2 for the actual simulated trait values at all the nodes.

# Examples
```jldoctest
julia> phy = readTopology(joinpath(dirname(pathof(PhyloNetworks)), "..", "examples", "carnivores_tree.txt"));

julia> par = ParamsBM(1, 0.1) # BM with expectation 1 and variance 0.1.
ParamsBM:
Parameters of a BM with fixed root:
mu: 1
Sigma2: 0.1


julia> using Random; Random.seed!(17920921); # for reproducibility

julia> sim = simulate(phy, par) # Simulate on the tree.
TraitSimulation:
Trait simulation results on a network with 16 tips, using a BM model, with parameters:
mu: 1
Sigma2: 0.1


julia> traits = sim[:Tips] # Extract simulated values at the tips.
16-element Array{Float64,1}:
  2.17618427971927   
  1.0330846124205684 
  3.048979175536912  
  3.0379560744947876 
  2.189704751299587  
  4.031588898597555  
  4.647725850651446  
 -0.8772851731182523 
  4.625121065244063  
 -0.5111667949991542 
  1.3560351170535228 
 -0.10311152349323893
 -2.088472913751017  
  2.6399137689702723 
  2.8051193818084057 
  3.1910928691142915 
```
"""
function simulate(net::HybridNetwork,
                  params::ParamsProcess,
                  checkPreorder=true::Bool)
    if isa(params, ParamsBM)
        model = "BM"
    else
        error("The 'simulate' function only works for a BM process (for now).")
    end
    !ismissing(params.shift) || (params.shift = ShiftNet(net))
    M = recursionPreOrder(net,
                          checkPreorder,
                          initSimulateBM,
                          updateRootSimulateBM!,
                          updateTreeSimulateBM!,
                          updateHybridSimulateBM!,
                          "c",
                          params)
    TraitSimulation(M, params, model)
end

# Initialization of the structure
function initSimulateBM(nodes::Vector{Node}, params::Tuple{ParamsBM})
    return(zeros(2, length(nodes)))
end

# Initialization of the root
function updateRootSimulateBM!(M::Matrix, i::Int, params::Tuple{ParamsBM})
    params = params[1]
    if (params.randomRoot)
        M[1, i] = params.mu # expectation
        M[2, i] = params.mu + sqrt(params.varRoot) * randn() # random value
    else
        M[1, i] = params.mu # expectation
        M[2, i] = params.mu # random value (root fixed)
    end
end

# Going down to a tree node
function updateTreeSimulateBM!(M::Matrix,
                               i::Int,
                               parentIndex::Int,
                               edge::Edge,
                               params::Tuple{ParamsBM})
    params = params[1]
    M[1, i] = M[1, parentIndex] + params.shift.shift[i] # expectation
    M[2, i] = M[2, parentIndex] + params.shift.shift[i] + sqrt(params.sigma2 * edge.length) * randn() # random value
end

# Going down to an hybrid node
function updateHybridSimulateBM!(M::Matrix,
                                 i::Int,
                                 parentIndex1::Int,
                                 parentIndex2::Int,
                                 edge1::Edge,
                                 edge2::Edge,
                                 params::Tuple{ParamsBM})
    params = params[1]
    M[1, i] =  edge1.gamma * M[1, parentIndex1] + edge2.gamma * M[1, parentIndex2] # expectation
    M[2, i] =  edge1.gamma * (M[2, parentIndex1] + sqrt(params.sigma2 * edge1.length) * randn()) + edge2.gamma * (M[2, parentIndex2] + sqrt(params.sigma2 * edge2.length) * randn()) # random value
end


# function updateSimulateBM!(i::Int, nodes::Vector{Node}, M::Matrix, params::Tuple{ParamsBM})
#     params = params[1]
#     parent = getParents(nodes[i]) #array of nodes (empty, size 1 or 2)
#     if(isempty(parent)) #nodes[i] is root
#         if (params.randomRoot)
#       M[1, i] = params.mu # expectation
#       M[2, i] = params.mu + sqrt(params.varRoot) * randn() # random value
#   else
#       M[1, i] = params.mu # expectation
#       M[2, i] = params.mu # random value (root fixed)
#   end
#
#     elseif(length(parent) == 1) #nodes[i] is tree
#         parentIndex = getIndex(parent[1],nodes)
#   l = getConnectingEdge(nodes[i],parent[1]).length
#   M[1, i] = params.mu  # expectation
#   M[2, i] = M[2, parentIndex] + sqrt(params.sigma2 * l) * randn() # random value
#
#     elseif(length(parent) == 2) #nodes[i] is hybrid
#         parentIndex1 = getIndex(parent[1],nodes)
#         parentIndex2 = getIndex(parent[2],nodes)
#         edge1 = getConnectingEdge(nodes[i],parent[1])
#         edge2 = getConnectingEdge(nodes[i],parent[2])
#         edge1.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[1].number) should be a hybrid egde")
#         edge2.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[2].number) should be a hybrid egde")
#   M[1, i] = params.mu  # expectation
#   M[2, i] =  edge1.gamma * (M[2, parentIndex1] + sqrt(params.sigma2 * edge1.length) * randn()) + edge2.gamma * (M[2, parentIndex2] + sqrt(params.sigma2 * edge2.length) * randn()) # random value
#     end
# end

# Extract the vector of simulated values at the tips
"""
    getindex(obj, d)

Getting submatrices of an object of type [`TraitSimulation`](@ref).

# Arguments
* `obj::TraitSimulation`: the matrix from which to extract.
* `d::Symbol`: a symbol precising which sub-matrix to extract. Can be:
  * `:Tips` columns and/or rows corresponding to the tips
  * `:InternalNodes` columns and/or rows corresponding to the internal nodes
"""
function Base.getindex(obj::TraitSimulation, d::Symbol, w=:Sim::Symbol)
     if w == :Exp
        return(getindex(obj.M, d)[1, :])
     end
    getindex(obj.M, d)[2, :]
end

# function extractSimulateTips(sim::Matrix, net::HybridNetwork)
#   mask = getTipsIndexes(net)
#   return(squeeze(sim[2, mask], 1))
# end

###############################################################################
###############################################################################
## Functions for Phylgenetic Network regression
###############################################################################
###############################################################################

# New type for phyloNetwork regression
"""
    PhyloNetworkLinearModel<:LinPredModel

Regression object for a phylogenetic regression. Result of fitting function [`phyloNetworklm`](@ref).
Dominated by the `LinPredModel` class, from package `GLM`.

The following StatsBase functions can be applied to it:
`coef`, `nobs`, `vcov`, `stderror`, `confint`, `coeftable`, `dof_residual`, `dof`, `deviance`,
`residuals`, `response`, `predict`, `loglikelihood`, `nulldeviance`, `nullloglikelihood`,
`r2`, `adjr2`, `aic`, `aicc`, `bic`.

The following StatsModels functions can also be applied to it:
`ModelFrame`, `ModelMatrix`, `Formula`.

Estimated variance and mean of the BM process used can be retrieved with
functions [`sigma2_estim`](@ref) and [`mu_estim`](@ref).

If a Pagel's lambda model is fitted, the parameter can be retrieved with function
[`lambda_estim`](@ref).

An ancestral state reconstruction can be performed from this fitted object using function:
[`ancestralStateReconstruction`](@ref).

The `PhyloNetworkLinearModel` object has fields: `lm`, `V`, `Vy`, `RL`, `Y`, `X`, `logdetVy`, `ind`, `nonmissing`, `model`, `lambda`.
Type in "?PhyloNetworkLinearModel.field" to get help on a specific field.
"""
mutable struct PhyloNetworkLinearModel <: LinPredModel
    "lm: a GLM.LinearModel object, fitted on the cholesky-tranformend problem"
    lm::GLM.LinearModel # result of a lm on a matrix
    "V: a MatrixTopologicalOrder object of the network-induced correlations"
    V::MatrixTopologicalOrder
    "Vy: the sub matrix corresponding to the tips and actually used for the correction"
    Vy::Matrix
    "RL: a LowerTriangular matrix, Cholesky transform of Vy=RL*RL'"
    RL::LowerTriangular
    "Y: the vector of data"
    Y::Vector
    "X: the matrix of regressors"
    X::Matrix
    "logdetVy: the log-determinent of Vy"
    logdetVy::Float64
    "ind: vector matching the tips of the network against the names of the dataframe provided. 0 if the match could not be performed."
    ind::Vector{Int}
    "nonmissing: vector indicating which tips have non-missing data"
    nonmissing::BitArray{1}
    "model: the model used for the fit"
    model::String
    "If applicable, value of lambda (default to 1)."
    lambda::Float64
end

PhyloNetworkLinearModel(lm_fit, V, Vy, RL, Y, X, logdetVy, ind, nonmissing, model) =
  PhyloNetworkLinearModel(lm_fit,V,Vy, RL, Y, X, logdetVy, ind, nonmissing, model, 1.0)

# Function for lm with net residuals
function phyloNetworklm(X::Matrix,
                        Y::Vector,
                        net::HybridNetwork;
                        nonmissing=trues(length(Y))::BitArray{1},
                        model="BM"::AbstractString,
                        ind=[0]::Vector{Int},
                        startingValue=0.5::Real,
                        fixedValue=missing::Union{Real,Missing})
    ## Choose Model
    if model == "BM"
        # Geting variance covariance
        V = sharedPathMatrix(net)
        # Fit
        return phyloNetworklm_BM(X, Y, V;
                                 nonmissing=nonmissing, ind=ind)
    end
    if model == "lambda"
        # Geting variance covariance
        V = sharedPathMatrix(net)
        # Get gammas and heights
        gammas = getGammas(net)
        times = getHeights(net)
        # Fit
        return phyloNetworklm_lambda(X, Y, V, gammas, times;
                                     nonmissing=nonmissing, ind=ind,
                                     startingValue=startingValue, fixedValue=fixedValue)
    end
    if (model == "scalingHybrid")
        # Get gammas
        preorder!(net)
        gammas = getGammas(net)
        # Fit
        return phyloNetworklm_scalingHybrid(X, Y, net, gammas;
                                            nonmissing=nonmissing, ind=ind,
                                            startingValue=startingValue, fixedValue=fixedValue)
    end
end

###############################################################################
## Fit BM

function phyloNetworklm_BM(X::Matrix,
                           Y::Vector,
                           V::MatrixTopologicalOrder;
                           nonmissing=trues(length(Y))::BitArray{1}, # Which tips are not missing ?
                           ind=[0]::Vector{Int},
                           model="BM"::AbstractString,
                           lambda=1.0::Real)
    # Extract tips matrix
    Vy = V[:Tips]
    # Re-order if necessary
    if (ind != [0]) Vy = Vy[ind, ind] end
    # Keep only not missing values
    Vy = Vy[nonmissing, nonmissing]
    # Cholesky decomposition
    R = cholesky(Vy)
    RL = R.L
    # Fit
    PhyloNetworkLinearModel(lm(RL\X, RL\Y), V, Vy, RL, Y, X, LinearAlgebra.logdet(Vy), ind, nonmissing, model, lambda)
end

###############################################################################
## Fit Pagel's Lambda

"""
    getGammas(net)

Get inheritance γ's of major hybrid edges. Assume pre-order calculated already
(with up-to-date field `nodes_changed`). See [`setGammas!`](@ref)
"""
function getGammas(net::HybridNetwork)
    isHybrid = [n.hybrid for n in net.nodes_changed]
    gammas = ones(size(isHybrid))
    for i in 1:size(isHybrid, 1)
        if isHybrid[i]
            majorHybrid = [n.hybrid & n.isMajor for n in net.nodes_changed[i].edge]
            gammas[i] = net.nodes_changed[i].edge[majorHybrid][1].gamma
        end
    end
    return gammas
end

"""
    setGammas!(net, γ vector)

Set inheritance γ's of hybrid edges, using input vector for *major* edges.
Assume pre-order calculated already, with up-to-date field `nodes_changed`.
See [`getGammas`](@ref).

Very different from [`setGamma!`](@ref), which focuses on a single hybrid event,
updates the field `isMajor` according to the new γ, and is not used here.

May assume a tree-child network.
"""
function setGammas!(net::HybridNetwork, gammas::Vector)
    isHybrid = [n.hybrid for n in net.nodes_changed]
    for i in 1:size(isHybrid, 1)
        if isHybrid[i]
            nod = net.nodes_changed[i]
            majorHybrid = [edg.hybrid &  edg.isMajor for edg in nod.edge]
            # worry: assume tree-child network? getMajorParent and getMinorParent would be safer
            minorHybrid = [edg.hybrid & !edg.isMajor for edg in nod.edge]
            nod.edge[majorHybrid][1].gamma = gammas[i]
            if any(minorHybrid) # case where gamma = 0.5 exactly
                nod.edge[minorHybrid][1].gamma = 1 - gammas[i]
            else
                nod.edge[majorHybrid][2].gamma = 1 - gammas[i]
            end
        end
    end
    return nothing
end

"""
    getHeights(net)

Return the height (distance to the root) of all nodes, assuming a time-consistent network
(where all paths from the root to a given hybrid node have the same length).
Also assumes that the network has been preordered, because it uses
[`getGammas`](@ref) and [`setGammas!`](@ref)).
"""
function getHeights(net::HybridNetwork)
    gammas = getGammas(net)
    setGammas!(net, ones(net.numNodes))
    V = sharedPathMatrix(net)
    setGammas!(net, gammas)
    return(LinearAlgebra.diag(V[:All]))
end

function maxLambda(times::Vector, V::MatrixTopologicalOrder)
    maskTips = indexin(V.tipNumbers, V.nodeNumbersTopOrder)
    maskNodes = indexin(V.internalNodeNumbers, V.nodeNumbersTopOrder)
    return minimum(times[maskTips]) / maximum(times[maskNodes])
    # res = minimum(times[maskTips]) / maximum(times[maskNodes])
    # res = res * (1 - 1/5/maximum(times[maskTips]))
end

function transform_matrix_lambda!(V::MatrixTopologicalOrder, lam::AbstractFloat,
                                  gammas::Vector, times::Vector)
    for i in 1:size(V.V, 1)
        for j in 1:size(V.V, 2)
            V.V[i,j] *= lam
        end
    end
    maskTips = indexin(V.tipNumbers, V.nodeNumbersTopOrder)
    for i in maskTips
        V.V[i, i] += (1-lam) * (gammas[i]^2 + (1-gammas[i])^2) * times[i]
    end
    #   V_diag = Matrix(Diagonal(diag(V.V)))
    #   V.V = lam * V.V .+ (1 - lam) .* V_diag
end

function logLik_lam(lam::AbstractFloat,
                    X::Matrix, Y::Vector,
                    V::MatrixTopologicalOrder,
                    gammas::Vector, times::Vector;
                    nonmissing=trues(length(Y))::BitArray{1}, # Which tips are not missing ?
                    ind=[0]::Vector{Int})
    # Transform V according to lambda
    Vp = deepcopy(V)
    transform_matrix_lambda!(Vp, lam, gammas, times)
    # Fit and take likelihood
    fit_lam = phyloNetworklm_BM(X, Y, Vp; nonmissing=nonmissing, ind=ind)
    res = - loglikelihood(fit_lam)
    # Go back to original V
    # transform_matrix_lambda!(V, 1/lam, gammas, times)
    return res
end

# Code for optim taken from PhyloNetworks.jl/src/optimization.jl, lines 276 - 331
const fAbsTr = 1e-10
const fRelTr = 1e-10
const xAbsTr = 1e-10
const xRelTr = 1e-10

function phyloNetworklm_lambda(X::Matrix,
                               Y::Vector,
                               V::MatrixTopologicalOrder,
                               gammas::Vector,
                               times::Vector;
                               nonmissing=trues(length(Y))::BitArray{1}, # Which tips are not missing ?
                               ind=[0]::Vector{Int},
                               ftolRel=fRelTr::AbstractFloat,
                               xtolRel=xRelTr::AbstractFloat,
                               ftolAbs=fAbsTr::AbstractFloat,
                               xtolAbs=xAbsTr::AbstractFloat,
                               startingValue=0.5::Real,
                               fixedValue=missing::Union{Real,Missing})
    if ismissing(fixedValue)
        # Find Best lambda using optimize from package NLopt
        opt = NLopt.Opt(:LN_BOBYQA, 1)
        NLopt.ftol_rel!(opt, ftolRel) # relative criterion
        NLopt.ftol_abs!(opt, ftolAbs) # absolute critetion
        NLopt.xtol_rel!(opt, xtolRel) # criterion on parameter value changes
        NLopt.xtol_abs!(opt, xtolAbs) # criterion on parameter value changes
        NLopt.maxeval!(opt, 1000) # max number of iterations
        NLopt.lower_bounds!(opt, 1e-100) # Lower bound
        # Upper Bound
        up = maxLambda(times, V)
        up = up-up/1000
        NLopt.upper_bounds!(opt, up)
        @info "Maximum lambda value to maintain positive branch lengths: " * @sprintf("%.6g", up)
        count = 0
        function fun(x::Vector{Float64}, g::Vector{Float64})
            x = convert(AbstractFloat, x[1])
            res = logLik_lam(x, X, Y, V, gammas, times; nonmissing=nonmissing, ind=ind)
            count =+ 1
            #println("f_$count: $(round(res, digits=5)), x: $(x)")
            return res
        end
        NLopt.min_objective!(opt, fun)
        fmin, xmin, ret = NLopt.optimize(opt, [startingValue])
        # Best value dans result
        res_lam = xmin[1]
    else
        res_lam = fixedValue
    end
    transform_matrix_lambda!(V, res_lam, gammas, times)
    res = phyloNetworklm_BM(X, Y, V; nonmissing=nonmissing, ind=ind, model="lambda", lambda=res_lam)
#    res.lambda = res_lam
#    res.model = "lambda"
    return res
end

###############################################################################
## Fit scaling hybrid

function matrix_scalingHybrid(net::HybridNetwork, lam::AbstractFloat,
                              gammas::Vector)
    setGammas!(net, 1.0 .- lam .* (1. .- gammas))
    V = sharedPathMatrix(net)
    setGammas!(net, gammas)
    return V
end

function logLik_lam_hyb(lam::AbstractFloat,
                        X::Matrix, Y::Vector,
                        net::HybridNetwork, gammas::Vector;
                        nonmissing=trues(length(Y))::BitArray{1}, # Which tips are not missing ?
                        ind=[0]::Vector{Int})
    # Transform V according to lambda
    V = matrix_scalingHybrid(net, lam, gammas)
    # Fit and take likelihood
    fit_lam = phyloNetworklm_BM(X, Y, V; nonmissing=nonmissing, ind=ind)
    return -loglikelihood(fit_lam)
end

function phyloNetworklm_scalingHybrid(X::Matrix,
                                      Y::Vector,
                                      net::HybridNetwork,
                                      gammas::Vector;
                                      nonmissing=trues(length(Y))::BitArray{1}, # Which tips are not missing ?
                                      ind=[0]::Vector{Int},
                                      ftolRel=fRelTr::AbstractFloat,
                                      xtolRel=xRelTr::AbstractFloat,
                                      ftolAbs=fAbsTr::AbstractFloat,
                                      xtolAbs=xAbsTr::AbstractFloat,
                                      startingValue=0.5::Real,
                                      fixedValue=missing::Union{Real,Missing})
    if ismissing(fixedValue)
        # Find Best lambda using optimize from package NLopt
        opt = NLopt.Opt(:LN_BOBYQA, 1)
        NLopt.ftol_rel!(opt, ftolRel) # relative criterion
        NLopt.ftol_abs!(opt, ftolAbs) # absolute critetion
        NLopt.xtol_rel!(opt, xtolRel) # criterion on parameter value changes
        NLopt.xtol_abs!(opt, xtolAbs) # criterion on parameter value changes
        NLopt.maxeval!(opt, 1000) # max number of iterations
        #NLopt.lower_bounds!(opt, 1e-100) # Lower bound
        #NLopt.upper_bounds!(opt, 1.0)
        count = 0
        function fun(x::Vector{Float64}, g::Vector{Float64})
            x = convert(AbstractFloat, x[1])
            res = logLik_lam_hyb(x, X, Y, net, gammas; nonmissing=nonmissing, ind=ind)
            #count =+ 1
            #println("f_$count: $(round(res, digits=5)), x: $(x)")
            return res
        end
        NLopt.min_objective!(opt, fun)
        fmin, xmin, ret = NLopt.optimize(opt, [startingValue])
        # Best value dans result
        res_lam = xmin[1]
    else
        res_lam = fixedValue
    end
    V = matrix_scalingHybrid(net, res_lam, gammas)
    res = phyloNetworklm_BM(X, Y, V; nonmissing=nonmissing, ind=ind)
    res.lambda = res_lam
    res.model = "scalingHybrid"
    return res
end


"""
    phyloNetworklm(f, fr, net, model="BM",
        fTolRel=1e^-10, fTolAbs=1e^-10, xTolRel=1e^-10, xTolAbs=1e^-10,
        startingValue=0.5)`

Phylogenetic regression, using the correlation structure induced by the network.

Returns an object of class [`PhyloNetworkLinearModel`](@ref). See documentation for this type and
example to see all the functions that can be applied to it.

# Arguments
* `f::Formula`: formula to use for the regression (see the `DataFrame` package)
* `fr::AbstractDataFrame`: DataFrame containing the data and regressors at the tips. It should have an extra column labelled "tipNames", that gives the names of the taxa for each observation.
* `net::HybridNetwork`: phylogenetic network to use. Should have labelled tips.
* `model::AbstractString="BM"`: the model to use, "BM" (default) or "lambda" (for Pagel's lambda).
* `no_names::Bool=false`: if `true`, force the function to ignore the tips names. The data is then assumed to be in the same order as the tips of the network. Default to false, setting it to true is dangerous, and strongly discouraged.
If `model="lambda"`, these parameters control the optimization of lambda:
* `fTolRel::AbstractFloat=1e-10`: relative tolerance on the likelihood value for the optimization in lambda.
* `fTolAbs::AbstractFloat=1e-10`: absolute tolerance on the likelihood value for the optimization in lambda.
* `xTolRel::AbstractFloat=1e-10`: relative tolerance on the parameter value for the optimization in lambda.
* `xTolAbs::AbstractFloat=1e-10`: absolute tolerance on the parameter value for the optimization in lambda.
* `startingValue::Real=0.5`: the starting value for the parameter in the optimization in lambda.

# See also

Type [`PhyloNetworkLinearModel`](@ref), Function [`ancestralStateReconstruction`](@ref)

# Examples

```jldoctest
julia> phy = readTopology(joinpath(dirname(pathof(PhyloNetworks)), "..", "examples", "caudata_tree.txt"));

julia> using CSV # to read data file, next

julia> dat = CSV.read(joinpath(dirname(pathof(PhyloNetworks)), "..", "examples", "caudata_trait.txt"));

julia> using StatsModels # for stat model formulas

julia> fitBM = phyloNetworklm(@formula(trait ~ 1), dat, phy);

julia> fitBM # Shows a summary
StatsModels.DataFrameRegressionModel{PhyloNetworkLinearModel,Array{Float64,2}}

Formula: trait ~ +1

Model: BM

Parameter(s) Estimates:
Sigma2: 0.00294521

Coefficients:
───────────────────────────────────────────────────
             Estimate  Std.Error  t value  Pr(>|t|)
───────────────────────────────────────────────────
(Intercept)     4.679   0.330627  14.1519    <1e-31
───────────────────────────────────────────────────
Log Likelihood: -78.9611507833
AIC: 161.9223015666

julia> round(sigma2_estim(fitBM), digits=6) # rounding for jldoctest convenience
0.002945

julia> round(mu_estim(fitBM), digits=4)
4.679

julia> using StatsBase # for aic() stderror() loglikelihood() etc.

julia> round(loglikelihood(fitBM), digits=10)
-78.9611507833

julia> round(aic(fitBM), digits=10)
161.9223015666

julia> round(aicc(fitBM), digits=10)
161.9841572367

julia> round(bic(fitBM), digits=10)
168.4887090241

julia> round.(coef(fitBM), digits=4)
1-element Array{Float64,1}:
 4.679

julia> confint(fitBM)
1×2 Array{Float64,2}:
 4.02696  5.33104

julia> abs(round(r2(fitBM), digits=10)) # absolute value for jldoctest convenience
0.0

julia> abs(round(adjr2(fitBM), digits=10))
0.0

julia> round.(vcov(fitBM), digits=6)
1×1 Array{Float64,2}:
 0.109314

julia> round.(residuals(fitBM), digits=6)
197-element Array{Float64,1}:
 -0.237648
 -0.357937
 -0.159387
 -0.691868
 -0.323977
 -0.270452
 -0.673486
 -0.584654
 -0.279882
 -0.302175
  ⋮
 -0.777026
 -0.385121
 -0.443444
 -0.327303
 -0.525953
 -0.673486
 -0.603158
 -0.211712
 -0.439833

julia> round.(response(fitBM), digits=5)
197-element Array{Float64,1}:
 4.44135
 4.32106
 4.51961
 3.98713
 4.35502
 4.40855
 4.00551
 4.09434
 4.39912
 4.37682
 ⋮
 3.90197
 4.29388
 4.23555
 4.3517
 4.15305
 4.00551
 4.07584
 4.46729
 4.23917

julia> round.(predict(fitBM), digits=5)
197-element Array{Float64,1}:
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 ⋮
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679

```
""" #"
function phyloNetworklm(f::Formula,
                        fr::AbstractDataFrame,
                        net::HybridNetwork;
                        model="BM"::AbstractString,
                        no_names=false::Bool,
                        ftolRel=fRelTr::AbstractFloat,
                        xtolRel=xRelTr::AbstractFloat,
                        ftolAbs=fAbsTr::AbstractFloat,
                        xtolAbs=xAbsTr::AbstractFloat,
                        startingValue=0.5::Real,
                        fixedValue=missing::Union{Real,Missing})
    # Match the tips names: make sure that the data provided by the user will
    # be in the same order as the ordered tips in matrix V.
    preorder!(net)
    if no_names # The names should not be taken into account.
        ind = [0]
        @info """As requested (no_names=true), I am ignoring the tips names
             in the network and in the dataframe."""
    elseif (any(tipLabels(net) == "") || !any(DataFrames.names(fr) .== :tipNames))
        if (any(tipLabels(net) == "") && !any(DataFrames.names(fr) .== :tipNames))
            error("""The network provided has no tip names, and the input dataframe has
                  no column labelled tipNames, so I can't match the data on the network
                  unambiguously. If you are sure that the tips of the network are in the
                  same order as the values of the dataframe provided, then please re-run
                  this function with argument no_name=true.""")
        end
        if any(tipLabels(net) == "")
            error("""The network provided has no tip names, so I can't match the data
                  on the network unambiguously. If you are sure that the tips of the
                  network are in the same order as the values of the dataframe provided,
                  then please re-run this function with argument no_name=true.""")
        end
        if !any(DataFrames.names(fr) .== :tipNames)
            error("""The input dataframe has no column labelled tipNames, so I can't
                  match the data on the network unambiguously. If you are sure that the
                  tips of the network are in the same order as the values of the dataframe
                  provided, then please re-run this function with argument no_name=true.""")
        end
    else
        #        ind = indexin(V.tipNames, fr[:tipNames])
        ind = indexin(fr[:tipNames], tipLabels(net))
        if any(ind == 0) || length(unique(ind)) != length(ind)
            error("""Tips names of the network and names provided in column tipNames
                  of the dataframe do not match.""")
        end
        #   fr = fr[ind, :]
    end
    # Find the regression matrix and answer vector
    mf = ModelFrame(f,fr)
    if isequal(f.rhs, -1) # If there are no regressors
        mm = ModelMatrix(zeros(size(mf.df, 1), 0), [0])
    else
        mm = ModelMatrix(mf)
    end
    Y = convert(Vector{Float64}, StatsModels.model_response(mf))
    # Fit the model (Method copied from DataFrame/src/statsmodels/statsmodels.jl, lines 47-58)
    # (then StatsModels/src/statsmodels.jl lines 42-46)
    StatsModels.DataFrameRegressionModel(phyloNetworklm(mm.m, Y, net;
                                                       nonmissing=mf.nonmissing, model=model, ind=ind,
                                                       startingValue=startingValue, fixedValue=fixedValue), mf, mm)
end

### Methods on type phyloNetworkRegression

## Un-changed Quantities
# Coefficients of the regression
StatsBase.coef(m::PhyloNetworkLinearModel) = coef(m.lm)
# Number of observations
StatsBase.nobs(m::PhyloNetworkLinearModel) = nobs(m.lm)
# vcov matrix: sigma2_estim * inv(X' * X)
StatsBase.vcov(m::PhyloNetworkLinearModel) = vcov(m.lm)
# standard error of coefficients: sqrt(diag(vcov))
StatsBase.stderror(m::PhyloNetworkLinearModel) = stderror(m.lm)
# confidence Intervals for coefficients:
#  hcat(coef,coef) + stderror * quantile(TDist(dof_residual, (1.-level)/2.) * [1. -1.]
StatsBase.confint(m::PhyloNetworkLinearModel; level=0.95::Real) = confint(m.lm, level)
# coef table: t-values t=coef/se
#    CoefTable(hcat(coef,se,t,ccdf(FDist(1, dof_residual), abs2(t))),
#              ["Estimate","Std.Error","t value", "Pr(>|t|)"],
#              ["x$i" for i = 1:size(X, 2)], 4)
function StatsBase.coeftable(m::PhyloNetworkLinearModel)
    if size(m.lm.pp.X, 2) == 0
        return CoefTable([0], ["Fixed Value"], ["(Intercept)"])
    else
        coeftable(m.lm)
    end
end
# Degrees of freedom for residuals
StatsBase.dof_residual(m::PhyloNetworkLinearModel) =  nobs(m) - length(coef(m))
# Degrees of freedom consumed in the model
function StatsBase.dof(m::PhyloNetworkLinearModel)
    res = length(coef(m)) + 1 # (+1: dispersion parameter)
    if any(m.model .== ["lambda", "scalingHybrid"])
        res += 1 # lambda is one parameter
    end
    return res
end
# Deviance (sum of squared residuals with metric V)
StatsBase.deviance(m::PhyloNetworkLinearModel) = deviance(m.lm)

## Changed Quantities
# Compute the residuals
# (Rescaled by cholesky of variance between tips)
StatsBase.residuals(m::PhyloNetworkLinearModel) = m.RL * residuals(m.lm)
# Tip data
StatsBase.response(m::PhyloNetworkLinearModel) = m.Y
# Predicted values at the tips
# (rescaled by cholesky of tips variances)
StatsBase.predict(m::PhyloNetworkLinearModel) = m.RL * predict(m.lm)
#Log likelihood of the fitted linear model
StatsBase.loglikelihood(m::PhyloNetworkLinearModel) =  loglikelihood(m.lm) - 1/2 * m.logdetVy
# - 0.5*(nobs + nobs * log(2pi) + nobs * log(sigma2_estim) + logdetVy)
# Null  Deviance (sum of squared residuals with metric V)
# REMARK Not just the null deviance of the cholesky regression
# Might be something better to do than this, though.
function StatsBase.nulldeviance(m::PhyloNetworkLinearModel)
    vo = ones(length(m.Y), 1)
    vo = m.RL \ vo
    bo = inv(vo'*vo)*vo'*response(m.lm)
    ro = response(m.lm) - vo*bo
    return sum(ro.^2)
end
# Null Log likelihood (null model with only the intercept)
# Same remark
function StatsBase.nullloglikelihood(m::PhyloNetworkLinearModel)
    n = length(m.Y)
    return -n/2 * (log(2*pi * nulldeviance(m)/n) + 1) - 1/2 * m.logdetVy
end
# coefficient of determination (1 - SS_res/SS_null)
# Copied from GLM.jl/src/lm.jl, line 139
StatsBase.r2(m::PhyloNetworkLinearModel) = 1 - deviance(m)/nulldeviance(m)
# adjusted coefficient of determination
# Copied from GLM.jl/src/lm.jl, lines 141-146
function StatsBase.adjr2(obj::PhyloNetworkLinearModel)
    n = nobs(obj)
    # dof() includes the dispersion parameter
    p = dof(obj) - 1
    1 - (1 - r2(obj))*(n-1)/(n-p)
end

## REMARK
# As PhyloNetworkLinearModel <: LinPredModel, the following functions are automatically defined:
# aic, aicc, bic

## New quantities
# ML estimate for variance of the BM
"""
    sigma2_estim(m::PhyloNetworkLinearModel)

Estimated variance for a fitted object.
"""
sigma2_estim(m::PhyloNetworkLinearModel) = deviance(m.lm) / nobs(m)
# Need to be adapted manually to DataFrameRegressionModel beacouse it's a new function
sigma2_estim(m::StatsModels.DataFrameRegressionModel{PhyloNetworkLinearModel,T} where T) =
  sigma2_estim(m.model)
# ML estimate for ancestral state of the BM
"""
    mu_estim(m::PhyloNetworkLinearModel)

Estimated root value for a fitted object.
"""
function mu_estim(m::PhyloNetworkLinearModel)
    @warn """You fitted the data against a custom matrix, so I have no way
         to know which column is your intercept (column of ones).
         I am using the first coefficient for ancestral mean mu by convention,
         but that might not be what you are looking for."""
    if size(m.lm.pp.X,2) == 0
        return 0
    else
        return coef(m)[1]
    end
end
# Need to be adapted manually to DataFrameRegressionModel beacouse it's a new function
function mu_estim(m::StatsModels.DataFrameRegressionModel{PhyloNetworkLinearModel,T} where T)
    if (!m.mf.terms.intercept)
        error("The fit was done without intercept, so I cannot estimate mu")
    end
    return coef(m)[1]
end
# Lambda estim
"""
    lambda_estim(m::PhyloNetworkLinearModel)

Estimated lambda parameter for a fitted object.
"""
lambda_estim(m::PhyloNetworkLinearModel) = m.lambda
lambda_estim(m::StatsModels.DataFrameRegressionModel{PhyloNetworkLinearModel,T} where T) = lambda_estim(m.model)

### Functions specific to DataFrameRegressionModel
StatsModels.ModelFrame(m::StatsModels.DataFrameRegressionModel) = m.mf
StatsModels.ModelMatrix(m::StatsModels.DataFrameRegressionModel) = m.mm
StatsModels.Formula(m::StatsModels.DataFrameRegressionModel) = Formula(m.mf.terms)
StatsModels.response(m::StatsModels.DataFrameRegressionModel) = response(m.model)

### Print the results
# Variance
function paramstable(m::PhyloNetworkLinearModel)
    Sig = sigma2_estim(m)
    res = "Sigma2: " * @sprintf("%.6g", Sig)
    if any(m.model .== ["lambda", "scalingHybrid"])
        Lamb = lambda_estim(m)
        res = res*"\nLambda: " * @sprintf("%.6g", Lamb)
    end
    return(res)
end
function Base.show(io::IO, obj::PhyloNetworkLinearModel)
    println(io, "$(typeof(obj)):\n\nParameter(s) Estimates:\n", paramstable(obj), "\n\nCoefficients:\n", coeftable(obj))
end
# For DataFrameModel. see also Base.show in
# https://github.com/JuliaStats/StatsModels.jl/blob/master/src/statsmodel.jl
function Base.show(io::IO, model::StatsModels.DataFrameRegressionModel{PhyloNetworkLinearModel,T} where T)
    ct = coeftable(model)
    println(io, "$(typeof(model))")
    println(io)
    println(io, Formula(model.mf.terms))
    println(io)
    println(io, "Model: $(model.model.model)")
    println(io)
    println(io,"Parameter(s) Estimates:")
    println(io, paramstable(model.model))
    println(io)
    println(io,"Coefficients:")
    show(io, ct)
    println(io)
    println(io, "Log Likelihood: "*"$(round(loglikelihood(model), digits=10))")
    println(io, "AIC: "*"$(round(aic(model), digits=10))")
end

###############################################################################
###############################################################################
## Anova - using ftest from GLM - Need version 0.8.1
###############################################################################
###############################################################################

function GLM.ftest(objs::StatsModels.DataFrameRegressionModel{PhyloNetworkLinearModel,T}...)  where T
    objsModels = [obj.model for obj in objs]
    return ftest(objsModels...)
end

function GLM.ftest(objs::PhyloNetworkLinearModel...)
    objslm = [obj.lm for obj in objs]
    return ftest(objslm...)
end

###############################################################################
###############################################################################
## Anova - old version - kept for tests purposes - do not export
###############################################################################
###############################################################################

"""
    anova(objs::PhyloNetworkLinearModel...)

Takes several nested fits of the same data, and computes the F statistic for each
pair of models.

The fits must be results of function [`phyloNetworklm`](@ref) called on the same
data, for models that have more and more effects.

Returns a DataFrame object with the anova table.
"""
function anova(objs::StatsModels.DataFrameRegressionModel{PhyloNetworkLinearModel,T}...) where T
    objsModels = [obj.model for obj in objs]
    return(anova(objsModels...))
end

function anova(objs::PhyloNetworkLinearModel...)
    anovaTable = Array{Any}(undef, length(objs)-1, 6)
    ## Compute binary statistics
    for i in 1:(length(objs) - 1)
      anovaTable[i, :] = anovaBin(objs[i], objs[i+1])
    end
    ## Transform into a DataFrame
    anovaTable = DataFrame(anovaTable)
    names!(anovaTable, [:dof_res, :RSS, :dof, :SS, :F, Symbol("Pr(>F)")])
    return(anovaTable)
end

function anovaBin(obj1::PhyloNetworkLinearModel, obj2::PhyloNetworkLinearModel)
    length(coef(obj1)) < length(coef(obj2)) || error("Models must be nested, from the smallest to the largest.")
    ## residuals
    dof2 = dof_residual(obj2)
    dev2 = deviance(obj2)
    ## reducted residuals
    dof1 = dof_residual(obj1) - dof2
    dev1 = deviance(obj1) - dev2
    ## Compute statistic
    F = (dev1 / dof1) / (dev2 / dof2)
    pval = GLM.ccdf.(GLM.FDist(dof1, dof2), F) # ccdf and FDist from Distributions, used by GLM
    return([dof2, dev2, dof1, dev1, F, pval])
end

###############################################################################
###############################################################################
## Ancestral State Reconstruction
###############################################################################
###############################################################################
# Class for reconstructed states on a network
"""
    ReconstructedStates

Type containing the inferred information about the law of the ancestral states
given the observed tips values. The missing tips are considered as ancestral states.

The following functions can be applied to it:
[`expectations`](@ref) (vector of expectations at all nodes), `stderror` (the standard error),
`predint` (the prediction interval).

The `ReconstructedStates` object has fields: `traits_nodes`, `variances_nodes`, `NodeNumbers`, `traits_tips`, `tipNumbers`, `model`.
Type in "?ReconstructedStates.field" to get help on a specific field.
"""
struct ReconstructedStates
    "traits_nodes: the infered expectation of 'missing' values (ancestral nodes and missing tips)"
    traits_nodes::Vector # Nodes are actually "missing" data (including tips)
    "variances_nodes: the variance covariance matrix between all the 'missing' nodes"
    variances_nodes::Matrix
    "NodeNumbers: vector of the nodes numbers, in the same order as `traits_nodes`"
    NodeNumbers::Vector{Int}
    "traits_tips: the observed traits values at the tips"
    traits_tips::Vector # Observed values at tips
    "TipNumbers: vector of tips numbers, in the same order as `traits_tips`"
    TipNumbers::Vector # Observed tips only
    "model: if not missing, the `PhyloNetworkLinearModel` used for the computations."
    model::Union{PhyloNetworkLinearModel, Missing} # if empirical, corresponding fitted object
end

"""
    expectations(obj::ReconstructedStates)

Estimated reconstructed states at the nodes and tips.
"""
function expectations(obj::ReconstructedStates)
    return DataFrame(nodeNumber = [obj.NodeNumbers; obj.TipNumbers], condExpectation = [obj.traits_nodes; obj.traits_tips])
end

"""
    expectationsPlot(obj::ReconstructedStates)

Compute and format the expected reconstructed states for the plotting function.
The resulting dataframe can be readily used as a `nodeLabel` argument to
`plot` from package [`PhyloPlots`](https://github.com/cecileane/PhyloPlots.jl).
Keyword argument `markMissing` is a string that is appended to predicted
tip values, so that they can be distinguished from the actual datapoints. Default to
"*". Set to "" to remove any visual cue.
"""
function expectationsPlot(obj::ReconstructedStates; markMissing="*"::AbstractString)
    # Retrieve values
    expe = expectations(obj)
    # Format values for plot
    expetxt = Array{AbstractString}(undef, size(expe, 1))
    for i=1:size(expe, 1)
        expetxt[i] = string(round(expe[i, 2], digits=2))
    end
    # Find missing values
    if !ismissing(obj.model)
        nonmissing = obj.model.nonmissing
        ind = obj.model.ind
        missingTipNumbers = obj.model.V.tipNumbers[ind][.!nonmissing]
        indexMissing = indexin(missingTipNumbers, expe[:nodeNumber])
        expetxt[indexMissing] .*= markMissing
    end
    return DataFrame(nodeNumber = [obj.NodeNumbers; obj.TipNumbers], PredInt = expetxt)
end

StatsBase.stderror(obj::ReconstructedStates) = sqrt.(LinearAlgebra.diag(obj.variances_nodes))

"""
    predint(obj::ReconstructedStates; level=0.95::Real)

Prediction intervals with level `level` for internal nodes and missing tips.
"""
function predint(obj::ReconstructedStates; level=0.95::Real)
    if ismissing(obj.model)
        qq = quantile(Normal(), (1. - level)/2.)
    else
        qq = quantile(GLM.TDist(dof_residual(obj.model)), (1. - level)/2.) # TDist from Distributions
        # @warn "As the variance is estimated, the predictions intervals are not exact, and should probably be larger."
    end
    tmpnode = hcat(obj.traits_nodes, obj.traits_nodes) .+ (stderror(obj) * qq) .* [1. -1.]
    return vcat(tmpnode, hcat(obj.traits_tips, obj.traits_tips))
end

function Base.show(io::IO, obj::ReconstructedStates)
    println(io, "$(typeof(obj)):\n",
            CoefTable(hcat(vcat(obj.NodeNumbers, obj.TipNumbers), vcat(obj.traits_nodes, obj.traits_tips), predint(obj)),
                      ["Node index", "Pred.", "Min.", "Max. (95%)"],
                      fill("", length(obj.NodeNumbers)+length(obj.TipNumbers))))
end

"""
    predintPlot(obj::ReconstructedStates; level=0.95::Real, withExp=false::Bool)

Compute and format the prediction intervals for the plotting function.
The resulting dataframe can be readily used as a `nodeLabel` argument to
`plot` from package [`PhyloPlots`](https://github.com/cecileane/PhyloPlots.jl).
Keyworks argument `level` control the confidence level of the
prediction interval. If `withExp` is set to true, then the best
predicted value is also shown along with the interval.
"""
function predintPlot(obj::ReconstructedStates; level=0.95::Real, withExp=false::Bool)
    # predInt
    pri = predint(obj; level=level)
    pritxt = Array{AbstractString}(undef, size(pri, 1))
    # Exp
    withExp ? exptxt = expectationsPlot(obj, markMissing="") : exptxt = ""
    for i=1:length(obj.NodeNumbers)
        !withExp ? sep = ", " : sep = "; " * exptxt[i, 2] * "; "
        pritxt[i] = "[" * string(round(pri[i, 1], digits=2)) * sep * string(round(pri[i, 2], digits=2)) * "]"
    end
    for i=(length(obj.NodeNumbers)+1):size(pri, 1)
        pritxt[i] = string(round(pri[i, 1], digits=2))
    end
    return DataFrame(nodeNumber = [obj.NodeNumbers; obj.TipNumbers], PredInt = pritxt)
end

# """
# 'plot(net::HybridNetwork, obj::ReconstructedStates; kwargs...)
#
# Plot the reconstructed states computed by function `ancestralStateReconstruction`
# on a network.
#
# # Arguments
# * `net::HybridNetwork`: a phylogenetic network.
# * `obj::ReconstructedStates`: the reconstructed states on the network.
# * `kwargs...`: further arguments to be passed to the netwotk `plot` function.
#
# See documentation for function `ancestralStateReconstruction(obj::PhyloNetworkLinearModel[, X_n::Matrix])` for examples.
#
# """
# function Gadfly.plot(net::HybridNetwork, obj::ReconstructedStates; kwargs...)
#   plot(net, nodeLabel = predintPlot(obj); kwargs...)
# end

"""
    ancestralStateReconstruction(net::HybridNetwork, Y::Vector, params::ParamsBM)

Compute the conditional expectations and variances of the ancestral (un-observed)
traits values at the internal nodes of the phylogenetic network (`net`),
given the values of the traits at the tips of the network (`Y`) and some
known parameters of the process used for trait evolution (`params`, only BM with fixed root
works for now).

This function assumes that the parameters of the process are known. For a more general
function, see `ancestralStateReconstruction(obj::PhyloNetworkLinearModel[, X_n::Matrix])`.

"""
function ancestralStateReconstruction(net::HybridNetwork,
                                      Y::Vector,
                                      params::ParamsBM)
    V = sharedPathMatrix(net)
    ancestralStateReconstruction(V, Y, params)
end

function ancestralStateReconstruction(V::MatrixTopologicalOrder,
                                      Y::Vector,
                                      params::ParamsBM)
    # Variances matrices
    Vy = V[:Tips]
    Vz = V[:InternalNodes]
    Vyz = V[:TipsNodes]
    R = cholesky(Vy)
    RL = R.L
    temp = RL \ Vyz
    # Vectors of means
    m_y = ones(size(Vy)[1]) .* params.mu # !! correct only if no predictor.
    m_z = ones(size(Vz)[1]) .* params.mu # !! works if BM no shift.
    # Actual computation
    ancestralStateReconstruction(Vz, temp, RL,
                                 Y, m_y, m_z,
                                 V.internalNodeNumbers,
                                 V.tipNumbers,
                                 params.sigma2)
end

# Reconstruction from all the needed quantities
function ancestralStateReconstruction(Vz::Matrix,
                                      VyzVyinvchol::Matrix,
                                      RL::LowerTriangular,
                                      Y::Vector, m_y::Vector, m_z::Vector,
                                      NodeNumbers::Vector,
                                      TipNumbers::Vector,
                                      sigma2::Real,
                                      add_var=zeros(size(Vz))::Matrix, # Additional variance for BLUP
                                      model=missing::Union{PhyloNetworkLinearModel,Missing})
    m_z_cond_y = m_z + VyzVyinvchol' * (RL \ (Y - m_y))
    V_z_cond_y = sigma2 .* (Vz - VyzVyinvchol' * VyzVyinvchol)
    ReconstructedStates(m_z_cond_y, V_z_cond_y + add_var, NodeNumbers, Y, TipNumbers, model)
end

# """
# `ancestralStateReconstruction(obj::PhyloNetworkLinearModel, X_n::Matrix)`
# Function to find the ancestral traits reconstruction on a network, given an
# object fitted by function phyloNetworklm, and some predictors expressed at all the nodes of the network.
#
# - obj: a PhyloNetworkLinearModel object, or a
# DataFrameRegressionModel{PhyloNetworkLinearModel}, if data frames were used.
# - X_n a matrix with as many columns as the number of predictors used, and as
# many lines as the number of unknown nodes or tips.
#
# Returns an object of type ancestralStateReconstruction.
# """

# Empirical reconstruciton from a fitted object
# TO DO: Handle the order of internal nodes for matrix X_n
function ancestralStateReconstruction(obj::PhyloNetworkLinearModel, X_n::Matrix)
    if (size(X_n)[2] != length(coef(obj)))
        error("""The number of predictors for the ancestral states (number of columns of X_n)
              does not match the number of predictors at the tips.""")
    end
    if size(X_n)[1] != length(obj.V.internalNodeNumbers) + sum(.!obj.nonmissing)
        error("""The number of lines of the predictors does not match
              the number of nodes plus the number of missing tips.""")
    end
    m_y = predict(obj)
    m_z = X_n * coef(obj)
    # If the tips were re-organized, do the same for Vyz
    if (obj.ind != [0])
#       iii = indexin(1:length(obj.nonmissing), obj.ind[obj.nonmissing])
#       iii = iii[iii .> 0]
#       jjj = [1:length(obj.V.internalNodeNumbers); indexin(1:length(obj.nonmissing), obj.ind[!obj.nonmissing])]
#       jjj = jjj[jjj .> 0]
#       Vyz = Vyz[iii, jjj]
        Vyz = obj.V[:TipsNodes, obj.ind, obj.nonmissing]
        missingTipNumbers = obj.V.tipNumbers[obj.ind][.!obj.nonmissing]
        nmTipNumbers = obj.V.tipNumbers[obj.ind][obj.nonmissing]
    else
        @warn """There were no indication for the position of the tips on the network.
             I am assuming that they are given in the same order.
             Please check that this is what you intended."""
        Vyz = obj.V[:TipsNodes, collect(1:length(obj.V.tipNumbers)), obj.nonmissing]
        missingTipNumbers = obj.V.tipNumbers[.!obj.nonmissing]
        nmTipNumbers = obj.V.tipNumbers[obj.nonmissing]
    end
    temp = obj.RL \ Vyz
    U = X_n - temp' * (obj.RL \ obj.X)
    add_var = U * vcov(obj) * U'
    # Warn about the prediction intervals
    @warn """These prediction intervals show uncertainty in ancestral values,
         assuming that the estimated variance rate of evolution is correct.
         Additional uncertainty in the estimation of this variance rate is
         ignored, so prediction intervals should be larger."""
    # Actual reconstruction
    ancestralStateReconstruction(obj.V[:InternalNodes, obj.ind, obj.nonmissing],
                                 temp,
                                 obj.RL,
                                 obj.Y,
                                 m_y,
                                 m_z,
                                 [obj.V.internalNodeNumbers; missingTipNumbers],
                                 nmTipNumbers,
                                 sigma2_estim(obj),
                                 add_var,
                                 obj)
end

"""
    ancestralStateReconstruction(obj::PhyloNetworkLinearModel[, X_n::Matrix])

Function to find the ancestral traits reconstruction on a network, given an
object fitted by function [`phyloNetworklm`](@ref). By default, the function assumes
that the regressor is just an intercept. If the value of the regressor for
all the ancestral states is known, it can be entered in X_n, a matrix with as
many columns as the number of predictors used, and as many lines as the number
of unknown nodes or tips.

Returns an object of type [`ReconstructedStates`](@ref).
See documentation for this type and examples for functions that can be applied to it.

# Examples

```jldoctest
julia> using CSV # to read data file

julia> phy = readTopology(joinpath(dirname(pathof(PhyloNetworks)), "..", "examples", "carnivores_tree.txt"));

julia> dat = CSV.read(joinpath(dirname(pathof(PhyloNetworks)), "..", "examples", "carnivores_trait.txt"));

julia> using StatsModels # for statistical model formulas

julia> fitBM = phyloNetworklm(@formula(trait ~ 1), dat, phy);

julia> ancStates = ancestralStateReconstruction(fitBM) # Should produce a warning, as variance is unknown.
┌ Warning: These prediction intervals show uncertainty in ancestral values,
│ assuming that the estimated variance rate of evolution is correct.
│ Additional uncertainty in the estimation of this variance rate is
│ ignored, so prediction intervals should be larger.
└ @ PhyloNetworks ~/build/crsl4/PhyloNetworks.jl/src/traits.jl:2163
ReconstructedStates:
───────────────────────────────────────────────
  Node index      Pred.        Min.  Max. (95%)
───────────────────────────────────────────────
        -5.0   1.32139   -0.288423     2.9312
        -8.0   1.03258   -0.539072     2.60423
        -7.0   1.41575   -0.0934395    2.92495
        -6.0   1.39417   -0.0643135    2.85265
        -4.0   1.39961   -0.0603343    2.85955
        -3.0   1.51341   -0.179626     3.20644
       -13.0   5.3192     3.96695      6.67145
       -12.0   4.51176    2.94268      6.08085
       -16.0   1.50947    0.0290151    2.98992
       -15.0   1.67425    0.241696     3.10679
       -14.0   1.80309    0.355568     3.2506
       -11.0   2.7351     1.21896      4.25123
       -10.0   2.73217    1.16545      4.29889
        -9.0   2.41132    0.639075     4.18357
        -2.0   2.04138   -0.0340955    4.11686
        14.0   1.64289    1.64289      1.64289
         8.0   1.67724    1.67724      1.67724
         5.0   0.331568   0.331568     0.331568
         2.0   2.27395    2.27395      2.27395
         4.0   0.275237   0.275237     0.275237
         6.0   3.39094    3.39094      3.39094
        13.0   0.355799   0.355799     0.355799
        15.0   0.542565   0.542565     0.542565
         7.0   0.773436   0.773436     0.773436
        10.0   6.94985    6.94985      6.94985
        11.0   4.78323    4.78323      4.78323
        12.0   5.33016    5.33016      5.33016
         1.0  -0.122604  -0.122604    -0.122604
        16.0   0.73989    0.73989      0.73989
         9.0   4.84236    4.84236      4.84236
         3.0   1.0695     1.0695       1.0695
───────────────────────────────────────────────

julia> expectations(ancStates)
31×2 DataFrames.DataFrame
│ Row │ nodeNumber │ condExpectation │
│     │ Int64      │ Float64         │
├─────┼────────────┼─────────────────┤
│ 1   │ -5         │ 1.32139         │
│ 2   │ -8         │ 1.03258         │
│ 3   │ -7         │ 1.41575         │
│ 4   │ -6         │ 1.39417         │
│ 5   │ -4         │ 1.39961         │
│ 6   │ -3         │ 1.51341         │
│ 7   │ -13        │ 5.3192          │
⋮
│ 24  │ 7          │ 0.773436        │
│ 25  │ 10         │ 6.94985         │
│ 26  │ 11         │ 4.78323         │
│ 27  │ 12         │ 5.33016         │
│ 28  │ 1          │ -0.122604       │
│ 29  │ 16         │ 0.73989         │
│ 30  │ 9          │ 4.84236         │
│ 31  │ 3          │ 1.0695          │

julia> predint(ancStates)
31×2 Array{Float64,2}:
 -0.288423    2.9312
 -0.539072    2.60423
 -0.0934395   2.92495
 -0.0643135   2.85265
 -0.0603343   2.85955
 -0.179626    3.20644
  3.96695     6.67145
  2.94268     6.08085
  0.0290151   2.98992
  0.241696    3.10679
  ⋮
  0.542565    0.542565
  0.773436    0.773436
  6.94985     6.94985
  4.78323     4.78323
  5.33016     5.33016
 -0.122604   -0.122604
  0.73989     0.73989
  4.84236     4.84236
  1.0695      1.0695

julia> expectationsPlot(ancStates) # format the ancestral states
31×2 DataFrames.DataFrame
│ Row │ nodeNumber │ PredInt   │
│     │ Int64      │ Abstract… │
├─────┼────────────┼───────────┤
│ 1   │ -5         │ 1.32      │
│ 2   │ -8         │ 1.03      │
│ 3   │ -7         │ 1.42      │
│ 4   │ -6         │ 1.39      │
│ 5   │ -4         │ 1.4       │
│ 6   │ -3         │ 1.51      │
│ 7   │ -13        │ 5.32      │
⋮
│ 24  │ 7          │ 0.77      │
│ 25  │ 10         │ 6.95      │
│ 26  │ 11         │ 4.78      │
│ 27  │ 12         │ 5.33      │
│ 28  │ 1          │ -0.12     │
│ 29  │ 16         │ 0.74      │
│ 30  │ 9          │ 4.84      │
│ 31  │ 3          │ 1.07      │

julia> using PhyloPlots # next: plot ancestral states on the tree

julia> plot(phy, :RCall, nodeLabel = expectationsPlot(ancStates));

julia> predintPlot(ancStates)
31×2 DataFrames.DataFrame
│ Row │ nodeNumber │ PredInt       │
│     │ Int64      │ Abstract…     │
├─────┼────────────┼───────────────┤
│ 1   │ -5         │ [-0.29, 2.93] │
│ 2   │ -8         │ [-0.54, 2.6]  │
│ 3   │ -7         │ [-0.09, 2.92] │
│ 4   │ -6         │ [-0.06, 2.85] │
│ 5   │ -4         │ [-0.06, 2.86] │
│ 6   │ -3         │ [-0.18, 3.21] │
│ 7   │ -13        │ [3.97, 6.67]  │
⋮
│ 24  │ 7          │ 0.77          │
│ 25  │ 10         │ 6.95          │
│ 26  │ 11         │ 4.78          │
│ 27  │ 12         │ 5.33          │
│ 28  │ 1          │ -0.12         │
│ 29  │ 16         │ 0.74          │
│ 30  │ 9          │ 4.84          │
│ 31  │ 3          │ 1.07          │

julia> plot(phy, :RCall, nodeLabel = predintPlot(ancStates));

julia> using DataFrames # to use allowmissing!

julia> allowmissing!(dat, :trait);

julia> dat[[2, 5], :trait] = missing; # missing values allowed to fit model

julia> fitBM = phyloNetworklm(@formula(trait ~ 1), dat, phy);

julia> ancStates = ancestralStateReconstruction(fitBM);
┌ Warning: These prediction intervals show uncertainty in ancestral values,
│ assuming that the estimated variance rate of evolution is correct.
│ Additional uncertainty in the estimation of this variance rate is
│ ignored, so prediction intervals should be larger.
└ @ PhyloNetworks ~/build/crsl4/PhyloNetworks.jl/src/traits.jl:2163

julia> expectations(ancStates)
31×2 DataFrames.DataFrame
│ Row │ nodeNumber │ condExpectation │
│     │ Int64      │ Float64         │
├─────┼────────────┼─────────────────┤
│ 1   │ -5         │ 1.42724         │
│ 2   │ -8         │ 1.35185         │
│ 3   │ -7         │ 1.61993         │
│ 4   │ -6         │ 1.54198         │
│ 5   │ -4         │ 1.53916         │
│ 6   │ -3         │ 1.64984         │
│ 7   │ -13        │ 5.33508         │
⋮
│ 24  │ 7          │ 0.773436        │
│ 25  │ 10         │ 6.94985         │
│ 26  │ 11         │ 4.78323         │
│ 27  │ 12         │ 5.33016         │
│ 28  │ 1          │ -0.122604       │
│ 29  │ 16         │ 0.73989         │
│ 30  │ 9          │ 4.84236         │
│ 31  │ 3          │ 1.0695          │

julia> predint(ancStates)
31×2 Array{Float64,2}:
 -0.31245     3.16694
 -0.625798    3.3295
 -0.110165    3.35002
 -0.0710391   3.15501
 -0.0675924   3.14591
 -0.197236    3.49692
  3.89644     6.77373
  2.8741      6.22808
 -0.0358627   3.12834
  0.182594    3.2534
  ⋮
  0.542565    0.542565
  0.773436    0.773436
  6.94985     6.94985
  4.78323     4.78323
  5.33016     5.33016
 -0.122604   -0.122604
  0.73989     0.73989
  4.84236     4.84236
  1.0695      1.0695

julia> expectationsPlot(ancStates) # format node <-> ancestral state
31×2 DataFrames.DataFrame
│ Row │ nodeNumber │ PredInt   │
│     │ Int64      │ Abstract… │
├─────┼────────────┼───────────┤
│ 1   │ -5         │ 1.43      │
│ 2   │ -8         │ 1.35      │
│ 3   │ -7         │ 1.62      │
│ 4   │ -6         │ 1.54      │
│ 5   │ -4         │ 1.54      │
│ 6   │ -3         │ 1.65      │
│ 7   │ -13        │ 5.34      │
⋮
│ 24  │ 7          │ 0.77      │
│ 25  │ 10         │ 6.95      │
│ 26  │ 11         │ 4.78      │
│ 27  │ 12         │ 5.33      │
│ 28  │ 1          │ -0.12     │
│ 29  │ 16         │ 0.74      │
│ 30  │ 9          │ 4.84      │
│ 31  │ 3          │ 1.07      │

julia> plot(phy, :RCall, nodeLabel = expectationsPlot(ancStates));

julia> predintPlot(ancStates) # prediction intervals, in data frame, useful to plot
31×2 DataFrames.DataFrame
│ Row │ nodeNumber │ PredInt       │
│     │ Int64      │ Abstract…     │
├─────┼────────────┼───────────────┤
│ 1   │ -5         │ [-0.31, 3.17] │
│ 2   │ -8         │ [-0.63, 3.33] │
│ 3   │ -7         │ [-0.11, 3.35] │
│ 4   │ -6         │ [-0.07, 3.16] │
│ 5   │ -4         │ [-0.07, 3.15] │
│ 6   │ -3         │ [-0.2, 3.5]   │
│ 7   │ -13        │ [3.9, 6.77]   │
⋮
│ 24  │ 7          │ 0.77          │
│ 25  │ 10         │ 6.95          │
│ 26  │ 11         │ 4.78          │
│ 27  │ 12         │ 5.33          │
│ 28  │ 1          │ -0.12         │
│ 29  │ 16         │ 0.74          │
│ 30  │ 9          │ 4.84          │
│ 31  │ 3          │ 1.07          │

julia> plot(phy, :RCall, nodeLabel = predintPlot(ancStates));
```
"""
function ancestralStateReconstruction(obj::PhyloNetworkLinearModel)
    # default reconstruction for known predictors
    if ((size(obj.X)[2] != 1) || !any(obj.X .== 1)) # Test if the regressor is just an intercept.
        error("""Predictor(s) other than a plain intercept are used in this `PhyloNetworkLinearModel` object.
    These predictors are unobserved at ancestral nodes, so they cannot be used
    for the ancestral state reconstruction. If these ancestral predictor values
    are known, please provide them as a matrix argument to the function.
    Otherwise, you might consider doing a multivariate linear regression (not implemented yet).""")
    end
  X_n = ones((length(obj.V.internalNodeNumbers) + sum(.!obj.nonmissing), 1))
    ancestralStateReconstruction(obj, X_n)
end
# For a DataFrameRegressionModel
function ancestralStateReconstruction(obj::StatsModels.DataFrameRegressionModel{PhyloNetworkLinearModel,T} where T)
    ancestralStateReconstruction(obj.model)
end
function ancestralStateReconstruction(obj::StatsModels.DataFrameRegressionModel{PhyloNetworkLinearModel,T} where T, X_n::Matrix)
    ancestralStateReconstruction(obj.model, X_n)
end

"""
    ancestralStateReconstruction(fr::AbstractDataFrame, net::HybridNetwork; kwargs...)

Function to find the ancestral traits reconstruction on a network, given some data at the tips.
Uses function [`phyloNetworklm`](@ref) to perform a phylogenetic regression of the data against an
intercept (amounts to fitting an evolutionary model on the network, BM being the only option
available for now).

See documentation on [`phyloNetworklm`](@ref) and `ancestralStateReconstruction(obj::PhyloNetworkLinearModel[, X_n::Matrix])`
for further details.

Returns an object of type [`ReconstructedStates`](@ref).
"""
function ancestralStateReconstruction(fr::AbstractDataFrame,
                                      net::HybridNetwork;
                                      kwargs...)
    nn = names(fr)
    datpos = nn .!= :tipNames
    if sum(datpos) > 1
        error("""Besides one column labelled 'tipNames', the dataframe fr should have
              only one column, corresponding to the data at the tips of the network.""")
    end
    f = @eval(@formula($(nn[datpos][1]) ~ 1))
    reg = phyloNetworklm(f, fr, net; kwargs...)
    return ancestralStateReconstruction(reg)
end

#################################################
## Old version of phyloNetworklm (naive)
#################################################

# function phyloNetworklmNaive(X::Matrix, Y::Vector, net::HybridNetwork, model="BM"::AbstractString)
#   # Geting variance covariance
#   V = sharedPathMatrix(net)
#   Vy = extractVarianceTips(V, net)
#   # Needed quantities (naive)
#   ntaxa = length(Y)
#   Vyinv = inv(Vy)
#   XtVyinv = X' * Vyinv
#   logdetVy = logdet(Vy)
#        # beta hat
#   betahat = inv(XtVyinv * X) * XtVyinv * Y
#        # sigma2 hat
#   fittedValues =  X * betahat
#   residuals = Y - fittedValues
#   sigma2hat = 1/ntaxa * (residuals' * Vyinv * residuals)
#        # log likelihood
#   loglik = - 1 / 2 * (ntaxa + ntaxa * log(2 * pi) + ntaxa * log(sigma2hat) + logdetVy)
#   # Result
# # res = phyloNetworkRegression(betahat, sigma2hat[1], loglik[1], V, Vy, fittedValues, residuals)
#   return((betahat, sigma2hat[1], loglik[1], V, Vy, logdetVy, fittedValues, residuals))
# end
