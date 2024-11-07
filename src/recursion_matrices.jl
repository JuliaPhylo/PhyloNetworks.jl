"""
    MatrixTopologicalOrder

Matrix associated to a [`HybridNetwork`](@ref) in which rows and/or columns
correspond to nodes in the network. In this matrix, nodes are indexed by
their topological order.
For example, if rows list nodes in network `net`, then `V[i,:]` corresponds to
node `net.nodes_changed[i]`.

The following functions and extractors can be applied to it:
[`tipLabels`](@ref),
`obj[:Tips]`, `obj[:InternalNodes]`, `obj[:TipsNodes]`.
See documentation for [`getindex(::MatrixTopologicalOrder, ::Symbol)`](@ref)).

A `MatrixTopologicalOrder` object has field
`V` (the matrix of interest)
and fields for mapping indices in `V` to node numbers
`nodeNumbersTopOrder`,
`internalNodeNumbers`,
`tipNumbers`,
`tipNames`,
`indexation`.
Type in "?MatrixTopologicalOrder.field" to get documentation on a specific field.
"""
struct MatrixTopologicalOrder
    "V: the matrix per se. V[i,:] and/or V[:,i] has data for the i-th node in topological order."
    V::Matrix
    "vector of node numbers, with nodes listed in topological order as in V"
    nodeNumbersTopOrder::Vector{Int} # for ordering of the matrix
    "vector of internal node numbers, as listed in net.node"
    internalNodeNumbers::Vector{Int}
    "vector of tip numbers, listed in the same order as in net.node (should be same as in net.leaf...)"
    tipNumbers::Vector{Int}
    "vector of tip names, listed in the same order as in net.node"
    tipNames::Vector
    """
    indexation: a string giving the type of matrix `V`:
    - `:r`: rows only are indexed by the network nodes
    - `:c`: columns only are indexed by the network nodes
    - `:b`: both rows and columns are indexed by nodes.
    """
    indexation::Symbol
end

function Base.show(io::IO, obj::MatrixTopologicalOrder)
    println(io, "$(typeof(obj)):\n$(obj.V)")
end

tipLabels(obj::MatrixTopologicalOrder) = obj.tipNames

function MatrixTopologicalOrder(V::Matrix, net::HybridNetwork, indexation)
    length(net.nodes_changed) == length(net.node) ||
        error("run preorder! on the network first")
    nodenumbers_internal = Int[]
    nodenumbers_leaves = Int[]
    for n in net.node
        push!((n.leaf ? nodenumbers_leaves : nodenumbers_internal), n.number)
    end
    return MatrixTopologicalOrder(
        V,
        [n.number for n in net.nodes_changed],
        nodenumbers_internal,
        nodenumbers_leaves,
        [n.name for n in net.leaf],
        indexation
    )
end

"""
    getindex(obj::MatrixTopologicalOrder,
             d::Symbol,[ indTips, nonmissing])

Get submatrices of a [`MatrixTopologicalOrder`](@ref) object.
In `obj.V`, row and/column `i` corresponds to the `i`th node in topological order.
In contrast, in matrix `obj[:tips]` for example,
row and/or column `i` corresponds to the `i`th tip when tips are listed in the
same order as in the network's original `.node` vector.

Arguments:

- `obj::MatrixTopologicalOrder`: the matrix from which to extract.
- `d::Symbol`: a symbol specifying which sub-matrix to extract. Can be:
  * `:Tips` columns and/or rows corresponding to the tips
  * `:InternalNodes` columns and/or rows corresponding to the internal nodes
    Includes tips not listed in `indTips` or missing data according to `nonmissing`.
  * `:TipsNodes` columns corresponding to internal nodes, and row to tips (works only is indexation="b")
- `indTips::Vector{Int}`: optional argument precising a specific order for the tips (internal use).
- `nonmissing::BitArray{1}`: optional argument saying which tips have data (internal use).
   Tips with missing data are treated as internal nodes.
"""
function Base.getindex(
    obj::MatrixTopologicalOrder,
    d::Symbol,
    indTips::Vector{Int}=collect(1:length(obj.tipNumbers)),
    nonmissing::BitArray{1}=trues(length(obj.tipNumbers))
)
    tipnums = obj.tipNumbers[indTips][nonmissing]
    maskTips = indexin(tipnums, obj.nodeNumbersTopOrder)
    if d == :Tips # Extract rows and/or columns corresponding to the tips with data
        obj.indexation == :b && return obj.V[maskTips, maskTips] # both columns and rows are indexed by nodes
        obj.indexation == :c && return obj.V[:, maskTips] # Only the columns
        obj.indexation == :r && return obj.V[maskTips, :] # Only the rows
    end
    intnodenums = [obj.internalNodeNumbers ; setdiff(obj.tipNumbers, tipnums)]
    maskNodes = indexin(intnodenums, obj.nodeNumbersTopOrder)
    #= indices in obj.nodeNumbersTopOrder, in this order:
    1. internal nodes, in the same order as in obj.internalNodeNumbers,
       that is, same order as in net.node (excluding leaves)
    2. tips absent from indTips or missing data according to nonmissing,
       in the same order as in obj.tipNumbers.
    =#
    if d == :InternalNodes # Idem, for internal nodes
        obj.indexation == :b && return obj.V[maskNodes, maskNodes]
        obj.indexation == :c && return obj.V[:, maskNodes]
        obj.indexation == :r && return obj.V[maskNodes, :]
    end
    if d == :TipsNodes
        obj.indexation == :b && return obj.V[maskTips, maskNodes]
        obj.indexation == :c && error("""Both rows and columns must be net
                                       ordered to take the submatrix tips vs internal nodes.""")
        obj.indexation == :r && error("""Both rows and columns must be net
                                       ordered to take the submatrix tips vs internal nodes.""")
    end
    d == :All && return obj.V
end

"""
    vcv(net::HybridNetwork; model::AbstractString="BM",
                            corr::Bool=false,
                            checkpreorder::Bool=true)

This function computes the variance covariance matrix between the tips of the
network, assuming a Brownian model of trait evolution (with unit variance).
If optional argument `corr` is set to `true`, then the correlation matrix is returned instead.

The function returns a `DataFrame` object, with columns named by the tips of the network.

The calculation of the covariance matrix requires a pre-ordering of nodes to be fast.
If `checkpreorder` is true (default), then [`preorder!`](@ref) is run on the network beforehand.
Otherwise, the network is assumed to be already in pre-order.

This function internally calls [`sharedPathMatrix`](@ref), which computes the variance
matrix between all the nodes of the network.

# Examples

```jldoctest
julia> tree_str = "(((t2:0.14,t4:0.33):0.59,t3:0.96):0.14,(t5:0.70,t1:0.18):0.90);";

julia> tree = readnewick(tree_str);

julia> C = vcv(tree)
5×5 DataFrame
 Row │ t2       t4       t3       t5       t1      
     │ Float64  Float64  Float64  Float64  Float64 
─────┼─────────────────────────────────────────────
   1 │    0.87     0.73     0.14      0.0     0.0
   2 │    0.73     1.06     0.14      0.0     0.0
   3 │    0.14     0.14     1.1       0.0     0.0
   4 │    0.0      0.0      0.0       1.6     0.9
   5 │    0.0      0.0      0.0       0.9     1.08

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
(for the model, see Bastide et al. 2018)
```jldoctest
julia> net = readnewick("((t1:1.0,#H1:0.1::0.30):0.5,((t2:0.9)#H1:0.2::0.70,t3:1.1):0.4);");

julia> C = vcv(net)
3×3 DataFrame
 Row │ t1       t2       t3      
     │ Float64  Float64  Float64 
─────┼───────────────────────────
   1 │    1.5     0.15      0.0
   2 │    0.15    1.248     0.28
   3 │    0.0     0.28      1.5
```
"""
function vcv(
    net::HybridNetwork;
    model::AbstractString="BM",
    corr::Bool=false,
    checkpreorder::Bool=true
)
    @assert (model == "BM") "The 'vcv' function only works for a BM process (for now)."
    V = sharedPathMatrix(net; checkpreorder=checkpreorder)
    C = V[:Tips]
    corr && StatsBase.cov2cor!(C, sqrt.(diag(C)))
    Cd = DataFrame(C, map(Symbol, V.tipNames))
    return(Cd)
end


"""
    sharedPathMatrix(net::HybridNetwork; checkpreorder::Bool=true)

This function computes the shared path matrix between all the nodes of a
network. It assumes that the network is in the pre-order. If checkpreorder is
true (default), then it runs function `preorder!` on the network beforehand.

Returns an object of type [`MatrixTopologicalOrder`](@ref).

"""
function sharedPathMatrix(net::HybridNetwork; checkpreorder::Bool=true)
    check_nonmissing_nonnegative_edgelengths(net,
        """The variance-covariance matrix of the network is not defined.
           A phylogenetic regression cannot be done.""")
    checkpreorder && preorder!(net)
    V = traversal_preorder(
            net.nodes_changed,
            initsharedPathMatrix,
            traversalupdate_default!,
            updateTreeSharedPathMatrix!,
            updateHybridSharedPathMatrix!,
        )
    M = MatrixTopologicalOrder(V, net, :b) # nodes in both columns & rows
    return M
end

function initsharedPathMatrix(nodes::Vector{Node},)
    n = length(nodes)
    return(zeros(Float64,n,n))
end

function updateTreeSharedPathMatrix!(
    V::Matrix,
    i::Int,
    parentIndex::Int,
    edge::Edge,
)
    for j in 1:(i-1)
        V[i,j] = V[j,parentIndex]
        V[j,i] = V[j,parentIndex]
    end
    V[i,i] = V[parentIndex,parentIndex] + edge.length
    return true
end

function updateHybridSharedPathMatrix!(
    V::Matrix,
    i::Int,
    parindx::AbstractVector{Int},
    paredge::AbstractVector{Edge},
)
    for j in 1:(i-1)
        for (pi,pe) in zip(parindx, paredge)
            V[i,j] += pe.gamma * V[pi,j]
        end
        V[j,i] = V[i,j]
    end
    for k1 in eachindex(paredge)
        p1i = parindx[k1]
        p1e = paredge[k1]
        V[i,i] +=  p1e.gamma^2 * (V[p1i,p1i] + p1e.length)
        for k2 in (k1+1):length(paredge)
            V[i,i] += 2 * p1e.gamma * paredge[k2].gamma * V[p1i,parindx[k2]]
        end
    end
    return true
end

"""
    descendenceMatrix(net::HybridNetwork; checkpreorder::Bool=true)

Descendence matrix between all the nodes of a network:
object `D` of type [`MatrixTopologicalOrder`](@ref) in which
`D[i,j]` is the proportion of genetic material in node `i` that can be traced
back to node `j`. If `D[i,j]>0` then `j` is a descendent of `i` (and `j` is
an ancestor of `i`).
The network is assumed to be pre-ordered if `checkpreorder` is false.
If `checkpreorder` is true (default), `preorder!` is run on the network beforehand.
"""
function descendenceMatrix(net::HybridNetwork; checkpreorder::Bool=true)
    checkpreorder && preorder!(net)
    V = traversal_postorder(
        net.nodes_changed,
        initDescendenceMatrix,
        traversalupdate_default!, # does nothing
        updateNodeDescendenceMatrix!,
    )
    M = MatrixTopologicalOrder(V, net, :r) # nodes in rows
    return M
end

function updateNodeDescendenceMatrix!(
    V::Matrix,
    i::Int,
    childrenIndex::Vector{Int},
    edges::Vector{Edge},
)
    for j in 1:length(edges)
        V[:,i] .+= edges[j].gamma .* V[:,childrenIndex[j]]
    end
    return true
end

function initDescendenceMatrix(nodes::Vector{Node},)
    n = length(nodes)
    return(Matrix{Float64}(I, n, n)) # identity matrix
end
