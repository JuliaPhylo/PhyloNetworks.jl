"""
    MatrixTopologicalOrder

Matrix associated to a [`HybridNetwork`](@ref) in which rows and/or columns
correspond to nodes in the network. In this matrix, nodes are indexed by
their topological order.
For example, if rows list nodes in network `net`, then `V[i,:]` corresponds to
node `net.vec_node[i]`.

The following functions and extractors can be applied to it:
[`tiplabels`](@ref),
`obj[:tips]`, `obj[:internalnodes]`, `obj[:tipsnodes]`.
See documentation for [`getindex(::MatrixTopologicalOrder, ::Symbol)`](@ref)).

A `MatrixTopologicalOrder` object has field
`V` (the matrix of interest)
and fields for mapping indices in `V` to node numbers
`nodenumbers_toporder`,
`internalnodenumbers`,
`tipnumbers`,
`tipnames`,
`indexation`.
Type in "?MatrixTopologicalOrder.field" to get documentation on a specific field.
"""
struct MatrixTopologicalOrder
    "V: the matrix per se. V[i,:] and/or V[:,i] has data for the i-th node in topological order."
    V::Matrix
    "vector of node numbers, with nodes listed in topological order `net.vec_node` as in V"
    nodenumbers_toporder::Vector{Int} # for ordering of the matrix
    "vector of internal node numbers, as listed in net.node"
    internalnodenumbers::Vector{Int}
    "vector of tip numbers, listed in the same order as in net.node (should be same as in net.leaf...)"
    tipnumbers::Vector{Int}
    "vector of tip names, listed in the same order as in net.node"
    tipnames::Vector
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

tiplabels(obj::MatrixTopologicalOrder) = obj.tipnames

function MatrixTopologicalOrder(V::Matrix, net::HybridNetwork, indexation)
    length(net.vec_node) == length(net.node) ||
        error("run preorder! on the network first")
    nodenumbers_internal = Int[]
    nodenumbers_leaves = Int[]
    for n in net.node
        push!((n.leaf ? nodenumbers_leaves : nodenumbers_internal), n.number)
    end
    return MatrixTopologicalOrder(
        V,
        [n.number for n in net.vec_node],
        nodenumbers_internal,
        nodenumbers_leaves,
        [n.name for n in net.leaf],
        indexation
    )
end

"""
    getindex(obj::MatrixTopologicalOrder,
             d::Symbol,[ indtips, nonmissing])

Get submatrices of a [`MatrixTopologicalOrder`](@ref) object.
In `obj.V`, row and/column `i` corresponds to the `i`th node in topological order.
In contrast, in matrix `obj[:tips]` for example,
row and/or column `i` corresponds to the `i`th tip when tips are listed in the
same order as in the network's original `.node` vector.

Arguments:

- `obj::MatrixTopologicalOrder`: the matrix from which to extract.
- `d`: symbol specifying which sub-matrix to extract. Can be:
  * `:tips` columns and/or rows corresponding to the tips
  * `:internalnodes` columns and/or rows corresponding to the internal nodes
    Includes tips not listed in `indtips` or missing data according to `nonmissing`.
  * `:tipsnodes` columns corresponding to internal nodes, and row to tips (works only is indexation="b")
  * `:all` nodes, listed in topological order
- `indtips::Vector{Int}`: optional argument precising a specific order for the tips (internal use).
- `nonmissing::BitVector`: optional argument saying which tips have data (internal use).
   Tips with missing data are treated as internal nodes.
"""
function Base.getindex(
    obj::MatrixTopologicalOrder,
    d::Symbol,
    indtips::Vector{Int}=collect(1:length(obj.tipnumbers)),
    nonmissing::BitVector=trues(length(obj.tipnumbers))
)
    d == :all && return obj.V
    # otherwise, do extra work
    tipnums = obj.tipnumbers[indtips][nonmissing]
    maskTips = indexin(tipnums, obj.nodenumbers_toporder)
    if d == :tips # Extract rows and/or columns corresponding to the tips with data
        obj.indexation == :b && return obj.V[maskTips, maskTips] # both columns and rows are indexed by nodes
        obj.indexation == :c && return obj.V[:, maskTips] # Only the columns
        obj.indexation == :r && return obj.V[maskTips, :] # Only the rows
    end
    intnodenums = [obj.internalnodenumbers ; setdiff(obj.tipnumbers, tipnums)]
    maskNodes = indexin(intnodenums, obj.nodenumbers_toporder)
    #= indices in obj.nodenumbers_toporder, in this order:
    1. internal nodes, in the same order as in obj.internalnodenumbers,
       that is, same order as in net.node (excluding leaves)
    2. tips absent from indtips or missing data according to nonmissing,
       in the same order as in obj.tipnumbers.
    =#
    if d == :internalnodes # Idem, for internal nodes
        obj.indexation == :b && return obj.V[maskNodes, maskNodes]
        obj.indexation == :c && return obj.V[:, maskNodes]
        obj.indexation == :r && return obj.V[maskNodes, :]
    end
    if d == :tipsnodes
        obj.indexation == :b && return obj.V[maskTips, maskNodes]
        obj.indexation == :c && error("""Both rows and columns must be net
                                       ordered to take the submatrix tips vs internal nodes.""")
        obj.indexation == :r && error("""Both rows and columns must be net
                                       ordered to take the submatrix tips vs internal nodes.""")
    end
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

This function internally calls [`sharedpathmatrix`](@ref), which computes the variance
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
    V = sharedpathmatrix(net; checkpreorder=checkpreorder)
    C = V[:tips]
    corr && StatsBase.cov2cor!(C, sqrt.(diag(C)))
    Cd = DataFrame(C, map(Symbol, V.tipnames))
    return(Cd)
end


"""
    sharedpathmatrix(net::HybridNetwork; checkpreorder::Bool=true)

Matrix Ω of shared path lengths from the root to a pair of nodes, for all pairs
of nodes in network `net`. This is the covariance under a Brownian motion model,
returned as a data frame by [`vcv`](@ref).
The network is assumed to be pre-ordered if `checkpreorder` is false.
If `checkpreorder` is true (default), `preorder!` is run on the network beforehand.

Returns an object of type [`MatrixTopologicalOrder`](@ref).

See also [`descendenceweight`](@ref), which returns a matrix P independent
of edge lengths. P only depends on the network topology and γ.
Then Ω = PDP' where D is a diagonal matrix with the network's edge lengths
on its diagonal.
"""
function sharedpathmatrix(net::HybridNetwork; checkpreorder::Bool=true)
    check_nonmissing_nonnegative_edgelengths(net,
        """The variance-covariance matrix of the network is not defined.
           A phylogenetic regression cannot be done.""")
    check_valid_gammas(net, "The variance-covariance matrix is not defined.")
    checkpreorder && preorder!(net)
    V = traversal_preorder(
            net.vec_node,
            init_sharedpathmatrix,
            traversalupdate_default!,
            updatetree_sharedpathmatrix!,
            updatehybrid_sharedpathmatrix!,
        )
    M = MatrixTopologicalOrder(V, net, :b) # nodes in both columns & rows
    return M
end

function init_sharedpathmatrix(nodes::Vector{Node},)
    n = length(nodes)
    return(zeros(Float64,n,n))
end

function updatetree_sharedpathmatrix!(
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

function updatehybrid_sharedpathmatrix!(
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
    descendencematrix(net::HybridNetwork; checkpreorder::Bool=true)

Descendence matrix between all the nodes of a network:
object `D` of type [`MatrixTopologicalOrder`](@ref) in which
`D[i,j]` is the proportion of genetic material in node `i` that can be traced
back to node `j`. If `D[i,j]>0` then `i` is a descendent of `j` (and `j` is
an ancestor of `i`).
The network is assumed to be pre-ordered if `checkpreorder` is false.
If `checkpreorder` is true (default), `preorder!` is run on the network beforehand.
"""
function descendencematrix(net::HybridNetwork; checkpreorder::Bool=true)
    checkpreorder && preorder!(net)
    check_valid_gammas(net, "The descendence matrix is not defined.")
    V = traversal_postorder(
        net.vec_node,
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

"""
    descendenceweight(net::HybridNetwork; checkpreorder::Bool=true)

Matrix of weights of all paths starting from an edge and ending to a node,
in network `net`. By default, `preorder!` is called first, to calculate
a pre-ordering of nodes. If `checkpreorder` is `false` instead, then
this pre-ordering is assumed to have been done already.

output: object `M` of type [`MatrixTopologicalOrder`](@ref), in which
nodes are in rows and edges are in columns, with edge numbers used as indices.
In other words, `M[i,j]` is the sum of weights of all paths from edge number `j`
to node `net.vec_node[i]`.

warning: `net`'s edge numbers are reset to be consecutive numbers from 1 to the
total number of edges. If any edge number needs to be modified for this,
a warning is sent, via [`resetedgenumbers!`](@ref).

See also [`descendencematrix`](@ref), which returns the matrix of
path weights from each *node* to any other node.

# example

```jldoctest descendenceweight
julia> net = readnewick("((t4:1.5,((t3:1.27,t1:1.27):1.06,#H8:0.0::0.4):1.18):1.16,(t2:1.32)#H8:1.34::0.6);");

julia> m = PhyloNetworks.descendenceweight(net);

julia> m[:tips]
4×9 Matrix{Float64}:
 1.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
 0.0  1.0  0.0  1.0  0.0  1.0  1.0  0.0  0.0
 0.0  0.0  1.0  1.0  0.0  1.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.4  0.4  0.4  1.0  0.6

julia> tiplabels(net) # order of tips along rows
4-element Vector{String}:
 "t4"
 "t3"
 "t1"
 "t2"
```

This matrix, which contains information about the network topology and
inheritance weights γ, can be combined with the vector of branch lengths
to obtain the variance-covariance between tips under a Brownian motion model:

```jldoctest descendenceweight
julia> d = Dict(e.number => e.length for e in net.edge);

julia> edgelengths = [d[i] for i in 1:length(net.edge)];

julia> using LinearAlgebra

julia> Ω = m[:tips] * Diagonal(edgelengths) * transpose(m[:tips]) # covariance
4×4 Matrix{Float64}:
 2.66   1.16   1.16   0.464
 1.16   4.67   3.4    0.936
 1.16   3.4    4.67   0.936
 0.464  0.936  0.936  2.1768

julia> Matrix(vcv(net)) ≈ Ω # same as from `vcv` function
true
```
"""
function descendenceweight(net::HybridNetwork; checkpreorder::Bool=true)
    check_valid_gammas(net, "The descendence weight matrix is not defined.")
    resetedgenumbers!(net)
    checkpreorder && preorder!(net)
    V = traversal_preorder(
            net.vec_node,
            init_descendenceweight,
            traversalupdate_default!,
            updatetree_descendenceweight!,
            updatehybrid_descendenceweight!,
            length(net.edge)
        )
    M = MatrixTopologicalOrder(V, net, :r) # nodes in rows
    return M
end

function init_descendenceweight(nodes::Vector{Node}, ne::Integer)
    return zeros(Float64, length(nodes), ne) # nodes in rows, edge in columns
end

function updatetree_descendenceweight!(
    V::Matrix,
    i::Int,
    parentIndex::Int,
    edge::Edge,
    args...
)
    V[i,:] .= V[parentIndex,:] # useless for later edges V[i,j] = 0 already, but unknown j's
    V[i,edge.number] = 1 # = edge.gamma
    return true
end

function updatehybrid_descendenceweight!(
    V::Matrix,
    i::Int,
    parindx::AbstractVector{Int},
    paredge::AbstractVector{Edge},
    args...
)
    for (pi,pe) in zip(parindx, paredge)
        V[i,:] += pe.gamma * V[pi,:]
    end
    for pe in paredge
        V[i,pe.number] =  pe.gamma
    end
    return true
end

"""
    numberpathsmatrix(net::HybridNetwork; preorder::Bool=true)

Matrix `M` containing the number of up-down paths between each pair of nodes
in `net`, as a [`MatrixTopologicalOrder`](@ref) object, in which nodes are in
both rows and columns. By default, `preorder!` is called first, to calculate
a pre-ordering of nodes, unless `preorder` is `false`.

An up-down path between node n1 and n2 is a path that starts from n1, goes up
to a common ancestor then down to n2 (and never using the same edge twice).
This set of paths remain the same if the network is rerooted, so it is
well-defined on semidirected networks
(see [Xu & Ané 2023](https://doi.org/10.1007/s00285-022-01847-8)).
If the input network is a tree, then all path numbers are 1.

# example

```jldoctest
julia> net = readnewick("((t4,((t3,#H0),#H1)),(((t1)#H0,t2))#H1);");

julia> m = PhyloNetworks.numberpathsmatrix(net);

julia> m[:tips]
4×4 Matrix{Int64}:
 1  1  3  2
 1  1  3  2
 3  3  1  3
 2  2  3  1

julia> print(tiplabels(net)) # order of tips along rows
["t4", "t3", "t1", "t2"]
```
"""
function numberpathsmatrix(net::HybridNetwork; preorder::Bool=true)
    preorder && preorder!(net)
    V = traversal_preorder(
            net.vec_node,
            init_intmatrix, # matrix of Int, undef
            fillwith_ones!,
            updatetree_numberpathsmatrix!,
            updatehybrid_numberpathsmatrix!,
        )
    M = MatrixTopologicalOrder(V, net, :b) # nodes in both columns & rows
    return M
end
function updatetree_numberpathsmatrix!(V::Matrix,i::Int,pari::Int,::Edge)
    for j in 1:(i-1)
        V[i,j] = V[j,pari]
        V[j,i] = V[j,pari]
    end
    return true
end
function updatehybrid_numberpathsmatrix!(
    V::Matrix,
    i::Int,
    pari::AbstractVector{Int},
    ::AbstractVector{Edge},
)
    for j in 1:(i-1)
        V[i,j] = sum(ip -> V[ip, j], pari)
        V[j,i] = V[i,j]
    end
    return true
end
