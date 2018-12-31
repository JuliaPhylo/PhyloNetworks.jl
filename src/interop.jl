# helper functions used by PhyloPlot to export (then plot)
# a network as a phylo R object.

"""
    makemissing!(x::AbstractVector)

Turn to `missing` any element of `x` exactly equal to -1.0.
Used for branch lengths and γs. `x` needs to accept missing values.
If not, this can be done with `allowmissing(x)`.
"""
@inline function makemissing!(x::AbstractVector)
    for i in 1:length(x)
        if x[i] == -1.0
            x[i] = missing
        end
    end
end

"""
    majoredgematrix(net::HybridNetwork)

Matrix of major edges from `net` where edge[i,1] is the number of the
parent node of edge i and edge[i,2] is the number of the child node of edge i.
Assume `nodes_changed` was updated, to list nodes in pre-order.

# Examples

```jldoctest
julia> net = readTopology("(A,(B,(C,D)));");

julia> PhyloNetworks.resetNodeNumbers!(net);

julia> PhyloNetworks.majoredgematrix(net)
6×2 Array{Int64,2}:
 5  1
 5  6
 6  2
 6  7
 7  3
 7  4
```
"""
function majoredgematrix(net::HybridNetwork)
    edge = Matrix{Int}(undef, length(net.edge)-length(net.hybrid), 2) # major edges
    i = 1 #row index for edge matrix
    for n in net.nodes_changed # topological pre-order
        !n.leaf || continue # skip leaves: associate node with children edges
        for e in n.edge
            if e.node[e.isChild1 ? 2 : 1] == n && e.isMajor
                #exclude parent edge and minor hybrid edges
                edge[i,1] = n.number # parent
                edge[i,2] = e.node[e.isChild1 ? 1 : 2].number # child
                i += 1
            end
        end
    end
    return edge
end

"""
    majoredgelength(net::HybridNetwork)

Generate vector of edge lengths of major `net` edges organized in the same order
as the `edge` matrix created via `majoredgematrix`. Considers values of `-1.0` as
missing values, recognized as NA in `R`.
Output: vector allowing for missing values.

Assume `nodes_changed` was updated, to list nodes in pre-order.

# Examples

```jldoctest
julia> net = readTopology("(((A:3.1,(B:0.2)#H1:0.3::0.9),(C,#H1:0.3::0.1):1.1),D:0.7);");

julia> directEdges!(net); preorder!(net);

julia> PhyloNetworks.majoredgelength(net)
8-element Array{Union{Missing, Float64},1}:
  missing
 0.7     
  missing
 1.1     
  missing
 3.1     
 0.3     
 0.2     
```
""" #"
function majoredgelength(net::HybridNetwork)
    edgeLength = Array{Union{Float64,Missing}}(undef, length(net.edge)-length(net.hybrid))
    i=1
    for n in net.nodes_changed # topological pre-order
        if !n.leaf
            for e in n.edge # for parent major edge below
                if e.node[e.isChild1 ? 2 : 1] == n && e.isMajor
                    edgeLength[i] = e.length
                    i=i+1
                end
            end
        end
    end
    #-1.0 interpreted as missing
    makemissing!(edgeLength)
    return edgeLength
end

"""
    minorreticulationmatrix(net::HybridNetwork)

Matrix of integers, representing the minor hybrid edges in `net`.
edge[i,1] is the number of the parent node of the ith minor hybrid edge,
and edge[i,2] is the number of its child node.
Node numbers may be negative, unless they were modified by `resetNodeNumbers!`.
Assumes correct `isChild1` fields.

# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);");
julia> PhyloNetworks.minorreticulationmatrix(net)
1×2 Array{Int64,2}:
 -6  3
```
""" #"
function minorreticulationmatrix(net::HybridNetwork)
    reticulation = Matrix{Int}(undef, length(net.hybrid), 2) # initialize
    j = 1 # row index, row = reticulate edge
    for e in net.edge
        if !e.isMajor # minor (hybrid) edges only
            reticulation[j,1] = getParent(e).number
            reticulation[j,2] = getChild(e).number
            j += 1
        end
    end
    return reticulation
end

"""
    minorreticulationlength(net::HybridNetwork)

Vector of lengths for the minor hybrid edges, organized in the same order
as in the matrix created via `minorreticulationmatrix`.
Replace values of `-1.0` with missing values recognized by `R`.
Output: vector allowing for missing values.

# Examples

```jldoctest
julia> net = readTopology("(((A:3.1,(B:0.2)#H1:0.4::0.9),(C,#H1:0.3::0.1):1.1),D:0.7);");

julia> PhyloNetworks.minorreticulationlength(net)
1-element Array{Union{Missing, Float64},1}:
 0.3
```
""" #"
function minorreticulationlength(net::HybridNetwork)
    reticulationLength = Vector{Union{Float64,Missing}}(undef, 0) # initialize
    for e in net.edge
        if !e.isMajor #find minor hybrid edge
            push!(reticulationLength, e.length)
        end
    end
    makemissing!(reticulationLength)
    return reticulationLength
end

"""
    minorreticulationgamma(net::HybridNetwork)

Vector of minor edge gammas (inheritance probabilities) organized in the
same order as in the matrix created via `minorreticulationmatrix`.
Considers values of `-1.0` as missing values, recognized as NA in `R`.
Output: vector allowing for missing values.

# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);");

julia> PhyloNetworks.minorreticulationgamma(net)
1-element Array{Union{Float64, Missings.Missing},1}:
 0.1
```
 """ #"
function minorreticulationgamma(net::HybridNetwork)
    reticulationGamma = Vector{Union{Float64,Missing}}(undef, 0) #initialize
    for e in net.edge
        if !e.isMajor # minor hybrid edges only
            push!(reticulationGamma, e.gamma)
        end
    end
    makemissing!(reticulationGamma)
    return reticulationGamma
end
