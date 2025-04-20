"""
    expectedf2matrix(net::HybridNetwork; checkpreorder::Bool=true)

Matrix of f2 statistics expected from `net`. The rows and columns
correspond to the taxa in the network, in the same order as listed
by `tiplabels(net)`.

See for example [Lipson 2020](https://doi.org/10.1111/1755-0998.13230).

# example

```jldoctest
julia> net2 = readnewick("(O:5.5,(((E:1.5)#H1:2.5::0.7,((#H1:0,D:1.5):1.5,((C:1,B:1):1)#H2:1::0.6):1.0):1.0,(#H2:0,A:2):3):0.5);");

julia> f2 = PhyloNetworks.expectedf2matrix(net2)
6×6 Matrix{Float64}:
  0.0   9.95  11.0   9.56  9.56  11.0
  9.95  0.0    5.45  5.95  5.95   8.95
 11.0   5.45   0.0   6.16  6.16  10.0
  9.56  5.95   6.16  0.0   2.0    6.16
  9.56  5.95   6.16  2.0   0.0    6.16
 11.0   8.95  10.0   6.16  6.16   0.0

julia> tiplabels(net2) # order of taxa in rows and columns of f2 matrix above
6-element Vector{String}:
 "O"
 "E"
 "D"
 "C"
 "B"
 "A"
```
"""
function expectedf2matrix(
    net::HybridNetwork;
    checkpreorder::Bool=true,
)
    m = descendenceweight(net; checkpreorder=checkpreorder)
    P = m[:tips]
    # change to :all to return a MatrixTopologicalOrder instead
    #= idea: f2[i,j] = Ω[i,i] + Ω[j,j] -2 Ω[i,j] where
    Ω = m[:tips] * Diagonal(edgelengths) * transpose(m[:tips])
    is a linear function of edge lengths.
    express: f2[i,j] = (nn × nn × ne) edgeweight array * edgelength vector
    edgeweight[i,j,e] = (m[i,e] - m[j,e])^2
    =#
    nn,ne = size(P) # number of nodes, number of edges
    edgeweight = zeros(eltype(P), nn,nn,ne)
    for (I, p_ik) in pairs(P)
        i = I[1]; k = I[2] # k = index of edge e, also edge number
        for j in 1:(i-1)
            w = (P[j,k] - p_ik)^2
            edgeweight[j,i,k] = w
            edgeweight[i,j,k] = w
        end
    end
    f2mat = zeros(eltype(P), nn,nn)
    for e in net.edge
        f2mat .+= e.length .* edgeweight[:,:,e.number]
    end
    # M = MatrixTopologicalOrder(f2mat, net, :b) # nodes in both columns & rows
    return f2mat
end

# similarly: f3[x;i,j] = Ω[x,x] + Ω[i,j] - Ω[x,i] -2 Ω[x,j]
# and f4[i1,i2; i3,i4] =
