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

"""
    expectedf4table(net::HybridNetwork;
                    showprogressbar::Bool=true,
                    checkpreorder::Bool=true)

Calculate the f4 statistics expected from `net`. Output: `(q,t)` where
`t` is a list of taxa and `q` is a list of 4-taxon set objects of type
[`PhyloNetworks.QuartetT{datatype}`](@ref).
In each element of `q`, `taxonnumber` gives the indices in `taxa` of the 4 taxa
of interest; and `data` contains the 3 primary expected f4-statistics, for the
following 3 ordering of the 4 taxa:
`t1,t2|t3,t4`, `t1,t3|t2,t4` and `t1,t4|t2,t3`.
This output is similar to that of `PhyloNetworks.countquartetsintrees`,
with 4-taxon sets listed in the same order
(same output `t`, then same order of 4-taxon sets in `q`).

For background of f-statistics, see for example
[Lipson 2020](https://doi.org/10.1111/1755-0998.13230).
f4-statistics are linear combination of f2-statistics:

`f4[t1,t2|t3,t4] = (f2[t1,t4] + f2[t2,t3] - f2[t1,t3] - f2[t2,t4])/2`

Given a set of 4 taxa, there are 12 ways to order them, but there are only
2 "degrees of freedom" in the associated 12 f4-statistics, thanks to symmetries:
* f4[t2,t1|t3,t4] = - f4[t1,t2|t3,t4]
* f4[t3,t4|t1,t2] =   f4[t1,t2|t3,t4]
* f4[t1,t2|t3,t4] + f4[t1,t3|t2,t4] + f4[t1,t4|t2,t3] = 0

# example

The first example is a tree: on which some f4 values are expected to be 0.

Next, we use a network with 2 reticulations, which contains a "32 cycle".

```jldoctest
julia> net = readTopology("(D:1,((C:1,#H25:0):0.1,((((B1:10,B2:1):1.5,#H1:0):10.8,((A1:1,A2:1):0.001)#H1:0::0.5):0.5)#H25:0::0.501):1);");

julia> # using PhyloPlots; plot(net, showedgelength=true);

julia> q,t = expectedf4table(net);
```
"""
function expectedf4table(
    net::HybridNetwork;
    showprogressbar::Bool=true,
    checkpreorder::Bool=true
)
    f2mat = expectedf2matrix(net; checkpreorder=checkpreorder)
    # f4s are linear combinations of edge lengths: from f2s, or from
    # f4[i1,i2; i3,i4] = Ω[i1,i3] + Ω[i2,i4] - Ω[i1,i4] - Ω[i2,i3]
    taxa = sort!(tipLabels(net))
    # todo: also get permutation to sort
    taxonnumber = Dict(taxa[i] => i for i in eachindex(taxa))
    ntax = length(taxa)
    nCk = nchoose1234(ntax) # matrix to rank 4-taxon sets
    qtype = MVector{3,Float64} # 3 floats: f4_12_34, f4_13_24, f4_14_23
    numq = nCk[ntax+1,4]
    quartet = Vector{QuartetT{qtype}}(undef, numq)
    ts = [1,2,3,4]
    for qi in 1:numq
        quartet[qi] = QuartetT(qi, SVector{4}(ts), MVector(0.,0.,0.))
        # next: find the 4-taxon set with the next rank,
        #       faster than using the direct mapping function
        ind = findfirst(x -> x>1, diff(ts))
        if ind === nothing ind = 4; end
        ts[ind] += 1
        for j in 1:(ind-1)
            ts[j] = j
        end
    end
    if showprogressbar
        nstars = (numq < 50 ? numq : 50)
        nquarnets_perstar = (numq/nstars)
        println("Calculation of expected f4 statistics for $numq quartets...")
        print("0+" * "-"^nstars * "+100%\n  ")
        stars = 0
        nextstar = Integer(ceil(nquarnets_perstar))
    end
    for qi in 1:numq
        # todo: modify quartet[qi] with its 3 f4 statistics
        if showprogressbar && qi >= nextstar
            print("*")
            stars += 1
            nextstar = Integer(ceil((stars+1) * nquarnets_perstar))
        end
    end
    showprogressbar && print("\n")
    return quartet, taxa
end

# similarly: f3[x;i,j] = Ω[x,x] + Ω[i,j] - Ω[x,i] - Ω[x,j]
