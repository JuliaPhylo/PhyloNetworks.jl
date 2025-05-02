"""
    expectedf2matrix(net::HybridNetwork; preorder::Bool=true)

Matrix of f2 statistics expected from `net`, assuming that branch lengths in
`net` represent "f2 distance". When data are based on allele frequencies,
edge lengths measure drift units scaled by a variance factor p(1-p) depending
on the allele frequency p at the root.

The rows and columns correspond to the taxa in the network, in the same order
as listed by `tiplabels(net)`.

For background on f-statistics, see for example
[Patterson et al. 2012](https://doi.org/10.1534/genetics.112.145037) or
[Lipson 2020](https://doi.org/10.1111/1755-0998.13230).

# example

```jldoctest
julia> net2 = readnewick("(O:5.5,(((E:1.5)#H1:2.5::0.7,((#H1:0,D:1.5):1.5,((C:1,B:1):1)#H2:1::0.6):1.0):1.0,(#H2:0,A:2):3):0.5);");

julia> f2 = expectedf2matrix(net2)
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
    preorder::Bool=true,
)
    m = descendenceweight(net; checkpreorder=preorder)
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
                    preorder::Bool=true)

Calculate the f4 statistics expected from `net`, assuming that branch lengths in
`net` represent "f2 distance".
Output: `(q,t)` where `t` is a list of taxa and `q` is a list of 4-taxon set
objects of type [`PhyloNetworks.QuartetT{datatype}`](@ref).
In each element of `q`, `taxonnumber` gives the indices in `taxa` of the 4 taxa
of interest; and `data` contains the 3 primary expected f4-statistics,
for the following 3 ordering of the 4 taxa:

    t1,t2|t3,t4   t1,t3|t4,t2   t1,t4|t2,t3.

This output is similar to that of [`countquartetsintrees`](@ref),
with 4-taxon sets listed in the same alphabetical order
(same output `t`, then same order of 4-taxon sets in `q`).

For background on f-statistics, see for example
[Patterson et al. 2012](https://doi.org/10.1534/genetics.112.145037) and
[Lipson 2020](https://doi.org/10.1111/1755-0998.13230).

f4-statistics are linear combination of f2-statistics:

`f4[t1,t2|t3,t4] = (f2[t1,t4] + f2[t2,t3] - f2[t1,t3] - f2[t2,t4])/2`

Given a set of 4 taxa, there are 12 ways to order them, but there are only
2 "degrees of freedom" in the associated 12 f4-statistics, thanks to symmetries:
* f4[t2,t1|t3,t4] = - f4[t1,t2|t3,t4]
* f4[t3,t4|t1,t2] =   f4[t1,t2|t3,t4]
* f4[t1,t2|t3,t4] + f4[t1,t3|t4,t2] + f4[t1,t4|t2,t3] = 0

# example

The first example is a tree: on which some f4 values are expected to be 0.
```jldoctest
julia> net0 = readnewick("((D:0.6,((a1:.1,a2:.1):0.1,B:0.2):0.3),C:0.4);");

julia> # using PhyloPlots; plot(net0, showedgelength=true);

julia> f4,t = expectedf4table(net0);
Calculation of expected f4 for 5 4-taxon sets...
0+-----+100%
  *****

julia> show(t)
["B", "C", "D", "a1", "a2"]
julia> show(f4[1].taxonnumber) # taxa numbered 1-4 are: B,C,D,a1
[1, 2, 3, 4]
julia> for q in f4
         println(join(t[q.taxonnumber],",") * ": " * string(round.(q.data, digits=3)))
       end
B,C,D,a1: [-0.3, 0.3, 0.0]
B,C,D,a2: [-0.3, 0.3, 0.0]
B,C,a1,a2: [0.0, -0.1, 0.1]
B,D,a1,a2: [0.0, -0.1, 0.1]
C,D,a1,a2: [0.0, -0.4, 0.4]
```
The zeros correspond to splits in the tree:
Ba1|CD, Ba2|CD, BC|a1a2, BD|a1a2, CD|a1a2.
The other values correspond to the internal path length for each split.

Next, we use a network with 2 reticulations, each time between sister species
(resulting in 3-cycle blobs).

```jldoctest
julia> net = readnewick("(D:1,((C:1,#H25:0):0.1,
        ((((B1:10,B2:1):1.5,#H1:0):10.8,
        ((A1:1,A2:1):0.001)#H1:0::0.5):0.5)#H25:0::0.501):1);");

julia> # plot(net, showedgelength=true);

julia> f4,t = expectedf4table(net, showprogressbar=false);

julia> using DataFrames

julia> df = tablequartetf4(f4, t) |> DataFrame
15×8 DataFrame
 Row │ qind   t1      t2      t3      t4      f4_12_34  f4_13_42  f4_14_23 
     │ Int64  String  String  String  String  Float64   Float64   Float64  
─────┼─────────────────────────────────────────────────────────────────────
   1 │     1  A1      A2      B1      B2           0.0    -4.201     4.201
   2 │     2  A1      A2      B1      C            0.0     2.699    -2.699
   3 │     3  A1      A2      B2      C            0.0     2.699    -2.699
   4 │     4  A1      B1      B2      C           -6.9     6.9       0.0
   5 │     5  A2      B1      B2      C           -6.9     6.9       0.0
   6 │     6  A1      A2      B1      D            0.0     2.699    -2.699
   7 │     7  A1      A2      B2      D            0.0     2.699    -2.699
   8 │     8  A1      B1      B2      D           -6.9     6.9       0.0
   9 │     9  A2      B1      B2      D           -6.9     6.9       0.0
  10 │    10  A1      A2      C       D            0.0    -3.176     3.176
  11 │    11  A1      B1      C       D            0.0    -5.875     5.875
  12 │    12  A2      B1      C       D            0.0    -5.875     5.875
  13 │    13  A1      B2      C       D            0.0    -5.875     5.875
  14 │    14  A2      B2      C       D            0.0    -5.875     5.875
  15 │    15  B1      B2      C       D            0.0   -12.775    12.775
```
"""
function expectedf4table(
    net::HybridNetwork;
    showprogressbar::Bool=true,
    preorder::Bool=true
)
    f2div2 = expectedf2matrix(net; preorder=preorder) ./2
    # f4s are linear combinations of edge lengths: from f2s, or from
    # f4[i1,i2; i3,i4] = Ω[i1,i3] + Ω[i2,i4] - Ω[i1,i4] - Ω[i2,i3]
    taxa = tiplabels(net) # order in f2 matrix
    o = sortperm(taxa)
    taxa .= taxa[o]   # order to construct rows, for f4
    # todo: also get permutation to sort
    taxonnumber = Dict(taxa[i] => i for i in eachindex(taxa))
    ntax = length(taxa)
    nCk = nchoose1234(ntax) # matrix to rank 4-taxon sets
    qtype = MVector{3,Float64} # 3 floats: f4_12_34, f4_13_42, f4_14_23
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
        println("Calculation of expected f4 for $numq 4-taxon sets...")
        print("0+" * "-"^nstars * "+100%\n  ")
        stars = 0
        nextstar = Integer(ceil(nquarnets_perstar))
    end
    for qi in 1:numq # modify quartet[qi].data to contain the 3 f4-stat
        qu = quartet[qi]
        ts .= o[qu.taxonnumber] # indices in f2 matrix. reuse memory from ts
        # @debug "4-taxon set: $(taxa[qu.taxonnumber]), index in f2: $ts"
        qf4 = qu.data
        # f4[t1,t2|t3,t4] = (f2[t1,t4] + f2[t2,t3] - f2[t1,t3] - f2[t2,t4])/2
        s14_23 = f2div2[ts[1],ts[4]] + f2div2[ts[2],ts[3]]
        s13_24 = f2div2[ts[1],ts[3]] + f2div2[ts[2],ts[4]]
        s12_34 = f2div2[ts[1],ts[2]] + f2div2[ts[3],ts[4]]
        qf4[1] = s14_23 - s13_24
        qf4[2] = s12_34 - s14_23
        qf4[3] = s13_24 - s12_34
        if showprogressbar && qi >= nextstar
            print("*")
            stars += 1
            nextstar = Integer(ceil((stars+1) * nquarnets_perstar))
        end
    end
    showprogressbar && print("\n")
    return quartet, taxa
end

"""
    expectedf3matrix(net::HybridNetwork, reftaxon; preorder::Bool=true)

Matrix of f3 statistics using reference taxon `reftaxon`, expected from `net`
assuming that its branch lengths represent "f2 distance".
f3-statistics are linear combination of f2-statistics, using `o` to denote
the reference taxon (for "outgtroup"):

`f3[o|t1,t2] = (f2[o,t1] + f2[o,t2] - f2[o,o] - f2[t1,t2])/2`

See [`expectedf3matrix`](@ref) and [expectedf4table](@ref).
"""
function expectedf3matrix(
    net::HybridNetwork,
    reftaxon::AbstractString;
    preorder::Bool=true
)
    taxa = tiplabels(net)
    iref = findall(isequal(reftaxon), taxa)
    length(iref) == 1 ||
        error("reference taxon $(reftaxon) found $(length(iref)) times in network")
    i0 = iref[1]
    f2 = - expectedf2matrix(net; preorder=preorder) ./2
    # f3[x;i,j] = Ω[x,x] +  Ω[i,j] -  Ω[x,i] -  Ω[x,j] = f4[x,i;x,j]
    #       = (- f2[x,x] - f2[i,j] + f2[x,i] + f2[x,j])/2    and f2[x,x]=0
    # modify f2 in place, but do *not* touch f2[i0,:] or f2[:,i0]
    for i in axes(f2,1)
        i == i0 && continue
        for j in 1:(i-1)
            j == i0 && continue
            f2[i,j] -= f2[j,i0] + f2[i,i0]
            f2[j,i] = f2[i,j]
        end
    end
    for i in axes(f2,1)
        f2[i,i0] = 0.0
        f2[i0,i] = 0.0
        f2[i,i]  = 0.0 # instead of -0 from taking -f2/2
    end
    return f2
end
