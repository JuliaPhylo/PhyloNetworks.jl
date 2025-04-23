
"""
    quartetrank(t1,t2,t3,t4, nCk::Matrix)
    quartetrank([t1,t2,t3,t4], nCk)

Return the rank of a four-taxon set with taxon numbers `t1,t2,t3,t4`,
assuming that `ti`s are positive integers such that t1<t2, t2<t3 and t3<t4
(assumptions not checked!).
`nCk` should be a matrix of "n choose k" binomial coefficients:
see [`nchoose1234`](@ref).

# examples

```jldoctest
julia> nCk = PhyloNetworks.nchoose1234(5)
6×4 Matrix{Int64}:
 0   0   0  0
 1   0   0  0
 2   1   0  0
 3   3   1  0
 4   6   4  1
 5  10  10  5

julia> PhyloNetworks.quartetrank([1,2,3,4], nCk)
1

julia> PhyloNetworks.quartetrank([3,4,5,6], nCk)
15
```
"""
@inline function quartetrank(tnum::AbstractVector, nCk::Matrix)
    quartetrank(tnum..., nCk)
end
@inline function quartetrank(t1::Int, t2::Int, t3::Int, t4::Int, nCk::Matrix)
    # rank-1 = t1-1 choose 1 + t2-1 choose 2 + t3-1 choose 3 + t4-1 choose 4
    return nCk[t1,1] + nCk[t2,2] + nCk[t3,3] + nCk[t4,4] + 1
end


"""
    nchoose1234(nmax)

`nmax+1 x 4` matrix containing the binomial coefficient
"n choose k" in row `n+1` and column `k`. In other words,
`M[i,k]` gives "i-1 choose k". It is useful to store these
values and look them up to rank (a large number of) 4-taxon sets:
see [`quartetrank`](@ref).
"""
function nchoose1234(nmax::Int)
    # compute nC1, nC2, nC3, nC4 for n in [0, nmax]: used for ranking quartets
    M = Matrix{Int}(undef, nmax+1, 4)
    for i in 1:(nmax+1)
        M[i,1] = i-1 # n choose 1 = n. row i is for n=i-1
    end
    M[1,2:4] .= 0 # 0 choose 2,3,4 = 0
    for i in 2:(nmax+1)
        for k in 2:4 # to choose k items in 1..n: the largest could be n, else <= n-1
            M[i,k] = M[i-1,k-1] + M[i-1,k]
        end
    end
    return M
end


"""
    tablequartetCF(quartetlist::Vector{QuartetT} [, taxonnames];
                   keepQwithoutgenes=true, colnames=nothing)

Convert a vector of [`QuartetT`](@ref) objects to a table with 1 row for
each four-taxon set in the list. Each four-taxon set contains quartet data of
some type `T`, which determines the number of columns in the table.
This data type `T` should be a vector of length 3 or 4, or a 3×n matrix.
The output is a `NamedTuple`, to which we can apply common table operations
as a [`Table.jl`](https://github.com/JuliaData/Tables.jl)-compatible source.
It can easily be converted to other table formats (e.g. these packages
[integrate with Tables.jl](https://github.com/JuliaData/Tables.jl/blob/master/INTEGRATIONS.md))
such as a `DataFrame`.

In the output table, the columns are, in this order:
- `qind`: contains the quartet's `number`
- `t1, t2, t3, t4`: contain the quartet's `taxonnumber`s if no `taxonnames`
  are given, or the taxon names otherwise. The name of taxon number `i` is
  taken to be `taxonnames[i]`.
- 3 or more columns for the quartet's `data`: see [`quartetdata_columnnames`](@ref).

In short, the first 3 columns of data are named `CF12_34, CF13_24, CF14_23`.

If quartets have 4 data entries, then the 4th column is named `ngenes`,
and a quartet with a value `ngenes` of 0 is skipped (excluded) from the table,
unless `keepQwithoutgenes=true` (which is the default: 4-taxon sets are kept
by default even without any informative genes).

If quartets have a data matrix with 3 rows and `d` columns, then the
columns 4,5,6 in the table are named `V2_12_34, V2_13_24, V2_14_23`
and contain the data in the second column of the quartet's data matrix.
And so on.

For the table to have non-default column names, provide the desired
3, 4, or 3×d names as a vector via the optional argument `colnames`.

See [`countquartetsintrees`](@ref) for examples.
"""
function tablequartetCF(
    quartets::Vector{QuartetT{T}},
    taxa::AbstractVector{<:AbstractString}=Vector{String}();
    keepQwithoutgenes::Bool=true,
    colnames=nothing,
) where T <: Union{StaticVector{3}, StaticVector{4}, StaticMatrix{3,N} where N}
    return tablequartetdata(quartets, taxa; keepQwithoutgenes=keepQwithoutgenes,
            colnames=colnames, prefix="CF")
end

"""
    tablequartetdata(quartetlist::Vector{QuartetT} [, taxonnames];
                     keepQwithoutgenes=true, colnames=nothing, prefix="")

Convert a vector of [`QuartetT`](@ref) objects to a table with 1 row for
each four-taxon set in the list. Each four-taxon set contains quartet data of
some type `T`, which determines the number of columns in the table.
This data type `T` should be a vector of length 3 or 4, or a 3×n matrix.
The output is a `NamedTuple`, to which we can apply common table operations
as a [`Table.jl`](https://github.com/JuliaData/Tables.jl)-compatible source.
It can easily be converted to other table formats (e.g. these packages
[integrate with Tables.jl](https://github.com/JuliaData/Tables.jl/blob/master/INTEGRATIONS.md))
such as a `DataFrame`.

In the output table, the columns are, in this order:
- `qind`: contains the quartet's `number`
- `t1, t2, t3, t4`: contain the quartet's `taxonnumber`s if no `taxonnames`
  are given, or the taxon names otherwise. The name of taxon number `i` is
  taken to be `taxonnames[i]`.
- 3 or more columns for the quartet's `data`: see [`quartetdata_columnnames`](@ref).

In short, the first 3 columns of data are named `12_34, 13_24, 14_23`.

If quartets have 4 data entries, then the 4th column is named `ngenes`,
and a quartet with a value `ngenes` of 0 is skipped (excluded) from the table,
unless `keepQwithoutgenes=true` (which is the default: 4-taxon sets are kept
by default even without any informative genes).

If quartets have a data matrix with 3 rows and `d` columns, then the
columns 4,5,6 in the table are named `V2_12_34, V2_13_24, V2_14_23`
and contain the data in the second column of the quartet's data matrix.
And so on.

For the table to have non-default column names, provide the desired
3, 4, or 3×d names as a vector via the optional argument `colnames`.

See [`countquartetsintrees`](@ref) for examples.
"""
function tablequartetdata(
    quartets::Vector{QuartetT{T}},
    taxa::AbstractVector{<:AbstractString}=Vector{String}();
    keepQwithoutgenes::Bool=true,
    colnames=nothing,
    prefix="",
) where T <: Union{StaticVector{3}, StaticVector{4}, StaticMatrix{3,N} where N}
    V = eltype(T)
    colnames_data = quartetdata_columnnames(T,prefix)
    if !isnothing(colnames)
        if length(colnames) == length(colnames_data)
            colnames_data = colnames
        else
          @error "'colnames' needs to be of length $(length(colnames_data)).\nwill use default column names."
        end
    end
    translate = !isempty(taxa)
    tnT = (translate ? eltype(taxa) : Int) # Type for taxon names
    if translate
        taxstring = x -> taxa[x] # will error if a taxonnumber > length(taxa): okay
    else
        taxstring = x -> x
    end
    tupletxn = [:t1, :t2, :t3, :t4]   # txn = taxon names
    tupledat = Symbol.(colnames_data) # dat = data names
    tuplenames = [:qind, tupletxn..., tupledat...]
    tuplevecs = [Int[], tnT[], tnT[], tnT[], tnT[]]
    for _ in eachindex(colnames_data) push!(tuplevecs, V[]); end
    nt = (; zip(tuplenames,tuplevecs)...) # no copy
    indngenes = findfirst(isequal(:ngenes), tupledat)
    hasngenes = !keepQwithoutgenes && !isnothing(indngenes)
    for q in quartets
        hasngenes && q.data[indngenes] == 0 && continue # skip quartets with 0 genes
        push!(nt[:qind], q.number)
        for (i,ti) in enumerate(tupletxn)
            push!(nt[ti], taxstring(q.taxonnumber[i]))
        end
        for (j,dj) in enumerate(tupledat)
            push!(nt[dj], q.data[j])
        end
    end
    return nt
end


"""
    quartetdata_columnnames(T, prefix) where T <: StaticArray

Vector of column names to hold the quartet data of type `T` in a table.
If T is a length-3 vector type, they are
"`prefix`12_34","`prefix`13_24","`prefix`14_23".
For example, with `prefix` "CF", these first 3 column names are
"CF12_34","CF13_24","CF14_23".
If T is a length-4 vector type, the 4th name is "ngenes".
If T is a 3×n matrix type, the output vector contains 3×n names.
Each group of 3 is as desribed above, with prefix `prefix` for the first 3,
and prefixes "V2_", "V3_", ... "Vn_" the following 3(n-1) names.

Used by [`tablequartetdata`](@ref) and [`tablequartetCF`](@ref)
to build a table from a vector of [`QuartetT`](@ref) objects.
"""
function quartetdata_columnnames(::Type{T}, prefix) where
    T <: StaticArray{Tuple{3},S,1} where S
    return ["$(prefix)12_34","$(prefix)13_24","$(prefix)14_23"]
end
function quartetdata_columnnames(::Type{T}, prefix) where
    T <: StaticArray{Tuple{4},S,1} where S
    return ["$(prefix)12_34","$(prefix)13_24","$(prefix)14_23","ngenes"]
end
function quartetdata_columnnames(::Type{T}, prefix) where
    T <: StaticArray{Tuple{3,N},S,2} where {N,S}
    # for a 3×N matrix: N names
    N > 0 || error("expected at least 1 column of data")
    colnames_q = ["12_34","13_24","14_23"]
    colnames = prefix .* colnames_q
    for i in 2:N append!(colnames, "V$(i)_" .* colnames_q); end
    return colnames
end



"""
    countquartetsintrees(trees [, taxonmap]; which=:all, weight_byallele=true)

Calculate the quartet concordance factors (CF) observed in the `trees` vector.
If present, `taxonmap` should be a dictionary that maps each allele name to it's species name.
To save to a file, first convert to a data frame using [`tablequartetCF`](@ref).
When `which=:all`, quartet CFs are calculated for all 4-taxon sets.
(Other options are not implemented yet.)

The algorithm runs in O(mn⁴) where m is the number of trees and n is the number
of tips in the trees.

CFs are calculated at the species level only, that is, considering 4-taxon sets
made of 4 distinct species, even if the gene trees may have multiple alleles
from the same species. For 4 distinct species `a,b,c,d`, all alleles from
each species (`a` etc.) will be considered to calculate the quartet CF.

By default, each gene has a weight of 1. So if there are `n_a` alleles from `a`,
`n_b` alleles from `b` etc. in a given gene, then each set of 4 alleles has a
weight of `1/(n_a n_b b_c n_c)` in the calculation of the CF for `a,b,c,d`.
With option `weight_byallele=true`, then each set of 4 alleles is given a
weight of 1 instead. This inflates the total number of sets used to calculate
the quartet CFs (to something larger than the number of genes). This may also
affect the CF values if the number of alleles varies across genes: genes with
more alleles will be given more weight.

# examples

```jldoctest quartet
julia> tree1 = readnewick("(E,(A,B),(C,D),O);"); tree2 = readnewick("(((A,B),(C,D)),E);");

julia> q,t = countquartetsintrees([tree1, tree2]);
Reading in trees, looking at 15 quartets in each...
0+--+100%
  **

julia> t # taxon order: t[i] = name of taxon number i
6-element Vector{String}:
 "A"
 "B"
 "C"
 "D"
 "E"
 "O"

julia> length(q) # 15 four-taxon sets on 6 taxa
15

julia> q[1] # both trees agree on AB|CD: resolution 1
4-taxon set number 1; taxon numbers: 1,2,3,4
data: [1.0, 0.0, 0.0, 2.0]

julia> q[8] # tree 2 is missing O (taxon 6), tree 1 wants resolution 3: AO|CD
4-taxon set number 8; taxon numbers: 1,3,4,6
data: [0.0, 0.0, 1.0, 1.0]

julia> q[11] # tree 1 has ACEO unresolved, and tree 2 is missing O: no data for this quartet
4-taxon set number 11; taxon numbers: 1,3,5,6
data: [0.0, 0.0, 0.0, 0.0]
```

In the next example, each tree has 2 individuals from population A.

```jldoctest quartet
julia> tree1 = readnewick("(E,(a1,B),(a2,D),O);"); tree2 = readnewick("(((a1,a2),(B,D)),E);");

julia> q,t = countquartetsintrees([tree1, tree2], Dict("a1"=>"A", "a2"=>"A"); showprogressbar=false);

julia> t
5-element Vector{String}:
 "A"
 "B"
 "D"
 "E"
 "O"

julia> q[1] # tree 1 has discordance: a1B|DE and a2D|BE. tree 2 has AE|BD for both alleles of A
4-taxon set number 1; taxon numbers: 1,2,3,4
data: [0.25, 0.25, 0.5, 2.0]

julia> q[3] # tree 2 is missing O (taxon 5), and a2 is unresolved in tree 1. There's only a1B|EO
4-taxon set number 3; taxon numbers: 1,2,4,5
data: [1.0, 0.0, 0.0, 0.5]
```

Next we show how to convert these objects to a table using [`tablequartetCF`](@ref).
The output is a `NamedTuple`. It can be saved later to a `DataFrame` for example,
using option `copycols=false` to avoid copying the columns (which can be very
large if there are many 4-taxon sets). Data frames are easier to visualize,
filter etc., but performance can be better on named tuples.

```jldoctest quartet
julia> nt = tablequartetCF(q,t); # named tuple

julia> using DataFrames

julia> df = DataFrame(nt, copycols=false); # convert to a data frame, without copying the column data

julia> show(df, allcols=true) # data frames are displayed much more nicely than named tuples
5×9 DataFrame
 Row │ qind   t1      t2      t3      t4      CF12_34  CF13_24  CF14_23  ngenes  
     │ Int64  String  String  String  String  Float64  Float64  Float64  Float64 
─────┼───────────────────────────────────────────────────────────────────────────
   1 │     1  A       B       D       E          0.25     0.25      0.5      2.0
   2 │     2  A       B       D       O          0.5      0.5       0.0      1.0
   3 │     3  A       B       E       O          1.0      0.0       0.0      0.5
   4 │     4  A       D       E       O          1.0      0.0       0.0      0.5
   5 │     5  B       D       E       O          0.0      0.0       0.0      0.0

julia> # using CSV; CSV.write(nt, "filename.csv");
```

Note that `CSV.write` can take a data frame or a named tuple as input, to write
the table to a file.

Finally, the example below shows the effect of using `weight_byallele=true` when
caculating quartet concordance factors from gene trees with multiple alleles per
population, and of filtering out 4-taxon sets with no informative genes with
`keepQwithoutgenes=false` when converting to a table.

```jldoctest quartet
julia> tree2 = readnewick("((A,(B,D)),E);");

julia> q,t = countquartetsintrees([tree1, tree2], Dict("a1"=>"A", "a2"=>"A"); weight_byallele=true);
Reading in trees, looking at 5 quartets in each...
0+--+100%
  **

julia> nt = tablequartetCF(q,t; keepQwithoutgenes=false); # qind=5 excluded: 0 genes

julia> show(DataFrame(nt, copycols=false), allcols=true)
4×9 DataFrame
 Row │ qind   t1      t2      t3      t4      CF12_34   CF13_24   CF14_23   ngenes  
     │ Int64  String  String  String  String  Float64   Float64   Float64   Float64 
─────┼──────────────────────────────────────────────────────────────────────────────
   1 │     1  A       B       D       E       0.333333  0.333333  0.333333      3.0
   2 │     2  A       B       D       O       0.5       0.5       0.0           2.0
   3 │     3  A       B       E       O       1.0       0.0       0.0           1.0
   4 │     4  A       D       E       O       1.0       0.0       0.0           1.0
```
"""
function countquartetsintrees(
    tree::Vector{HybridNetwork},
    taxonmap::Dict=Dict{String,String}();
    whichQ::Symbol=:all, weight_byallele::Bool=false,
    showprogressbar::Bool=true
)
    whichQ in [:all, :intrees] || error("whichQ must be either :all or :intrees, but got $whichQ")
    if isempty(taxonmap)
        taxa = tiplabels(tree)
    else
        taxa = sort!(collect(Set(haskey(taxonmap, l.name) ? taxonmap[l.name] : l.name
                                 for t in tree for l in t.leaf)))
    end
    taxonnumber = Dict(taxa[i] => i for i in eachindex(taxa))
    ntax = length(taxa)
    nCk = nchoose1234(ntax) # matrix used to ranks 4-taxon sets
    qtype = MVector{4,Float64} # 4 floats: CF12_34, CF13_24, CF14_23, ngenes; initialized at 0.0
    if whichQ == :all
        numq = nCk[ntax+1,4]
        quartet = Vector{QuartetT{qtype}}(undef, numq)
        ts = [1,2,3,4]
        for qi in 1:numq
            quartet[qi] = QuartetT(qi, SVector{4}(ts), MVector(0.,0.,0.,0.))
            # next: find the 4-taxon set with the next rank,
            #       faster than using the direct mapping function
            ind = findfirst(x -> x>1, diff(ts))
            if ind === nothing ind = 4; end
            ts[ind] += 1
            for j in 1:(ind-1)
                ts[j] = j
            end
        end
    else
        error("whichQ = :intrees not implemented yet")
        # fixit: read all the trees, go through each edge, check if the current quartet list covers the edge, if not: add a quartet
    end
    totalt = length(tree)
    if showprogressbar
        nstars = (totalt < 50 ? totalt : 50)
        ntrees_perstar = (totalt/nstars)
        println("Reading in trees, looking at $numq quartets in each...")
        print("0+" * "-"^nstars * "+100%\n  ")
        stars = 0
        nextstar = Integer(ceil(ntrees_perstar))
    end
    for i in 1:totalt # number of times each quartet resolution is seen in each tree
        countquartetsintrees!(quartet, tree[i], whichQ, weight_byallele, nCk, taxonnumber, taxonmap)
        if showprogressbar && i >= nextstar
            print("*")
            stars += 1
            nextstar = Integer(ceil((stars+1) * ntrees_perstar))
        end
    end
    showprogressbar && print("\n")
    # normalize counts to frequencies & number of genes
    for q in quartet
        d = q.data
        d[4] = d[1]+d[2]+d[3] # number of genes
        if d[4] > 0.0
            d[1:3] /= d[4]
        end # otherwise: no genes with data on this quartet (missing taxa or polytomy): leave all zeros (NaN if /0)
    end
    return quartet, taxa
end
function countquartetsintrees!(quartet::Vector{QuartetT{MVector{4,Float64}}},
            tree::HybridNetwork, whichQ::Symbol, weight_byallele::Bool, nCk::Matrix,
            taxonnumber::Dict{<:AbstractString,Int}, taxonmap::Dict{<:AbstractString,<:AbstractString})
    tree.numhybrids == 0 || error("input phylogenies must be trees")
    # next: reset node & edge numbers so that they can be used as indices: 1,2,3,...
    resetnodenumbers!(tree; checkpreorder=true, type=:postorder) # leaves first & post-order
    resetedgenumbers!(tree)
    # next: build list leaf number -> species ID, using the node name then taxon map
    nleaves = length(tree.leaf)
    nnodes  = length(tree.node)
    taxID = Vector{Int}(undef, nleaves)
    for n in tree.leaf
        taxID[n.number] = haskey(taxonmap, n.name) ? taxonnumber[taxonmap[n.name]] : taxonnumber[n.name]
    end
    # number of individuals from each species: needed to weigh the quartets at the individual level
    # weight of t1,t2,t3,t4: 1/(taxcount[t1]*taxcount[t2]*taxcount[t3]*taxcount[t4])
    taxcount = zeros(Int, length(taxonnumber))
    for ti in taxID taxcount[ti] += 1; end
    # next: build data structure to get descendant / ancestor clades
    below,above = ladderpartition(tree) # re-checks that node numbers can be used as indices, with leaves first
    # below[n][1:2 ]: left & clades below node number n
    # above[n][1:end]: grade of clades above node n
    for n in (nleaves+1):nnodes # loop through internal nodes indices only
        bn = below[n]
        an = above[n]
        for c1 in 2:length(bn) # c = child clade, loop over all pairs of child clades
          for t1 in bn[c1]     # pick 1 tip from each child clade
            s1 = taxID[t1]
            for c2 in 1:(c1-1)
              for t2 in bn[c2]
                s2 = taxID[t2]
                s1 != s2 || continue # skip quartets that have repeated species
                t12max = max(t1,t2)
                leftweight = 1/(taxcount[s1]*taxcount[s2])
                for p1 in 1:length(an) # p = parent clade
                  for t3i in 1:length(an[p1])
                    t3 = an[p1][t3i]
                    s3 = taxID[t3]
                    (s3 != s1 && s3 != s2) || continue
                    for t4i in 1:(t3i-1) # pick 2 distinct tips from the same parent clade
                        t4 = an[p1][t4i]
                        t3 > t12max || t4 > t12max || continue   # skip: would be counted twice otherwise
                        s4 = taxID[t4]
                        (s4 != s1 && s4 != s2 && s4 != s3) || continue
                        rank,res = quartetRankResolution(s1, s2, s3, s4, nCk)
                        weight = ( weight_byallele ? 1.0 : leftweight / (taxcount[s3]*taxcount[s4]) )
                        quartet[rank].data[res] += weight
                    end
                    for p2 in 1:(p1-1) # distinct parent clade: no risk of counting twice
                      for t4 in an[p2]
                        s4 = taxID[t4]
                        (s4 != s1 && s4 != s2 && s4 != s3) || continue
                        rank,res = quartetRankResolution(s1, s2, s3, s4, nCk)
                        weight = ( weight_byallele ? 1.0 : leftweight / (taxcount[s3]*taxcount[s4]) )
                        quartet[rank].data[res] += weight
                      end
                    end
                  end
                end
              end
            end
          end
        end
    end
end


function quartetRankResolution(t1::Int, t2::Int, t3::Int, t4::Int, nCk::Matrix)
    # quartet: t1 t2 | t3 t4, but indices have not yet been ordered: t1<t2, t3<t4, t1<min(t3,t4)
    if t3 > t4 # make t3 smallest of t3, t4
        (t3,t4) = (t4,t3)
    end
    if t1 > t2 # make t1 smallest of t1, t2
        (t1,t2) = (t2,t1)
    end
    if t1 > t3 # swap t1 with t3, t2 with t4 - makes t1 smallest
        (t1,t3) = (t3,t1)
        (t2,t4) = (t4,t2)
    end
    if t2 < t3 # t2 2nd smallest: order t1 < t2 < t3 < t4
        resolution = 1; # 12|34 after ordering indices
        rank = quartetrank(t1, t2, t3, t4, nCk)
    else # t3 2nd smallest
        if t2 < t4 # order t1 < t3 < t2 < t4
            resolution = 2; # 13|24 after ordering
            rank = quartetrank(t1, t3, t2, t4, nCk);
        else # order t1 < t3 < t4 < t2
            resolution = 3; # 14|23 after ordering
            rank = quartetrank(t1, t3, t4, t2, nCk);
        end
    end
    return rank, resolution
end
