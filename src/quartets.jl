
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
some type (vector of length 3 or 4, or 3×n matrix), which determines the
number of columns in the table (3, 4, or 3n).
The first 3 columns are named "CF12_34", "CF13_24", "CF14_23",
unless the full list of column names is provided in `colnames`.
See [`tablequartetdata`](@ref) for details.

The output is a `NamedTuple`, to which we can apply common table operations
as a [`Table.jl`](https://github.com/JuliaData/Tables.jl)-compatible source.
It can easily be converted to other table formats such as a `DataFrame`.

See [`countquartetsintrees`](@ref) for examples.
"""
function tablequartetCF(
    quartets::Vector{QuartetT{T}},
    taxa::AbstractVector{<:AbstractString}=Vector{String}();
    keepQwithoutgenes::Bool=true,
    colnames=nothing,
) where T <: Union{StaticVector{3}, StaticVector{4}, StaticMatrix{3,N} where N}
    return tablequartetdata(quartets, taxa; keepQwithoutgenes=keepQwithoutgenes,
            colnames=colnames, prefix="CF", basenames=["12_34","13_24","14_23"])
end

"""
    tablequartetf4(quartetlist::Vector{QuartetT} [, taxonnames];
                   keepQwithoutgenes=true, colnames=nothing)

Convert a vector of [`QuartetT`](@ref) objects to a table with 1 row for
each four-taxon set in the list. Each four-taxon set contains quartet data of
some type (vector of length 3 or 4, or 3×n matrix), which determines the
number of columns in the table (3, 4, or 3n).
The first 3 columns are named "f4_12_34", "f4_13_42" and "f4_14_23",
unless the full list of column names is provided in `colnames`.
See [`tablequartetdata`](@ref) for details.

The output is a `NamedTuple`, to which we can apply common table operations
as a [`Table.jl`](https://github.com/JuliaData/Tables.jl)-compatible source.
It can easily be converted to other table formats such as a `DataFrame`.

See [`countquartetsintrees`](@ref) for examples.
"""
function tablequartetf4(
    quartets::Vector{QuartetT{T}},
    taxa::AbstractVector{<:AbstractString}=Vector{String}();
    keepQwithoutgenes::Bool=true,
    colnames=nothing,
) where T <: Union{StaticVector{3}, StaticVector{4}, StaticMatrix{3,N} where N}
    return tablequartetdata(quartets, taxa; keepQwithoutgenes=keepQwithoutgenes,
            colnames=colnames, prefix="f4_", basenames=["12_34","13_42","14_23"])
end

"""
    tablequartetdata(quartetlist::Vector{QuartetT} [, taxonnames];
                     keepQwithoutgenes=true,
                     prefix="X", basenames=["1","2","3"], colnames=nothing)

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

In short, the first 3 columns of data are named "X1", "X2", "X3" by default.
If a `prefix` and/or `basenames` is provided, this prefix and/or base names
are used for these first 3 column names,
for example: "`CF12_34`", "`CF13_24`", "`CF14_23`" if `prefix="CF"` and
`basename=["12_34","13_24","14_23"]`.

If quartets have 4 data entries, then the 4th column is named `ngenes`,
and a quartet with a value `ngenes` of 0 is skipped (excluded) from the table,
unless `keepQwithoutgenes=true` (which is the default: 4-taxon sets are kept
by default even without any informative genes).

If quartets have a data matrix with 3 rows and `d` columns, then the
columns 4,5,6 in the table are named `V2_1, V2_2, V2_3`
and contain the data in the second column of the quartet's data matrix.
And so on.

For the table to have non-default column names, provide the desired
3, 4, or 3×d names as a vector via the optional argument `colnames`.
This argument takes precedence over `prefix`, that is, a prefix is ignored if
the full list of column names is given.

Used by [`tablequartetCF`](@ref) and [`tablequartetf4`](@ref).
See [`countquartetsintrees`](@ref) for examples.
"""
function tablequartetdata(
    quartets::Vector{QuartetT{T}},
    taxa::AbstractVector{<:AbstractString}=Vector{String}();
    keepQwithoutgenes::Bool=true,
    colnames=nothing,
    prefix="X",
    basenames=["1","2","3"],
) where T <: Union{StaticVector{3}, StaticVector{4}, StaticMatrix{3,N} where N}
    V = eltype(T)
    colnames_data = quartetdata_columnnames(T,prefix,basenames)
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
    quartetdata_columnnames(T, prefix, basenames) where T <: StaticArray

Vector of column names to hold the quartet data of type `T` in a table.
If `T` is a length-3 vector type, they are
"`prefix`1","`prefix`2","`prefix`3"
if `basenames=["1","2","3"]`, say, the default used in [`tablequartetdata`](@ref)
to build a table from a vector of [`QuartetT`](@ref) objects.

For example, [`tablequartetCF`](@ref) uses prefix "CF" and base names
"12_34","13_24","14_23", so the first 3 column names are:
"CF12_34","CF13_24","CF14_23".
[`tablequartetf4`](@ref) uses different base names.

If `T` is a length-4 vector type, the 4th name is "ngenes".

If `T` is a 3×n matrix type, the output vector contains 3×n names.
Each group of 3 is as desribed above, with prefix `prefix` for the first 3,
and prefixes "V2_", "V3_", ... "Vn_" the following 3(n-1) names.
"""
function quartetdata_columnnames(
    ::Type{T},
    prefix,
    basenames,
) where T <: StaticArray{Tuple{3},S,1} where S
    return prefix .* basenames
end
function quartetdata_columnnames(
    ::Type{T},
    prefix,
    basenames,
) where T <: StaticArray{Tuple{4},S,1} where S
    return [(prefix .* basenames)..., "ngenes"]
end
function quartetdata_columnnames(
    ::Type{T},
    prefix,
    basenames,
) where T <: StaticArray{Tuple{3,N},S,2} where {N,S}
    # for a 3×N matrix: 3N names
    N > 0 || error("expected at least 1 column of data")
    colnames = prefix .* basenames
    for i in 2:N append!(colnames, "V$(i)_" .* basenames); end
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

julia> df = DataFrame(nt, copycols=false); # convert to data frame, without copying column data

julia> show(df, allcols=true) # data frames displayed more nicely than named tuples
5×9 DataFrame
 Row │ qind   t1      t2      t3      t4      CF12_34  CF13_24  CF14_23  ngenes  
     │ Int64  String  String  String  String  Float64  Float64  Float64  Float64 
─────┼───────────────────────────────────────────────────────────────────────────
   1 │     1  A       B       D       E          0.25     0.25      0.5      2.0
   2 │     2  A       B       D       O          0.5      0.5       0.0      1.0
   3 │     3  A       B       E       O          1.0      0.0       0.0      0.5
   4 │     4  A       D       E       O          1.0      0.0       0.0      0.5
   5 │     5  B       D       E       O          0.0      0.0       0.0      0.0

julia> using CSV; CSV.write("quartetCFs_fromgenetrees.csv", nt); # save to file if needed
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

"""
    quartetdisplayprobability(net::HybridNetwork; showprogressbar=true)

Display probabilities of quartet trees, under a displayed-tree model
(see [Xu & Ané (2024)](https://doi.org/10.1007/s00285-022-01847-8) or
[Rhodes et al. (2025)](https://doi.org/10.1016/j.aam.2024.102804)).

A tree `T` is displayed in the network if it can be obtained by keeping all but
1 hybrid edge above each hybrid node. The display probability `γ(T)` of `T` is
the probability of choosing the edges that are kept to obtain `T`, that is:
the product of `γ(e)` over all edges `e` in `T`.
For a quartet tree `q = ab|cd` of 4 taxa, its display probability `γ(q)` is
the sum of `γ(T)` over all trees T that have `ab|cd`.

Polytomies in `net` can give positive probability to unresolved quartet trees,
in which case the probabilities of the 3 resolved quartet trees on a given
4-taxon set would sum to less than 1.

Output: `(q,t)` where `t` is a list of taxa,
and `q` is a list of 4-taxon set objects of type `PhyloNetworks.QuartetT{datatype}`.
In each element of `q`, `taxonnumber` gives the indices in `taxa`
of the 4 taxa of interest; and `data` contains the 3 quartet γs, for the
3 unrooted (resolved) topologies in the following order:
`t1,t2|t3,t4`, `t1,t3|t2,t4` and `t1,t4|t2,t3`.
This output is similar to that of `PhyloNetworks.countquartetsintrees` when
1 individual = 1 taxon, with 4-taxon sets listed in the same order
(same output `t`, then same order of 4-taxon sets in `q`).

Assumption: the network should have non-missing and valid inheritance γ values.

# example

```jldoctest
julia> net = readnewick("(((C,#H2),((((B1,B2,B3),#H1))#H2:::0.6,((A)#H3:::.8)#H1:::0.5)),(#H3,D));");

julia> # using PhyloPlots; plot(net, showgamma=true);

julia> q,t = PhyloNetworks.quartetdisplayprobability(net,showprogressbar=false);

julia> for qi in q
         println(join(t[qi.taxonnumber],",") * ": " *
            string(round.(qi.data, sigdigits=3)))
       end
A,B1,B2,B3: [0.0, 0.0, 0.0]
A,B1,B2,C: [0.0, 0.0, 1.0]
A,B1,B3,C: [0.0, 0.0, 1.0]
A,B2,B3,C: [0.0, 0.0, 1.0]
B1,B2,B3,C: [0.0, 0.0, 0.0]
A,B1,B2,D: [0.0, 0.0, 1.0]
A,B1,B3,D: [0.0, 0.0, 1.0]
A,B2,B3,D: [0.0, 0.0, 1.0]
B1,B2,B3,D: [0.0, 0.0, 0.0]
A,B1,C,D: [0.64, 0.0, 0.36]
A,B2,C,D: [0.64, 0.0, 0.36]
B1,B2,C,D: [1.0, 0.0, 0.0]
A,B3,C,D: [0.64, 0.0, 0.36]
B1,B3,C,D: [1.0, 0.0, 0.0]
B2,B3,C,D: [1.0, 0.0, 0.0]

julia> nt = PhyloNetworks.tablequartetdata(q,t, prefix="p", basenames=["12_34","13_24","14_23"]);

julia> using DataFrames; df = DataFrame(nt, copycols=false)
15×8 DataFrame
 Row │ qind   t1      t2      t3      t4      p12_34   p13_24   p14_23  
     │ Int64  String  String  String  String  Float64  Float64  Float64 
─────┼──────────────────────────────────────────────────────────────────
   1 │     1  A       B1      B2      B3         0.0       0.0     0.0
   2 │     2  A       B1      B2      C          0.0       0.0     1.0
   3 │     3  A       B1      B3      C          0.0       0.0     1.0
   4 │     4  A       B2      B3      C          0.0       0.0     1.0
   5 │     5  B1      B2      B3      C          0.0       0.0     0.0
   6 │     6  A       B1      B2      D          0.0       0.0     1.0
   7 │     7  A       B1      B3      D          0.0       0.0     1.0
   8 │     8  A       B2      B3      D          0.0       0.0     1.0
   9 │     9  B1      B2      B3      D          0.0       0.0     0.0
  10 │    10  A       B1      C       D          0.64      0.0     0.36
  11 │    11  A       B2      C       D          0.64      0.0     0.36
  12 │    12  B1      B2      C       D          1.0       0.0     0.0
  13 │    13  A       B3      C       D          0.64      0.0     0.36
  14 │    14  B1      B3      C       D          1.0       0.0     0.0
  15 │    15  B2      B3      C       D          1.0       0.0     0.0

```

From this output, we see that Bi,Bj are sister in all quartets. We also see
that (B1,B2,B3) is unresolved, from γ's not summing up to 1 for taxon sets
having all 3 Bs and another taxon.
The other 4-taxon sets are of the form {A,Bi,C,D} (rows 10,11,13).
For them only 2 trees are displayed: ABi|CD (with probability 0.64) and
AD|CBi (with probability 0.36). AC|BiD is *not* displayed.
These quarnets have circular order (A,Bi,C,D).
"""
function quartetdisplayprobability(
    net::HybridNetwork;
    showprogressbar=true,
)
    getroot(net).leaf && error("The root can't be a leaf.")
    check_valid_gammas(net,
        "inheritance probabilities (γ) are needed to calculate quartet display probabilities.")
    taxa = sort!(tiplabels(net))
    taxonnumber = Dict(taxa[i] => i for i in eachindex(taxa))
    ntax = length(taxa)
    nCk = nchoose1234(ntax) # matrix to rank 4-taxon sets
    qtype = MVector{3,Float64} # 3 floats: CF12_34, CF13_24, CF14_23; initialized at 0.0
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
        println("Calculation quartet CFs for $numq quartets...")
        print("0+" * "-"^nstars * "+100%\n  ")
        stars = 0
        nextstar = Integer(ceil(nquarnets_perstar))
    end
    for qi in 1:numq
        quartetdisplayprobability!(quartet[qi], net, taxa, taxonnumber)
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
    quartetdisplayprobability!(quartet::QuartetT, net::HybridNetwork, taxa, taxonnumber)

Update `quartet.data` to contain the quartet display probabilities
(from inheritance γ) of 4-taxon quartet trees from `net` for the 4-taxon set
`taxa[quartet.taxonnumber]`.
`taxa` should contain the tip labels in `net`.
`quartet.taxonnumber` gives the indices in `taxa` of the 4 taxa of interest.
`taxonnumber` should be a dictionary mapping taxon labels in to their indices
in `taxa`, for easier lookup.

`net` is not modified.
"""
function quartetdisplayprobability!(
    quartet::QuartetT{MVector{3,Float64}},
    net::HybridNetwork,
    taxa,
    taxonnumber,
)
    net = deepcopy(net)
    removedegree2nodes!(net)
    # delete all taxa except for the 4 in the quartet
    for taxon in taxa
        taxonnumber[taxon] in quartet.taxonnumber && continue
        deleteleaf!(net, taxon, simplify=true, unroot=false)
        # would like unroot=true but deleteleaf! throws an error when the root is connected to 2 outgoing hybrid edges
    end
    quartet.data .= quartetdisplayprobability!(net, taxa[quartet.taxonnumber])
    return quartet
end

"""
    quartetdisplayprobability!(net::HybridNetwork, fourtaxa)

Display probabilities, from inheritance γ in `net`, of the 3 four-taxon
resolved quartet trees. These 3 topologies are ordered following the
ordering of taxon names in `fourtaxa`, that is: if `fourtaxa` is a,b,c,d,
then the quartets are listed in this order:

    (γ(ab|cd), γ(ac|bd), γ(ad,bc))

A quartet not displayed by `net` has γ=0. The display probability of the
star tree (unresolved) is 1 minus the sum of these 3 γ probabilities.

Assumptions about `net`:
- has 4 taxa, and those are the same as `fourtaxa`
- no degree-2 nodes, except perhaps for the root
- hybrid edge γ's are non-missing

The network is modified as follows: what's above the LSA is removed,
the 2 edges incident to the root are fused (if the root is of degree 2),
external degree-2 blobs are removed, and 2-cycles are removed (retaining
their major edge). 3-cycles are *not* shrunk, because shrinking a 3-cycle
could remove (and hide) an unresolved displayed tree in non-binary networks.
These operations do not modify the set of displayed
unrooted quartets and their probabilities (from γ's).
Then `net` is simplified recursively
by removing hybrid edges for extracting displayed trees.
"""
function quartetdisplayprobability!(
    net::HybridNetwork,
    fourtaxa,
)
    deleteaboveLSA!(net)
    length(getroot(net).edge) <= 2 && fuseedgesat!(net.rooti, net)
    deleteexternal2blobs!(net)
    delete2cycles!(net, true) # unroot=true
    ndes = 4 # number of taxa descendant from lowest hybrid node
    if net.numhybrids > 0
        preorder!(net)
        # find a lowest hybrid node and # of taxa below it
        hyb = net.vec_node[findlast(n -> n.hybrid, net.vec_node)]
        funneledge = [e for e in hyb.edge if getparent(e) === hyb]
        ispolytomy = length(funneledge) > 1
        funneldescendants = union([descendants(e) for e in funneledge]...)
        ndes = length(funneldescendants)
        n2 = (ispolytomy ? hyb : getchild(funneledge[1]))
        ndes > 2 && n2.leaf && error("2+ descendants below the lowest hybrid, yet n2 is a leaf. taxa: $(fourtaxa)")
    end
    if ndes >= 2 # then 0 or 1 quartet only. find cut edge
        # pool of cut edges below. contains NO external edge, bc n2 not leaf (if reticulation), nice tree ow
        cutpool = (ndes == 2 ? (ispolytomy ? Edge[] : funneledge) :
            # otherwise ndes = 3 or 4: children edges of n2
            (net.numhybrids == 0 ? net.edge :
                [e for e in n2.edge if getparent(e) === n2]))
        filter!(e -> !getchild(e).leaf, cutpool)
        net.numhybrids > 0 || length(cutpool) <= 1 ||
            error("2+ cut edges, yet 4-taxon tree, degree-3 root and no degree-2 nodes. taxa: $(fourtaxa)")
        sistertofirst = 0 # initialize as if 3-way polytomy (no cut edge)
        for e in cutpool
            hwc = hardwiredcluster(e, fourtaxa)
            sistertofirst = findnext(x -> x == hwc[1], hwc, 2)
        end
        qγ = (sistertofirst == 0 ? MVector{3,Float64}(0.,0.,0.) :
             (sistertofirst == 2 ? MVector{3,Float64}(1.,0.,0.) :
             (sistertofirst == 3 ? MVector{3,Float64}(0.,1.,0.) :
                                   MVector{3,Float64}(0.,0.,1.)   )))
        return qγ
    end
    ndes > 0 || error("weird: hybrid node has no descendant taxa")
    # by now, there are 1 or 2 taxa below the lowest hybrid
    qγ = MVector{3,Float64}(0.,0.,0.) # mutated later
    parenthedge = [e for e in hyb.edge if getchild(e) === hyb]
    all(h.hybrid for h in parenthedge) || error("hybrid $(hyb.number) has a parent edge that's a tree edge")
    parenthnumber = [p.number for p in parenthedge]
    nhe = length(parenthedge)
    if ndes == 1 # weighted qγs average of the nhe (often = 2) displayed networks
        for i in 1:nhe # keep parenthedge[i], remove all others
            gamma = parenthedge[i].gamma
            simplernet = ( i < nhe ? deepcopy(net) : net ) # last case: to save memory allocation
            for j in 1:nhe
                j == i && continue # don't delete hybrid edge i!
                pe_index = findfirst(e -> e.number == parenthnumber[j], simplernet.edge)
                deletehybridedge!(simplernet, simplernet.edge[pe_index],
                    false,true,false,true,false) # ., unroot=true, ., simplify=true,.
            end
            qγ .+= gamma .* quartetdisplayprobability!(simplernet, fourtaxa)
        end
        return qγ
    end
    return qγ
end

"""
quarnetdistancematrix(
    net::HybridNetwork;
    showprogressbar=false,
    cost=:nanuqplus
)

Distance (or dissimilarity) matrix `d` between pairs of taxa, based on the
quarnets in `net`: the subnetworks induced by subsets of 4 taxa.
For taxa `x` and `y`, `d(x,y)` is the sum of `cost(x,y | q)` over all quarnets
`q` on 4-taxon subsets `{x,y,w,z}` that include both `x` and `y`.
The penalty `cost(x,y | q)` depends on the quartet trees *displayed* by quarnet `q`,
and possibly on their display probabilities (depending on the cost).
See [`quartetdisplayprobability`](@ref).

Output: n×n matrix where n is the number of taxa in the network, listed in the
same order as in `tiplabels(net)`.

Cost: scheme to penalize pairwise relationships within a quarnet.

- With `cost=:nanuqplus`, the "modified NANUQ" cost `ρ` is used (see
  [Allman et al. 2025](https://doi.org/10.1186/s13015-025-00274-w)):
  `ρc = 0.5` for **c**herries,
  `ρs = 1` for **s**plits (2 taxa **s**eparated by a tree-split),
  `ρa = 0.5` for **a**djacent pairs (in a quarnet that admits 1 circular order)
  `ρo = 1` for **o**pposite taxa (diagonal from each other in circular quarnets),
    also if the quarnet displays all 3 resolved quartets.
- With `cost=:nanuq` we use the original NANUQ costs: the same as above except
  that `ρc = 0`.
- With `cost=:mgamma`, the cost depends on the taxon pair relationship within
  the quarnet, and also on the inheritance probabilities:
  `cost(ab in q={a,b,c,d}) = 1-γ(ab|cd)`.
  This scheme is similar to NANUQ's, as it simplifies to `ρc = 0` for cherries,
  `ρs = 1` and `ρo = 1` for split and opposite taxa. But the cost of
  adjacent taxa depends on the "strength" of adjacency.
- To use custom costs, we can provide a named tuple with desired cost values:
  `cost = (cherry=0.5, split=2, adjacent=0.5, opposite=1)` 

**Note**: the (modified) NANUQ distance is defined as `2d(x,y) + 2n-4`,
where `n` is the number of taxa.
`d(x,y)` calculated here does **not** include the factor 2
nor the constant off-diagonal term `2n-4`. With these extra terms,
the NANUQ distance on a tree was proved to be additive on that tree
[(Rhodes 2019)](https://doi.org/10.1109/TCBB.2019.2917204).
It corresponds to edge lengths `w(e) = |X_1| |X_2| + |Y_1| |Y_2|`
for an internal edge `e` with quadripartition `X1,X2|Y_1,Y2`, and
`w(e) = |Y_1| |Y_2|` for an external edge `e` with tripartition `{x},Y_1,Y_2`
(see p.4 of [Rhodes (2019)](https://doi.org/10.1109/TCBB.2019.2917204)
for the general case with polytomies).

**Polytomies**: some quarnets with polytomies can display the star quartet.
The (modified) NANUQ distance definition in
[Allman et al. 2025](https://doi.org/10.1186/s13015-025-00274-w))
is for binary networks, which only display resolved quartets.
This implementation handles displayed star quartets `(abcd)` in the following way:
- For `cost=:mgamma`, we still have `cost(ab in q={a,b,c,d}) = 1-γ(ab|cd)`.
- When using `:nanuq`, `:nanuplus` or custom, the cost is a weighted average:
  `cost(ab in q={a,b,c,d}) = (1-γ(star))ρx + γ(star) * ρo`, where x=c,s,a, or o
  based on whether a and b are a cherry, star, adjacent, or opposite in the quarnet on {a,b,c,d},
  and γ(star) is the probability that the star quartet is displayed under the quarnet. 
  This is consistent with the metric for trees with polytomies in 
  [Rhodes (2019)](https://doi.org/10.1109/TCBB.2019.2917204).

```jldoctest nqd
julia> net = readnewick("(O:5.5,(((E:1.5)#H1:2.5::0.7,((#H1:0,D:1.5):1.5,((C:1,B:1):1)#H2:1::0.6):1.0):1.0,(#H2:0,A:2):3):0.5);");

julia> # using PhyloPlots; plot(net, showgamma=true);

julia> const PN = PhyloNetworks; # to write less later

julia> PN.quarnetdistancematrix(net; cost=:nanuqplus)
6×6 Matrix{Float64}:
 0.0  3.5  4.5  5.5  5.5  3.0
 3.5  0.0  3.0  5.5  5.5  4.5
 4.5  3.0  0.0  4.5  4.5  5.5
 5.5  5.5  4.5  0.0  3.0  4.5
 5.5  5.5  4.5  3.0  0.0  4.5
 3.0  4.5  5.5  4.5  4.5  0.0

julia> print(tiplabels(net))
["O", "E", "D", "C", "B", "A"]

julia> PN.quarnetdistancematrix(net; cost=:mgamma)
6×6 Matrix{Float64}:
 0.0   3.36  4.2   5.42  5.42  1.6
 3.36  0.0   1.68  5.4   5.4   4.16
 4.2   1.68  0.0   4.56  4.56  5.0
 5.42  5.4   4.56  0.0   0.0   4.62
 5.42  5.4   4.56  0.0   0.0   4.62
 1.6   4.16  5.0   4.62  4.62  0.0

julia> PN.quarnetdistancematrix(net;
        cost = (cherry=0.001, split=2, adjacent=0.5, opposite=1))
6×6 Matrix{Float64}:
 0.0    4.001  5.001  8.5    8.5    2.002
 4.001  0.0    2.002  8.5    8.5    5.001
 5.001  2.002  0.0    7.5    7.5    6.001
 8.5    8.5    7.5    0.0    0.006  7.5
 8.5    8.5    7.5    0.006  0.0    7.5
 2.002  5.001  6.001  7.5    7.5    0.0
```

The next example uses a tree. The distance we get is additive on that tree
(unrooted), so we can re-build the unrooted tree from distances,
using neighbor-joining for example. We get branch lengths that reflect the
quartet distance: **not** any distance from the original phylogeny (if any).

```jldoctest nqd
julia> caterpillar = readnewick("(((((a1,a2),a3),a4),a5),a6);");

julia> PN.quarnetdistancematrix(caterpillar; cost=:mgamma) # same as nanuq
6×6 Matrix{Float64}:
 0.0  0.0  3.0  5.0  6.0  6.0
 0.0  0.0  3.0  5.0  6.0  6.0
 3.0  3.0  0.0  4.0  5.0  5.0
 5.0  5.0  4.0  0.0  3.0  3.0
 6.0  6.0  5.0  3.0  0.0  0.0
 6.0  6.0  5.0  3.0  0.0  0.0

julia> d = PN.quarnetdistancematrix(caterpillar; cost=:nanuqplus)
6×6 Matrix{Float64}:
 0.0  3.0  4.5  5.5  6.0  6.0
 3.0  0.0  4.5  5.5  6.0  6.0
 4.5  4.5  0.0  5.0  5.5  5.5
 5.5  5.5  5.0  0.0  4.5  4.5
 6.0  6.0  5.5  4.5  0.0  3.0
 6.0  6.0  5.5  4.5  3.0  0.0

julia> tre = PN.nj!(copy(d), tiplabels(caterpillar)); writenewick(tre)
"((((a1:1.5,a2:1.5):1.0,a3:2.0):1.0,a4:2.0):1.0,a6:1.5,a5:1.5);"

julia> dtre = pairwisetaxondistancematrix(tre); dtre ≈ d
true
```
"""
function quarnetdistancematrix(
    net::HybridNetwork;
    showprogressbar=false,
    cost=:nanuqplus
)
    isa(cost, NamedTuple) || cost ∈ (:nanuq, :nanuqplus, :mgamma) ||
        error("invalid cost specification: $cost")
    addcost! =
        (cost == :nanuqplus ? addQDcost_nanuqplus! :
        (cost == :nanuq     ? addQDcost_nanuq! :
        (cost == :mgamma    ? addQDcost_gamma! :
        (d, qi, qγ) -> addQDcost_rho!(d, qi, qγ, cost[:cherry], cost[:split], cost[:adjacent], cost[:opposite])
        )))
    taxa = tiplabels(net)
    nn = length(taxa)
    nanuqd = zeros(Float64, nn,nn)
    quartet,t = quartetdisplayprobability(net; showprogressbar=showprogressbar)
    nn == length(t) || error("quartetdisplayprobability got a different number of taxa.")
    o = indexin(t, taxa) # taxa[o[i]] = t[i]:
    # place results for tip t[i] on nanuqd's row & column o[i]
    for q in quartet
        addcost!(nanuqd,  o[q.taxonnumber], q.data)
    end
    return nanuqd
end

function addQDcost_rho!(d, qi, qγ, chry, splt, adjt, opps)
    qdisp = map(x -> x ≉ 0.0, qγ)
    ndisp = sum(qdisp)
    # ndisp=1: 1 cherry, 2 splits. ndisp=2: 2 adjacent, 1 opposite
    # ndisp=0 or ndisp=3: all 3 considered opposite
    sumγ = sum(qγ) # P(resolved quartet)
    # common cost to all 6 pairs: star, or all 3 displayed
    c3 = (ndisp == 0 || ndisp == 3 ? opps : opps * (1-sumγ))
    (c1,c2) = (ndisp == 0 || ndisp == 3 ? (c3, c3) :
              (ndisp == 1 ? (sumγ*chry + c3, sumγ*splt + c3) :
                            (sumγ*adjt + c3, sumγ*opps + c3)))
    for jrow in 1:3 # resolution 1,i2|i3,i4
        c = (qdisp[jrow] ? c1 : c2)
        (i2,i3,i4) = (jrow == 1 ? (2,3,4) : (jrow == 2 ? (3,4,2) : (4,2,3)))
        d[qi[ 1],qi[i2]] += c; d[qi[i2],qi[ 1]] = d[qi[ 1],qi[i2]]
        d[qi[i3],qi[i4]] += c; d[qi[i4],qi[i3]] = d[qi[i3],qi[i4]]
    end
    return nothing
end
addQDcost_nanuq!(d,qi,qγ) = addQDcost_rho!(d,qi,qγ, 0. , 1., 0.5, 1.)
addQDcost_nanuqplus!(d,qi,qγ) = addQDcost_rho!(d,qi,qγ, 0.5, 1., 0.5, 1.)
function addQDcost_gamma!(d, qi, qγ)
    # like nanuq: chry=0, splt=opps=1. but adjt=1-γ instead of 1/2
    for jrow in 1:3 # resolution 1,i2|i3,i4
        γ = qγ[jrow]
        c = 1-γ
        (i2,i3,i4) = (jrow == 1 ? (2,3,4) : (jrow == 2 ? (3,4,2) : (4,2,3)))
        d[qi[ 1],qi[i2]] += c; d[qi[i2],qi[ 1]] = d[qi[ 1],qi[i2]]
        d[qi[i3],qi[i4]] += c; d[qi[i4],qi[i3]] = d[qi[i3],qi[i4]]
    end
    return nothing
end
