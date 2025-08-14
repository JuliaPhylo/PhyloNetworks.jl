```@setup nj
using PhyloNetworks
using DataFrames, CSV
```
# Neighbor joining

A tree can be inferred from pairwise distances using the neighbor
joining algorithm ([Satou & Nei
1987](https://doi.org/10.1093/oxfordjournals.molbev.a040454))

## calculating distances

To read a DNA sequence alignment, see [`PhyloNetworks.readphylip`](@ref) and
[`PhyloNetworks.readfastatoarray`](@ref).
Then, to estimate distances from aligned sequences,
see [`PhyloNetworks.hammingdistancematrix`](@ref) and consider
[`PhyloNetworks.distancecorrection_JC!`](@ref), with an example given
in the [hamming distance](@ref) section.

Next, we use a pre-computed distance matrix.

## getting the joining algorithm tree

The [`nj`](@ref) function takes a data frame of pairwise distances as input
and constructs a tree using the neighbor joining method.  The column
names (headers) are used as taxon names.  Rows are assumed to
correspond to taxa in the same order as they do in columns.

```@repl nj
D = DataFrame(CSV.File(joinpath(dirname(pathof(PhyloNetworks)), "..","examples","caudata_dist.txt")); copycols=false);
tree = nj(D)
```

There is also a method [`PhyloNetworks.nj!`](@ref), which takes a distance
matrix and a vector of the names as argument.  This function, however,
would modify `D`.  One also has to make sure the vector of names match
the columns/rows of the distance matrix.
