```@setup nj
using PhyloNetworks
using CSV
```
# Neighbor joining

A tree can be inferred from pairwise distances using the neighbor
joining algorithm ([Satou & Nei
1987](https://doi.org/10.1093/oxfordjournals.molbev.a040454))

The [`nj`](@ref) function takes a data frame of pairwise distances as input
and constructs a tree using the neighbor joining method.  The column
names (headers) are used as taxon names.  Rows are assumed to
correspond to taxa in the same order as they do in columns.

```@repl nj
D = CSV.read(joinpath(dirname(pathof(PhyloNetworks)), "..","examples","caudata_dist.txt"));
tree = nj(D)
```

There is also a method [`PhyloNetworks.nj!`](@ref), which takes a distance
matrix and a vector of the names as argument.  This function, however,
would modify `D`.  One also has to make sure the vector of names match
the columns/rows of the distance matrix.
