# Neighbor joining

Neighbor joining is also implemented in the package, and we can using
pairwise distances to construct a tree.

The `nj` function takes a data frame of pairwise distances as input,
and construct a tree using the neighbor joining method.  The column
names (headers) are used as taxon names.

```@example
using PhyloNetworks
D = CSV.read(joinpath(dirname(pathof(PhyloNetworks)), "..","examples","caudata_dist.txt"); types=[Float64 for i in 1:197]);
```

```@repl
nj(D)
```

There is also a method `PhyloNetworks.nj!`, which takes a distance
matrix and a vector of the names as argument.  This function, however,
would modify `D`.  One also has to make sure the vector of names match
the columns/rows of the distance matrix.
