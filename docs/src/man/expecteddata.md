```@setup edata
using PhyloNetworks
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
```

# Pairwise and quartet data

We show here some functionalities to calculate data expected
from a given network, or observed in data.
To calculate expectations under a network, this network's edges
need to have branch lengths and γ inheritance values.

We use 2 example networks in this section:
`net0` without reticulations (a tree) and
`net2` with 2 reticulations.
`net0` is in fact `net2`'s major tree: obtained by
deleting every minor hybrid edge.

```@example edata
net0 = readnewick("(O:5.5,((E:4.0,(D:3.0,(C:1.0,B:1.0):2.0):1.0):1.0,A:5.0):0.5);");
net2 = readnewick("(O:5.5,(((E:1.5)#H1:2.5::0.7,((#H1:0,D:1.5):1.5,
    ((C:1,B:1):1)#H2:1::0.6):1.0):1.0,(#H2:0,A:2):3):0.5);")
nothing # hide
```

```@setup edata
using RCall, PhyloPlots
R"svg"(figname("expectedata_fig_net02.svg"), width=7, height=3);
R"layout"([1 2]);
R"par"(mar=[0,0,0.5,0]);
plot(net0, useedgelength=true, showedgelength=true, tipoffset=0.1);
R"mtext"("net0, showing edge lengths", line=-1);
plot(net2, showgamma=true, useedgelength=true, style=:majortree, arrowlen=0.1, tipoffset=0.1);
R"mtext"("net2, showing γ's", line=-1);
R"dev.off"();
```
![net02-distquartet](../assets/figures/expectedata_fig_net02.svg)

## average pairwise distances

One distance between pairs of taxa, say between t1 and t2, is the
average length of all "up-down" paths in the network to go from t1 to t2
(see [Xu & Ané 2023](https://doi.org/10.1007/s00285-022-01847-8) for example).
In a tree, there is a single such path: going up from t1 to their
most recent common ancestor, then down to t2.
In a general network, there can be multiple paths.
Each path has an inheritance weight: the product of the inheritance γ's of
all edges in the path.
Under a simple model, this is the proportion of genetic material that took
this path.
The average distance betwen t1 and t2 is the weighted average of these
paths' lengths, weighted by the paths' inheritance.

It can be calculated with [`pairwisetaxondistancematrix`](@ref).
On our tree `net0`, we get this:

```@repl edata
aveD_net0 = pairwisetaxondistancematrix(net0)
```
The order of rows and columns in this matrix is the same as
given by `tiplabels`. For example, we can convert the distance matrix
`aveD_net0` to a data frame with a column named for each taxon like this.

```@repl edata
taxonlist0 = tiplabels(net0)
using DataFrames
DataFrame(aveD_net0, taxonlist0)
```

Using the network `net2`, we see that its reticulations bring
E closer to D; and bring B & C closer to A (and away from D & E
since A is distant from them):

```@repl edata
aveD_net2 = pairwisetaxondistancematrix(net2);
taxonlist2 = tiplabels(net2);
DataFrame(aveD_net2, taxonlist2)
```

To calculate pairwise distances observed in data, such as
along a DNA alignment across multiple sites, we can use
[`PhyloNetworks.hammingdistancematrix`](@ref).

```@repl edata
dna_data = [
    ['T','G','T','A','G'], # taxon 1
    ['T','G','A','A','G'],
    ['T','G','A','A','C'],
    ['T','G','A',missing,missing],
    ['T','G','A','T','C'], # taxon 5
];
d = PhyloNetworks.hammingdistancematrix(dna_data)
```

Perhaps we want to give weights to each trait, such as if an alignment is
summarized by keeping each site pattern once, weighted by the number of sites
having this pattern:

```@repl edata
site_count = [3,2,1,1,1]; # invariable site patterns have higher counts
d = PhyloNetworks.hammingdistancematrix(dna_data, site_count)
```

A Jukes-Cantor correction can be applied with
[`PhyloNetworks.distancecorrection_JC!`](@ref):

```@repl edata
PhyloNetworks.distancecorrection_JC!(d, 4) # 4 states
```

## expected f2-statistics

The f2-statistic gives another measure of dissimilarity between pairs of taxa.
The expected value of f2 between taxa t₁ and t₂ is
```math
f_2(t_1, t_2) = E(X(t_1) - X(t_2))^2
```
under a Brownian motion model for trait X evolving along the network,
where X(t₁) and X(t₂) are the values of X for taxa t₁ and t₂.

In the network, branch lengths measures units of drifts when the data X
are allele frequencies, and with a variance factor dependent on the allele
frequency at the root.
For background, see for example
[Patterson et al. 2012](https://doi.org/10.1534/genetics.112.145037) or
[Lipson 2020](https://doi.org/10.1111/1755-0998.13230).

If the network is a tree, then this is exactly the average distance
(or simply, the length of the unique path) between t₁ and t₂.

It can be calculated with [`expectedf2matrix`](@ref).
On our tree `net0`, we get an f2 matrix equal to the average distance matrix:

```@repl edata
f2D_net0 = expectedf2matrix(net0)
f2D_net0 == aveD_net0
```

Again, taxa are listed along rows and along columns in the same
order as listed by `tiplabels()`.  
On our network `net2`, the f2 and average distances differ:

```@repl edata
f2D_net2 = expectedf2matrix(net2);
DataFrame(f2D_net2, taxonlist2)
```

## expected f4-statistics

```@repl edata
f4,t = expectedf4table(net0);
t # taxa, but ordered alphabetically: not as in tiplabels(net0)
first(f4,2) # first 2 4-taxon sets, each with 3 quartet f4s
```

The taxon numbers above are indices in the taxon list `t`.
Here is a way to print all f4 statistics expected from our tree `net0`:
```@repl edata
for q in f4
    println(join(t[q.taxonnumber],",") * ": " * string(round.(q.data, sigdigits=2)))
end
```

For each set of 4 taxa, the three f4s sum up to 0: as it should be.
On a tree with a split `t1,t2|t4,t5`, the corresponding f4 value should
be 0, and the other 2 should give ± the length of the internal path separating
the 2 groups of 2 taxa.
(Again, the branch lengths unit depends on the data being considered.)

Here for example, the first 4-taxon set is `A,B,C,D`. In our tree,
`BC` is a clade, separated from `AD` by a branch of length 2.
Accordingly, the third f4 value, corresponding to `AD|BC`, is 0;
and the other two f4s are 2 or -2.

![net02-distquartet-again](../assets/figures/expectedata_fig_net02.svg)

We can see how adding reticulations to our tree affects expected f4s.
Note the colums names below: they correspond to the order of taxa
for each f4.

```@repl edata
f4,t = expectedf4table(net2, showprogressbar=false); # same t: alphabetically
nt = tablequartetdata(f4, t; colnames="f4_" .* ["12_34", "13_42", "14_23"]);
df = DataFrame(nt) # convert table to data frame
```

We still have f4=0 on the 3rd column for first taxon set `A,B,C,D`,
because `BC` are still sister in the network.
But the last taxon set for example, `C,D,E,O`, has no 0 values of f4
due to the reticulation between ancestors of `D` and `E`:
in the subnetwork for `C,D,E,O`, there is a cycle of 4 edges.

We may also calculate f3 statistics using [`PhyloNetworks.expectedf3matrix`](@ref).
For this, we need a reference taxon. Here we use the outgroup O:

```@repl edata
f3_net2 = PhyloNetworks.expectedf3matrix(net2, "O");
DataFrame(f3_net2, taxonlist2)
```

## quartet concordance factors

The concordance factor of a quartet `ab|cd` is the proportion of the
genome whose genealogy has this unrooted topology.

Tools to calculate quartet concordance factors expected from
a network under the coalescent model are provided in package
[QGoF](https://github.com/JuliaPhylo/QuartetNetworkGoodnessFit.jl):
see its documentation about [expected concordance factors](https://juliaphylo.github.io/QuartetNetworkGoodnessFit.jl/stable/man/expected_qCFs/#expected-concordance-factors).

To calculate quartet concordance factors observed in data, one option is
to count the number of gene trees that display each quartet, using
[`countquartetsintrees`](@ref) and [`tablequartetCF`](@ref).

In the example below, the 4th gene tree is missing taxon A, and
the 6th gene tree has a polytomies (unresolved ABE clade), such as if
a branch of low support was collapsed.
The number of genes underlying each quartet is captured in the table below.

```@repl edata
sixgenetrees_nwk = [
  "(E,((A,B),(C,D)),O);","(((A,B),(C,D)),(E,O));","(A,B,((C,D),(E,O)));",
  "(B,((C,D),(E,O)));","((C,D),(A,(B,E)),O);","((C,D),(A,B,E),O);"];
genetrees = readnewick.(sixgenetrees_nwk);
q,t = countquartetsintrees(genetrees, showprogressbar=false);
df = tablequartetCF(q,t) |> DataFrame
```

If low-support branches are not collapsed, this counting method does not
account for gene tree estimation error.
(It biases concordance factors towards 1/3 for each resolution of a
4-taxon tree: lack of knowledge is mistaken as lack of concordance).
Estimation error can be accounted for in a Bayesian framework:
see [PhyloUtilities](https://juliaphylo.github.io/PhyloUtilities/)
for a pipeline.
