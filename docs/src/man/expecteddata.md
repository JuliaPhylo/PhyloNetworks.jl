```@setup edata
using PhyloNetworks
using DataFrames, CSV
using RCall, PhyloPlots
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
```

# Pairwise and quartet data

We show here some functionalities to calculate data expected
from a given network, or observed in data.
To calculate expectations under a given network, this network
needs to have branch lengths and γ inheritances at hybrids.

We use 2 example networks in this section:
`net0` without reticulations (a tree) and
`net2` with 2 reticulations.
`net0` is in fact `net2`'s major tree: obtained by
deleting every minor hybrid edge.

```@example edata
net0 = readnewick("(O:5.5,((E:4.0,(D:3.0,(C:1.0,B:1.0):2.0):1.0):1.0,A:5.0):0.5);");
net2 = readnewick("(O:5.5,(((E:1.5)#H1:2.5::0.7,((#H1:0,D:1.5):1.5,((C:1,B:1):1)#H2:1::0.6):1.0):1.0,(#H2:0,A:2):3):0.5);")
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
![net0 and net2](../assets/figures/expectedata_fig_net02.svg)

## average pairwise distances

One distance between pairs of taxa, say between t1 and t2, is the
average length of all "up-down" paths in the network to go from t1 to t2,
see [Xu & Ané 2023](https://doi.org/10.1007/s00285-022-01847-8) for example.
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

```@example edata
aveD_net2 = pairwisetaxondistancematrix(net2);
taxonlist2 = tiplabels(net2)
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

The f2-statistic gives another measure of dissimilarity between pairs of taxa
(see [Lipson 2020](https://doi.org/10.1111/1755-0998.13230) for example).
The expected value of f2 between taxa t₁ and t₂ is
```math
f_2(t_1, t_2) = E(X(t_1) - X(t_2))^2
```
under a Brownian motion model for trait X evolving along the network,
where X(t₁) and X(t₂) are the values of X for taxa t₁ and t₂.

If the network is a tree, then this is exactly the average distance
(or simply, the length of the unique path) between t₁ and t₂.

It can be calculated with [`PhyloNetworks.expectedf2matrix`](@ref).
On our tree `net0`, we an f2 matrix equal to the average distance matrix:

```@example edata
const PN = PhyloNetworks; # for lazy typing below!
f2D_net0 = PN.expectedf2matrix(net0)
f2D_net0 == aveD_net0
```

Again, taxa are listed along rows and along columns in the same
order as listed by `tiplabels()`.  
On our network `net2`, the f2 and average distances differ:

```@example edata
f2D_net2 = PN.expectedf2matrix(net2)
DataFrame(f2D_net2, taxonlist2)
```

## expected f4-statistics

coming next: example to use
- a new function to calculate f3, and
- a new function to calculate f4 expected from a network

## quartet concordance factors

Tools to calculate quartet concordance factors expected from
a network are provided in package
[QGoF](https://github.com/JuliaPhylo/QuartetNetworkGoodnessFit.jl):
see its documentation about [expected concordance factors](@extref QGoF).

coming next: example to use
[`countquartetsintrees`](@ref) and [`tablequartetCF`](@ref)
Refer to QGoF for expected qCFs.
