```@setup edata
using PhyloNetworks
using DataFrames, CSV
using RCall, PhyloPlots
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
```

```@eval edata
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

coming next: example to use
[`pairwisetaxondistancematrix`](@ref)

perhaps: extract from [`startingBL!`](@ref) the code to calculate
pairwise distances from data, either Hamming or JC-corrected,
and show how to use it here.

## expected f2-statistics

coming next: example to use
[`expectedf2matrix`](@ref)

## expected f4-statistics

coming next: example to use
- a new function to calculate f3, and
- a new function to calculate f4 expected from a network

## quartet concordance factors

coming next: example to use
[`countquartetsintrees`](@ref) and [`tablequartetCF`](@ref)
Refer to QGoF for expected qCFs.
