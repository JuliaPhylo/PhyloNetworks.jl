# PhyloNetworks.jl

[PhyloNetworks](https://github.com/juliaphylo/PhyloNetworks.jl)
is a [Julia](http://julialang.org) package with core utilities for
phylogenetic networks.
Phylogenetic networks represent the evolutionary relationships between a set
of organisms, e.g. populations, species, languages, cultures etc.
They provide an explicit representation of splitting events (when populations
diverge from one another) and of merging events (when populations mix together
due to migration of individuals, hybridization, polyploidization, recombination,
etc.).
Phylogenetic networks are called *admixture graphs* in population genetics,
when merging events are interpreted as admixture between populations.
They can also summarize ancestral recombination graphs.

`PhyloNetworks` is a core package that supports reading, writing, manipulating
phylogenetic networks, and other standard tools.
It is used by other packages for more specialized tasks, such as
- [`PhyloPlots`](https://github.com/juliaphylo/PhyloPlots.jl)
  to visualize phylogenetic networks,
- [`SNaQ`](https://github.com/juliaphylo/SNaQ.jl)
  for the inference of phylogenetic networks from quartet concordance factors (qCF),
- [`QuartetNetworkGoodnessFit`](https://github.com/JuliaPhylo/QuartetNetworkGoodnessFit.jl)
  to calculate quartet concordance factors expected from a general network
  and test the adequacy of a network to qCF data,
- [`PhyloTraits`](https://github.com/juliaphylo/PhyloTraits.jl)
  for the inference of trait evolution along a phylogenetic network,
- [`PhyloCoalSimulations`](https://github.com/juliaphylo/PhyloCoalSimulations.jl)
  to simulate gene trees under a coalescent process along a phylogenetic networks.

---

**How to get help**

- the package [wiki](https://github.com/juliaphylo/PhyloNetworks.jl/wiki) has a
  step-by-step tutorial, done for the MBL workshop (last revised 2022),
  with background on networks and explanations.
- [tutorial](https://cecileane.github.io/networkPCM-workshop/) for
  comparative methods, including network calibration (2023 workshop)
- the [google group](https://groups.google.com/forum/#!forum/phylonetworks-users)
  has answers to common questions.
- the [Manual](@ref) and below has a quick tutorial (navigation on the left).
- the [Library](@ref) further below has the full list of documented functions.

## References

See their [bibtex format](https://github.com/juliaphylo/PhyloNetworks.jl/blob/master/CITATION.bib).

for the package in particular, please cite:
- Claudia Solís-Lemus, Paul Bastide and Cécile Ané (2017).
  PhyloNetworks: a package for phylogenetic networks.
  [Molecular Biology and Evolution](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msx235/4103410/PhyloNetworks-a-package-for-phylogenetic-networks?guestAccessKey=230afceb-df28-4160-832d-aa7c73f86369)
  34(12):3292–3298.
  [doi:10.1093/molbev/msx235](https://doi.org/10.1093/molbev/msx235)

## Manual

```@contents
Pages = [
    "man/installation.md",
    "man/netmanipulation.md",
    "man/net_plot.md",
    "man/dist_reroot.md",
    "man/network_support.md",
    "man/parsimony.md",
    "man/nj.md"
]
Depth = 3
```

## Library

For help on individual functions, see the library:

```@contents
Pages = [
    "lib/public.md",
    "lib/internals.md",
]
Depth = 3
```
