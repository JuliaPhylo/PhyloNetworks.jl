# PhyloNetworks.jl

[PhyloNetworks](https://github.com/juliaphylo/PhyloNetworks.jl)
is a [Julia](http://julialang.org) package for the
manipulation, visualization, inference of phylogenetic networks,
and their use for trait evolution.

---

**How to get help**

- the package [wiki](https://github.com/juliaphylo/PhyloNetworks.jl/wiki) has a step-by-step
  tutorial, done for the MBL workshop (last revised 2022), with background on networks and
  explanations.
- [tutorial](https://cecileane.github.io/networkPCM-workshop/) for
  comparative methods, including network calibration (2023 workshop)
- the [google group](https://groups.google.com/forum/#!forum/phylonetworks-users)
  has answers to common questions.
- the [Manual](@ref) below has a quick tutorial (navigation on the left).
- the [Library](@ref) and
  and [Index](@ref main-index) further below has the full list of documented functions.

## References

See their [bibtex format](https://github.com/juliaphylo/PhyloNetworks.jl/blob/master/CITATION.bib).

for the package:
- Claudia Solís-Lemus, Paul Bastide and Cécile Ané (2017).
  PhyloNetworks: a package for phylogenetic networks.
  [Molecular Biology and Evolution](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msx235/4103410/PhyloNetworks-a-package-for-phylogenetic-networks?guestAccessKey=230afceb-df28-4160-832d-aa7c73f86369)
  34(12):3292–3298.
  [doi:10.1093/molbev/msx235](https://doi.org/10.1093/molbev/msx235)

for trait evolution:
- Teo, Rose, Bastide & Ané (2022).
  Accounting for intraspecific variation in continuous trait evolution
  on a reticulate phylogeny.
  [bioRxiv](https://doi.org/10.1101/2022.05.12.490814)
- Karimi, Grover, Gallagher, Wendel, Ané & Baum (2020). Reticulate evolution
  helps explain apparent homoplasy in floral biology and pollination in baobabs
  (*Adansonia*; Bombacoideae; Malvaceae).
  Systematic Biology, 69(3):462-478.
  [doi:10.1093/sysbio/syz073](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syz073/5613901?guestAccessKey=a32e7dd3-27fd-4a13-b171-7ff5d6da0e01).
- Bastide, Solís-Lemus, Kriebel, Sparks, Ané (2018).
  Phylogenetic Comparative Methods for Phylogenetic Networks with Reticulations.
  Systematic Biology, 67(5):800–820.
  [doi:10.1093/sysbio/syy033](https://doi.org/10.1093/sysbio/syy033).

for network inference:
- Claudia Solís-Lemus and Cécile Ané (2016).
  Inferring Phylogenetic Networks with Maximum Pseudolikelihood under Incomplete Lineage Sorting.
  [PLoS Genet](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896)
  12(3):e1005896. [doi:10.1371/journal.pgen.1005896](https://doi.org/10.1371/journal.pgen.1005896)

## Manual

```@contents
Pages = [
    "man/installation.md",
    "man/netmanipulation.md",
    "man/inputdata.md",
    "man/ticr_howtogetQuartetCFs.md",
    "man/snaq_plot.md",
    "man/dist_reroot.md",
    "man/fixednetworkoptim.md",
    "man/expectedCFs.md",
    "man/bootstrap.md",
    "man/multiplealleles.md",
    "man/trait_tree.md",
    "man/parsimony.md",
    "man/fitdiscrete.md",
    "man/nj.md"
]
Depth = 3
```

## Library

```@contents
Pages = ["lib/public.md", "lib/internals.md"]
Depth = 2
```

## [Index](@id main-index)

### Functions

```@index
Pages = ["lib/public.md", "lib/internals.md"]
Order = [:function]
```

### Types

```@index
Pages = ["lib/public.md", "lib/internals.md"]
Order = [:type]
```
