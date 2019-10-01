# PhyloNetworks.jl

[PhyloNetworks](https://github.com/crsl4/PhyloNetworks.jl)
is a [Julia](http://julialang.org) package for the
manipulation, visualization, inference of phylogenetic networks,
and their use for trait evolution.

---

**How to get help**

- the package [wiki](https://github.com/crsl4/PhyloNetworks.jl/wiki) has a step-by-step
  tutorial, done for the 2019 MBL workshop, with background on networks and
  explanations.
- the [google group](https://groups.google.com/forum/#!forum/phylonetworks-users)
  has answers to common questions.
- the Manual below has a quick tutorial (navigation on the left).
- the [Index](@ref main-index) further below has the full list of documented functions.

## References

- Claudia Solís-Lemus, Paul Bastide and Cécile Ané (2017).
  PhyloNetworks: a package for phylogenetic networks.
  [Molecular Biology and Evolution](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msx235/4103410/PhyloNetworks-a-package-for-phylogenetic-networks?guestAccessKey=230afceb-df28-4160-832d-aa7c73f86369)
  34(12):3292–3298.
  [doi:10.1093/molbev/msx235](https://doi.org/10.1093/molbev/msx235)
- Bastide, Solís-Lemus, Kriebel, Sparks, Ané (2018).
  Phylogenetic Comparative Methods for Phylogenetic Networks with Reticulations.
  Systematic Biology, 67(5):800–820.
  [doi:10.1093/sysbio/syy033](https://doi.org/10.1093/sysbio/syy033).
- Claudia Solís-Lemus and Cécile Ané (2016).
  Inferring Phylogenetic Networks with Maximum Pseudolikelihood under Incomplete Lineage Sorting.
  [PLoS Genet](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896)
  12(3):e1005896. [doi:10.1371/journal.pgen.1005896](https://doi.org/10.1371/journal.pgen.1005896)

## Manual Outline

```@contents
Pages = [
    "man/installation.md",
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

## Library Outline

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
