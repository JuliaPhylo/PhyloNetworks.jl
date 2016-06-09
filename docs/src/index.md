# `PhyloNetworks`

PhyloNetworks is a Julia package for the manipulation, visualization
and inference of phylogenetic networks.  SNaQ implements the
statistical inference method in
[Sol&iacute;s-Lemus and An&eacute; 2016](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896).
The procedure involves a numerical optimization of branch lengths and inheritance
probabilities and a heuristic search in the space of phylogenetic networks.

## Manual Outline

This manual is a quick tutorial on SNaQ and PhyloNetworks,
and [here](http://pages.stat.wisc.edu/~claudia/smallTutorial.pdf)
are slides with background on networks and explanations.
But check out the
[google group](https://groups.google.com/forum/#!forum/phylonetworks-users) for common questions.

```@contents
Pages = [
    "man/installation.md",
    "man/simpleJulia.md",
    "man/inputdata.md",
    "man/ticr_howtogetQuartetCFs.md",
    "man/snaq_plot.md",
    "man/fixednetworkoptim.md",
    "man/expectedCFs.md",
    "man/bootstrap.md",
]
Depth = 2
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
