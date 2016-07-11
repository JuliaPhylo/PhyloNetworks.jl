# PhyloNetworks.jl

PhyloNetworks is a Julia package for the manipulation, visualization
and inference of phylogenetic networks.  SNaQ implements the
statistical inference method in
[Sol&iacute;s-Lemus and An&eacute; 2016](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896).
The procedure involves a numerical optimization of branch lengths and inheritance
probabilities and a heuristic search in the space of phylogenetic networks.

## How to get help

- the package [wiki](https://github.com/crsl4/PhyloNetworks.jl/wiki) has a step-by-step
  tutorial, done for the Evolution 2016 conference, with background on networks and
  explanations.
- the [google group](https://groups.google.com/forum/#!forum/phylonetworks-users)
  has answers to common questions.
- the Manual below has a quick tutorial on SNaQ and PhyloNetworks.
- the [Index](@ref main-index) further below has the full list of documented functions.

## Manual Outline

```@contents
Pages = [
    "man/installation.md",
    "man/inputdata.md",
    "man/ticr_howtogetQuartetCFs.md",
    "man/snaq_plot.md",
    "man/fixednetworkoptim.md",
    "man/expectedCFs.md",
    "man/bootstrap.md",
    "man/multiplealleles.md"
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
