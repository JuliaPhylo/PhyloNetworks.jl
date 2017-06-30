# PhyloNetworks.jl

PhyloNetworks is a [Julia](http://julialang.org) package for the
manipulation, visualization, inference of phylogenetic networks,
and their use for trait evolution.

## How to get help

- the package [wiki](https://github.com/crsl4/PhyloNetworks.jl/wiki) has a step-by-step
  tutorial, done for the Evolution 2016 conference, with background on networks and
  explanations.
- the [google group](https://groups.google.com/forum/#!forum/phylonetworks-users)
  has answers to common questions.
- the Manual below has a quick tutorial (navigation on the left).
- the [Index](@ref main-index) further below has the full list of documented functions.

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
    "man/trait_tree.md"
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
