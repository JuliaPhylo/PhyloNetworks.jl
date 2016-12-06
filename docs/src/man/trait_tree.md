# Trait Evolution

Once the network is inferred, we can use it for downstream analysis, including
for Phylogenetic Comparative Methods (PCM). In the following, we will show how
to use the functions related to trait evolution, limiting ourselves to the
special case where the studied network is in fact a phylogenetic tree. We will
assume that it is rooted, fixed, and with branch lengths.

In the examples below, we will use a dataset described in Mahler et al.
(2013)[^fn1], and available on Dryad [^fn2].  We can read it in from the file
"examples/lizard_tree.txt", that is in Newick format:
```julia
phy = readTopology(joinpath(Pkg.dir("PhyloNetworks"), "examples", "lizard_tree.txt"));
```
This example file can be accessed
[here](https://github.com/crsl4/PhyloNetworks/master/examples/lizard_tree.txt).

[^fn1]: Mahler DL, Ingram T, Revell LJ, Losos JB (2013). Exceptional
convergence on the macroevolutionary landscape in island lizard radiations.
Science 341(6143): 292-295. http://dx.doi.org/10.1126/science.1232392

[^fn2]: Mahler DL, Ingram T, Revell LJ, Losos JB (2013). Data from: Exceptional
convergence on the macroevolutionary landscape in island lizard radiations.
Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.9g182

## Simulation

First, we might want to simulate some trait data for a given set of species. For
now, only simulations according to a simple Brownian Motion (BM) are implemented.

We first need to create an object of class `paramsBM <: paramsProcess`:
```julia
params_simu = paramsBM(2, 1) # BM with mean 2 and variance 1
```
We can then simulate according to these parameters on the phylogeny, using
function `simulate`.
```julia
sim = simulate(phy, params_simu) # simulate a BM on phy
```
This creates an object of class `traitSimulation`, from which we can extract 
the data at the tips:
```julia
pred = sim[:Tips]
```
We can also extract the simulated values at the internal nodes of the
phylogeny:
```julia
predNodes = sim[:InternalNodes]
```
Both extractors create an `Array` with one columns, and as many lines as the
number of tips (respectively, internal nodes) there are in the phylogeny.  It
is sorted in the same order as the tips (respectively, internal nodes) of the
phylogeny used to simulate it.

## Phylogenetic Regression

The main function that can be used to do phylogenetic regression is the
function `phyloNetworklm`. It is based on function `lm` from package
[GLM](https://github.com/JuliaStats/GLM.jl), and inherit from a lot of its
features.

We first need to get some data to analyze. Here, we can use the predictor we
just simulated to create a new trait that depends linearly on the predictor,
with a noise that has a phylogenetic structure:
```julia
noise = simulate(phy, paramsBM(0, 0.1)) # Phylogenetic residuals
trait = 10 + 2 * pred + noise[:Tips] # Trait to study
```
We then need to create a data frame that contains the trait, the predictor, and
the tips names.  In order to avoid confusion, it is important that the data
frame contains an extra column specifying the tips names, and labeled
`tipsNames`.
```julia
using DataFrames
dat = DataFrame(trait = trait, pred = pred, tipsNames = tipLabels(sim))
```

The fit can then be done in a fashion that is very similar to the
[GLM](https://github.com/JuliaStats/GLM.jl) package, using formulas. We just
need to specify the phylogeny we are using:
```julia
fitBM = phyloNetworklm(trait ~ pred, dat, phy)
``` 
This returns an object of type `phyloNetworkLinearModel<:LinPredModel`. It is
dominated by the GLM type `LinPredModel`, which means that all base functions
from the [StatsBase]https://github.com/JuliaStats/StatsBase.jl) package can be
applied to it. See the documentation for this type for a list of all functions
that can be used.  Some functions allow the user to retrieve directly the
estimated parameters of the BM, and are specific to this object.
```julia
@doc phyloNetworkLinearModel # List all the base functions
sigma2_estim(fitBM) # The estimated variance of the BM
mu_estim(fitBM) # The estimated root value of the BM
```

## Ancestral State Reconstruction

## Phylogenetic Anova

## Pagel's Lambda

