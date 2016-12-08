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
params_simu = paramsBM(2, 0.5) # BM with mean 2 and variance 0.5
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
Both extractors create an `Array` with one column, and as many lines as the
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
fitTrait = phyloNetworklm(trait ~ pred, dat, phy)
``` 
This returns an object of type `phyloNetworkLinearModel<:LinPredModel`. It is
dominated by the GLM type `LinPredModel`, which means that all base functions
from the [StatsBase](https://github.com/JuliaStats/StatsBase.jl) package can be
applied to it. See the documentation for this type for a list of all functions
that can be used. Some functions allow the user to retrieve directly the
estimated parameters of the BM, and are specific to this object.
```julia
@doc phyloNetworkLinearModel # List all the base functions
sigma2_estim(fitTrait) # The estimated variance of the BM
mu_estim(fitTrait) # The estimated root value of the BM
```

## Ancestral State Reconstruction

When we know the model of evolution that the traits followed, we can do some
ancestral state reconstruction, finding the Best Linear Unbiased Predictor
(BLUP) and some prediction intervals for the values of the traits at the
internal nodes of the phylogeny.

### From known parameters

Here, we simulated the predictor trait ourselves, so we know the exact model of
evolution that generated the data. In this favorable case, an ancestral state
reconstruction can be obtained using function `ancestralStateReconstruction`:
```julia
ancPred = ancestralStateReconstruction(phy, pred, params_simu)
```
The object created has type `reconstructedStates`. Several extractors, as well
as a plot function, can be applied to it:
```julia
plot(phy, ancPred)
expectations(ancPred) # The predictors
stderr(ancPred) # The standard errors associated
predint(ancPred) # The prediction interval (default to 95%)
```
As we know the true ancestral states here, we can compare them to our
estimation. In particular, we can count the number of times the true value lies
in the prediction interval. This should be true about 95% of the times.
```julia
sum((predNodes .< predint(ancPred)[1:99, 2]) & (predNodes .> predint(ancPred)[1:99, 1]))
```

### From estimated parameters

In real applications though, we do not have access to the true parameters of
the process that generated the data. If we make an assumption on the type of
process (assuming here that the trait evolves as a BM), we can however compute
an approximation of it, using the previous function. To fit a regular BM, we
just need to do a regression of the trait against a simple intercept:
```julia
fitPred = phyloNetworklm(pred ~ 1, dat, phy)
```
and then, we can apply the `ancestralStateReconstruction` function directly to
the fitted object:
```julia
ancPredApprox = ancestralStateReconstruction(fitPred)
```
As we estimated the parameters, the prediction interval we get are less
accurate, and the function reminds us that by throwing a warning. Apart from
that, it creates an object of the same `reconstructedStates` type as earlier,
and the same extractors can be applied to it:
```julia
plot(phy, ancPredApprox)
```
For convenience reasons, the two steps described above (fitting against the
intercept, and then do ancestral state reconstruction) can be done all at once
with a single call of the function `ancestralStateReconstruction` on a
DataFrame with the trait to reconstruct, and the tips labels:
```julia
datPred = DataFrame(pred = pred, tipsNames = tipLabels(sim))
ancPredApprox = ancestralStateReconstruction(datPred, phy)
plot(phy, ancPredApprox)
``` 
This produces the exact same results.

### With known predictors

At this point, it might be tempting to apply this function to the trait we
simulated earlier as a linear combination of the predictor and a phylogenetic
noise.  However, this can be done directly:
```julia
ancTrait = ancestralStateReconstruction(fitTrait) # Throws an error !
```
This is because the model we used to fit the trait (a regression with one
predictor and an intercept) is not compatible with the simple model of Brownian
evolution that we assumed for the ancestral state reconstruction. As the
predictor used is not known for ancestral states, it is not possible to
reconstruct the trait for this particular model.

The only option we have is to provide the function with some ancestral states,
if they are known. Thankfully, in this toy example, as we generated the
predictors ourselves, it is indeed the case, and we can reconstruct our trait
doing the following:
```julia
ancTrait = ancestralStateReconstruction(fitTrait, [ones(99, 1) predNodes])
plot(phy, ancTrait)
```
Where we provided the ancestral predictors as a matrix, containing the
intercept, and the known predictor at the nodes. The user must be very careful
with this function, as no check is done for the order of the predictors, that
must be in the same order as the internal nodes of the phylogeny. As ancestral
predictors are often unknown, the use of this functionality is discouraged.

## Data Imputation

There is no theoretical difference between an internal node, for which we could
not measure the value of the trait, and a missing value at a tip of the
phylogeny. Consequently, the previous `ancestralStateReconstruction` function
can be used to do data imputation. To see this, let's add some missing values
in the `pred` array we used earlier:
```julia
datPred[[2, 6, 15, 38, 54, 89], :pred] = NA;
ancPredApprox = ancestralStateReconstruction(datPred, phy)
plot(phy, ancPredApprox)
```
In the plotting function, a prediction interval is shown for the missing
values.

## Phylogenetic Anova

As mentioned above, the `phyloNetworklm`function is based on the `lm` function
from package [GLM](https://github.com/JuliaStats/GLM.jl). This means that it
inherits from most of its features, and in particular, it can handle formulas
with factors or interactions. To see this, let's load a subset of the lizard
dataset [^fn1] we already used:
```julia
dat = readtable(joinpath(Pkg.dir("PhyloNetworks"), "examples", "lizard_trait.txt"));
```
This example file can be accessed
[here](https://github.com/crsl4/PhyloNetworks/master/examples/lizard_trait.txt).
This DataFrame has four columns: `tipsNames` contains the tips names of the
species, matching the one on the phylogeny we already loaded earlier; `AVG_SVL`
is the average body size ("Snout to Vent Length"); `AVG_ltoe_IV` and
`AVG_lfing_IV` are the average length of the fourth toe and finger; and
`region` is a factor with 4 levels coding for the region where each species is
found. See the description of the dataset on Dryad[^fn2] for further
informations, and other measured traits.

For the sake of the example, let's assume that we want to do a regression of
the SVL length against the toe and finger lengths, taking the region into
account, and with an interaction between the finger length and the region.
(This model has no biological basis, and is just used to show the possibilities
of the function). First, we need to make sure that the region information is
indeed considered as a factor. This is done by transforming the region column
into a `PooledDataArray`:
```julia
dat[:region] = PooledDataArray(dat[:region]); # Pool by region
```
Then, we just need to apply `phyloNerworklm` with the right formula:
```julia
fitAnova = phyloNetworklm(AVG_SVL ~ AVG_ltoe_IV + AVG_lfing_IV * region, dat, phy)
```
As before, a summary is shown by default, and all the functions that can be applied
to a `lm` result still apply here.

## Pagel's Lambda

One classical question in Phylogenetic Comparative Method is to assess the
"Phylogenetic Signal" of a dataset, that is, the importance of the tree
structure to explain the observed traits. One way of doing that is to use
Pagel's lambda[^fn3] transformation of the branch lengths. This model assumes a
BM on a tree where the internal branches are multiplied by a factor $\lambda$,
while the external branches are modified so that the total hight of the tree is
constant. Hence, $\lambda$ varies between $0$ (the tree has no influence on
the data) and $1$ (the tree is unchanged).

We can illustrate this with the predictor trait we used earlier. We use the
same function as before, only indicating the model we want to use:
```julia
fitPagel = phyloNetworklm(pred ~ 1, datPred, phy, model = "lambda")
```
As it is indeed generated according to a plain BM on the phylogeny, the
estimated $\lambda$ should be close to $1$.

Note that we took a version of the dataset with missing values, that can
be readily handed by function `phyloNetworklm`.


[^fn3]: Pagel M (1999). Inferring the historical patterns of biological
        evolution. Nature. 401: 877–884. doi:10.1038/44766

