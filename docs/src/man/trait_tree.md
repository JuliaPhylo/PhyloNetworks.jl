```@setup tree_trait
using PhyloNetworks
mkpath("../assets/figures")
```
# Continuous Trait Evolution

Once the network is inferred, we can take
these species relationships into account when studying the distribution of quantitative
traits measured for extant species.
This is the goal of phylogenetic comparative methods (PCM).
More details can be found on the developments below in Bastide et al. 2018 [^B18]

We assume a fixed network, correctly rooted, with branch lengths
proportional to calendar time. Here, we consider the true network that was
used in the previous sections, and which is ultrametric (all the tips are contemporary).
```@example tree_trait
truenet = readTopology("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
```
As previously, we can plot the network thanks to the `RCall` package.
The `name` function is only instrumental here, to ensure that the figure is
saved in the correct directory when the documentation is built.
We only show the commands to actually save the plot in this first example for
the interested reader, but we will hide those in the rest of the chapter, for
the sake of clarity.
```@example tree_trait
using PhyloPlots, RCall
R"name <- function(x) file.path('..', 'assets', 'figures', x)"
R"svg(name('truenet.svg'), width=8, height=4)"
R"par"(mar=[0,0,0,0])
plot(truenet, :R, useEdgeLength=true, showGamma=true);
R"dev.off()"
nothing # hide
```
![truenet](../assets/figures/truenet.svg)

## Model and Variance Matrix

Assuming that the network is known and that the continuous traits evolve like a
Brownian Motion (BM) in time, it is possible to compute the expected variance
covariance matrix between tip measurements. This can be done using function
[`vcv`](@ref), whose syntax is inspired from the well known corresponding
[`ape`](https://CRAN.R-project.org/package=ape) function.
```@repl tree_trait
C = vcv(truenet)
```
The matrix is returned as a `DataFrame`, with columns named by the
tips of the network to allow for easy identification.
Each row also corresponds to a tip in the network, and rows are
ordered in the same way as columns.

The computation of this matrix is based on the more general function
[`sharedPathMatrix`](@ref). It is at the core of all the Phylogenetic
Comparative Methods described below.


## Trait simulation

We start by generating continuous traits to study. We simulate three
traits on the network (two independent, one dependent),
using a Brownian Motion (BM) model of trait evolution on the network. We start
by choosing the parameters of the BM (ancestral mean and variance), by creating
objects of class [`ParamsBM`](@ref)`<:ParamsProcess`.
```@example tree_trait
params_trait1 = ParamsBM( 2, 0.5) # BM with mean  2 and variance 0.5
params_trait2 = ParamsBM(-2, 1)   # BM with mean -2 and variance 1.0
nothing # hide
```
We then simulate the independent traits according to these parameters, using
function [`simulate`](@ref) (fixing the seed, for reproducibility).
```@example tree_trait
using Random
Random.seed!(18480224);
sim1 = simulate(truenet, params_trait1) # simulate a BM on truenet
sim2 = simulate(truenet, params_trait2)
nothing # hide
```
This creates objects of class [`TraitSimulation`](@ref), from which we can
extract the data at the tips, thanks to the method
[`getindex(::TraitSimulation, ::Symbol)`](@ref).
```@example tree_trait
trait1 = sim1[:Tips] # trait 1 at the tips (data)
trait2 = sim2[:Tips]
nothing # hide
```
This extractor creates an `Array` with one column, and as many lines as the
number of tips there are in the phylogeny.  It is sorted in the same order as
the tips of the phylogeny used to simulate it.  
If needed, we could also extract the simulated values at the internal nodes
in the network:
```@example tree_trait
sim1[:InternalNodes]
nothing # hide
```

Finally, we generate the last trait correlated with trait 1
(but not trait 2), with phylogenetic noise.
```@example tree_trait
Random.seed!(18700904);
noise = simulate(truenet, ParamsBM(0, 0.1)) # phylogenetic residuals
trait3 = 10 .+ 2 * trait1 .+ noise[:Tips] # trait to study. independent of trait2
nothing # hide
```

## Phylogenetic regression

Assume that we measured the three traits above, and that we wanted to study the
impact of traits 1 and 2 on trait 3. To do that, we can perform a phylogenetic
regression.

In order to avoid confusion, the function takes in a `DataFrame`, that has an
extra column with the names of the tips of the network, labeled `tipNames`.
Here, we generated the traits ourselves, so they are all in the same order.
```@repl tree_trait
using DataFrames
dat = DataFrame(trait1 = trait1, trait2 = trait2, trait3 = trait3,
                tipNames = tipLabels(sim1))
```

Phylogenetic regression / ANOVA is based on the
[GLM](https://github.com/JuliaStats/GLM.jl) package, with the network as an
extra argument, using function [`phyloNetworklm`](@ref).
```@repl tree_trait
using StatsModels # for statistical model formulas
fitTrait3 = phyloNetworklm(@formula(trait3 ~ trait1 + trait2), dat, truenet)
```
From this, we can see that the intercept, the coefficient for trait 1
and the variance of the noise are correctly estimated
(given that there are only 6 taxa).
In addition, the Student T test for the coefficient
associated with trait 2 has a high p-value, which means that this coefficient
is not significantly different from 0. This is consistent with the
way we simulated trait 3.

The function returns an object of type [`PhyloNetworkLinearModel`](@ref)`<:GLM.LinPredModel`.
It is a subtype of the GLM type `LinPredModel`, which means that all base
functions from Julia [StatsBase](https://github.com/JuliaStats/StatsBase.jl) can
be applied to it. See the documentation for this type for a list of all
functions that can be used. Some functions allow the user to retrieve directly
the estimated parameters of the BM, and are specific to this object.
```@repl tree_trait
sigma2_estim(fitTrait3) # estimated variance of the BM
mu_estim(fitTrait3) # estimated root value of the BM
```

## Ancestral State Reconstruction


### From known parameters

If we assume that we know the exact model of evolution that generated the
traits, we can do ancestral trait reconstruction. Here, we simulated trait 1
ourselves, so we can use the true process, with the true parameters.
In other words, we can reconstruct the state at the internal nodes,
given the values at the tips, the known value at the root
and the known BM variance.
```@example tree_trait
ancTrait1 = ancestralStateReconstruction(truenet, trait1, params_trait1)
nothing # hide
```
Function [`ancestralStateReconstruction`](@ref) creates an object with type
[`ReconstructedStates`](@ref). Several extractors can be applied to it:
```@repl tree_trait
expectations(ancTrait1) # predictions
using StatsBase # for stderror(), aic(), likelihood() etc.
stderror(ancTrait1) # associated standard errors
predint(ancTrait1, level=0.9) # prediction interval (with level 90%)
```
We can plot the ancestral states or prediction intervals on the tree, using the
`nodeLabel` argument of the `plot` function.
```@example tree_trait
ancExpe = expectationsPlot(ancTrait1); # format expected ancestral states for the plot
R"svg(name('ancestral_expe.svg'), width=8, height=4)" # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, :R, nodeLabel = ancExpe);
R"dev.off()" # hide
nothing # hide
```
![ancestral_expe](../assets/figures/ancestral_expe.svg)

```@example tree_trait
ancInt = predintPlot(ancTrait1) # format the prediction intervals for the plot
R"svg(name('ancestral_predint.svg'), width=8, height=4)" # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet,:R, nodeLabel = ancInt);
R"dev.off()" # hide
nothing # hide
```
![ancestral_predint](../assets/figures/ancestral_predint.svg)

The `predint` and `predintPlot` functions have an optional argument to state
the `level` of the prediction interval. If not given, the default value is
0.95.

It is also possible to plot both the reconstructed state and the predicted value
on the same plot, using the optional keyword argument `withExp`.
As shown below, we could also use the `RCall` method from the
[`plot`](https://cecileane.github.io/PhyloPlots.jl/stable/lib/public/) function.
```@example tree_trait
plot(truenet, :R, nodeLabel = predintPlot(ancTrait1, withExp=true));
nothing # hide
```
These plots tend to be quite busy, even for small networks.

As we know the true ancestral states here, we can compare them to our
estimation.
```@repl tree_trait
predictions = DataFrame(infPred=predint(ancTrait1)[1:7, 1],
                        trueValue=sim1[:InternalNodes],
                        supPred=predint(ancTrait1)[1:7, 2])
```

### From estimated parameters

In real applications though, we do not have access to the true parameters of the
process that generated the data. We can estimate it using the previous function.
To fit a regular BM, we just need to do a regression of trait 1 against a simple
intercept:
```@example tree_trait
fitTrait1 = phyloNetworklm(@formula(trait1 ~ 1), dat, truenet)
nothing # hide
```
We can then apply the [`ancestralStateReconstruction`](@ref) function directly
to the fitted object:
```@example tree_trait
ancTrait1Approx = ancestralStateReconstruction(fitTrait1)
nothing # hide
```
The prediction intervals ignore the fact that we estimated the process
parameters, so they are less accurate and the function throws a warning.
The output is an object of the same [`ReconstructedStates`](@ref) type as earlier,
and the same extractors can be applied to it:
```@example tree_trait
R"svg(name('ancestral1.svg'), width=8, height=4)" # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, :R, nodeLabel = expectationsPlot(ancTrait1Approx));
R"dev.off()" # hide
nothing # hide
```
![ancestral1](../assets/figures/ancestral1.svg)

For convenience, the two steps described above (fitting against the
intercept, and then do ancestral state reconstruction) can be done all at once
with a single call of the function [`ancestralStateReconstruction`](@ref) on a
DataFrame with the trait to reconstruct, and the tip labels:
```@example tree_trait
datTrait1 = DataFrame(trait1 = trait1, tipNames = tipLabels(sim1))
ancTrait1Approx = ancestralStateReconstruction(datTrait1, truenet)
nothing # hide
```
```@example tree_trait
R"svg(name('ancestral2.svg'), width=8, height=4)" # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, :R, nodeLabel = predintPlot(ancTrait1Approx, level=0.9));
R"dev.off()" # hide
nothing # hide
```
![ancestral2](../assets/figures/ancestral2.svg)

This produces the exact same results. Here, we chose a `level` of 90% for the
plotted prediction intervals.

### Data imputation

Note that there is no theoretical difference between an internal node, for which
we could not measure the value of the trait, and a missing value at a tip of the
network. Consequently, the previous [`ancestralStateReconstruction`](@ref)
function can be used to do data imputation. To see this, let's add some missing
values in trait 1.
```@example tree_trait
allowmissing!(datTrait1, :trait1)
datTrait1[2, :trait1] = missing; # second row: for taxon C
ancTrait1Approx = ancestralStateReconstruction(datTrait1, truenet)
nothing # hide
```
```@example tree_trait
R"svg(name('ancestral3.svg'), width=8, height=4)" # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, :R, nodeLabel = predintPlot(ancTrait1Approx));
R"dev.off()" # hide
nothing # hide
```
![ancestral3](../assets/figures/ancestral3.svg)

A prediction interval is shown for the missing values.

### With known predictors

At this point, it might be tempting to apply this function to trait 3 we
simulated earlier as a linear combination of trait 1 and a phylogenetic
noise. However, this cannot be done directly:
```julia
ancTrait3 = ancestralStateReconstruction(fitTrait3) # Throws an error !
```
This is because the model we used to fit the trait (a regression with one
predictor and an intercept) is not compatible with the simple model of Brownian
evolution that we assumed for the ancestral state reconstruction. As the
predictor used is not known for ancestral states, it is not possible to
reconstruct the trait for this particular model.

The only option we have is to provide the function with the predictor's
ancestral states, if they are known. They are known indeed in this
toy example that we generated ourselves, so we can reconstruct our trait
doing the following:
```@example tree_trait
ancTrait3 = ancestralStateReconstruction(fitTrait3,
              [ones(7, 1) sim1[:InternalNodes] sim2[:InternalNodes]])
nothing # hide
```
```@example tree_trait
R"svg(name('ancestral4.svg'), width=8, height=4)" # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, :R, nodeLabel = predintPlot(ancTrait3));
R"dev.off()" # hide
nothing # hide
```
![ancestral4](../assets/figures/ancestral4.svg)

where we provided the ancestral predictors as a matrix, containing the
intercept, and the known predictor at the nodes. The user must be very careful
with this function, as no check is done for the order of the predictors, that
must be in the same order as the internal nodes of the phylogeny. As ancestral
predictors are often unknown, the use of this functionality is discouraged.


## Phylogenetic ANOVA

The [`phyloNetworklm`](@ref) function is based on the `lm` function
from [GLM](https://github.com/JuliaStats/GLM.jl). This means that it
inherits from most of its features, and in particular, it can handle formulas
with factors or interactions.
For example, in lizards, we might want to do a regression of toe length against
body length and the region where each species is found, where this region is coded
into 4 categories (say). We might also want to include an interaction effect
between body length and region.
(This model has no biological basis. It is just meant to show the possibilities
of the function).

To illustrate the use of categorical predictors of particular interest
in a network with reticulations, let's assume that some transgressive evolution took place
after the hybridization event, so that tips "A" and "B" have larger mean
compared to the others
(see [^B18] for transgressive evolution after a reticulation event).
```@example tree_trait
delta = 5.0; # value of heterosis
underHyb = [(n == "A" || n == "B") for n in tipLabels(sim1)] # tips under hybrid
underHyb
for i in 1:length(trait3)
    underHyb[i] && (trait3[i]+=delta) # add delta to tips A and B
end
nothing # hide
```
```@repl tree_trait
trait3 # changed: +5 was added by the previous loop to A and B
```
The categorical variable `underHyb` separates tips "A" and "B" from the others.
We need to mark it as a categorical variable, not a numerical variable,
i.e. as a `PooledDataArray`.
```@example tree_trait
dat = DataFrame(trait1 = trait1, trait2 = trait2, trait3 = trait3,
                underHyb = underHyb,
                tipNames = tipLabels(sim1))
categorical!(dat, :underHyb)
nothing # hide
```
```@repl tree_trait
dat
```
Now we can include this reticulation variable in the regression.
```@example tree_trait
fitTrait = phyloNetworklm(@formula(trait3 ~ trait1 + underHyb), dat, truenet)
```
In this case, the categorical variable indicating which tips are descendants
of the reticulation event is indeed relevant, and the transgressive evolution effect
is recovered.

This is a very simple example of how to include transgressive evolution,
but some general
functions to test for it, on networks with more than on hybrid, are also
available.


## Pagel's Lambda

One classical question about trait evolution is the amount of
"phylogenetic signal" in a dataset, that is, the importance of the tree
structure to explain variation in the observed traits.
One way of doing measuring that is to use
Pagel's lambda transformation of the branch lengths [^P99].
This model assumes a
BM on a tree where the internal branches are multiplied by a factor λ,
while the external branches are modified so that the total height of the tree is
constant. Hence, λ varies between 0 (the tree has no influence on
the data) and 1 (the tree is unchanged).
Using the same branch length transformations, this model can
be straightforwardly extended to phylogenetic networks.

We can illustrate this with the predictor trait we used earlier. We use the
same function as before, only indicating the model we want to use:
```@example tree_trait
fitPagel = phyloNetworklm(@formula(trait1 ~ 1), dat, truenet, model="lambda")
```
As it is indeed generated according to a plain BM on the phylogeny, the
estimated λ should be close to 1. It can be extracted with function
`lambda_estim`:
```@repl tree_trait
lambda_estim(fitPagel)
```

## Shifts and transgressive evolution

In the ANOVA section above, we showed how to include transgressive evolution
in a simple case.
In general, transgressive evolution can be seen as a particular example
of a *shifted BM* on the phylogenetic network.

### Simulation of a Shifted BM

In a shifted BM, the trait evolves as a BM on the network most of
the time, but *shifts* on some of the branches.
The positions and values of the shifts can be stored in a [`ShiftNet`](@ref)
object. For identifiability reasons, shifts are only allowed on tree-like
branches. The position of the shifts can be given using vector of edges.
To see this, let's first plot the network with its associated edges and node
numbers.
```@example tree_trait
R"svg(name('truenet_with_numbers.svg'), width=8, height=4)" # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, :R, useEdgeLength=true, showEdgeNumber=true);
R"dev.off()" # hide
nothing # hide
```
![truenet_with_numbers](../assets/figures/truenet_with_numbers.svg)

Let's say that we want to add a shift with value 5.0 on the branch directly
following the hybridization event, in order to model transgressive evolution.
We can see on the
plot above that this branch is number 6, so we define the following object:
```@example tree_trait
shift = ShiftNet(truenet.edge[6], 5.0,  truenet)
nothing # hide
```
Note that the edge numbers and values of a `ShiftNet` object can be retrieved
thanks to functions [`getShiftEdgeNumber`](@ref) and [`getShiftValue`](@ref).
The constructor can take a single edge and associated value, like here,
or two vectors of edges and matching values.

Because we often need to put shifts only on edges right after hybrids,
there is a special function [`shiftHybrid`](@ref) to do that, so that 
we do not have to find out their edges number. Here, the `shift` object
could hence have been defined as:
```@example tree_trait
shift = shiftHybrid(5.0,  truenet)
```

The parameters for the simulation are then defined as above, just adding
the `ShiftNet` object as a parameter.

```@example tree_trait
params_sh = ParamsBM(2, 0.5, shift) # BM with mean 2, variance 0.5, and shifts.
nothing # hide
```
The traits are simulated using the same function [`simulate`](@ref), and
extracted at the tips as before.
```@example tree_trait
Random.seed!(18700904)
sim_sh = simulate(truenet, params_sh) # simulate a shifted BM on truenet
trait_sh = sim_sh[:Tips]              # trait at the tips (data)
nothing # hide
```

### Fit of a Shifted BM

Let's assume that we measured `trait_sh`, and that we want to test whether
there were some ancestral hybridizations. To do that, we can use the 
custom columns of the [`descendenceMatrix`](@ref), that can be directly
defined thanks to function [`regressorHybrid`](@ref).
```@example tree_trait
df_shift = regressorHybrid(truenet) # Regressors matching Hybrid Shifts
nothing # hide
```
This creates a dataframe, with as many columns as the number of hybrids
in the network, each named according to the number of the edge after the
hybrid.
We can use this dataframe as regressors in the `phyloNetworklm` function.

```@example tree_trait
dat = DataFrame(trait = trait_sh, tipNames = tipLabels(sim_sh))  # Data
dat = join(dat, df_shift, on=:tipNames)                          # join the two
fit_sh = phyloNetworklm(@formula(trait ~ shift_6), dat, truenet) # fit
```
Here, because there is only one hybrid in the network, we can directly
see whether the ancestral transgressive evolution is significant or not thanks to the
Student T test on the coefficient associated with `shift_6`. In more
complex cases, it is possible to do a Fisher F test, thanks to the `GLM`
function `ftest`.
```@example tree_trait
fit_null = phyloNetworklm(@formula(trait ~ 1), dat, truenet) # fit against the null (no shift)
ftest(fit_sh, fit_null)                                      # nested models, from more complex to most simple
```
Here, this test is equivalent to the Fisher F test, and gives the same p-value.

Note that, for conventional reasons, the `ftest` function always takes the
*most complex* model as the first one. This means that, in the table of
results, the models are actually named in a reverse order, so that "Model 2" is
actually our model under H₀ (null model), and "Model 1" the one under H₁
(model with shifts).

---

### References

[^B18]: Bastide, Solís-Lemus, Kriebel, Sparks, Ané (2018):
    Phylogenetic Comparative Methods for Phylogenetic Networks with Reticulations.
    Systematic Biology 67(5):800–820. doi:10.1093/sysbio/syy033

[^P99]: Pagel M (1999). Inferring the historical patterns of biological
    evolution. Nature. 401: 877–884. doi:10.1038/44766
