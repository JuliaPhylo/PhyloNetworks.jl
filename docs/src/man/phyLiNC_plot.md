```@setup snaqplot
using PhyloNetworks
mkpath("../assets/figures")
fastafile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","h1_net.fasta")
startingtree = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","h1_net.treefile")
net = readTopology(joinpath(dirname(pathof(PhyloNetworks)), "..","examples","h1.net"))
```

# Estimate a Network from Concatenated DNA Sequences using PhyLiNC

Estimate a phylogenetic network from concatenated DNA data using
maximum likelihood, ignoring incomplete lineage sorting.
Starting at a given starting tree or network, PhyLiNC uses a hill-climbing algorithm
to search for the best network structure that best fits the given DNA sequences.
(phyLiNC: phylogenetic Likelihood Network from Concatenated data)

The hill-climbing search is made up of nearest neighbor interchange moves,
root change moves, add hybrid, and remove hybrid proposals.

## Getting Started
To estimate a phylogenetic network with PhyLiNC, you'll need concatenated DNA
sequence data and a starting network.

Data: The DNA sequences should be in fasta format.

Starting Network or Tree: The starting network can be built using a tree-building
tool like IQTree or inferred using the analyst's knowledge of the field. It should be
in Newick string format. It should then be read in using PhyloNetwork's ReadTopology()
function.

Usage:
The network is constrained to have `maxhybrid` reticulations at most,
but can be of any level.
The search starts at (or near) the network `net`,
using a local hill-climbing search to optimize the topology
(nearest-neighbor interchange moves, add hybridizations,
and remove hybridizations). Also optimized are evolutionary rates,
amount of rate variation across sites, branch lengths and inheritance γs.
This search strategy is run `nruns` times, and the best of the `nruns`
networks is returned.

Return a [`StatisticalSubstitutionModel`](@ref) object, say `obj`, which
contains the estimated network in `obj.net`.

Required arguments:
- `net`: a network or tree of type `HybridNetwork`, to serve as a starting point
  in the search for the best network.
  Newick strings can be converted to this format with [`readTopology`] (@ref).
- `fastafile`: file with the sequence data in FASTA format.
- `substitutionModel`: A symbol indicating which substitution model is used.
  Choose `:JC69` f [`JC69`] (@ref) for the Jukes-Cantor model or `:HKY85` for
  the Hasegawa, Kishino, Yano model [`HKY85`] (@ref).

## Network Estimation

By default, only three arguments are required for estimating a network with phyLiNC:
a starting tree, a file of concatenated DNA sequences, and a subsitution model.
- `net`: a network or tree of type `HybridNetwork`, to serve as a starting point
  in the search for the best network.
  Newick strings can be converted to this format with [`readTopology`] (@ref).
- `fastafile`: file with the sequence data in FASTA format.
- `substitutionModel`: A symbol indicating which substitution model is used.
  Choose `:JC69` f [`JC69`] (@ref) for the Jukes-Cantor model or `:HKY85` for
  the Hasegawa, Kishino, Yano model [`HKY85`] (@ref).

```julia
net = phyLiNC!(startingtree, fastafile, :JC69)
```
At the beginning of optimization, the screen will show this: TODO add

At the end, it will give a final topology, as well as a model of evolution over
the network: TODO add

## Beyond Defaults: Customize Optimization with Optional Arguments
In many cases, you may want to move beyond these defaults for a more customized
estimation.

PhyLiNC's optional arguments (default value in parenthesis):
-  symbol for the model of rate variation across sites
  (`:noRV` for no rate variation):
  use `:G` or `:Gamma` for Gamma-distributed rates,
  `:I` or `:Inv` for a proportion of invariable sites, and
  `:GI` or `:GammaInv` for a combination (not recommended).
- integer (4) for the number of categories to use in estimating
  evolutionary rates using a discretized gamma model. When allowing for rate
  variation, four categories is standard. With 1 category,
  no rate variation is assumed. See [`RateVariationAcrossSites`] (@ref).

Main optional keyword arguments (default value in parenthesis):
- `speciesfile` (""): path to a csv file with samples in rows and two columns:
  species (column 1), individual (column 2)
  Include this file to group individuals by species.
- `cladefile` (""): path to a csv file containing two columns:
  clades and individuals used to create one or more clade topology
  constraints to meet during the search.
  (NOTE: clade contraints not yet implemented.)
- `filename` ("phyLiNC"): root name for the output files (`.out`, `.err`).
  If empty (""), files are *not* created, progress log goes to the screen only
  (standard out).
- `maxhybrid` (1): maximum number of hybridizations allowed.
  `net` (starting network) must have `maxhybrid` or fewer reticulations.
- `nruns` (10): number of independent starting points for the search
- `no3cycle` (true): prevents 3-cycles, which are (almost) not
  identifiable
- `nohybridladder` (true): prevents hybrid ladder in network. If true,
  the input network must not have hybrid ladders.
  Warning: Setting this to true will avoid most hybrid ladders, but some can still
  occur in some cases when deleting hybrid edges. This might be replaced with an
  option to avoid all non-tree-child networks in the future.

Optional arguments controlling the search:
- `seed` (default 0 to get it from the clock): seed to replicate a given search
- `nreject` (75): maximum number of times that new topologies are
  proposed and rejected in a row. Lower values of `nreject` result in a less
  thorough but faster search. Controls when to stop proposing new
  network topologies.
- `probST` (0.5): probability to use `net` as the starting topology
  for each given run. If probST < 1, the starting topology is k NNI moves
  away from `net`, where k is drawn from a geometric distribution: p (1-p)ᵏ,
  with success probability p = `probST`.
- `maxmoves` (100): maximum number of topology moves before branch lengths,
  hybrid γ values, evolutionary rates, and rate variation parameters are
  reestimated.
- `verbose` (true): set to false to turn off screen output
- `alphamin` (0.02): minimum value for shape parameter alpha in rate variation
  across sites model.
- `alphamax` (50.0): maximum value for shape parameter alpha in rate variation
  across sites model.

The following optional arguments control when to stop the optimization of branch
lengths and gamma values on each individual candidate network. Defaults in
parentheses.
- `ftolRel` (1e-6) and `ftolAbs` (1e-6): relative and absolute differences of the
  network score between the current and proposed parameters
- `xtolRel` (1e-5) and `xtolAbs` (1e-5): relative and absolute differences
  between the current and proposed parameters.
Greater values will result in a less thorough but faster search. These parameters
are used when evaluating candidate networks only.
Regardless of these arguments, once a final topology is chosen, branch lenghts
are optimized using stricter tolerances (1e-10, 1e-12, 1e-10, 1e-10) for better
estimates.

## Estimating in Parallel

By default, phyLiNC estimates ten networks, then returns the network with the highest
likelihood. To speed this process, these ten runs can be done in parallel.
For example, if your machine has more than four processors (or cores), you can tell
julia to add 4 worker processors by starting julia with `julia -p 4`. This will
give Julia five processors: one main and four workers.
Alternatively, you can start Julia as usual (`julia`), then adding processors with:

```julia
using Distributed
addprocs(4)
```

If you're unsure how many you've added, youcan check how many processors Julia
currently has using `nprocs()`.

```julia
nprocs()
```

After adding processors, whether with `-p` or `addprocs()`, be sure to reload any
packages you might need using the Distributed package's macro `@everywhere`. This
ensures that all processors will have access to the packages.

After adding these processors, PhyLiNC will automatically run any estimation in
this session using all cores.

```julia
@everywhere using PhyloNetworks
net = PhyloNetworks.phyLiNC!(startingtree, fastafile, :JC69)
```

TODO show output?

## Network Vizualization
To plot your estimated networksm you can use the visualization package, PhyloPlots.

```@example snaqplot
using PhyloPlots
using RCall # hide
R"name <- function(x) file.path('..', 'assets', 'figures', x)" # hide
R"svg(name('snaqplot_net0_1.svg'), width=4, height=3)" # hide
R"par"(mar=[0,0,0,0]) # hide
plot(net0, :R);
R"dev.off()"; # hide
nothing # hide
```

## Branch Length Estimation and Identifiability
The length of the edge below a reticulation is not identifiable.
Therefore, phyLiNC estimates the canonical version of the network: with
reticulations **unzipped**: edges below reticulations are set to 0, and
hybrid edges (parental lineages) have estimated lengths that are
increased accordingly.

If any branch lengths are missing in the input network, phyLiNC estimates the
starting branch lengths using pairwise distances. Otherwise, it uses the input branch
lengths as starting branch lengths, only unzipping all reticulations, as
described above.

## Choosing a Substitution Model
PhyLiNC allows you to choose from two Markov Process-based Evolutionary Rate Substitution models and four options for adding rate variation. These can be combined for eight possible ways to model evolutionary rates.

### Markov Model Substitution Models
There are two possible nucleic acid substitution models to choose from: JC69 (1969) and HKY85 (Hasegawa et al. 1985). The relative version of the
model is used.

JC69: Jukes Cantor (1969) nucleic acid substitution model has a single rate parameter.
`rate` corresponds to the absolute diagonal elements, that is, the rate of change
(to any of the other 2 states). Individual rates are `rate`/3.
The transition matrix [`Q`](@ref) is normalized to an average of 1 transition per unit of time: in which case `rate` is set to 1.0.

HKY85: A nucleic acid substitution model based on Hasegawa et al. 1985 substitution model.
`rate` should be a vector of 1 or 2 rates, and `pi` a vector of 4 probabilities summing to 1.
If `relative` is false, the 2 rates represent the transition rate and the transversion rate,
α and β. If `relative` is true (default), only the first rate is used and represents the transition/transversion ratio: κ=α/β. The rate transition matrix Q is normalized to have 1 change / unit of time on average, i.e. the absolute version of Q is divided by
`2(piT*piC + piA*piG)α + 2(piY*piR)β`.

Withever of these models, you choose, PhyLiNC will optimize the rate paramters
throughout network estimation, then show the transition rate matrix `Q` for the
final network when optimization is complete. The model is stored in the object
that phyLiNC returns, which we've called `obj` here. To explore the model further,
use the following commands.

```julia
julia> model = obj.model
HKY85 Substitution Model base frequencies: [0.30095987201706437, 0.19137448340221302, 0.18764164778029596, 0.3200239968004266]
relative rate version with transition/tranversion ratio kappa = 4.03502,
 scaled so that there is one substitution per unit time
rate matrix Q:
               A       C       G       T
       A       *  0.1320  0.5223  0.2207
       C  0.2076       *  0.1294  0.8907
       G  0.8377  0.1320       *  0.2207
       T  0.2076  0.5327  0.1294       *

julia> nstates(model)
4
julia> getlabels(obj.model)
4-element Array{BioSymbols.DNA,1}:
 DNA_A
 DNA_C
 DNA_G
 DNA_T
```

### Rate Variation and Invariable Sites
Concatenated models are well-suited to data with high percentages of invariable
sites (TODO add citation). To this end, we provide functionality to estimate a
rate of invariant sites using the invariable-sites model
(+I, Hasegawa, Kishino & Yano 1985 J Mol Evol). To enable this, run PhyLiNC like
this:

```julia
net = phyLiNC!(startingtree, fastafile, :JC69, :Inv)
```

We can also estimate a variable rates across sites using
the discrete Gamma model (+G, Yang 1994, Journal of Molecular Evolution). To use
this model using the :G symbol in your PhyLiNC call:

```julia
net = phyLiNC!(startingtree, fastafile, :JC69, :G)
```

You can combine the two types of rate variation (+G+I, Gu, Fu & Li 1995, Mol Biol Evol)
but this is discouraged (Jia, Lo & Ho 2014 PLOS One). Using rate variation
increases the number of parameters by one (+G or +I) or by two (+G+I).
Because the mean of the desired distribution or rates is 1, we use a Gamma
distribution with shape α and scale θ=1/α (rate β=α) if no invariable sites,
or scale θ=1/(α(1-pinv)), that is rate β=α(1-pinv) with a proportion pinv
of invariable sites.

The shape parameter is referred to as alpha here.

The Gamma distribution is discretized into `ncat` categories.

In each category, the category's rate multiplier is a normalized quantile of the gamma distribution.

The rate model is also returned by PhyLiNC. To explore it further, use the
following commands.
For a network estimation with rate variation across sites:
```julia
juila> ratemodel = obj.ratemodel
Rate variation across sites: discretized Gamma
alpha: 2.0
categories for Gamma discretization: 4
rates: [0.319, 0.683, 1.109, 1.889]

julia> nparams(ratemodel)
1
```

For an invariable sites model:
```julia
juila> ratemodel = obj.ratemodel
Rate variation across sites: +I (invariable sites)
pinv: 0.3
rates: [0.0, 1.429]

julia> nparams(ratemodel)
1
```

If the network has both types of rate variation, the output will look like this:
```julia
juila> ratemodel = obj.ratemodel
Rate variation across sites: discretized Gamma+I
pinv: 0.3
alpha: 2.0
categories for Gamma discretization: 4
rates: [0.0, 0.456, 0.976, 1.584, 2.698]
probabilities: [0.3, 0.175, 0.175, 0.175, 0.175]

julia> nparams(ratemodel)
2
```

## Limitations
TODO
