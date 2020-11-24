```@setup phylinc
using PhyloNetworks
mkpath("../assets/figures")
fastafile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","h1_net.fasta")
startingtreefile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","h1_net.treefile")
```

# Introduction
The phylogenetic network estimation method, PhyLiNC (phylogenetic Likelihood Network from Concatenated data) is part of the `PhyloNetworks` Julia package.

Assumptions on the phylogenetic network. First, we assume that the evolutionary
history of the set of taxa is not heavily impacted by incomplete lineage sorting.
Second, assume that the data we use to infer the network are independent. That is, we assume that sites within the concatenated sequences are not grouped into loci.
Third, we assume that the network does not include parallel edges, also known as
two-cycles. These are not or nearly not-identifiable under our model and cannot
be inferred. Finally, a hybrid node can have only two parent edges, meaning that
each hybridization event can only occur between two species.

# Estimate a Network from Concatenated DNA Sequences using PhyLiNC

Starting at a given starting tree or network, PhyLiNC uses a hill-climbing algorithm
to search for the best network structure that best fits the given DNA sequences.
(phyLiNC: phylogenetic Likelihood Network from Concatenated data)

The hill-climbing search is made up of nearest neighbor interchange moves,
root change moves, add hybrid, and remove hybrid proposals.

## Getting Started
To estimate a phylogenetic network with PhyLiNC, you'll need concatenated DNA
sequence data, a starting network, and a substitution model.

## Data
The DNA sequences should be in fasta format. You'll need a path to the file.

## Starting Network or Tree
The starting network can be built using a tree-building tool like IQTree or
inferred using the analyst's knowledge of the field. It should be
in Newick string format. It should then be read in using PhyloNetwork's ReadTopology()
function. The network is constrained to have `maxhybrid` reticulations at most
(or one reticulation , if `maxhybrid` is not provided). It can be of any level.
## Estimate a Network
```julia
startingnet = readTopology(startingtreefile)
net = PhyloNetworks.phyLiNC(startingtree, fastafile, :HKY85; maxhybrid=1, nreject=5)
```
## Choose a Markov Substitiution Model
A symbol indicating which substitution model is used. Choose `:JC69` [`JC69`] (@ref)
for the Jukes-Cantor model or `:HKY85` for the Hasegawa, Kishino, Yano model [`HKY85`] (@ref).

## Additional Options
- `maxhybrid`: the maximum number of hybridizations to include in the estimated
  network. Allows one hybridization by default.
- `nreject`: a number controlling how precise the estimation will be. A higher
  value will lead to more precise estimation and higher computing times. For this
  example, for the sake of timing, we've kept it small, at 5.

```julia
net = PhyloNetworks.phyLiNC(startingtree, fastafile, :HKY85; maxhybrid=1, nreject=5)
```

The search starts at (or near) the network `net`,
using a local hill-climbing search to optimize the topology
(nearest-neighbor interchange moves, add hybridizations,
and remove hybridizations). Evolutionary rates,
amount of rate variation across sites, branch lengths and inheritance inheritance
weights are also optimized.
This search strategy is run `nruns` times, and the best of the `nruns`
networks is returned.

At the beginning of optimization, phyLiNC will write to the screen "PhyLiNC
network estimation starting." and provide a list of parameters.

Once estimation is complete, phyLiNC will write a final topology, as well as a
model of evolution over the network. For this example it looks like this:
Final network:
(2:1.0e-8,1:1.0e-8,((3:1.0e-8,4:0.0966237305741647):0.20868719647880038,(5:1.0e-8,6:0.09662373126075129):1.0e-8):0.20802463974006694);
Total time elapsed: 12.8 seconds (includes final branch length and gamma optimization)
Final time: 2020-11-18 16:56:7.316
---------------------
NOTE: This network should be interpreted as unrooted, since a likelihood model cannot identify
the root. The network can be rooted with an outgroup or other external information.
---------------------
PhyloNetworks.StatisticalSubstitutionModel:
HKY85 Substitution Model base frequencies: [0.15714285714285714, 0.34285714285714286, 0.32857142857142857, 0.17142857142857143]
  relative rate version with transition/tranversion ratio kappa = 2.35056,
   scaled so that there is one substitution per unit time
  rate matrix Q:
                 A       C       G       T
         A       *  0.3366  0.7582  0.1683
         C  0.1543       *  0.3226  0.3956
         G  0.3626  0.3366       *  0.1683
         T  0.1543  0.7912  0.3226       *
on a network with 0 reticulations
data:
  6 species
  11 sites
  9 distinct patterns
log-likelihood: -37.5133

## Beyond Defaults: Customize the Network's Optimization with Optional Arguments
In many cases, you may want to move beyond this simple example for a more customized
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

## Rate Variation and Invariable Sites
By default, the model does not allow rate variation across sites. To add rate
variation to your selection Markov model of evolution, add one of these optional
arguments:
  use `:G` or `:Gamma` for Gamma-distributed rates,
  `:I` or `:Inv` for a proportion of invariable sites, and
  `:GI` or `:GammaInv` for a combination (not recommended).

To specify the number of categories to use in estimating evolutionary rates
using a discretized gamma model. (When allowing for rate variation, four
categories is standard. With 1 category, no rate variation is assumed.)
For more details, see [`RateVariationAcrossSites`] (@ref).

To allow for variable rates across sites using
the discrete Gamma model (+G, Yang 1994, Journal of Molecular Evolution). To use
this model using the :G symbol in your PhyLiNC call:

```julia
obj = PhyloNetworks.phyLiNC(startingtree, fastafile, :JC69, :G; nreject = 2)
```

Our model is well-suited to data with a high percentages of invariable sites. To
this end, we provide functionality to estimate a rate of invariant sites using
the invariable-sites model (+I, Hasegawa, Kishino & Yano 1985 J Mol Evol). To
add the +I to your Markov model, run phyLiNC with the :Inv symbol.

```julia
obj_I = PhyloNetworks.phyLiNC(startingtree, fastafile, :JC69, :Inv; nreject = 2)
```

You can combine the two types of rate variation (+G+I, Gu, Fu & Li 1995, Mol Biol Evol)
but this is discouraged (Jia, Lo & Ho 2014 PLOS One). Using rate variation
increases the number of parameters by one (+G or +I) or by two (+G+I).
Because the mean of the desired distribution or rates is 1, we use a Gamma
distribution with shape α and scale θ=1/α (rate β=α) if no invariable sites,
or scale θ=1/(α(1-pinv)), that is rate β=α(1-pinv) with a proportion pinv
of invariable sites.
```julia
obj_GI = PhyloNetworks.phyLiNC(startingtree, fastafile, :JC69, :GI; nreject = 2)
```

After running phyLiNC, the rate model will be returned. For a network estimated
with gamma rate variation across sites, you could explore this model in this way.
```julia
juila> ratemodel = obj_G.ratemodel
Rate variation across sites: discretized Gamma
alpha: 2.0
categories for Gamma discretization: 4
rates: [0.319, 0.683, 1.109, 1.889]

julia> nparams(ratemodel)
1
```

For an invariable sites model, it might look like this:
```julia
juila> ratemodel = obj_I.ratemodel
Rate variation across sites: +I (invariable sites)
pinv: 0.3
rates: [0.0, 1.429]

julia> nparams(ratemodel)
1
```

# Multiple Alleles or Individuals within a Species
We allow analysts to force groupings of alleles or individuals within a species group.


To plot your estimated networks you can use PhyloNetworks' associated visualization
package, PhyloPlots.
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
nworkers()
```

After adding processors, whether with `-p` or `addprocs()`, be sure to reload any
packages you might need using the Distributed package's macro `@everywhere`. This
ensures that all processors will have access to the packages.

After adding these processors, PhyLiNC will automatically run any estimation in
this session using all cores.

```julia
@everywhere using PhyloNetworks
net = PhyloNetworks.phyLiNC(startingtree, fastafile, :JC69; nreject = 2)
```
