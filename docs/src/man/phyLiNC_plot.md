```@setup snaqplot
using PhyloNetworks
mkpath("../assets/figures")
fastafile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","h1_net.fasta")
startingtree = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","h1_net.treefile")
net = readTopology(joinpath(dirname(pathof(PhyloNetworks)), "..","examples","h1.net"))
```

# Estimate a Network from Concatenated DNA Sequences

## Network Estimation
After estimating a starting tree, we can estimate a network.

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

## Beyond Defaults
To see the defaults in phyLiNC, see the documentation for [`phyLiNC!`](@ref).
In many cases, you may want to move beyond these defaults for a more customized
estimation.

TODO

## Network Visualization
We can plot the resulting network with PhyloPlots like this.

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



## parallel computations

By default, phyLiNC estimates 10 networks, then provides the user with the best one.
To speed this process, you can run this in parallel. For example, if your machine has 4 or more processors (or cores),you can tell julia to use 4 processors by starting julia with `julia -p 4`,
or by starting julia the usual way (`julia`) and then adding processors with:

```julia
using Distributed
addprocs(4)
```

If we load a package (`using PhyloNetworks`) before adding processors,
then we need to re-load it again so that all processors have access to it. Then,
simply call phyLiNC as usual.

```julia
@everywhere using PhyloNetworks
net = PhyloNetworks.phyLiNC!(startingtree, fastafile, :JC69)
```