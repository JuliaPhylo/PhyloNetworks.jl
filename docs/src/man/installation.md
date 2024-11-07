# Installation

## Install Julia

Julia is a high-level and interactive programming language (like R or Matlab),
but it is also high-performance (like C).
To install Julia, follow instructions [here](http://julialang.org/downloads/).
For a quick & basic tutorial on Julia, see
[learn x in y minutes](http://learnxinyminutes.com/docs/julia/).

Editors:

- [Visual Studio Code](https://code.visualstudio.com) provides an editor
  and an integrated development environment (IDE) for Julia: highly recommended!
  [Positron](https://github.com/posit-dev/positron) is a great (and similar)
  alternative.
- Install the [Julia extension](https://code.visualstudio.com/docs/languages/julia)
  in VS Code or Positron.
- We can also run Julia in a Pluto notebook.
  [Pluto.jl](https://plutojl.org/) is a great to get started with Julia.

Julia code is just-in-time compiled. This means that the
first time we run a function, it will be compiled at that moment.
Future calls to the function will be much faster.
Trying out toy examples for the first calls is a good idea.

## Install PhyloNetworks

To install the package, type inside Julia:
```julia
using Pkg
Pkg.add("PhyloNetworks")
```
If you already installed the package and want
the latest registered version, do this to update all of your packages:
```julia
Pkg.update()
```
It is important to update the package regularly as it is
undergoing constant development. Join the google group for updates
[here](https://groups.google.com/forum/#!forum/phylonetworks-users/new).

`Pkg.update()` will install the latest registered version, but there
could be other improvements in the `master` branch of the
repository. If you want to update to the latest unregistered version
of the package, you can do
`Pkg.add(PackageSpec(name="PhyloNetworks", rev="master"))`
just beware that the latest changes could be not as robust.
If you want to go back to the registered package, you can do
`Pkg.free("PhyloNetworks")`.

Similarly, you can pin a version of the package
`Pkg.pin("PhyloNetworks")` so that `Pkg.update()` will not modify
it. You can always free a pinned package with
`Pkg.free("PhyloNetworks")`. More on package management
[here](https://docs.julialang.org/en/v1/stdlib/Pkg/).

The PhyloNetworks package has dependencies like
[NLopt](https://github.com/jump-dev/NLopt.jl) and
[DataFrames](https://dataframes.juliadata.org/stable/)
(see the `Project.toml` file for the full list), but everything is installed automatically.

The companion package [PhyloPlots](https://github.com/juliaphylo/PhyloPlots.jl)
has utilities to visualize networks, and for interoperability,
such as to export networks to R (which can then be plotted via R).
To install:

```julia
using Pkg
Pkg.add("PhyloPlots")
```

PhyloPlots depends on PhyloNetworks, and has further dependencies
like
[RCall](https://github.com/JuliaInterop/RCall.jl)

## Test example

To check that your installation worked, type this in Julia to load the package.
This is something to type every time you start a Julia session:
```@example install
using PhyloNetworks;
```
This step can also take a while, to pre-compile the code (after a package
update for instance).
Here is a very small test for the installation of PhyloNetworks.

```@repl install
net = readTopology("(A,(B,(C,D)));");
tipLabels(net)
```

You can see a list of all the functions with
```julia
varinfo(PhyloNetworks)
```
and press `?` inside Julia to switch to help mode,
followed by the name of a function (or type) to get more details about it.


## Julia types

Each object in Julia has a *type*. We show here small examples on how to get more
info on an object.
If we want to know the type of a particular object, use `typeof`
For example, let's read a list of gene trees.
First, we need the file name.


```@repl install
raxmltreefile = joinpath(dirname(dirname(pathof(PhyloNetworks))), "examples","raxmltrees.tre")
typeof(raxmltreefile)
```
The object `raxmltreefile` is a basic string (of letters).
Let's create our list of gene trees by reading this file.

```@repl install
genetrees = readmultitopology(raxmltreefile); # the semicolon suppresses info on the result
typeof(genetrees)
```
which shows us that `genetrees` is of type `Vector{HybridNetwork}`, that is,
a vector containing networks.
If we want to know about the attributes the object has, we can type `?` in Julia,
followed by `HybridNetwork` for a description.

Typing `varinfo()` will provide a list of objects and packages in memory,
including `raxmltreefile` and `genetrees` that we just created.

## Quick start

Here we could check the length of our list of gene trees, as a sanity check
to make sure we have all gene trees we expected, and check that the third tree
has whatever taxon names we expected:

```@repl install
length(genetrees)
tipLabels(genetrees[3])
```

We can also see some basic information on the third gene tree, say:
```@repl install
genetrees[3]
```
To visualize any of these gene trees, use the
[PhyloPlots](https://github.com/juliaphylo/PhyloPlots.jl) package:
```@example install
using PhyloPlots
using RCall # hide
mkpath("../assets/figures") # hide
R"name <- function(x) file.path('..', 'assets', 'figures', x)" # hide
R"svg(name('inputdata_gene3.svg'), width=4, height=3)" # hide
R"par"(mar=[0,0,0,0])                          # hide
plot(genetrees[3]); # tree for 3rd gene
R"dev.off()"                                   # hide
nothing # hide
```
![gene3](../assets/figures/inputdata_gene3.svg)


We can also read a network from a newick formatted string,
and, for example, print a list of its edges:

```@repl install
newickstring = "(A,((B,#H1),(C,(D)#H1)));";
net = readTopology(newickstring);
printEdges(net)
```

We see that the edges do not have branch lengths,
and the hybrid edges do not have gamma (inheritance) values.
We can set them with

```@repl install
setlength!(net.edge[1], 1.9)
setGamma!(net.edge[3],  0.8)
printEdges(net)
```
where 1 and 3 correspond to the position of the given edge to modify in the list of edges.
We can only change the γ value of hybrid edges,
not tree edges (for which γ=1 necessarily).
Such an attempt below will cause an error with a message to explain that
the edge was a tree edge:
```julia
setGamma!(net.edge[4], 0.7)
# should return this:
# ERROR: cannot change gamma in a tree edge
```
