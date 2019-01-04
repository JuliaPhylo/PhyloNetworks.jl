# Installation

## Installation of Julia

Julia is a high-level and interactive programming language (like R or Matlab),
but it is also high-performance (like C).
To install Julia, follow instructions [here](http://julialang.org/downloads/).
For a quick & basic tutorial on Julia, see
[learn x in y minutes](http://learnxinyminutes.com/docs/julia/).

Editors:

- [Visual Studio Code](https://code.visualstudio.com) provides an editor
  and an integrated development environment (IDE) for Julia: highly recommended!
- [Juno](http://junolab.org) provides an IDE for Julia,
  based on the [Atom](http://atom.io/) editor.
- you can also run Julia within a [Jupyter](http://jupyter.org) notebook
  (formerly IPython notebook).

IMPORTANT: Julia code is just-in-time compiled. This means that the
first time you run a function, it will be compiled at that moment. So,
please be patient! Future calls to the function will be much much
faster. Trying out toy examples for the first calls is a good idea.

## Installation of the package PhyloNetworks

To install the package, type inside Julia:
```julia
using Pkg
Pkg.add("PhyloNetworks")
```
The first step can take a few minutes, be patient. If you already installed the package and want
the latest registered version, just do this (which will update all of your packages):
```julia
Pkg.update()
```
Warning: It is important to update the package regularly as it is
undergoing constant development. Join the google group for updates
[here]
(https://groups.google.com/forum/#!forum/phylonetworks-users/new).

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
[NLopt](https://github.com/JuliaOpt/NLopt.jl) and
[DataFrames](http://juliadata.github.io/DataFrames.jl/stable/)
(see the REQUIRE file for the full list), but everything is installed automatically.

The companion package [PhyloPlots](https://github.com/cecileane/PhyloPlots.jl)
has utilities to visualize networks, and for interoperability,
such as to export networks to R (which can then be plotted via R).
To install:

```julia
using Pkg
Pkg.add("PhyloPlots")
```

PhyloPlots depends on PhyloNetworks, and has further dependencies
like [Gadfly](http://gadflyjl.org/stable/) and
[RCall](https://github.com/JuliaInterop/RCall.jl)

## Test example

To check that your installation worked, type this in Julia to load the package.
This is something to type every time you start a Julia session:
```@example install
using PhyloNetworks;
```
This step can also take a while, if Julia needs to pre-compile the code (after a package
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
info on an object, what's its type, and how to manipulate objects.
For example, let's take an object `raxmlCF` created from reading in some data
(see [Input for SNaQ](@ref)):

```@repl install
raxmltrees = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","raxmltrees.tre");
raxmlCF = readTrees2CF(raxmltrees);
```

Typing `varinfo()` will provide a list of objects and packages in memory,
including `raxmlCF` that we just created.
If we want to know the type of a particular object, we do:
```@repl install
typeof(raxmlCF)
```
which shows us that `raxmlCF` is of type `DataCF`.
If we want to know about the attributes the object has, we can type `?` in Julia,
followed by `DataCF` for a description.
We can also ask for a list of all its attributes with

```@repl install
fieldnames(typeof(raxmlCF))
```
For example, we see that one attribute is `numQuartets`: its the number of 4-taxon subsets
in the data. To see what this number is:
```@repl install
raxmlCF.numQuartets
```
We also noticed an attribute `quartet`. It is a vector of Quartet objects inside `raxmlCF`, so
```@repl install
raxmlCF.quartet[2].taxon
```
will provide the list of taxon names for the second 4-taxon subset in the data.
To see the observed CF, we can type
```@repl install
raxmlCF.quartet[2].obsCF
```
We can verify the type with
```@repl install
typeof(raxmlCF.quartet[2])
```
We can also read a simple network in Julia and print the list of edges
```@repl install
str = "(A,((B,#H1),(C,(D)#H1)));";
net = readTopology(str);
printEdges(net)
```
We see that the edges do not have branch lengths,
and the hybrid edges do not have gamma values. We can set them with
```@repl install
setLength!(net.edge[1],1.9)
setGamma!(net.edge[3],0.8)
printEdges(net)
```
where 1 and 3 correspond to the position of the given edge to modify in the list of edges.
We can only change the gamma value of hybrid edges (not tree edges).
Such an attempt below will cause an error with a message to explain that
the edge was a tree edge:
```julia
setGamma!(net.edge[4],0.7)
# should return this:
# ERROR: cannot change gamma in a tree edge
```
