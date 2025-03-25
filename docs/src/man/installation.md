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
[here](https://groups.google.com/forum/#!forum/juliaphylo-users/new).

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

To check that your installation worked, type this in Julia to load the package.
This is something to type every time you start a Julia session:
```@example install
using PhyloNetworks
```
This step can also take a while, to pre-compile the code (after a package
update for instance).
