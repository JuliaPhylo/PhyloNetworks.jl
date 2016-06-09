# Installation

## Installation of Julia

Julia is a high-level and interactive programming language (like R or Matlab),
but it is also high-performance (like C).
To install Julia, go to http://julialang.org/downloads/
For a basic tutorial on Julia, see http://learnxinyminutes.com/docs/julia/

IMPORTANT: Julia code is just-in-time compiled. This means that the
first time you run a function, it will be compiled at that moment. So,
please be patient! Future calls to the function will be much much
faster. Trying out toy examples for the first calls is a good idea.

Users should have Julia 0.4 or above.


## Installation of the package PhyloNetworks

To install the package, type inside Julia:
```julia
Pkg.add("PhyloNetworks")
```
The first step can take a few minutes, be patient. If you already installed the package and want
the latest registered version, just do this (which will update all of your packages):
```julia
Pkg.update()
```
WARNING: It is important to update the package regularly as it is
undergoing constant development. Join the google group for updates
[here]
(https://groups.google.com/forum/#!forum/phylonetworks-users/new).

`Pkg.update()` will install the latest registered version, but there
could be other improvements in the `master` branch of the
repository. If you want to update to the latest unregistered version
of the package, you can do `Pkg.checkout("PhyloNetworks")` just beware
that the latest changes could be not as robust. If you want to go back to the registered package, you can do `Pkg.free("PhyloNetworks")`.

Similarly, you can pin a version of the package
`Pkg.pin("PhyloNetworks")` so that `Pkg.update()` will not modify
it. You can always free a pinned package with
`Pkg.free("PhyloNetworks")`. More on checked and pinned packages [here]
(http://docs.julialang.org/en/release-0.4/manual/packages/).

<!--
If, for some reason, the *Pkg.update* function does not update to the latest version,
users can do the following through the terminal:
    cd HOME/.julia/v0.4/PhyloNetworks/
    git pull
where HOME is replaced by your home directory.
-->

The PhyloNetworks package has dependencies like NLopt and Gadfly
(see the REQUIRE file for the full list), but everything is installed automatically.


## Test example

To check that your installation worked, type this in Julia to load the package.
This is something to type every time you start a Julia session:
```@example 1
using PhyloNetworks;
```
This step can also take a while, if Julia needs to pre-compile the code (after a package
update for instance).
Here is a very small test for the installation of PhyloNetworks.

```@example 1
net = readTopology("(A,(B,(C,D)));");
tipLabels(net)
```

You can see a list of all the functions with
```julia
whos(PhyloNetworks)
```
and press `?` inside Julia to switch to help mode,
followed by the name of a function (or type) to get more details about it.
 
Here are small examples on how [Julia types](@ref) can be accessed.
