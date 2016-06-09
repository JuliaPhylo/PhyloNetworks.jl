# Getting a Network

## Network Estimation

After [Input for SNaQ](@ref), we can estimate the network using the
input data `d` and starting from tree (or network) `T`:

```julia
net1=snaq!(T,d,filename="net1_snaq");
less("net1_snaq.err")
less("net1_snaq.out")
less("net1_snaq.networks")
net2=snaq!(T,d,hmax=2, filename="net2_snaq");
```
when viewing the result files "net1_snaq.err" and "net1_snaq.out" with `less`
within Julia, use arrows to scroll down and type `q` to quit viewing the files.
The "net1_snaq.networks" file contains a list of networks obtained from moving
the placement of the hybrid node to another node inside the cycle,
along with its pseudolikelihood score.

The option `hmax` corresponds to the maximum number of hybridizations allowed,
1 by default.
The function name `snaq!` ends with ! because it modifies the argument `d`
by including the expected CF. Type `?` then `snaq!` to get help on that function.

The estimation function creates a `.out` file (`snaq.out` by default) with the estimated
network in parenthetical format, which you can also print directly to the screen like this:
```julia
net1
writeTopology(net1)                   # topology to screen, full precision for branch lengths and γ
writeTopology(net1,di=true)           # γ omitted: for dendroscope
writeTopology(net1, "bestnet_h1.tre") # topology to file 'bestnet_h1.tre': creates or overwrites file
less("bestnet_h1.tre")                # just view the file
```
The option `di=true` is for the parenthetical format used by
[Dendroscope](http://dendroscope.org/) (without reticulation heritabilities).
Copy this parenthetical description and paste it into Dendroscope,
or use the plotting function described below.

### SNaQ error

Please report any bugs and errors to *claudia@stat.wisc.edu*, so we can debug it.
The easiest way to do it is by checking the `.err` file which will show the number of runs that
failed by a bug and the corresponding seed to replicate the run.
This is an example of what the `.err` file looks like:
`Total errors: 1 in seeds [4545]`.
You need to run the following function with the same settings that caused the error:

```julia
snaqDebug(T,d,hmax=2,seed=4545)
```

This will create two files:
*snaqDebug.log* and *debug.log* which you can then send to
*claudia@stat.wisc.edu* with subject "SNaQ bug found" or something
similar. I will not have access to any part of your data, the files
simply print out the steps to retrace the bug, and hopefully fix it.

## Network Visualization

To visualize the network:
```julia
p = plot(net1, showGamma=true)
```
This function will open a browser where the plot will appear. To get a pdf version of the plot:
```julia
using Gadfly
draw(PDF("bestnet_h1.pdf", 4inch, 4inch),p)
```
The plot function has many options. Type `?` to switch to the help mode
of Julia, then type the name of the function, here `plot`.
Edge colors can be modified, for instance.
```julia
# using Gadfly # if not done earlier
plot(net1, showEdgeLength=true, minorHybridEdgeColor=colorant"tan")
```

## Re-rooting networks

SNaQ infers an unrooted semi-directed network, in the sense
that the direction of tree edges cannot be inferred, but the direction
of hybrid edges can be inferred. To obtain a representative visualization,
it is best to root the network first, using one or more outgroup.
If there is a single outgroup, the network can be rooted with this outgroup,
if compatible, like this:
```julia
rootatnode!(net1, "4")
plot(net1, showGamma=true)
```
Here we used taxon "4" as outgroup. More options are available, to root the
network either at a give node or along a given edge. Use the help mode (type `?`)
to get help on the functions `rootatnode!` and `rootonedge!` to get more info.

If the network is plotted with crossing edges, you may identify
ways to rotate the children edges at some nodes to untangle some crossing edges.
This can be done using the function `rotate!`. Type `?` then `rotate!` to get
help and examples.

## Candidate Network Evaluation

From a set of candidate networks, one might simply need to score of each network
to pick the best. Here, the score is the negative log pseudo-likelihood, and the
lower the better. See the section to get the score of a [Fixed Network](@ref).

Else: go next to [Extract Expected CFs](@ref) to see how your network fits your data,
or go for a [Bootstrap](@ref) analysis to quantify support for tree edges and hybrid edges.
