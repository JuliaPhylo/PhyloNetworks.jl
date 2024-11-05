```@setup snaqplot
using PhyloNetworks
mkpath("../assets/figures")
raxmltrees = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","raxmltrees.tre")
raxmlCF = readTrees2CF(raxmltrees, writeTab=false, writeSummary=false)
astralfile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","astral.tre")
astraltree = readMultiTopology(astralfile)[102] # 102th tree = last tree here
net0 = readTopology(joinpath(dirname(pathof(PhyloNetworks)), "..","examples","net0.out"))
net1 = readTopology(joinpath(dirname(pathof(PhyloNetworks)), "..","examples","net1.out"))
rotate!(net1, -6)
net2 = readTopology(joinpath(dirname(pathof(PhyloNetworks)), "..","examples","net2.out"))
net3 = readTopology(joinpath(dirname(pathof(PhyloNetworks)), "..","examples","net3.out"))
net0.loglik = 53.53150526187732
net1.loglik = 28.31506721890958
net2.loglik = 28.31506721890957
net3.loglik = 28.315067218909626
```
# Network Visualization


```julia
net0 = snaq!(astraltree,raxmlCF, hmax=0, filename="net0", seed=1234)
```

```@example snaqplot
using PhyloPlots
using RCall # hide
R"name <- function(x) file.path('..', 'assets', 'figures', x)" # hide
R"svg(name('snaqplot_net0_1.svg'), width=4, height=3)" # hide
R"par"(mar=[0,0,0,0]) # hide
plot(net0);
R"dev.off()"; # hide
nothing # hide
```
![net0_1](../assets/figures/snaqplot_net0_1.svg)

We can visualize the estimated network and its inheritance values γ, which
measure the proportion of genes inherited via each parent at a reticulation event
(e.g. proportion of genes inherited via gene flow).
```@example snaqplot
R"svg(name('snaqplot_net1_1.svg'), width=4, height=3)"; # hide
R"par"(mar=[0,0,0,0]); # hide
plot(net1, showgamma=true);
R"dev.off()"; # hide
nothing # hide
```
![net1_1](../assets/figures/snaqplot_net1_1.svg)

This network has A as a hybrid, 80.4% sister to B,
and 19.6% sister to E (which is otherwise sister to O).
C & D are sister to each other.

The main output file, here `net1.out` (or `snaq.out` by default) has the estimated
network in parenthetical format, but we can also print it directly to the screen:
```@repl snaqplot
net1
writeTopology(net1)  # writes to screen, full precision for branch lengths and γ
writeTopology(net1, round=true, digits=2)
writeTopology(net1, di=true) # γ omitted: for dendroscope
writeTopology(net1, "bestnet_h1.tre") # writes to file: creates or overwrites file
rm("bestnet_h1.tre") # hide
```
The option `di=true` is for the parenthetical format used by
[Dendroscope](http://dendroscope.org/) (without reticulation heritabilities).
Copy this parenthetical description and paste it into Dendroscope,
or use the plotting function described below.

We can go on and let the network have up to 2 or 3 hybrid nodes:
```julia
net2 = snaq!(net1,raxmlCF, hmax=2, filename="net2", seed=3456)
net3 = snaq!(net0,raxmlCF, hmax=3, filename="net3", seed=4567)
```
and plot them (they are identical and they both have a single reticulation):
```@example snaqplot
R"svg(name('snaqplot_net23.svg'), width=7, height=3)" # hide
using RCall                  # to be able to tweak our plot within R
R"layout(matrix(1:2, 1, 2))" # to get 2 plots into a single figure: 1 row, 2 columns
R"par"(mar=[0,0,1,0])        # for smaller margins
plot(net2, showgamma=true);
R"mtext"("hmax=2")           # add text annotation: title here
plot(net3, showgamma=true);
R"mtext"("hmax=3")
R"dev.off()"; # hide
nothing # hide
```
![net23](../assets/figures/snaqplot_net23.svg)

## Network Visualization

To visualize a network, we can use the companion package
[PhyloPlots](https://github.com/juliaphylo/PhyloPlots.jl).
In the example below, julia creates and sends the plot to R
via [RCall](https://github.com/JuliaInterop/RCall.jl),
so we can tweak the plot in various ways via commands sent to R.
To save the plot in a file: we first tell R to create an image file,
then we send the plot of the network,
then we tell R to wrap up and save its image file.

```@example snaqplot
using PhyloPlots # to visualize networks
using RCall      # to send additional commands to R like this: R"..."
imagefilename = "../assets/figures/snaqplot_net1_2.svg"
R"svg"(imagefilename, width=4, height=3) # starts image file
R"par"(mar=[0,0,0,0]) # to reduce margins (no margins at all here)
plot(net1, showgamma=true, showedgenumber=true); # network is plotted & sent to file
R"dev.off()"; # wrap up and save image file
nothing # hide
```
![net1_2](../assets/figures/snaqplot_net1_2.svg)

The plot function has many options, to annotate nodes and edges. In the
example above, hybrid edges were annotated with their γ inheritance values
(in blue: light blue for the minor edge with γ<0.5, and dark blue for the
major edge with γ>0.5), and edges were annotated with their internal numbers.

Type `?` to switch to the help mode
of Julia, then type the name of the function, here `plot`.
Below are two visualizations.
The first uses the default style (`:fulltree`) and modified edge colors.
The second uses the `:majortree` style.
That style doesn't have an arrow by default for minor hybrid edges,
but we can ask for one by specifying a positive arrow length.
```@example snaqplot
R"svg(name('snaqplot_net1_3.svg'), width=7, height=3)" # hide
R"par"(mar=[0,0,0,0]) # hide
R"layout(matrix(1:2,1,2))";
plot(net1, showedgelength=true, minorhybridedgecolor="tan");
plot(net1, style=:majortree, arrowlen=0.07);
R"dev.off()"; # hide
nothing # hide
```
![net1_3](../assets/figures/snaqplot_net1_3.svg)

Edge lengths are shown, too. They were estimated in coalescent units:
number of generations / effective population size.
Some edge lengths are not identifiable, hence not shown.

Below is another example, where space was added between the network and
the taxon names via the `tipoffset` option.
Also, edge colors were changed, and the nodes numbers are shown (used internally)

```@example snaqplot
R"svg(name('snaqplot_net1_4.svg'), width=4, height=3)" # hide
R"par"(mar=[0,0,0,0]) # hide
plot(net1, tipoffset=0.5, shownodenumber=true, edgecolor="tomato4",
     minorhybridedgecolor="skyblue", majorhybridedgecolor="tan");
R"dev.off()"; # hide
nothing # hide
```
![net1_4](../assets/figures/snaqplot_net1_4.svg)

## Re-rooting networks

SNaQ and other methods infer a semidirected network: in which the root position
is unknown and the direction of tree edges is unknown, but the direction
of hybrid edges is known (and so the identify of hybrid nodes is known).

To obtain a representative visualization,
it is best to root the network first, using one or more outgroup.
Go to [Re-rooting trees and networks](@ref) for this.
If the outgroup conflicts with the direction of reticulations
in the estimated network, see section
[Candidate networks compatible with a known outgroup](@ref).
