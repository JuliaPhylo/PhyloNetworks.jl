# Extract Expected CFs

A good way to visualize the "goodness-of-fit" of a given estimated network to the data
is to plot the observed CF versus the expected CF. If the network is a good fit, then the dots
in the plot will be close to the diagonal (x=y line).
The following function will create a dataframe with the observed and expected CFs,
which are all saved in the DataCF object after running snaq:
```julia
df_wide = fittedQuartetCF(d) # same as fittedQuartetCF(d, :wide)
df_long = fittedQuartetCF(d, :long)
```
It is important to have run `snaq!`, `topologyQPseudolik!` or `topologyMaxQPseudolik!`
before making this plot, or the result would be meaningless.
These functions update the fitted concordance factors (those expected under the network)
inside the dataCF object `d`.

Now, we can plot them with any of the Julia packages for plotting. For example:
```julia
using Gadfly
p = plot(layer(df_long, x="obsCF", y="expCF", Geom.point, Stat.x_jitter(range=0.4)),
         layer(x=0:1,y=0:1, Geom.line))
```
This will pop up a browser window with the plot.
The plot can be saved as PDF (or many other formats, see
[Gadfly tutorial](http://dcjones.github.io/Gadfly.jl/)) with
```julia
draw(PDF("plot.pdf", 4inch, 3inch), p)
```
To highlight quartets that include taxon "6", say,
if we suspect that it is an unrecognized hybrid, one may do this.
```
using DataFramesMeta # install this package with Pkg.add("DataFramesMeta") if not done before
mycolor = @with(df_long, (:tx1 .== "6") | (:tx2 .== "6") | (:tx3 .== "6") | (:tx4 .== "6"));
# 'mycolor' is true for quartets having taxon "6", false for others
p = plot(layer(df_long, x="obsCF", y="expCF", color=mycolor, Geom.point),
         layer(x=0:1,y=0:1, Geom.line),
         Guide.colorkey("has taxon 6?"))
```

To export this table and explore the fit of the network with other tools:
```julia
writetable("fittedCF_net1_long.csv", df_long)
```

