# Bootstrap

## Running a bootstrap analysis

We can run a bootstrap analysis if we estimated quartet concordance factors (CFs)
with credibility intervals, such as if we used the [TICR pipeline](@ref).
This pipeline provides a CF table with extra columns for
credibility intervals.
```julia
using DataFrames
df = readtable("tableCFCI.csv")
bootnet = bootsnaq(T, df, hmax=1, nrep=10, runs=3, filename="bootstrap")
```
This example file can be accessed
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/tableCFCI.csv)
or [here](https://raw.githubusercontent.com/crsl4/PhyloNetworks/master/examples/tableCFCI.csv).
Make sure that your downloaded file has the name `tableCFCI.csv` (to match the name in the
Julia command), with no an extra extension `.txt`. Rename the file after download if necessary.

This example uses a number of replicates (10) that is definitely too small, to
make the example run faster. You might also increase the number of optimization
runs (`runs`) done for each bootstrap replicate.

The bootstrap networks are saved in the `boostrap.out` file, so they
can be read in a new session with `bootnet =
readMultiTopology("boostrap.out")`. To save the bootstrap networks to
a different file (perhaps after having re-rooted them with an
outgroup), you can do this:

```julia
writeMultiTopology(bootnet, "bootstrapNets_h1.tre")
length(bootnet) # 10 networks in the array 'bootnet'
less("bootstrapNets_h1.tre")
```
If you close your session and re-open it later, you can re-read the bootstrap networks
back in memory like this:
```julia
bootnet = readMultiTopology("bootstrapNets_h1.tre");
```

We can run bootstrap analysis from gene trees also. We need a text
file with a list of files containg the bootstrap trees (one file per
gene). This is the same format used by ASTRAL (see the wiki page for more details.).
```julia
bootTrees = readBootstrapTrees("BSlistfiles")
bootnet = bootsnaq(T, bootTrees, hmax=1, nrep=100, filename="bootsnaq1_raxmlboot")
```
If you used a specified list of quartets on the original data, you
should use that same list for the bootstrap analysis through the
option `quartetfile`.

## Summarizing bootstrap on the main tree

The `bootsnaq` function, as used before, returned a list of 10 networks (`nrep=10`).
Before summarizing them on the best network, it is best to re-read this network
to get a reproducible internal numbering of its nodes and edges, used later for mapping
bootstrap support to edges. We can then summarize bootstrap networks:
```julia
net1 = readTopology("net1_snaq.out") # reading from output file of snaq
BStable, tree1 = treeEdgesBootstrap(bootnet,net1)
```
This will calculate the major tree `tree1` displayed in `net1`, that is,
the tree obtained by following the major parent (γ>0.5) of each hybrid node.
This tree can be visualized like this, with edge numbers shown for later use.
```julia
rootatnode!(net1, "4")
rootatnode!(tree1, "4")
plot(tree1, showEdgeNumber=true)
```
Next, we can look at bootstrap table `BStable`, which has one row for
each tree edge in `net1`. One column contains the edge number
(same as shown in the plot) and another column contains the edge
bootstrap support: the proportion of bootstrap replicates in which this edge was
found in the major tree of the inferred network.
We can see the full bootstrap table and see
which tree edges have bootstrap support lower than 100% with
```julia
showall(BStable)
BStable[BStable[:proportion] .< 1.0, :]
```
Finally, we can map the bootstrap proportions onto the network or its main tree
by passing the bootstrap table to the `edgeLabel` option of `plot`:
```julia
plot(tree1, edgeLabel=BStable)
plot(net1,  edgeLabel=BStable)
```
(Here, it is important that the numbers assigned to edges when building the boostrap
table --those in `net1` at the time-- correspond to the current edge numbers
in `tree1` and `net1`. That was the purpose of reading the network from the
output file of snaq! earlier, for consistency across different Julia sessions.)

If we wanted to plot only certain bootstrap values, like those below 100% (1.0),
we could do this:
```julia
plot(net1, edgeLabel=BStable[BStable[:proportion] .< 1.0, :])
```

## Summarizing bootstrap on reticulations

Summarizing the placement of hybridization edges is not standard.
The function `hybridBootstrapSupport` attempts to do so.
The descendants of a given hybrid node form the "recipient" or "hybrid" clade,
and is obtained after removing all other reticulations.
If reticulation is due to gene flow or introgression, the minor hybrid edge (with γ<0.5)
represents this event. The descendants of the lineage from which gene flow originated
is then a second "sister" of the hybrid clade. Because of the reticulation event,
the hybrid clade has 2 sister clades, not 1: the major sister (through the major hybrid edge
with γ>0.5) and the minor sister (through the minor hybrid edge with γ<0.5).
Note that the network says *nothing* about the process: its shows the *relationships* only.
We can calculate the frequency that each clade is a hybrid clade, or a major or minor sister
for some other hybrid, in the bootstrap networks:
```julia
BSn, BSe, BSc, BSgam, BSedgenum = hybridBootstrapSupport(bootnet, net1);
```
Let's look at the results.
We can list all the clades and the percentage of bootstrap networks (bootstrap support)
in which each clade is a hybrid or sister to a hybrid:
```julia
BSn
```
If a clade contains a single taxon, it is listed with its taxon name.
The clade found in the best network is listed with its tag, starting with H (e.g. "H7").
The name of other clades start with "c_" followed by their number in the best network, if they
do appear in the best network.
The node numbers, as used internally in the best network, are listed in a separate column.
They can be used later to display the bootstrap support values onto the network.
Various columns give the bootstrap support that each clade is a hybrid, or a (major/minor) sister
to a hybrid. The last column gives the bootstrap support for the full relationship in the
best network: same hybrid with same two sisters.
These bootstrap values are associated with nodes (or possibly, their parent edges).

To see what is the clade named "H7", for instance:
```julia
BSc # this might be too big
showall(BSc)
BSc[:taxa][BSc[:H7]]
```
We can also get bootstrap values associated with edges, to describe the support that a given
hybrid clade has a given sister clade.
```julia
BSe
```
Here, each row describes a pair of 2 clades: one being the hybrid, the other being its sister,
connected by a hybrid edge. The first rows corresponds to hybrid edges in the best network. Other
rows correspond to edges seen in bootstrap networks but not in the best network.
```julia
BSedgenum
```
lists all the hybrid edges in the best network, two for each hybrid node:
the major parent edge and then the minor parent edge.
In our case, there is only one reticulation, so only 2 hybrid edges.

We can plot the bootstrap values of the 2 hybrid edges in the best network:
```julia
plot(net1, edgeLabel=BSe[[:edge,:BS_hybrid_edge]])
```
This is showing the bootstrap support each hybrid edge: percentage of bootstrap trees with an
edge from the same sister clade to the same hybrid clade.
Alternatively, we could show the bootstrap support for the full reticulation relationships in
the network, one at each hybrid node (support for same hybrid with same sister clades):
```julia
plot(net1, nodeLabel=BSn[[:hybridnode,:BS_hybrid_samesisters]])
```
On a different plot, we can show the bootstrap support for hybrid clades, shown on the parent
edge of each node with positive hybrid support:
```julia
plot(net1, edgeLabel=BSn[BSn[:BS_hybrid].>0, [:edge,:BS_hybrid]])
```
or plot the support for minor sister clades (shown along parent edge, and filtered to those with
sister support > 5%):
```julia
plot(net1, edgeLabel=BSn[BSn[:BS_minor_sister].>5, [:edge,:BS_minor_sister]])
```
To plot the same values next to the associated nodes, rather than edges, we use the option
`nodeLabel` and we select the column `node` or `hybridnode`:
```julia
plot(net1, nodeLabel=BSn[BSn[:BS_hybrid].>0, [:hybridnode,:BS_hybrid]])
plot(net1, nodeLabel=BSn[BSn[:BS_minor_sister].>5, [:node,:BS_minor_sister]])
plot(net1, nodeLabel=BSn[BSn[:BS_major_sister].>5, [:node,:BS_major_sister]])
```

The estimated heritability γ on hybrid edges in the reference network, when present in a
bootstrap network, was also extracted:
```julia
BSgam
```
γ=0 values are for bootstrap replicates that did not have the edge in their network.
Basic summaries on γ values for a given edge, say the minor parent,
could be obtained like this:
```julia
minimum(BSgam[:,2])
maximum(BSgam[:,2])
mean(BSgam[:,2])
std(BSgam[:,2])
```
