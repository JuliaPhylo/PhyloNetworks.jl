#PhyloNetworks: analysis for phylogenetic networks in [Julia](http://julialang.org)
## Maximum pseudolikelihood estimation of species network: SNaQ <img src="http://pages.stat.wisc.edu/~claudia/Images/snaq.png" align=right title="SNaQ logo" width=262.5 height=111>
<!-- ![SNaQ logo](http://pages.stat.wisc.edu/~claudia/Images/snaq.png)
original size: 525px × 222px-->

[![Build Status](https://travis-ci.org/crsl4/PhyloNetworks.jl.svg)](https://travis-ci.org/crsl4/PhyloNetworks.jl)

<!--
[![Coverage Status](https://coveralls.io/repos/crsl4/PhyloNetworks/badge.svg?branch=master&service=github)](https://coveralls.io/github/crsl4/PhyloNetworks?branch=master)
-->

SNaQ implements the statistical inference method in Sol&iacute;s-Lemus and An&eacute;
[(2016, PLoS Genetics)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896).
The procedure involves a
numerical optimization of branch lengths and inheritance probabilities
and a heuristic search in the space of phylogenetic
networks.

Below is a quick tutorial on SNaQ and PhyloNetworks,
and [here](http://pages.stat.wisc.edu/~claudia/smallTutorial.pdf)
are slides with background on networks and explanations.
But check out the
[google group](https://groups.google.com/forum/#!forum/phylonetworks-users) for common questions.

### Input for SNaQ

Two alternatives:

1. A list of estimated gene trees for each locus, which can be obtained using
  [MrBayes](http://mrbayes.sourceforge.net) or [RAxML](http://sco.h-its.org/exelixis/software.html). Or:
2. A table of concordance factors (CF) for each 4-taxon subset which can be
  obtained from [BUCKy](http://www.stat.wisc.edu/~ane/bucky/),
  to account for gene tree uncertainty

### Pipeline from sequence alignments

This [pipeline](https://github.com/nstenz/TICR) (with additional details
[here](https://github.com/crsl4/PhyloNetworks/blob/master/docs/howto.md))
can be used to obtain the table of CF needed as input
for SNaQ. The pipeline starts with input the sequence alignments, and
it runs MrBayes and then BUCKy, producing the
table of estimated CFs and their credibility intervals.


### Installation of Julia

Julia is a high-level and interactive programming language (like R or Matlab),
but it is also high-performance (like C).
To install Julia, go to http://julialang.org/downloads/
For a basic tutorial on Julia, see http://learnxinyminutes.com/docs/julia/

IMPORTANT: Julia code is just-in-time compiled. This means that the
first time you run a function, it will be compiled at that moment. So,
please be patient! Future calls to the function will be much much
faster. Trying out toy examples for the first calls is a good idea.

Users should have Julia 0.4 or above.


### Installation of the package PhyloNetworks

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
There is a known bug for Mac users where the *Pkg.update* function does not update to the latest version. We recommend Mac users to do the following through the terminal:
```shell
cd HOME/.julia/v0.4/PhyloNetworks/
git pull
```
where HOME is replaced by your home directory.
-->

The PhyloNetworks package has dependencies like NLopt and Gadfly
(see the REQUIRE file for the full list), but everything is installed automatically.

<!--
- GraphViz (version 0.0.3)
- NLopt (version 0.2.0)

The version in parenthesis correspond to the ones used when
implementing PhyloNetworks.
-->

### Small examples
Everytime you start a session in Julia, you should type:
```julia
using PhyloNetworks
```
This step can also take a while, because Julia compiles all the code at this moment.
This might change in future versions of Julia.
Here is a very small test for the installation of PhyloNetworks.

```julia
net = readTopology("(A,(B,(C,D)));");
tipLabels(net)
```

You can see a list of all the functions with
```julia
whos(PhyloNetworks)
```
and press `?` inside Julia to switch to help mode,
followed by the name of a function (or type) to get more details about it.

#### Input data

<!--The examples files for this section can be found within the
PhyloNetworks folder, typically in your
*HOME/.julia/v0.4/PhyloNetworks/examples/*. However, links to the
files are also included below.
-->
We suggest that you create a special directory for running these examples,
where input files can be downloaded and where output files will be
created (with estimated networks for instance). Enter this directory
and run Julia from there.

Suppose you have a file with a list of gene trees in parenthetical
format called *treefile.txt*.
You can access the example file of input trees
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/treefile.txt)
or
[here](https://raw.githubusercontent.com/crsl4/PhyloNetworks/master/examples/treefile.txt)
for easier download.

Do not copy-paste into a "smart" text-editor. Instead, save the file
directly into your working directory using "save link as" or "download linked file as".
This file contains 10 trees, each in parenthetical format on 6 taxa
like this:

(6:2.728,((3:0.655,5:0.655):1.202,(1:0.881,(2:0.783,4:0.783):0.098):0.976):0.871);

If 'treefile.txt' is in your working directory, you can view its content
within Julia:
```julia
less("treefile.txt")
```
Just type `q` to quit viewing this file.
You could read in these 10 trees and visualize the third one (say) like this:
```julia
tentrees = readMultiTopology("treefile.txt")
tentrees[3]
plot(tentrees[3])
```
To read in all gene trees and directly summarize them by a list
of quartet CFs (proportion of input trees with a given quartet):
```julia
d=readTrees2CF("treefile.txt", CFfile="tableCFall.txt")
less("tableCFall.txt")
```
`less("tableCFall.txt")` lets you see the content of the newly created
file "tableCFall.txt", within Julia. Again, type `q` to quit viewing this file.

If instead of all the 4-taxon subsets, you just want to use a random
sample of 10 4-taxon subsets:
```julia
d=readTrees2CF("treefile.txt", whichQ="rand", numQ=10, CFfile="tableCF10.txt")
```
Be careful to use a numQ value smaller than the total number of possible
4-taxon subsets, which is *n choose 4* on *n* taxa (e.g. 15 on 6 taxa).
To get a predictable random sample, you may set the seed with `srand(12321)`
(for instance) prior to sampling the quartets as above.

If you already have a table of CF values in a file *tableCF.txt*
in this format

|Taxon1 | Taxon2 | Taxon3 | Taxon4 | CF12_34 | CF13_24 | CF14_23 |
|-------|:-------|:-------|:-------|:--------|:--------|:--------|

you would read it like this:
```julia
d=readTableCF("tableCF.txt") # one step only: read file and convert to 'DataCF' object
# or in 2 steps:
using DataFrames
dat = readtable("tableCF.txt") # read file into a 'DataFrame' object
d = readTableCF(dat)           # convert DataFrame into DataCF
```
You can access this example file
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/tableCF.txt)
or
[here](https://raw.githubusercontent.com/crsl4/PhyloNetworks/master/examples/tableCF.txt).

Columns need to be in the right order. If you have the information on the number of genes used for each 4-taxon subset, you can add this information in a column named "ngenes".

If you have a tree *startTree.txt* in parenthetical format to
use as starting point for the optimization, you can read it with
```julia
T=readTopology("startTree.txt")
```
You can access this example file
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/startTree.txt)
(raw file
[here](https://raw.githubusercontent.com/crsl4/PhyloNetworks/master/examples/startTree.txt)).

If the topology `T` is to be used as a starting tree for SNaQ, it needs
to be of level 1. To make sure that it is, you may read it in with
`readTopologyLevel1` (which may also unroot the tree, resolve polytomies,
replace missing branch lengths by 1 for starting values etc.):
```julia
T=readTopologyLevel1("startTree.txt")
```
Note that all trees and all networks with 1 hybridization are of level 1.

#### Network Estimation

To estimate the network using the input data
*d* and starting from tree (or network) *T*, do this:

```julia
net1=snaq!(T,d,filename="net1_snaq");
less("net1_snaq.err")
less("net1_snaq.out")
net2=snaq!(T,d,hmax=2, filename="net2_snaq");
```
when viewing the result files "net1_snaq.err" and "net1_snaq.out" with `less`
within Julia, use arrows to scroll down and type `q` to quit viewing the files.
<!---
The "net1_snaq.networks" file contains a list of networks obtained from moving
the placement of the hybrid node to another node inside the cycle,
along with its pseudolikelihood score.
-->
The option *hmax* corresponds to the maximum number of hybridizations allowed,
1 by default.
The function name `snaq!` ends with ! because it modifies the argument d
by including the expected CF. Type `?` then `snaq!` to get help on that function.

The estimation function creates a .out file (snaq.out by default) with the estimated
network in parenthetical format, which you can also print directly to the screen like this:
```julia
net1
writeTopology(net1)                   # topology to screen, full precision for branch lengths and γ
writeTopology(net1,di=true)           # γ omitted: for dendroscope
writeTopology(net1, "bestnet_h1.tre") # topology to file 'bestnet_h1.tre': creates or overwrites file
less("bestnet_h1.tre")                # just view the file
```
The option *di=true* is for the parenthetical format used by
[Dendroscope](http://dendroscope.org/) (without reticulation heritabilities).
Copy this parenthetical description and paste it into Dendroscope, or use the plotting function described below.

#### Network Visualization
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


### Bootstrap

You can run a bootstrap analysis if you estimated CFs with credibility intervals,
such as if you used the TICR
pipeline (see above). The TICR pipeline provides a CF table with extra columns for
credibility intervals.
```julia
using DataFrames
df = readtable("tableCFCI.csv")
bootnet = bootsnaq(T, df, hmax=1, nrep=10, runs=3)
```
You can access this example file
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/tableCFCI.csv)
or [here](https://raw.githubusercontent.com/crsl4/PhyloNetworks/master/examples/tableCFCI.csv).
Make sure that your downloaded file has the name `tableCFCI.csv` (to match the name in the
Julia command), with no an extra extension `.txt`. Rename the file after download if necessary.

This example uses a number of replicates (10) that is definitely too small, to
make the example run faster. You might also increase the number of optimization
runs (`runs`) done for each bootstrap replicate.
To save the bootstrap networks to a file (especially if it took a while to get them!)
and to check the content of the created file, do this:

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

#### Summarizing bootstrap on the main tree

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

#### Summarizing bootstrap on reticulations

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

<!---
To summarize the hybridizations, we need an outgroup to root all the networks.
```julia
outgroup = "4"
HFmat,discTrees = hybridDetection(net_bs,net1,outgroup)
```

The function `hybridDetection` will provide a matrix `HFmat` that will
have one row per bootstrap network, and number of columns depending on
the number of hybrids in `net1`. If `net1` had 2 hybrids, `HFmat` will
have 4 columns:

- the first 2 columns indicate the presence (1) or absence (0) of each
  hybrid (column) for each bootstrap network (row)
- the last 2 columns indicate the estimated gamma in the bootstrap
  network if the hybrid was found (and 0.0 if it was not found)

Hybrid comparison between explicit networks only makes sense if the
underlying trees are the same, so `discTrees` has the list of trees
that do not match the underlying tree in `net1` (`tree1`).

Finally, you can summarize the information in `HFmat`
```julia
df_hyb = summarizeHFdf(HFmat)
```
`df_hyb` has one row per hybrid, and 5 columns:

- hybrid index
- number of trees that match the underlying tree in `net1` (same for all hybrids)
- number of networks with that hybrid
- mean estimated gamma among networks with the hybrid
- sd estimated gamma among networks with the hybrid

The last row contains in 3er column the number of networks that have
all same hybrids as `net1` (hybrid index, mean gamma and sd gamma are
meaningless in this row).

You can save all the information with
```julia
HFdf=convert(DataFrame,HFmat) #convert to dataframe to save
writetable("HFdf.csv",HFdf)
writetable("summaryHybridDetection.csv",df_hyb)
s=open("discrepantTrees.out","w")
for(t in discTrees)
    write(s,"$(writeTopology(t))\n")
end
close(s)
```
--->

### Simple use of Julia objects

For a small example on how Julia objects can be accessed, see
[here](https://github.com/crsl4/PhyloNetworks/blob/master/docs/simpleJulia.md)

### Multiple alleles

The usual settings for SNaQ consider each allele in the gene trees to
be its own tip in the network. If instead each allele can be mapped confidently
to a species, and if only the species-level network needs
to be estimated, this can be done with the following functions:
```julia
new_df = mapAllelesCFtable(mappingFile, CFtable);
new_d = readTableCF(new_df);
```
where the mapping file can be a text file (or csv) with two columns
named *allele* and *species*, mapping each allele name to a species
name. The CF table is the original table with allele names for each
4-taxon subset. This function will create a new CF data frame with the
species names instead of allele names, and will modify *new_df* by
removing rows like *sp1,sp1,sp1,sp1*, which contain no information about
between-species relationships.

Estimation will work the same way:
```julia
new_net = snaq!(new_T,new_d);
```
where *new_T* should be a starting topology with one tip per species, labelled with the species names.
<!--
WARNING: the current function works best if all alleles from the same
individual are given the same name (the individual's 'name') across
all genes for which that individual was sequenced.
-->

### Optimizing branch lengths and inheritance probabilities for a given network

For a given network topology, we can optimize the branch lengths and
inheritance probabilities (γ) with the pseudolikelihood.
This is useful if we have a few candidate networks to compare.
Each network can be optimized individually, and the network with the best
pseudolikelihood can be chosen.

The score being optimized is the pseudo-deviance, i.e.
the negative log pseudo-likelihood up to an additive constant,
such that a perfect fit corresponds to a deviance of 0.0 (the lower the better).
```julia
net1topo = readTopology("(2,(4,(3,(5,(6,#H1)))),(1)#H1);");
net1par = topologyMaxQPseudolik!(net1topo,d)
net1par.loglik # pseudo deviance, actually
```
For a more thorough optimization, we may increase the requirements before
the search stops:
```julia
net1par = topologyMaxQPseudolik!(net1topo,d, xtolRel=1e-10, xtolAbs=1e-10)
net1par.loglik
```

For a network with given branch lengths and γ heritabilies,
we can compute the pseudolikelihood with:
```julia
net1withBL = readTopology("(2,(4,(3,(5,(6,#H6:1.0::0.3):5.006):0.518):0.491):1.533,(1)#H6:1.0::0.7);");
topologyQPseudolik!(net1withBL,d)
net1withBL.loglik
```
This function is not maximizing the pseudolikelihood, it is simply computing the
pseudolikelihood (or deviance) for the given branch lengths and probabilities of
inheritance. At the moment, both of these functions require that the
given network is of level 1 (cycles don't overlap).

### Plot observed CF vs expected CF
A good way to visualize the "goodness-of-fit" of a given estimated network to the data is to plot the observed CF to the expected CF. If the network is a good fit, then the dots in the plot will be close to the y=x line.
The following function will create a dataframe with the observed and expected CF which are all saved in the DataCF object after running snaq:
```julia
df = dfObsExpCF(d)
```
It is important to have run snaq, topologyQPseudoLik or topologyMaxQPseudolik before making this plot or the expected CF would be meaningless.

Now, we can plot them with any of the Julia packages for plotting. In particular:
```julia
using Gadfly
p = plot(df,layer(x="obsCF12",y="expCF12",Geom.point,Theme(default_color=colorant"orange")),
            layer(x="obsCF13",y="expCF13",Geom.point,Theme(default_color=colorant"purple")),
            layer(x="obsCF14",y="expCF14",Geom.point,Theme(default_color=colorant"blue")),
            layer(x=0:1,y=0:1),Geom.line,Theme(default_color=colorant"black"))
```
This will pop out a browser window with the plot. The plot can be saved as PDF (or many other formats, see [Gadfly tutorial](http://dcjones.github.io/Gadfly.jl/)) with
```julia
draw(PDF("plot.pdf", 4inch, 3inch), p)
```

### Debugging: the .err file
Please report any bugs and errors to *claudia@stat.wisc.edu*. The easiest way to do it is by checking the .err file which will show the number of runs that failed by a bug and the corresponding seed to replicate the run. This is an example of what the .err file looks like:

Total errors: 1 in seeds [4545]

You need to run the following function with the same settings that caused the error:

```julia
snaqDebug(T,d,hmax=2,seed=4545)
```

This will create two files:
*snaqDebug.log* and *debug.log* which you can then send to
*claudia@stat.wisc.edu* with subject "SNaQ bug found" or something
similar. I will not have access to any part of your data, the files
simply print out the steps to retrace the bug, and hopefully fix it.

### Getting help

The easiest way to get help is to post (or email) a question to the
PhyloNetworks-users google group [here]
(https://groups.google.com/forum/#!forum/phylonetworks-users/new).  It
is a good idea to join the group to receive information on new
versions, bugs fixed, etc.
