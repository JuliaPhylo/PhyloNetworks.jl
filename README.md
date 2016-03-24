#PhyloNetworks: analysis for phylogenetic networks in [Julia](http://julialang.org)
## Maximum pseudolikelihood estimation of species network: SNaQ <img src="http://pages.stat.wisc.edu/~claudia/Images/snaq.png" align=right title="SNaQ logo" width=262.5 height=111>
<!-- ![SNaQ logo](http://pages.stat.wisc.edu/~claudia/Images/snaq.png)
original size: 525px × 222px-->

[![Build Status](https://travis-ci.org/crsl4/PhyloNetworks.svg)](https://travis-ci.org/crsl4/PhyloNetworks)

<!--
[![Coverage Status](https://coveralls.io/repos/crsl4/PhyloNetworks/badge.svg?branch=master&service=github)](https://coveralls.io/github/crsl4/PhyloNetworks?branch=master)
-->

SNaQ implements the statistical inference method in Sol&iacute;s-Lemus and An&eacute;
[(2016, PLoS Genetics)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896).
The procedure involves a
numerical optimization of branch lengths and inheritance probabilities
and a heuristic search in the space of phylogenetic
networks.

Below there is a quick tutorial on SNaQ and PhyloNetworks, but refer
to the [PDF documentation](https://github.com/crsl4/PhyloNetworks/blob/master/docs/PhyloNetworks.pdf) for more details, and check out the [google group](https://groups.google.com/forum/#!forum/phylonetworks-users) for common questions.

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
Pkg.clone("https://github.com/crsl4/PhyloNetworks.git")
Pkg.build("PhyloNetworks")
```
The first step can take a few minutes, be patient. If you already installed the package and want
the latest version, just do this (which will update all of your packages):
```julia
Pkg.update()
```
WARNING: It is important to update the package regularly as it is
undergoing constant development. Join the google group for updates
[here]
(https://groups.google.com/forum/#!forum/phylonetworks-users/new).

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
and press ? inside Julia, followed by the name of a functions to get more details about it.

#### Input data

The examples files for this section can be found within the
PhyloNetworks folder, typically in your
*HOME/.julia/v0.4/PhyloNetworks/examples/*. However, links to the
files are also included below.

Suppose you have a file with a list of gene trees in parenthetical
format called *treefile.txt*.
If 'treefile.txt' is in your directory, do this to read in all gene trees
and to summarize them with a list
of quartet CFs:
```julia
d=readTrees2CF("treefile.txt")
```
You can access this example file
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/treefile.txt)
or
[here](https://raw.githubusercontent.com/crsl4/PhyloNetworks/master/examples/treefile.txt)
for easier download.
This file contains 10 trees, each in parenthetical format on 6 taxa
like this:

(6:2.728,((3:0.655,5:0.655):1.202,(1:0.881,(2:0.783,4:0.783):0.098):0.976):0.871);


If instead of all the 4-taxon subsets, you just want to use a random
sample of 10 4-taxon subsets:
```julia
d=readTrees2CF("treefile.txt",whichQ="rand",numQ=10)
```
Be careful to use a numQ value smaller than the total number of possible
4-taxon subsets, which is *n choose 4* on *n* taxa (e.g. 15 on 6 taxa).

If you have already a table of CF values in a file *tableCF.txt*
in this format

|Taxon1 | Taxon2 | Taxon3 | Taxon4 | CF12_34 | CF13_24 | CF14_23 |
|-------|:-------|:-------|:-------|:--------|:--------|:--------|

you would read it like this:
```julia
d=readTableCF("tableCF.txt");
```
You can access this example file
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/tableCF.txt)
or
[here](https://raw.githubusercontent.com/crsl4/PhyloNetworks/master/examples/tableCF.txt).

Columns need to be in the right order. If you have the information on the number of genes used for each 4-taxon subset, you can add this information in a column named "ngenes".

If you have a tree *startTree.txt* in parenthetical format to
use as starting point for the optimization, you can read it with
```julia
T=readTopologyLevel1("startTree.txt")
writeTopology(T)
```
You can access this example file
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/startTree.txt)
(raw file
[here](https://raw.githubusercontent.com/crsl4/PhyloNetworks/master/examples/startTree.txt)).

#### Network Estimation

To estimate the network using the input data
*d* and starting from tree (or network) *T*, do this:

```julia
net1=snaq!(T,d,filename="net1_snaq");
net2=snaq!(T,d,hmax=2, filename="net2_snaq");
```
The option *hmax* corresponds to the maximum number of hybridizations allowed,
1 by default.
The function ends with ! because it modifies the argument d by including the expected CF.

The estimation function creates a .out file (snaq.out by default) with the estimated
network in parenthetical format, which you can also print directly to the screen like this:
```julia
net1
writeTopology(net1)
writeTopology(net1,di=true)
```
The option *di=true* is for the parenthetical format used by
[Dendroscope](http://dendroscope.org/) (without reticulation heritabilities).
Copy this parenthetical description and paste it into Dendroscope, or use the plotting function described below.

#### Network Visualization
To visualize the network:
```julia
p = plot(net1)
```

This function will open a browser where the plot will appear. To get a pdf version of the plot:
```julia
using Gadfly
draw(PDF("mynetwork.pdf", 4inch, 4inch),p)
```
The plot function has many options, type `?plot` to get a list inside Julia.

For a list of all the functions in the PhyloNetworks package, and all
the options on the SNaQ function, refer to the [PDF
documentation](https://github.com/crsl4/PhyloNetworks/blob/master/docs/PhyloNetworks.pdf).

### Bootstrap

You can run a bootstrap analysis if you estimated CFs with credibility intervals,
such as if you used the TICR
pipeline (see above). The TICR pipeline provides a CF table with extra columns for
credibility intervals.
```julia
using DataFrames
df = readtable("tableCFCI.txt", separator=';')
net_bs = bootsnaq(T,df,hmax=1,nrep=10, bestNet=net1, runs=3)
```
You can access this example file
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/tableCFCI.txt)
or [here](https://raw.githubusercontent.com/crsl4/PhyloNetworks/master/examples/tableCFCI.txt).

#### Summarizing bootstrap results

The `bootsnaq` function will return a list of 10 networks (`nrep=10`) which you can then summarize with
```julia
df_bs,tree1 = treeEdgesBootstrap(net_bs,net1)
```
This will calculate the major tree `tree1` displayed in `net1`, that is,
the tree obtained by following the major parent (gamma>0.5) of each hybrid node.
This tree can be visualized like this, with edge numbers shown for later use.
```julia
using Gadfly
plot(tree1, showEdgeNumber=true)
```
Next, we can look at table `df_bs`, which has row for
each tree edge in `net1`. One column contains the edge number
(same as shown in the plot) and another column contains the edge's
bootstrap support: the proportion of bootstrap replicates in which this edge was
found in the major tree of the inferred network.
We can see which tree edges have bootstrap support lower than 100% with
```julia
df_bs[df_bs[:bs] .< 1.0, :]
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
For a given network topology, you can optimize the branch lengths and
inheritance probabilities with the pseudolikelihood. Minus the logarithm of the
pseudolikelihood value for the network will be printed to the screen (the lower the better).
```julia
net1topo = readTopologyLevel1("(2,(4,(3,(5,(6,#H1)))),(1)#H1);");
topologyMaxQPseudolik!(net1topo,d)
writeTopology(net1topo)
```
This is useful if the user has a few candidate networks to compare.
Each network can be optimized individually, and the network with the best
pseudolikelihood can be chosen.
For a network with given branch lengths and heritabilies, we can compute the pseudolikelihood with:
```julia
net1withBL = readTopologyLevel1("(2,(4,(3,(5,(6,#H6:1.0::0.288):5.006):0.518):0.491):1.533,(1)#H6:1.0::0.712);");
topologyQPseudolik!(net1withBL,d)
```
This function is not maximizing the pseudolikelihood, it is simply computing the
pseudolikelihood for the given branch lenghts and probabilities of
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
p = plot(df,layer(x="obsCF1",y="expCF1",Geom.point,Theme(default_color=colorant"orange")),
            layer(x="obsCF2",y="expCF2",Geom.point,Theme(default_color=colorant"purple")),
            layer(x="obsCF3",y="expCF3",Geom.point,Theme(default_color=colorant"blue")),
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
