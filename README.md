#PhyloNetworks: analysis for phylogenetic networks in [Julia](http://julialang.org)
## Maximum pseudolikelihood estimation of species network: SNaQ <img src="http://pages.stat.wisc.edu/~claudia/Images/snaq.png" align=right title="SNaQ logo" width=262.5 height=111>
<!-- ![SNaQ logo](http://pages.stat.wisc.edu/~claudia/Images/snaq.png)
original size: 525px × 222px-->

[![Build Status](https://travis-ci.org/crsl4/PhyloNetworks.svg)](https://travis-ci.org/crsl4/PhyloNetworks)

SNaQ implements the statistical inference method in [Sol&iacute;s-Lemus and
An&eacute;](http://arxiv.org/pdf/1509.06075.pdf). The procedure involves a
numerical optimization of branch lengths and inheritance probabilities
and a heuristic search in the space of phylogenetic
networks.

Below there is a quick tutorial on SNaQ and PhyloNetworks, but refer
to the [PDF documentation](https://github.com/crsl4/PhyloNetworks/blob/master/docs/PhyloNetworks.pdf) for more details.

### Input for SNaQ

Two alternatives:

1. A list of estimated gene trees for each locus, which can be obtained using
  [MrBayes](http://mrbayes.sourceforge.net) or [RAxML](http://sco.h-its.org/exelixis/software.html). Or:
2. A table of concordance factors (CF) for each 4-taxon subset which can be
  obtained from [BUCKy](http://www.stat.wisc.edu/~ane/bucky/),
  to account for gene tree uncertainty

### Pipeline from sequence alignments

This [pipeline](https://github.com/nstenz/TICR)
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
WARNING: It is important to update the package regularly as it is undergoing constant development. There is a known bug for Mac users where the *Pkg.update* function does not update to the latest version. We recommend Mac users to do the following through the terminal:
```shell
cd HOME/.julia/v0.4/PhyloNetworks/
git pull
```
where HOME is replaced by your home directory.

The PhyloNetworks package has the following dependencies, but everything is installed automatically.

- GraphViz (version 0.0.3)
- NLopt (version 0.2.0)

The version in parenthesis correspond to the ones used when
implementing PhyloNetworks.

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
*HOME/.julia/v0.3/PhyloNetworks/examples/*. However, links to the
files are also included below.

Suppose you have a file with a list of gene trees in parenthetical
format called *treefile.txt*.
If 'treefile.txt' is in your directory, do this to read in all gene trees
and to summarize them with a list
of quartet CFs:
```julia
d=readTrees2CF("treefile.txt");
```
Make sure to have the semicolon (;) at the end to avoid useless output to the screen!
You can access this example file
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/treefile.txt).
This file contains 10 trees, each in parenthetical format on 6 taxa
like this:

(6:2.728,((3:0.655,5:0.655):1.202,(1:0.881,(2:0.783,4:0.783):0.098):0.976):0.871);


If instead of all the 4-taxon subsets, you just want to use a random
sample of 10 4-taxon subsets:
```julia
d=readTrees2CF("treefile.txt",whichQ="rand",numQ=10);
```
Be careful to use a numQ value smaller than the total number of possible
4-taxon subsets, which is *n choose 4* on *n* taxa (e.g. 15 on 6 taxa).

If you have already a table of CF values in a file *tableCF.txt*
in this format

|Taxon1 | Taxon2 | Taxon3 | Taxon4 | CF12vs34 | CF13vs24 | CF14vs23 |
|-------|:-------|:-------|:-------|:---------|:---------|:---------|

you would read it like this:
```julia
d=readTableCF("tableCF.txt");
```
You can access this example file
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/tableCF.txt).

If you have a tree *startTree.txt* in parenthetical format to
use as starting point for the optimization and want to
update the branch lengths according to the CF already read in the data
structure *d*, do this:
```julia
T=readStartTop("startTree.txt",d);
writeTopology(T)
```
You can access this example file
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/startTree.txt).

#### Network Estimation

To estimate the network using the input data
*d* and starting from tree (or network) *T*, do this:

```julia
net1=snaq(T,d,filename="net1_snaq");
net2=snaq(T,d,hmax=2, filename="net2_snaq");
```
Make sure to have the semicolon (;) at the end, to avoid much useless output
to the screen!
The option *hmax* corresponds to the maximum number of hybridizations allowed,
1 by default.

The estimation function creates a .out file (snaq.out by default) with the estimated
network in parenthetical format, which you can also print directly to the screen like this:
```julia
writeTopology(net1)
writeTopology(net1,di=true)
```
The option *di=true* is for the parenthetical format used by
[Dendroscope](http://dendroscope.org/) (without reticulation heritabilities).
Copy this parenthetical description and paste it into Dendroscope, or use the plotting function described below.

<!--
#### Network Visualization
To visualize the network:
```julia
plotPhylonet(net1)
plotPhylonet(net1,unrooted=true)
```
For now, this function will create an .svg figure file (netImage.svg by default).

WARNING: There is a known bug in the plotPhylonet function,
see the issue in the PhyloNetworks Github repository for details.
The error can be sometimes fixed by changing the position of the root
with the root function.

WARNING: on Mac computers, sometimes the function cannot be called directly sometimes, but this will work:
```julia
PhyloNetworks.plotPhylonet(net)
```
-->

For a list of all the functions in the PhyloNetworks package, and all the options on the SNaQ function, refer to the [PDF documentation](https://github.com/crsl4/PhyloNetworks/blob/master/docs/PhyloNetworks.pdf).

#### Multiple alleles

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
species names instead of allele names.

Estimation will work the same way:
```julia
new_net = snaq(new_T,new_d);
```
where *new_T* should be a starting topology with one tip per species, labelled with the species names.

WARNING: the current function works best if all alleles from the same
individual are given the same name (the individual's 'name') across
all genes for which that individual was sequenced.

##### Optimizing branch lengths and inheritance probabilities for a given network
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

#### Debugging: the .err file
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
