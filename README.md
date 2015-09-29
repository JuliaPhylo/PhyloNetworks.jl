#PhyloNetworks: analysis for phylogenetic networks in [Julia](http://julialang.org)
## Maximum pseudolikelihood estimation of species network: SNaQ <img src="http://pages.stat.wisc.edu/~claudia/Images/snaq.png" align=right title="SNaQ logo" width=262.5 height=111>
<!-- ![SNaQ logo](http://pages.stat.wisc.edu/~claudia/Images/snaq.png)
original size: 525px × 222px-->

[![Build Status](https://travis-ci.org/crsl4/PhyloNetworks.svg)](https://travis-ci.org/crsl4/PhyloNetworks.jl)

SNaQ implements the statistical inference method in [Sol&iacute;s-Lemus and
An&eacute;](http://arxiv.org/pdf/1509.06075.pdf). The procedure involves a
numerical optimization of branch lengths and inheritance probabilities
and a heuristic search in the space of phylogenetic
networks.

Below there is a quick tutorial on SNaQ and PhyloNetworks, but refer
to the [PDF documentation](https://github.com/crsl4/PhyloNetworks/blob/master/docs/PhyloNetworks.pdf) for more details.

### Input for SNaQ

Two alternatives:

1. A table of concordance factors (CF) for each 4-taxon subset which can be
obtained from [BUCKy](http://www.stat.wisc.edu/~ane/bucky/).

2. The estimated gene trees for each locus which can be obtained by
[MrBayes](http://mrbayes.sourceforge.net) or [RAxML](http://sco.h-its.org/exelixis/software.html).

### Pipeline from sequence alignments

This [pipeline](https://github.com/nstenz/TICR)
can be used to obtain the table of CF needed as input
for SNaQ. The pipeline starts with input the sequence alignments, and
it runs MrBayes and then BUCKy, producing the
table of estimated CFs and their confidence intervals.


### Installation of Julia

Julia is a high-level and interactive programming language (like R or Matlab),
but it is also high-performance (like C).
To install Julia, go to http://julialang.org/downloads/

PhyloNetworks was developed under Julia
version 0.3.5, and has been tested on different versions of 0.3.X.
We have not tested its robustness on Julia version 0.4 or above yet.


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
cd HOME/.julia/v0.3/PhyloNetworks/
git pull
```
where HOME is replaced by your home directory and v0.3 could be replaced by v0.4 if you have version 0.4 of Julia.

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
d=readTrees2CF("treefile.txt",whichQ=:rand,numQ=10);
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

If you have a tree *startTree.tre* in parenthetical format to
use as starting point for the optimization and want to
update the branch lengths according to the CF already read in the data
structure *d*, do this:
```julia
T=readStartTop("startTree.tre",d);
writeTopology(T)
```
You can access this example file
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/startTree.tre).

#### Network Estimation

To estimate the network using the input data
*d* and starting from tree (or network) *T*, do this:

```julia
net1=snaq(T,d);
net2=snaq(T,d,hmax=2);
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
Copy this parenthetical description and paste it into Dendroscope.

<!--#### Network Visualization
To visualize the network:
```julia
plotPhylonet(net)
plotPhylonet(net,unrooted=true)
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

##### Optimizing branch lengths and inheritance probabilities for a given network
For a given network topology, you can optimize the branch lengths and
inheritance probabilities with the pseudolikelihood. Minus the logarithm of the
pseudolikelihood value for the network will be printed to the screen.
```julia
topologyMaxQPseudolik!(net1,d)
topologyMaxQPseudolik!(net1,d,verbose=true)
```
This is in particular useful if the user has some network alternatives to test out individually and choose the best in terms of the pseudolikelihood among them.

For a given network with branch lengths, we can compute the pseudolikelihood with:
```julia
topologyQPseudolik!(net1,d)
topologyQPseudolik!(net1,d,verbose=true)
```
This function is not maximizing, it is simply computing the
pseudolikelihood for the given branch lenghts and probabilities of
inheritance. It is advisable though to maximize the pseudolikelihood even
when branch lengths and probabilities of inheritance are given. This
is simply because the optimization procedure is quite fast (so it does
not represent any computational burden) and it will provide a more
accurate pseudolikelihood value.

#### Debugging: the .err file
Please report any bugs and errors to *claudia@stat.wisc.edu*. The easiest way to do it is by checking the .err file which will show the number of runs that failed by a bug and the corresponding seed to replicate the run. This is an example of what the .err file looks like:

Total errors: 1 in seeds [4545]

You need to run the following function (with the two global variable sets) with the same settings that caused the error:

```julia
const DEBUG = true
const REDIRECT = true
snaqDebug(T,d,hmax=2,seed=4545)
```




This will create two files:
*snaqDebug.log* and *debug.log* which you can then send to
*claudia@stat.wisc.edu* with subject "SNaQ bug found" or something
similar. I will not have access to any part of your data, the files
simply print out the steps to retrace the bug, and hopefully fix it.
