#PhyloNetworks: phylogenetic analysis for networks in [Julia](http://julialang.org)

## Maximum pseudolikelihood estimation of species network: SNaQ

SNaQ implements the statistical inference method in [Solis-Lemus and
Ane, submitted](www.stat.wisc.edu/~claudia). The procedure involves a
numerical optimization of branch lengths and probabilities of
inheritance and heuristic search in the space of phylogenetic
networks.  

### Input for SNaQ

Two alternatives: 

1. A table of concordance factors (CF) for each quartet which can be
obtained from [BUCKy](http://www.stat.wisc.edu/~ane/bucky/).

2. The estimated gene trees for each locus which can be obtained by
[MrBayes](http://mrbayes.sourceforge.net) or [RAxML](http://sco.h-its.org/exelixis/software.html).

### Pipeline from sequence alignments

We developed a pipeline that obtains the table of CF needed as input
for SNaQ. The pipeline starts with input the sequence alignments, and
it runs MrBayes and then BUCKy, producing the
table of estimated CF and its confidence intervals.


### Installation of Julia

To install Julia, go to http://julialang.org/downloads/

PhyloNetworks was developed under Julia
version 0.3.5, and has been tested on different versions of 0.3.X.
We have not tested its robustness on Julia version 0.4 or above yet.


### Installation of the package PhyloNetworks

To install the package:

```julia
Pkg.add("PhyloNetworks")
Pkg.build("PhyloNetworks")
```

The PhyloNetworks package has the following dependencies, but everything is installed automatically.

*GraphViz (version...)

*NLopt (version...)

The version in parenthesis correspond to the ones used when
implementing PhyloNetworks.

### Small examples

#### Input data

Suppose you have a file with a list of gene trees in parenthetical
format called *treefile.txt* and you want to use all the possible
quartets for the taxa in those trees to calculate the CF: ```julia
d=readTrees2CF("treefile.txt") ```

If instead you want to use a random sample of 100 quartets:
```julia
d=readTrees2CF("treefile.txt",whichQ=:rand,numQ=100)
```

On the contrary, if you have already the CF in a file *tableCF.txt*
in the format:

|Taxon1 | Taxon2 | Taxon3 | Taxon4 | CF12vs34 | CF13vs24 | CF14vs23 |
|-------|:-------|:-------|:-------|:---------|:---------|:---------|

You would read it like:
```julia
d=readTableCF("tableCF.txt")
```

If you have a tree *startTree.tre* in parenthetical format to
use as starting point for the optimization and want to
update the branch lengths according to the CF already read in the data
structure *d*:
```julia
T=readStartTop("startTree.tre",d);
```

#### Network Estimation To estimate the network using the input data
*d,T*: 
```julia net=snaq(T,d); net=snaq(T,d,hmax=2); ``` 
The
estimation function also creates a .out file with the estimated
network in parenthetical format.

#### Network Visualization
To visualize the network:
```julia
plotPhylonet(net)
plotPhylonet(net,unrooted=true)
```

For a list of all the functions in the PhyloNetworks package, and all the options on the SNaQ function, refer to the PDF documentation.