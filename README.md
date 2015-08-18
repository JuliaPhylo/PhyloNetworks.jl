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

To install version 0.3.X, go to http://julialang.org/downloads/

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