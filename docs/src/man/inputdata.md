# Input for SNaQ

SNaQ is a method implemented in the package to estimate a phylogenetic network
from multiple molecular sequence alignments. There are two alternatives for the input data:

1. A list of estimated gene trees for each locus, which can be obtained using
  [MrBayes](http://mrbayes.sourceforge.net) or [RAxML](http://sco.h-its.org/exelixis/software.html). Or:
2. A table of concordance factors (CF) for each 4-taxon subset which can be
  obtained from [BUCKy](http://www.stat.wisc.edu/~ane/bucky/),
  to account for gene tree uncertainty

This [pipeline](https://github.com/nstenz/TICR) can be used to obtain the table of 
quartet CF needed as input for SNaQ. The pipeline starts from the sequence alignments,
runs MrBayes and then BUCKy (both parallelized), producing the
table of estimated CFs and their credibility intervals.
Additional details on this [TICR pipeline](@ref)
# https://github.com/crsl4/PhyloNetworks/blob/master/docs/src/man/ticr_howtogetQuartetCFs.md
 describe how to insert data at various stages (e.g. after running MrBayes on each locus).
 
## Tutorial data

<!--The examples files for this section can be found within the
PhyloNetworks folder, typically in your
*HOME/.julia/v0.4/PhyloNetworks/examples/*. However, links to the
files are also included below.
-->
We suggest that you create a special directory for running these examples,
where input files can be downloaded and where output files will be
created (with estimated networks for instance). Enter this directory
and run Julia from there.

### gene trees

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

### Table of quartet concordance factors (CFs)

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

