# Input for SNaQ

SNaQ is a method implemented in the package to estimate a phylogenetic network
from multiple molecular sequence alignments. There are two alternatives for the input data:

1. A list of estimated gene trees for each locus, which can be obtained using [MrBayes](http://mrbayes.sourceforge.net) or [RAxML](http://sco.h-its.org/exelixis/software.html). Or:
2. A table of concordance factors (CF), i.e. gene tree frequencies, for each 4-taxon subset. This table can be obtained from [BUCKy](http://www.stat.wisc.edu/~ane/bucky/), to account for gene tree uncertainty

This [pipeline](https://github.com/nstenz/TICR) can be used to obtain the table of
quartet CF needed as input for SNaQ
(see also the [wiki](https://github.com/crsl4/PhyloNetworks.jl/wiki/TICR:-from-alignments-to-quartet-concordance-factors).)
It starts from the sequence alignments,
runs MrBayes and then BUCKy (both parallelized), producing the
table of estimated CFs and their credibility intervals.
Additional details on this [TICR pipeline](@ref)
describe how to insert data at various stages (e.g. after running MrBayes on each locus).

## Tutorial data: gene trees

We suggest that you create a special directory for running these examples,
where input files can be downloaded and where output files will be
created (with estimated networks for instance). Enter this directory
and run Julia from there.

Suppose you have a file with a list of gene trees in parenthetical
format called `raxmltrees.tre`.
You can access the example file of input trees
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/raxmltrees.tre)
or
[here](https://raw.githubusercontent.com/crsl4/PhyloNetworks/master/examples/raxmltrees.tre)
for easier download.

Do not copy-paste into a "smart" text-editor. Instead, save the file
directly into your working directory using "save link as" or "download linked file as".
This file contains 30 gene trees, each in parenthetical format on 6 taxa
like this (with rounded branch lengths):

`(E:0.038,((A:0.014,B:0.010):0.010,(C:0.008,D:0.002):0.010):0.025,O:0.078);`

If `raxmltrees.tre` is in your working directory, you can view its content
within Julia:
```julia
less("raxmltrees.tre")
```
or like this, to view the version downloaded with the package:
```julia
raxmltrees = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","raxmltrees.tre")
less(raxmltrees)
```
Just type `q` to quit viewing this file.
You could read in these 30 trees and visualize the third one (say) like this:
```@example qcf
using PhyloNetworks
raxmltrees = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","raxmltrees.tre");
nothing # hide
```
```@repl qcf
genetrees = readMultiTopology(raxmltrees);
genetrees[3]
```
To visualize any of these input trees, use the
[PhyloPlots](https://github.com/cecileane/PhyloPlots.jl) package:
```@example qcf
using PhyloPlots
using RCall # hide
mkpath("../assets/figures") # hide
R"name <- function(x) file.path('..', 'assets', 'figures', x)" # hide
R"svg(name('inputdata_gene3.svg'), width=4, height=3)" # hide
R"par"(mar=[0,0,0,0])                          # hide
plot(genetrees[3], :R); # tree for 3rd gene
R"dev.off()"                                   # hide
nothing # hide
```
![gene3](../assets/figures/inputdata_gene3.svg)

To read in all gene trees and directly summarize them by a list
of quartet CFs (proportion of input trees with a given quartet):
```@repl qcf
q,t = observedquartetCF(genetrees); # read in trees, calculate quartet CFs
df = writeTableCF(q,t)   # data frame with observed CFs: gene frequencies
using CSV
CSV.write("tableCF.csv", df); # to save the data frame to a file
raxmlCF = readTableCF("tableCF.csv") # read in the file and produces a "DataCF" object
rm("tableCF.csv") # hide
```
`less("tableCF.csv")` lets you see the content of the newly created
file "tableCF.csv", within Julia. Again, type `q` to quit viewing this file.

In this table, each 4-taxon set is listed in one row.
The 3 "CF" columns gives the proportion of genes that has
each of the 3 possible trees on these 4 taxa.

For more help on any function, type `?` to enter the help mode,
then type the name of the function. For example: type `?` then `observedquartetCF`
for information on the various options of that function.

When there are many more taxa, the number of quartets
might be very large and we might want to use a subset to speed things up.
Here, if we wanted to use a random sample of 10 quartets
instead of all quartets, we could do:

`raxmlCF = readTrees2CF(raxmltrees, whichQ="rand", numQ=10, CFfile="tableCF10.txt")`

Be careful to use a numQ value smaller than the total number of possible
4-taxon subsets, which is *n choose 4* on *n* taxa (e.g. 15 on 6 taxa).
To get a predictable random sample, you may set the seed with
`using Random; Random.seed!(12321)`
(for instance) prior to sampling the quartets as above.
The `readTrees2CF` is *much* slower than the function `observedquartetCF`
to read in trees and calculate the quartet CFs observed in the trees,
when we want to get *all* quartet CFs. But for a small sample of quartets,
then `readTrees2CF` is available.

## Tutorial data: quartet CFs

If we already have a table of quartet concordance factor (CF) values
in a file `buckyCF.csv` in this format

| Taxon1 | Taxon2 | Taxon3 | Taxon4 | CF12_34 | CF13_24 | CF14_23
|:-------|:-------|:-------|:-------|:--------|:--------|:-------
| D      | A| E | O|   0.565 |       0.0903 |       0.3447
| ...    |  |   |  |         |              |       ...

we would read it in one step like this: `readTableCF("buckyCF.csv")`.
An example file comes with the package, available
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/buckyCF.csv)
or
[here](https://raw.githubusercontent.com/crsl4/PhyloNetworks/master/examples/buckyCF.csv).

```@repl qcf
buckyCFfile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","buckyCF.csv");
buckyCF = readTableCF(buckyCFfile)
```
The same thing could be done in 2 steps:
first to read the file and convert it to a 'DataFrame' object,
and then to convert this DataFrame into a DataCF object.
```@repl qcf
using CSV, DataFrames
dat = CSV.read(buckyCFfile);
first(dat, 6) # to see the first 6 rows
buckyCF = readTableCF(dat)
writeTableCF(buckyCF)
```
In the input file, columns need to be in the right order:
with the first 4 columns giving the names of the taxa in each 4-taxon set.
The CF values are assumed to be in columns named "CF12_34", etc.,
or else in columns 5,6,7.
If available, a column named "ngenes" will be taken to have the
the number of genes for each 4-taxon subset.

## Tutorial data: starting tree

If we have a tree for the data set at hand,
it can be used as a starting point for the optimization.
From our gene trees, we estimated a species tree with
[ASTRAL](https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md).
This tree comes with the package in file `astral.tre`
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/astral.tre).
This file has 102 trees: 100 bootstrap species trees,
followed by their greedy consensus,
followed by the best tree on the original data.
It's this last tree that we are most interested in.
We can read it with
```@example qcf
astralfile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","astral.tre");
astraltree = readMultiTopology(astralfile)[102] # 102th tree: last tree here
R"svg(name('inputdata_astraltree.svg'), width=4, height=3)" # hide
R"par"(mar=[0,0,0,0]) # hide
plot(astraltree, :R, showEdgeLength=true);
R"dev.off()"; # hide
nothing # hide
```
![astraltree](../assets/figures/inputdata_astraltree.svg)

To start its search, SNaQ will need a network of "level 1".
All trees and all networks with 1 hybridization are of level 1.
To make sure that a network with 2 or more hybridizations is of level 1,
we can read it in with
`readTopologyLevel1` (which also unroots the tree, resolves polytomies,
replaces missing branch lengths by 1 for starting values etc.):
```julia
T=readTopologyLevel1("startNetwork.txt")
```
(here `startNetwork.txt` is a hypothetical file: replace this by
the name of a file that contains your network of interest.)

Next: [Getting a Network](@ref)
