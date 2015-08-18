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

To install version 0.3.*, go to http://julialang.org/downloads/
PhyloNetworks was developed under Julia.
version 0.3.5, and has been tested on different versions of 0.3.X.
We have not tested its robustness on Julia version 0.4 or above yet.


### Installation of the package

#### PhyloNetworks

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



#### Visualization

here John's readme file along with a small example


Dissertation work on julia

Claudia August 2014

### Steps to setup git:

git clone https://github.com/crsl4/CFnetworks (only once)

**in Mac:git_laptop**

1. git pull (origin master (or other branch)):

better to do: git fetch
              git status (to check if remote is ahead of local)
              git merge FETCH_HEAD

2. make changes

3. git add .

4. git commit -m "message"

5. git push (origin master (or other branch))

**in stat:git_work**

need to compare to version in ane/public/quartetNetwork

2. run compare_jl.sh (with argument "classes",...) to check if two versions are equal, or check date file

1. git pull (origin master (or other branch))

better to do: git fetch
              git status (to check if remote is ahead of local)
              git merge FETCH_HEAD


3. adapt changes from public, consider different branches

4. make new changes if needed

5. git add .

6. git commit -m "message"

7. git push (origin master (or other branch))

8. run sync_public.sh to have same version in ane/public, add date in date file

### Create branches

git branch bla

git checkout bla

## Implementation of pseudolikelihood in Julia

Need to log into desk00 for julia

```julia
include("types.jl")
include("functions.jl")
include("case_f_example.jl")
```

## Procedure to create hybrid network:

1. create edges defined as hybrid or not

2. create nodes defined as hybrid, leaves or tree with such edges

3. setNode! to add nodes into edges

4. create hybrid network

5. updateInCycle! updateGammaz! updateGamma2z! updateContainRoot!
