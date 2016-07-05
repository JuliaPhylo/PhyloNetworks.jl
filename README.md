# PhyloNetworks: analysis for phylogenetic networks in [Julia](http://julialang.org)

## Maximum pseudolikelihood estimation of species network: SNaQ <img src="http://pages.stat.wisc.edu/~claudia/Images/snaq.png" align=right title="SNaQ logo" width=262.5 height=111>
<!-- ![SNaQ logo](http://pages.stat.wisc.edu/~claudia/Images/snaq.png)
original size: 525px × 222px-->

[![Build Status](https://travis-ci.org/crsl4/PhyloNetworks.jl.svg)](https://travis-ci.org/crsl4/PhyloNetworks.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://crsl4.github.io/PhyloNetworks.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://crsl4.github.io/PhyloNetworks.jl/latest)
<!--
[![Coverage Status](https://coveralls.io/repos/crsl4/PhyloNetworks/badge.svg?branch=master&service=github)](https://coveralls.io/github/crsl4/PhyloNetworks?branch=master)
-->

SNaQ implements the statistical inference method in Sol&iacute;s-Lemus and An&eacute;
[(2016)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896).
The procedure involves a
numerical optimization of branch lengths and inheritance probabilities
and a heuristic search in the space of phylogenetic
networks.

If you use the package, and for details on SNaQ, please cite

- Claudia Sol&iacute;s-Lemus and C&eacute;cile An&eacute; (2016).
  Inferring Phylogenetic Networks with Maximum Pseudolikelihood under Incomplete Lineage Sorting.
  [PLoS Genet](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896)
  12(3):e1005896. doi: 10.1371/journal.pgen.1005896

To get help, check

- the [latest documentation](https://crsl4.github.io/PhyloNetworks.jl/latest)
- the [wiki](https://github.com/crsl4/PhyloNetworks.jl/wiki) for a step-by-step tutorial
  (June 2016) with background on networks
- the [google group](https://groups.google.com/forum/#!forum/phylonetworks-users)
  for common questions. Join the group to post/email your questions,
  or to receive information on new versions, bugs fixed, etc.
