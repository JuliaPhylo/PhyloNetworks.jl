# PhyloNetworks: analysis for phylogenetic networks <img src="docs/src/logo_text.png" align=right>

[![Build Status](https://travis-ci.org/crsl4/PhyloNetworks.jl.svg)](https://travis-ci.org/crsl4/PhyloNetworks.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://crsl4.github.io/PhyloNetworks.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://crsl4.github.io/PhyloNetworks.jl/dev)
[![codecov.io](http://codecov.io/github/crsl4/PhyloNetworks.jl/coverage.svg?branch=master)](http://codecov.io/github/crsl4/PhyloNetworks.jl?branch=master)
[![Coverage Status](https://coveralls.io/repos/crsl4/PhyloNetworks.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/crsl4/PhyloNetworks?branch=master)

## Overview

PhyloNetworks is a [Julia](http://julialang.org) package with utilities to:
- read / write phylogenetic trees and networks,
  in (extended) Newick format.
  Networks are considered explicit: nodes represent ancestral species.
  They can be rooted or unrooted.
- manipulate networks: re-root, prune taxa, remove hybrid edges,
  extract the major tree from a network, extract displayed networks / trees
- compare networks / trees with dissimilarity measures
  (Robinson-Foulds distance on trees)
- summarize samples of bootstrap networks (or trees)
  with edge and node support
- estimate species networks from multilocus data (see below)
- phylogenetic comparative methods for continuous trait evolution
  on species networks / trees
- plot networks (and trees), via the companion package
  [PhyloPlots](https://github.com/cecileane/PhyloPlots.jl)

To get help, check

- the [latest documentation](https://crsl4.github.io/PhyloNetworks.jl/dev)
- the [wiki](https://github.com/crsl4/PhyloNetworks.jl/wiki) for a step-by-step tutorial
  (July 2018) with background on networks
- the [google group](https://groups.google.com/forum/#!forum/phylonetworks-users)
  for common questions. Join the group to post/email your questions,
  or to receive information on new versions, bugs fixed, etc.

If you use the package, please cite

- Claudia Sol&iacute;s-Lemus, Paul Bastide and C&eacute;cile An&eacute; (2017).
  PhyloNetworks: a package for phylogenetic networks.
  [Molecular Biology and Evolution](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msx235/4103410/PhyloNetworks-a-package-for-phylogenetic-networks?guestAccessKey=230afceb-df28-4160-832d-aa7c73f86369)
  34(12):3292–3298.
  [doi:10.1093/molbev/msx235](https://doi.org/10.1093/molbev/msx235)

## Maximum pseudolikelihood estimation of species network: SNaQ <img src="docs/src/snaq.png" align=right title="SNaQ logo" width=262.5 height=111>
<!-- ![SNaQ logo](http://pages.stat.wisc.edu/~claudia/Images/snaq.png)
original size: 525px × 222px-->

SNaQ implements the statistical inference method in Sol&iacute;s-Lemus and An&eacute;
[(2016)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896).
The procedure involves a
numerical optimization of branch lengths and inheritance probabilities
and a heuristic search in the space of phylogenetic
networks.

If you use SNaQ, please cite

- Claudia Sol&iacute;s-Lemus and C&eacute;cile An&eacute; (2016).
  Inferring Phylogenetic Networks with Maximum Pseudolikelihood under Incomplete Lineage Sorting.
  [PLoS Genet](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896)
  12(3):e1005896.
  [doi:10.1371/journal.pgen.1005896](https://doi.org/10.1371/journal.pgen.1005896)

## Phylogenetic comparative methods for trait evolution

For continuous traits, study based on the Brownian motion process,
with or without transgressive evolution after reticulations:

- Bastide, Solís-Lemus, Kriebel, Sparks, Ané (2018).
  Phylogenetic Comparative Methods for Phylogenetic Networks with Reticulations.
  Systematic Biology, 67(5):800–820.
  [doi:10.1093/sysbio/syy033](https://doi.org/10.1093/sysbio/syy033).
  SI on [dryad](http://dx.doi.org/10.5061/dryad.nt2g6)
  including a [tutorial for trait evolution](https://datadryad.org/bitstream/handle/10255/dryad.177752/xiphophorus_PCM_analysis.html?sequence=1)
  and a [tutorial for network calibration](https://datadryad.org/bitstream/handle/10255/dryad.177754/xiphophorus_networks_calibration.html?sequence=1).
  