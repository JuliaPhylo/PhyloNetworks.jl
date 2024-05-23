# PhyloNetworks: analysis for phylogenetic networks <img src="docs/src/logo_text.png" align=right>

[![doc stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://crsl4.github.io/PhyloNetworks.jl/stable)
[![doc dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://crsl4.github.io/PhyloNetworks.jl/dev)
[![Build status](https://github.com/crsl4/PhyloNetworks.jl/workflows/CI/badge.svg?branch=master)](https://github.com/crsl4/PhyloNetworks.jl/actions/workflows/ci.yml)
[![coverage](https://codecov.io/gh/crsl4/PhyloNetworks.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/crsl4/PhyloNetworks.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/P/PhyloNetworks.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/report.html)

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
- estimate species networks from multilocus data: SNaQ
- phylogenetic comparative methods for continuous trait evolution
  on species networks / trees
- plot networks (and trees), via the companion package
  [PhyloPlots](https://github.com/cecileane/PhyloPlots.jl)

To get help, check

- the [latest documentation](https://crsl4.github.io/PhyloNetworks.jl/dev)
- the [wiki](https://github.com/crsl4/PhyloNetworks.jl/wiki) for a step-by-step tutorial
  with background on networks (last revised 2022)
- [tutorial](https://cecileane.github.io/networkPCM-workshop/) for
  comparative methods, including network calibration (2023 workshop)
- the [google group](https://groups.google.com/forum/#!forum/phylonetworks-users)
  for common questions. Join the group to post/email your questions,
  or to receive information on new versions, bugs fixed, etc.

If you use the package, please cite ([bibtex format here](CITATION.bib))

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
  including a tutorial for trait evolution
  and a tutorial for network calibration.

Continuous traits, accounting for within-species variation:

- Benjamin Teo, Jeffrey P. Rose, Paul Bastide & Cécile Ané (2022).
  Accounting for intraspecific variation in continuous trait evolution
  on a reticulate phylogeny.
  [bioRxiv](https://doi.org/10.1101/2022.05.12.490814)

For a discrete trait (influence of gene flow on the trait,
ancestral state reconstruction, rates):

- Karimi, Grover, Gallagher, Wendel, Ané & Baum (2020). Reticulate evolution
  helps explain apparent homoplasy in floral biology and pollination in baobabs
  (*Adansonia*; Bombacoideae; Malvaceae).
  Systematic Biology,
  69(3):462-478. doi: [10.1093/sysbio/syz073](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syz073/5613901?guestAccessKey=a32e7dd3-27fd-4a13-b171-7ff5d6da0e01).
