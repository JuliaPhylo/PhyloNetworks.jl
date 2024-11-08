# PhyloNetworks: utilities for phylogenetic networks <img src="docs/src/logo_text.png" align=right>

[![doc stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaphylo.github.io/PhyloNetworks.jl/stable)
[![doc dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaphylo.github.io/PhyloNetworks.jl/dev)
[![Build status](https://github.com/juliaphylo/PhyloNetworks.jl/workflows/CI/badge.svg?branch=master)](https://github.com/juliaphylo/PhyloNetworks.jl/actions/workflows/ci.yml)
[![coverage](https://codecov.io/gh/juliaphylo/PhyloNetworks.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaphylo/PhyloNetworks.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/P/PhyloNetworks.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/report.html)

## Overview

PhyloNetworks is a [Julia](http://julialang.org) package with utilities to
to handle phylogenetic trees and networks.
It serves as a core package that other packages can depend on, such as
[PhyloPlots](https://github.com/JuliaPhylo/PhyloPlots.jl)
to visualize phylogenies,
[SNaQ](https://github.com/JuliaPhylo/PhyloPlots.jl)
to infer phylogenies from genetic data,
[PhyloPlots](https://github.com/JuliaPhylo/PhyloPlots.jl)
to analyze the evolution of traits along phylogenies.

Phylogenetic networks represent the evolutionary relationships between a set
of organisms, e.g. populations, species, languages, cultures etc.
They are also called *admixture graphs* when their leaves are populations
of the same or closely related species.
They are *explicit* graph representations, in which nodes represent
ancestral populations / species, and edge lengths represent evolutionary time.

Utilities in this core package include:
- read / write phylogenies in (extended) Newick format
- manipulate networks: re-root, prune taxa, remove hybrid edges,
  transform a network with a semidirected nearest-neighbor interchange (sNNI),
  extract the major tree from a network, extract displayed networks / trees
- compare networks with dissimilarity measures
  (e.g. the Robinson-Foulds distance on trees)
- summarize samples of networks with support for local features
  (edges and clades of hybrid origin or sister to a hybrid clade)
- fit edges lengths from average pairwise distances between leaves,
  using least-squares
- network traversal routines

To get help, check

- the [latest documentation](https://juliaphylo.github.io/PhyloNetworks.jl/dev)
- the [wiki](https://github.com/juliaphylo/PhyloNetworks.jl/wiki) for a
  step-by-step tutorial with background on networks (last revised 2022)
- [tutorial](https://cecileane.github.io/networkPCM-workshop/) for
  comparative methods, including network calibration (2023 workshop)
- the [google group](https://groups.google.com/forum/#!forum/phylonetworks-users)
  for common questions. Join the group to post/email your questions,
  or to receive information on new versions, bugs fixed, etc.

If you use the package, please cite ([bibtex format here](CITATION.bib)).
For the PhyloNetworks package in particular, please cite:

- Claudia Sol&iacute;s-Lemus, Paul Bastide and C&eacute;cile An&eacute; (2017).
  PhyloNetworks: a package for phylogenetic networks.
  [Molecular Biology and Evolution](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msx235/4103410/PhyloNetworks-a-package-for-phylogenetic-networks?guestAccessKey=230afceb-df28-4160-832d-aa7c73f86369)
  34(12):3292–3298.
  [doi:10.1093/molbev/msx235](https://doi.org/10.1093/molbev/msx235)
