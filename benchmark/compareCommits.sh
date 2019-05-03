#!/bin/bash
# This benchmark compares the speed of your current version of PhyloNetworks to a 
# version in a previous commit.

# To use, enter PhyloNetworks' benchmark directory and run:
#     bash compareCommits.sh oldCommitNumber

# For example, in .julia/dev/PhyloNetworks/benchmark, run:
#     bash compareCommits.sh oldCommitNumber

# variables:
#    oldCommitNumber: a GitHub commit number

currBranch=$(git branch | sed -e '/^[^*]/d' -e 's/* \(.*\)/\1/')
git checkout $1
julia runBenchmark.jl
git checkout $currBranch
julia runBenchmark.jl