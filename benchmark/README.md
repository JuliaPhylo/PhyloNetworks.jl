# Using PkgBenchmark to Compare the Efficiency of Two Different Commits using Benchmark

PkgBenchmarks allows us to compare the performance of a package at different branches, commits, or tags. For full documentation, see the PkgBenchmark [documentation here] (https://juliaci.github.io/PkgBenchmark.jl/stable/)

# Comparing Two Commits on Speed using Benchmarks
This benchmark compares the speed of your current version of PhyloNetworks to the
version in a previous commit.

To use, enter PhyloNetworks' benchmark directory and run:
```bash
    bash compareCommits.sh oldCommitNumber
```
For example, in .julia/dev/PhyloNetworks/benchmark, run:
```bash
    bash compareCommits.sh oldCommitNumber
```
variables:
   oldCommitNumber: a GitHub commit number

# Adding New Benchmarks

To add new benchmarks, use the dictionary interface introduced by Benchmarks.jl. [docs here](https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/doc/manual.md#defining-benchmark-suites)

First, open <PKGROOT>/benchmark/benchmarks.jl. In this file, create a new suite. I've created a suite to test nucleic acid substitution models. It has two subparts, jc69 and hky85.
```julia
SUITE["nasm"] = BenchmarkGroup(["jc69", "hky85"])
```
We can then add to this suite:
```julia
SUITE["nasm"]["jc69"] = @benchmarkable JC69([0.5])
```

