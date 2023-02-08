using Pkg
Pkg.activate("/Users/cora/.julia/environments/net") #sets up development environment
using PkgBenchmark
using BenchmarkTools
using PhyloNetworks
benchmarkpkg("PhyloNetworks")