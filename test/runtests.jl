# see readme file in tests/ for description of tests
# Claudia July 2015
# modified to using PhyloNetworks always, all test files have commented out
# the include(...) or the using PhyloNetworks part
# Claudia May 2016

using Test
using PhyloNetworks
using CSV # for reading files
using DataFrames
using Distributed # parallel in test_correctLik.jl and test_bootstrap.jl
using GLM # for coef, nobs, residuals etc.
using LinearAlgebra: norm, diag, logdet, PosDefException # LinearAlgebra.rotate! not brought into scope
using Random
using StaticArrays # for rate substitution matrices
using Statistics
using StatsBase # for aic etc., stderr
using BioSymbols



    ## readTopology
    getIndexEdge = PhyloNetworks.getIndexEdge
    getIndexNode = PhyloNetworks.getIndexNode
    Edge = PhyloNetworks.Edge
    Node = PhyloNetworks.Node
    setNode! = PhyloNetworks.setNode!
    ## calculateExpCF
    approxEq = PhyloNetworks.approxEq
    searchHybridNode = PhyloNetworks.searchHybridNode
    ## add2hyb
    hybridEdges = PhyloNetworks.hybridEdges

    ## orderings_plot
    RootMismatch = PhyloNetworks.RootMismatch
    fuseedgesat! = PhyloNetworks.fuseedgesat!
    ## compareNetworks
    deletehybridedge! = PhyloNetworks.deletehybridedge!
    displayedNetworks! = PhyloNetworks.displayedNetworks!
 

tests = [
    "test_auxiliary.jl",
    "test_generatetopology.jl",
    "test_addHybrid.jl",
    "test_graph_components.jl",
    "test_manipulateNet.jl", "test_compareNetworks.jl",
    "test_bootstrap.jl",
    "test_moves_semidirected.jl",
    "test_parsimony.jl", ##fails. parsimony search broken
    "test_recursion_matrices.jl",
    "test_calibratePairwise.jl", "test_relaxed_reading.jl",
    "test_isMajor.jl", "test_interop.jl",
    "test_nj.jl"
]


anyerrors = false

for t in tests
    global anyerrors
    try
        @info "starting $t"
        include(t)
        println("\033[1m\033[32mPASSED\033[0m: $t")
    catch
        anyerrors = true
        println("\033[1m\033[31mFAILED\033[0m: $t")
    end
end
println("-------------------------------------")

if anyerrors
    throw("Tests failed")
else
    println("\033[1m\033[32mTests passed")
end
