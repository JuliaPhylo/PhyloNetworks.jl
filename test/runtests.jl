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
using LinearAlgebra
using Random
using StaticArrays # for rate substitution matrices
using Statistics
using StatsBase # for aic etc., stderr
using BioSymbols


PhyloNetworks.setCHECKNET(true)

    ## readTopology
    getIndexEdge = PhyloNetworks.getIndexEdge
    getIndexNode = PhyloNetworks.getIndexNode
    Edge = PhyloNetworks.Edge
    Node = PhyloNetworks.Node
    setNode! = PhyloNetworks.setNode!
    ## calculateExpCF
    approxEq = PhyloNetworks.approxEq
    Quartet = PhyloNetworks.Quartet
    extractQuartet! = PhyloNetworks.extractQuartet!
    identifyQuartet! = PhyloNetworks.identifyQuartet!
    eliminateHybridization! = PhyloNetworks.eliminateHybridization!
    updateSplit! = PhyloNetworks.updateSplit!
    updateFormula! = PhyloNetworks.updateFormula!
    calculateExpCF! = PhyloNetworks.calculateExpCF!
    parameters! = PhyloNetworks.parameters!
    searchHybridNode = PhyloNetworks.searchHybridNode
    updateInCycle! = PhyloNetworks.updateInCycle!
    updateContainRoot! = PhyloNetworks.updateContainRoot!
    updateGammaz! = PhyloNetworks.updateGammaz!
    ## correctLik
    calculateExpCFAll! = PhyloNetworks.calculateExpCFAll!
    logPseudoLik = PhyloNetworks.logPseudoLik
    optTopRun1! = PhyloNetworks.optTopRun1!
    ## partition
    addHybridizationUpdate! = PhyloNetworks.addHybridizationUpdate!
    deleteHybridizationUpdate! = PhyloNetworks.deleteHybridizationUpdate!
    ## partition2
    writeTopologyLevel1 = PhyloNetworks.writeTopologyLevel1
    printPartitions = PhyloNetworks.printPartitions
    cleanBL! = PhyloNetworks.cleanBL!
    cleanAfterRead! = PhyloNetworks.cleanAfterRead!
    identifyInCycle = PhyloNetworks.identifyInCycle
    updatePartition! = PhyloNetworks.updatePartition!
    ## deleteHybridizationUpdate
    checkNet = PhyloNetworks.checkNet
    ## add2hyb
    hybridEdges = PhyloNetworks.hybridEdges
    ## optBLparts
    update! = PhyloNetworks.update!
    ## orderings_plot
    RootMismatch = PhyloNetworks.RootMismatch
    fuseedgesat! = PhyloNetworks.fuseedgesat!
    ## compareNetworks
    deleteHybridEdge! = PhyloNetworks.deleteHybridEdge!
    displayedNetworks! = PhyloNetworks.displayedNetworks!
    ## perfect data
    writeExpCF = PhyloNetworks.writeExpCF
    optBL! = PhyloNetworks.optBL!
    ## traitLikDiscrete
    P = PhyloNetworks.P
    P! = PhyloNetworks.P!

tests = [
    "test_5taxon_readTopology.jl", "test_calculateExpCF.jl", "test_calculateExpCF2.jl",
    "test_hasEdge.jl", "test_parameters.jl", "test_correctLik.jl",
    "test_partition.jl", "test_partition2.jl", "test_deleteHybridizationUpdate.jl", "test_add2hyb.jl", "test_optBLparts.jl", "test_undirectedOtherNetworks.jl",
    "test_manipulateNet.jl", "test_compareNetworks.jl",
    "test_badDiamII.jl",
    "test_multipleAlleles.jl",
    "test_bootstrap.jl",
    "test_perfectData.jl", # "test_readme.jl"
    "test_moves_semidirected.jl",
    "test_lm.jl", "test_lm_tree.jl", "test_traits.jl", "test_simulate.jl",
    "test_parsimony.jl",
    "test_calibratePairwise.jl", "test_relaxed_reading.jl",
    "test_isMajor.jl", "test_interop.jl",
    "test_traitLikDiscrete.jl",
    "test_readInputData.jl",
    "test_nj.jl"
]

@show PhyloNetworks.CHECKNET

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
