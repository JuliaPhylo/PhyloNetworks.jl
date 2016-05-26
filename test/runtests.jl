# see readme file in tests/ for description of tests
# Claudia July 2015
# modified to using PhyloNetworks always, all test files have commented out
# the include(...) or the using PhyloNetworks part
# Claudia May 2016

using Base.Test

if !isdefined(:localtests) localtests = false; end

#localtests = true

if(!localtests)
    using PhyloNetworks
    using DataFrames
    PhyloNetworks.setCHECKNET(true)

    ## readTopology
    getIndexEdge = PhyloNetworks.getIndexEdge
    getIndexNode = PhyloNetworks.getIndexNode
    Edge = PhyloNetworks.Edge
    Node = PhyloNetworks.Node
    setNode! = PhyloNetworks.setNode!
    ## calculateExpCF
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
else
    const CHECKNET = true #for debugging only
    include("../src/types.jl")
    include("../src/functions.jl")
end

tests = ["test_5taxon_readTopology.jl", "test_calculateExpCF.jl", "test_calculateExpCF2.jl", "test_hasEdge.jl", "test_parameters.jl","test_correctLik.jl",
         "test_partition.jl", "test_partition2.jl","test_deleteHybridizationUpdate.jl", "test_add2hyb.jl", "test_optBLparts.jl",
         "test_orderings_plot.jl", "test_compareNetworks.jl"]#, "test_readme.jl"]

if isdefined(:PhyloNetworks)
    @show PhyloNetworks.CHECKNET
else
    @show CHECKNET
end

anyerrors = false

for t in tests
    try
        info("starting $t")
        include(t)
        println("\t\033[1m\033[32mPASSED\033[0m: $t")
    catch
        anyerrors = true
        println("\t\033[1m\033[31mFAILED\033[0m: $t")
    end
end
println("-------------------------------------")

if anyerrors
    throw("Tests failed")
else
    println("\t\033[1m\033[32mTests passed")
end
