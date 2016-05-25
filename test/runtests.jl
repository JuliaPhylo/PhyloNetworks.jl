# see readme file in tests/ for description of tests
# Claudia July 2015
# modified to using PhyloNetworks always, all test files have commented out
# the include(...) or the using PhyloNetworks part
# Claudia May 2016

using Base.Test

if !isdefined(:localtests) localtests = false; end

if(!localtests)
    using PhyloNetworks
    PhyloNetworks.setCHECKNET(true)
else
    const CHECKNET = true #for debugging only
    include("../src/types.jl")
    include("../src/functions.jl")
end

tests = ["test_5taxon_readTopology.jl", "test_calculateExpCF.jl", "test_calculateExpCF2.jl", "test_hasEdge.jl", "test_parameters.jl","test_correctLik.jl",
         "test_partition.jl", "test_partition2.jl","test_deleteHybridizationUpdate.jl", "test_add2hyb.jl", "test_optBLparts.jl",
         "test_orderings_plot.jl", "test_compareNetworks.jl"]#, "test_readme.jl"]


anyerrors = false

for t in tests
    try
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
