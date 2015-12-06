# see readme file in tests/ for description of tests
# Claudia July 2015

const CHECKNET = true #for debugging only
tests = ["test_5taxon_readTopology.jl", "test_calculateExpCF.jl", "test_calculateExpCF2.jl", "test_hasEdge.jl", "test_parameters.jl","test_correctLik.jl", "test_partition.jl", "test_partition2.jl","test_deleteHybridizationUpdate.jl"]
#tests = ["test_calculateExpCF.jl", "test_calculateExpCF2.jl", "test_hasEdge.jl", "test_parameters.jl","test_correctLik.jl"]

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
