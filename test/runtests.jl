using Test
using PhyloNetworks
using CSV # for reading files
using DataFrames
using Distributed # for parsimony search, currently broken
using LinearAlgebra: diag, Diagonal # LinearAlgebra.rotate! not brought into scope
using Random
using StableRNGs

const PN = PhyloNetworks

tests = [
    "test_deprecated.jl",
    "test_quartet.jl",
    "test_auxiliary.jl",
    "test_generatetopology.jl",
    "test_addHybrid.jl",
    "test_graph_components.jl",
    "test_manipulateNet.jl",
    "test_compareNetworks.jl",
    "test_bootstrap.jl",
    "test_moves_semidirected.jl",
    "test_parsimony.jl", # has broken tests: parsimony search broken
    "test_recursion_matrices.jl",
    "test_calibratePairwise.jl",
    "test_expectedfstat.jl",
    "test_relaxed_reading.jl",
    "test_isMajor.jl",
    "test_interop.jl",
    "test_nj.jl",
    "test_mu_rep.jl"
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
        @error "Error in $t" exception=(err, catch_backtrace())
    end
end
println("-------------------------------------")

if anyerrors
    
    throw("Tests failed")
else
    println("\033[1m\033[32mTests passed")
end
