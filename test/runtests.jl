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
    "test_murepresentation.jl"
]

@testset "PhyloNetworks.jl" begin
    for testfile in tests
        include(testfile)
    end
end

