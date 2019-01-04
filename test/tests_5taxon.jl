# tests for addHybridizationUpdate! for the 5 taxon network
# we have in the ipad summary
# Claudia September 2014
##########################################################

@warn "to run tests_5taxon.jl, you need to set updateGammaz to return true always. this is because we have bad triangles in here"
# types in "types.jl"
include("../src/types.jl")

# needed modules:
using Base.Collections # for updateInCycle with priority queue

# test functions
include("test_functions_5taxon.jl")

@warn "BUG IN CASE C: sometimes it shows errors, but if Julia is closed and reopened, no more error shown"

tests = ["F","G","H","J","I"];
wrong = String[];

for t in tests
    include("../src/functions.jl")
    include("tree_example.jl");
    tp = string("add_hybrid_case","$(t).jl");
    println("running $(tp)-----");
    try
        include(tp)
    catch
        println("error in $(tp)");
        push!(wrong,t);
    end
end

if(!isempty(wrong))
    for t in wrong
        include("../src/functions.jl")
        include("tree_example.jl");
        tp = string("add_hybrid_case","$(t).jl");
        println("running $(tp)");
        include(tp)
    end
else
    println("----------NO ERRORS!----------");
end


