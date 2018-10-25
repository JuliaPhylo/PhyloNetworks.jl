# tests for deleteHybridizationUpdate! for the 5 taxon network
# we have in the ipad summary
# Claudia September 2014
##########################################################

@warn "to run tests_5taxon_delete.jl, you need to set updateGammaz to return true always. this is because we have bad triangles in here"
# types in "types.jl"
include("../src/types.jl")

# needed modules:
using Base.Collections # for updateInCycle with priority queue

# test functions
include("test_functions_5taxon.jl")

tests = ["C","F","G","H","J","D","E","I"];
wrong = String[];
t="C"
for t in tests
    include("../src/functions.jl")
    include("tree_example.jl");
    tp = string("delete_case","$(t).jl");
    println("running $(tp)");
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
        tp = string("delete_case","$(t).jl");
        println("running $(tp)");
        include(tp)
    end
else
    println("----------NO ERRORS!----------");
end

