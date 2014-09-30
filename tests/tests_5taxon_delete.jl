# tests for deleteHybridizationUpdate! for the 5 taxon network
# we have in the ipad summary
# Claudia September 2014
##########################################################

# types in "types.jl"
include("../types.jl")

# needed modules:
using Base.Collections # for updateInCycle with priority queue

# test functions
include("test_functions_5taxon.jl")

tests = ["C","F","G","H","J","D","E","I"];
wrong = String[];

for t in tests
    include("../functions.jl")
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
        include("../functions.jl")
        include("tree_example.jl");
        tp = string("delete_case","$(t).jl");
        println("running $(tp)");
        include(tp)
    end
else
    println("----------NO ERRORS!----------");
end

