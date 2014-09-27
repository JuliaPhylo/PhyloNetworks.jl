# tests for addHybridizationUpdate! for the 5 taxon network
# we have in the ipad summary
# Claudia September 2014
##########################################################

# types in "types.jl"
include("types.jl")

# functions in "functions.jl"
include("functions.jl")

# needed modules:
using Base.Collections # for updateInCycle with priority queue

# examples
include("tree_example.jl");

# test functions
include("/Users/Clauberry/Documents/phylo/software/CFimplementation/julia/test_functions_5taxon.jl")

# we need a different chooseEdgesGamma to control the edges1,2 chosen
# change index1, index2 according to the specific case
# warning: not real chooseEdgesGamma, the real one is in functions.jl
function chooseEdgesGamma(net::HybridNetwork)
warn("function chooseEdgesGamma is deterministic")
    index1 = 1;
    index2 = 3;
    edge1 = net.edge[index1];
    edge2 = net.edge[index2];
    gamma = rand()*0.5;
    return edge1, edge2, gamma
end

# addHybridization! only to check that the createHybrid is working fine, before updating
#node = addHybridization!(net);

# index1=7, index2=6 => case F (bad diamond)
include("/Users/Clauberry/Documents/phylo/software/CFimplementation/julia/print_add.jl")
testCaseF(net)

# index1=3, index2=7 => case G
include("/Users/Clauberry/Documents/phylo/software/CFimplementation/julia/print_add.jl")
testCaseG(net)

# index1=1, index2=3 => case H
include("/Users/Clauberry/Documents/phylo/software/CFimplementation/julia/print_add.jl")
testCaseH(net)

# index1=5, index2=2 => case J
include("/Users/Clauberry/Documents/phylo/software/CFimplementation/julia/print_add.jl")
testCaseJ(net)

# index1=7, index2=1 => case D (bad triangle)
include("/Users/Clauberry/Documents/phylo/software/CFimplementation/julia/print_add.jl")
testCaseD(net)

# index1=1, index2=4 => case E
include("/Users/Clauberry/Documents/phylo/software/CFimplementation/julia/print_add.jl")
testCaseE(net)

# index1=6, index2=4 => case I
include("/Users/Clauberry/Documents/phylo/software/CFimplementation/julia/print_add.jl")
testCaseI(net)

