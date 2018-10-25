# function that chooses the right index1, index2
# to add hybridization for Case J to tree_example.jl
# warning: requires other functions not updated here
#          as this function will only run as part of
#          tests_5taxon.jl
# Claudia September 2015


# we need a different chooseEdgesGamma to control the edges1,2 chosen
# change index1, index2 according to the specific case
# warning: not real chooseEdgesGamma, the real one is in functions.jl
function chooseEdgesGamma(net::HybridNetwork)
@warn "function chooseEdgesGamma is deterministic"
    index1 = 5;
    index2 = 2;
    edge1 = net.edge[index1];
    edge2 = net.edge[index2];
    gamma = rand()*0.5;
    return edge1, edge2, gamma
end


# index1=5, index2=2 => case J
include("print_add.jl")
testCaseJ(net)
