# function that chooses the right index1, index2
# to add hybridization for Case E to tree_example.jl
# warning: requires other functions not updated here
#          as this function will only run as part of
#          tests_5taxon.jl
# it then deletes the hybridization and tests
# if the result is the original tree example
# Claudia September 2015


# we need a different chooseEdgesGamma to control the edges1,2 chosen
# change index1, index2 according to the specific case
# warning: not real chooseEdgesGamma, the real one is in functions.jl
function chooseEdgesGamma(net::HybridNetwork)
@warn "function chooseEdgesGamma is deterministic"
    index1 = 1;
    index2 = 4;
    edge1 = net.edge[index1];
    edge2 = net.edge[index2];
    gamma = rand()*0.5;
    return edge1, edge2, gamma
end


success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(net);
deleteHybridizationUpdate!(net,hybrid,false);
testTree(net)
