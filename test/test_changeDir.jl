# tests with change direction update
# if done twice, does it leave the same network?
# Claudia March 2015

# Case G
include("../examples/case_g_example.jl");
net0 = deepcopy(net);
printEdges(net)
printNodes(net)
node = searchHybridNode(net);
node[1].number
success,newnode = changeDirectionUpdate!(net,node[1]);
printEdges(net)

success,node = changeDirectionUpdate!(net,newnode);
printEdges(net)
printEdges(net0)
# not the same, bad diamond I found in between

# Case H
include("../examples/case_h_example.jl");
net0 = deepcopy(net);
printEdges(net)
printNodes(net)
node = searchHybridNode(net);
node[1].number
success,newnode = changeDirectionUpdate!(net,node[1]);
printEdges(net)

success,node = changeDirectionUpdate!(net,newnode);
printEdges(net)
printEdges(net0)
# not the same, bad diamond II found in between

# Case J
include("../examples/case_j_example.jl");
net0 = deepcopy(net);
printEdges(net)
printNodes(net)
node = searchHybridNode(net);
node[1].number
success,newnode = changeDirectionUpdate!(net,node[1]);
printEdges(net)

success,node = changeDirectionUpdate!(net,newnode);
printEdges(net)
printEdges(net0)
