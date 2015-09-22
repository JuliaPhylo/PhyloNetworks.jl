# test tree move NNI in a simple quartet
# Claudia February 2015


include("../src/types.jl")
include("../src/functions.jl")

using Base.Collections # for updateInCycle with priority queue

tree = "((1:0.1,2:0.2):0.5,3:0.3,4:0.4);"
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)

net = readTopology("prueba_tree.txt");
printEdges(net)
printNodes(net)

edge = chooseEdgeNNI(net);
println("$(edge.number)")
NNI!(net,edge)

for e in net.edge
    e.inCycle = -1
end

# Case 1
net.edge[1].inCycle = 2
net.edge[3].inCycle = 2
net.edge[4].inCycle = 2

flag,edge = chooseEdgeNNI(net,10);
NNI!(net,edge)

# Case 2
net.edge[1].inCycle = 2
net.edge[3].inCycle = 2
net.edge[4].inCycle = 2
net.edge[2].inCycle = 2

flag,edge = chooseEdgeNNI(net,10);
NNI!(net,edge)

# Case 3
net.edge[1].inCycle = 2
net.edge[2].inCycle = 3

flag,edge = chooseEdgeNNI(net,10);
NNI!(net,edge)

# Case 4
net.edge[1].inCycle = 2
net.edge[2].inCycle = 2

flag,edge = chooseEdgeNNI(net,10);
NNI!(net,edge)

# Case 5
net.edge[1].inCycle = 2
net.edge[2].inCycle = 2
net.edge[4].inCycle = 3
net.edge[5].inCycle = 3

flag,edge = chooseEdgeNNI(net,10);
NNI!(net,edge)

# Case 1
net.edge[2].inCycle = 2
net.edge[3].inCycle = 2
net.edge[4].inCycle = 2

net.node[3].inCycle = 2
net.node[6].inCycle = 2

flag,edge = chooseEdgeNNI(net,10)
NNI!(net,edge)
printEdges(net)
printNodes(net)
