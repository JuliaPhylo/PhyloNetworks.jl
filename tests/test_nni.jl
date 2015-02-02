# test tree move NNI in a simple quartet
# Claudia February 2015


include("types.jl")
include("functions.jl")

using Base.Collections # for updateInCycle with priority queue

tree = "((1,2),3,4);"
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)

net = readTopology("prueba_tree.txt");
printEdges(net)
printNodes(net)

