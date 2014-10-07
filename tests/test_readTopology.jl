# initial tests for readTopology
# reads only a tree and with taxon number, not taxon names
# Claudia October 2014
#########################################################

include("../types.jl")
include("../functions.jl")

using Base.Collections # for updateInCycle with priority queue

net = readTopology("prueba_tree3.txt");
printEdges(net)
printNodes(net)


net = readTopology("prueba_tree4.txt");
printEdges(net)
printNodes(net)
