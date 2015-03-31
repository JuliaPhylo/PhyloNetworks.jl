# Case C example
# Claudia December 2014


# types in "types.jl"
include("../types.jl")
include("../functions.jl")

# needed modules:
using Base.Collections # for updateInCycle with priority queue

tree = "((((6:0.1,4:1.5),(7:0.2)11#H1),11#H1),8:0.1,10:0.1);" # Case C: bad triangle II

f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");

printEdges(net)
printNodes(net)
