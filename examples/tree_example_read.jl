# Tree example
# Claudia February 2015


# types in "types.jl"
include("types.jl")
include("functions.jl")

tree = "(((6:0.1,4:1.5)1:0.2,7:0.2)5:0.1,8:0.1,10:0.1);" # normal tree
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");

#printEdges(net)
#printNodes(net)


