# Case I example: Bad diamond II
# Claudia January 2015


# types in "types.jl"
include("../src/types.jl")
include("../src/functions.jl")

tree = "((((8,10))#H1,7),6,(4,#H1));" # Case I Bad diamond II
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");

#printEdges(net)
#printNodes(net)
