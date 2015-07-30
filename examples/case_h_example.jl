# Case H example
# Claudia January 2015


# types in "types.jl"
include("../types.jl")
include("../functions.jl")

tree = "((((6:0.1,4:1.5):0.1,#H1),7:0.2):0.1,(8)#H1,10:0.1);" # Case H
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");


#printEdges(net)
#printNodes(net)
