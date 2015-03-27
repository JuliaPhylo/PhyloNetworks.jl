# Case J example
# Claudia January 2015


# types in "types.jl"
include("types.jl")
include("functions.jl")

tree = "((((6)#H1,4:1.5):0.2,7:0.2):0.1,8:0.1,(#H1,10));" # Case J
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");

#printEdges(net)
#printNodes(net)
