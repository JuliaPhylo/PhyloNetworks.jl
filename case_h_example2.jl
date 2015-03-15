# Case H example
# Claudia March 2015
# longer branches for it be more identifiable

# types in "types.jl"
include("types.jl")
include("functions.jl")

tree = "((((6:1.1,4:1.5):1.1,#H1),7:1.2):1.1,(8)#H1,10:1.1);" # Case H
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");

#printEdges(net)
#printNodes(net)
