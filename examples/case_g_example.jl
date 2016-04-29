# Case G example
# Claudia December 2014


# types in "types.jl"
include("../src/types.jl")
include("../src/functions.jl")

tree = "((((6:0.1,4:1.5)1:0.2,(7)11#H1)5:0.1,(11#H1,8)),10:0.1);" # Case G
#f = open("prueba_tree.txt","w")
#write(f,tree)
#close(f)
net = readTopologyUpdate(tree)

#printEdges(net)
#printNodes(net)
