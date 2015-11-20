# strange error in simple network
# Claudia November 2015

include("../src/types.jl")
include("../src/functions.jl")

tree= "(A,((B,#H1),(C,(D)#H1)));" #fails
tree= "(A,((C,(D)#H1),(B,#H1)));" #works
net = readTopology(tree)
net = readTopologyUpdate(tree)
printEdges(net)
printNodes(net)
