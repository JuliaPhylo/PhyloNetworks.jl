# To test readTopology we need many networks and trees to read
# INPUT: text file name with tree/network in parenthetical format
# Claudia (October 2014)

include("../types.jl")
include("../functions.jl")

using Base.Collections

file = "prueba_tree.txt" # change here

net = readTopology(file) # warning: run without ; to see errors
printEdges(net)
printNodes(net)
