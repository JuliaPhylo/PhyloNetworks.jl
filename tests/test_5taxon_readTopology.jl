# tests with the 5 taxon networks read from parenthetical format
# Claudia November 2014

include("../types.jl")
include("../functions.jl")

using Base.Collections # for updateInCycle with priority queue

tree = "(((6:0.1,4:1.5)1:0.2,7:0.2)5:0.1,8:0.1,10:0.1);" # normal tree
tree = "((((6:0.1,4:1.5),(7:0.2)11#H1),11#H1),8:0.1,10:0.1);" # Case C: bad triangle II
tree = "(((6:0.1,(4)11#H1)1:0.2,(11#H1,7))5:0.1,8:0.1,10:0.1);" # Case F: bad diamond I
tree = "((((6:0.1,4:1.5)1:0.2,(7)11#H1)5:0.1,(11#H1,8)),10:0.1);" # Case G
tree = "((((6,4),#H1),7),(8)#H1,10);" # Case H
tree = "((((6)#H1,4),7),8,(#H1,10));" # Case J
tree = "((((6,4))#H1,(#H1,7)),8,10);" # Case D Bad triangle I
tree = "(((((8,10))#H1,7),#H1),6,4);" # Case E Bad triangle I
tree = "((((8,10))#H1,7),6,(4,#H1));" # Case I Bad triangle II
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)

net = readTopology("prueba_tree.txt");
printEdges(net)
printNodes(net)
net = readTopologyUpdate("prueba_tree.txt");

