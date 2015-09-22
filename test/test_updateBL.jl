# test for updateBL function
# Claudia April 2015

include("../examples/tree_example.jl");
printEdges(net)

parts = edgesParts(net);
[isInternalEdge(e)?e.number:0 for e in net.edge]
i = 1
[n.number for n in parts[i].part1]
[n.number for n in parts[i].part2]
[n.number for n in parts[i].part3]
[n.number for n in parts[i].part4]


include("../src/types.jl")
include("../src/functions.jl")
net=readTopologyUpdate("1_astral.out");
printEdges(net)
net.names
df0 = readtable("HGT_truenet_expCF.csv")
d = readTableCF(df0); #expCF

parts = edgesParts(net);
df = makeTable(net,parts,d)
x = updateBL!(net,d)
