# test to extract quartet from a "good" 5taxon network (Case G)
# and calculate expCF
# Claudia November 2014


# types in "types.jl"
include("../types.jl")
include("../functions.jl")

# needed modules:
using Base.Collections # for updateInCycle with priority queue

tree = "((((6:0.1,4:1.5)1:0.2,(7)11#H1)5:0.1,(11#H1,8)),10:0.1);" # Case G
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");

printEdges(net)
printNodes(net)
net.names

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);

qnet = extractQuartet(net,q1)
printEdges(qnet)
printNodes(qnet)

calculateExpCFAll!(qnet)
# has errors: calculates t1 wrong, and qnet.leaf is not updated when extracting quartet (still has 5 leaves)
