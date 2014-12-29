# test of calculateExpCF function
# Claudia December 2014

# Case G -----------------

include("../case_g_example.jl")
net.names

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet(net,q1);
printEdges(qnet)
printNodes(qnet)

# en vez de calculateExpCFAll directo, hacer paso a paso a ver donde esta el error
calculateExpCFAll!(qnet)
# has errors: calculates t1 wrong, and qnet.leaf is not updated when extracting quartet (still has 5 leaves)



q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
qnet = extractQuartet(net,q2);
printEdges(qnet)
printNodes(qnet)

q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet(net,q3);
printEdges(qnet)
printNodes(qnet)

q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet(net,q4);
printEdges(qnet)
printNodes(qnet)

q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);
qnet = extractQuartet(net,q5);
printEdges(qnet)
printNodes(qnet)
