# test to extract quartet from a "good" 5taxon network (Case G)
# Claudia November 2014
# also test to extract all quartets from a Data object
# Claudia January 2015

# Case G ---------

include("../examples/case_g_example.jl");
net.names

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q1);
printEdges(qnet)
printNodes(qnet)

q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q2);
printEdges(qnet)
printNodes(qnet)

q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q3);
printEdges(qnet)
printNodes(qnet)

q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q4);
printEdges(qnet)
printNodes(qnet)

q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q5);
printEdges(qnet)
printNodes(qnet)

[n.number for n in qnet.leaf]
qnet.numTaxa

# Bad triangle ------------

## include("../bad_triangle_example.jl");
## printEdges(net)
## printNodes(net)
## net.names

## q1 = Quartet(1,["6","7","1","5"],[0.5,0.4,0.1]);
## qnet = extractQuartet!(net,q1)
## printEdges(qnet)
## printNodes(qnet)

## q2 = Quartet(2,["6","5","1","8"],[0.5,0.4,0.1]);
## qnet = extractQuartet!(net,q2);
## printEdges(qnet)
## printNodes(qnet)

## q3 = Quartet(3,["1","7","5","8"],[0.5,0.4,0.1]);
## qnet = extractQuartet!(net,q3);
## printEdges(qnet)
## printNodes(qnet)
## # correct, but does not update length of external edge 8 (should be 0.2, at is 0.1)

## q4 = Quartet(4,["6","1","7","8"],[0.5,0.4,0.1]);
## qnet = extractQuartet!(net,q4);
## printEdges(qnet)
## printNodes(qnet)
## # correct, but does not update length of external edge 9 (should be 0.2, at is 0.1)

## q5 = Quartet(5,["6","7","5","8"],[0.5,0.4,0.1]);
## qnet = extractQuartet!(net,q5);
## printEdges(qnet)
## printNodes(qnet)

## [n.number for n in qnet.leaf]
## qnet.numTaxa

# Bad diamond ------------

include("../examples/case_f_example.jl")
printEdges(net)
printNodes(net)
net.names
parameters!(net)
net.numht

parameters!(qnet,net)

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q1);
printEdges(qnet)
printNodes(qnet)
qnet.indexht
qnet.hasEdge

q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q2);
printEdges(qnet)
printNodes(qnet)

q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q3);
printEdges(qnet)
printNodes(qnet)

q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q4);
printEdges(qnet)
printNodes(qnet)

q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q5);
printEdges(qnet)
printNodes(qnet)

[n.number for n in qnet.leaf]
qnet.numTaxa

# ------------------------------------
# extract all quartets
include("../examples/case_g_example.jl");

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(net,d);
calculateExpCFAll!(d)
