# test of calculateExpCF function
# Claudia December 2014

# Case G -----------------

println("------ Case G ----------")
include("../case_g_example.jl");
#net.names

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet(net,q1);
#printEdges(qnet)
#printNodes(qnet)
try
    identifyQuartet!(qnet)
    qnet.which != 1 ? error("qnet which not correctly assigned") : nothing
    qnet.hybrid[1].k != 3 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
    qnet.hybrid[1].typeHyb != 2 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
    qnet.hybrid[1].prev.number != -6 ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 0 || qnet.numHybrids != 0 ? error("qnet should not have hybrid nodes anymore") : nothing
    qnet.t1 != 0.2-log(1-0.1*(1-exp(-1.1))) ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [1,1,2,2] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [2,1,2] ? error("qnet.formula not correctly assigned") :nothing

    calculateExpCF!(qnet)
    qnet.expCF != [1/3*exp(-qnet.t1),1-2/3*exp(-qnet.t1),1/3*exp(-qnet.t1)] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 1")
end



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
