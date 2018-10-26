# test function from a debugging case in simScript.jl
# in julia/debug/bigSimulation/n6/julia
# seed 4545
# Claudia August 2015

# to debug problem
tree="(3,(2,(((6,(5)#H9:0.91507):0.93066,(4,#H9:0.0):0.73688):0.0)#H7:1.79104::0.99498):0.11675,(1,#H7:0.04487::0.00502):0.4897);"
net0=readTopologyLevel1(tree);
#printEdges(net0)
net0.node[6].gammaz =1.0  #-5
net0.node[8].gammaz =0.067 #-7

q1 = Quartet(1,["6","5","4","1"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net0,q1);
#printEdges(qnet)

identifyQuartet!(qnet)
qnet.which != 2 ? error("qnet which not correctly assigned") : nothing
#printEdges(qnet)

eliminateHybridization!(qnet)
qnet.which != 2 ? error("qnet which not correctly assigned") : nothing
qnet.numHybrids == 1 || error("qnet not correctly eliminated hyb")
qnet.hybrid[1].k == 4 || error("qnet not correclty identified")
qnet.hybrid[1].typeHyb == 5 || error("qnet not correclty identified")
qnet.hybrid[1].isBadDiamondI || error("qnet forgot it is a bad diamondI")

updateSplit!(qnet)
qnet.split != [-1,-1,-1,-1] ? error("qnet.split not correctly assigned") : nothing

updateFormula!(qnet)
qnet.formula != [-1,-1,-1] ? error("qnet.formula not correctly assigned") : nothing

calculateExpCF!(qnet)
qnet.expCF[2] < 0. || error("expCF not correctly calculated")

