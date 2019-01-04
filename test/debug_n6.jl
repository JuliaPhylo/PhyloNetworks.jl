# tests for the usual functions: optBL, extractQuartet, moves
# for network with n=6 (1_astral.out, add hybridization; true network)
# messy code intended to find specific bugs, does not follow an order
# Claudia April 2015

include("../src/types.jl")
include("../src/functions.jl")

quartets = readListQuartets("allQuartets.txt");

truenet = readTopologyUpdate("truenetwork.txt");
printEdges(truenet)
printNodes(truenet)

tree = readTopologyUpdate("1_astral.out");
printEdges(tree)
writeTopologyLevel1(tree)
net = deepcopy(tree);
Random.seed!(1234)
addHybridizationUpdate!(net);
printEdges(net)
writeTopologyLevel1(net)


extractQuartet!(truenet, quartets) # no errors
extractQuartet!(net,quartets) # errors: after delete leaf!: fixed!

printNodes(net)
nettest = deepcopy(net);
extractQuartet!(nettest,quartets[1]) # error here

qnet = QuartetNetwork(net);
deleteLeaf!(qnet,qnet.node[3])
printEdges(qnet)

## qnet2 = deepcopy(qnet); # for testing
## printNodes(qnet2)

## deleteLeaf!(qnet2,qnet2.node[4]) # no error!!!

deleteLeaf!(qnet,qnet.node[4])
updateHasEdge!(qnet,net)
parameters!(qnet,net) # error here: fixed

# now, error in minor edge 12 set as identifiable AND node.k=-1 even after updateInCycle
tree = readTopologyUpdate("1_astral.out");
printEdges(tree)
writeTopologyLevel1(tree)
net = deepcopy(tree);
Random.seed!(1234)
hybrid = addHybridization!(net);
printNodes(net)
printEdges(net)
flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,hybrid);
net.node[11].k #4, so correct!
updateMajorHybrid!(net,hybrid);
flag2, edgesGammaz = updateGammaz!(net,hybrid,allow);

success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(net);


# ----------------
include("../src/types.jl")
include("../src/functions.jl")

quartets = readListQuartets("allQuartets.txt");
df = readtable("HGT_truenet_expCF.csv")
d = readTableCF(df); #expCF

tree = readTopologyUpdate("1_astral.out");
printEdges(tree)
net = deepcopy(tree);
Random.seed!(1234)
addHybridizationUpdate!(net);
printEdges(net)

currT = deepcopy(net);
optBL!(currT,d)
newT = deepcopy(currT);
count = 0
move = whichMove(currT)
move = :CHdir
move = :MVorigin
move = :MVtarget
move = :nni

movescount = zeros(Int,18)
movesfail = zeros(Int,6)
flag = proposedTop!(move,newT,true,count,1, movescount, movesfail)
printEdges(newT)
printNodes(newT)
optBL!(newT,d)

newT0 = deepcopy(newT);

qnet = QuartetNetwork(newT);
q=quartets[3];
q.taxon
newT.names
extractQuartet!(newT,quartets) # error

qnet = QuartetNetwork(newT);
printNodes(qnet)
deleteLeaf!(qnet,qnet.node[4])
printNodes(qnet)
printEdges(qnet)
qnet0 = deepcopy(qnet);
qnet = deepcopy(qnet0);

identifyQuartet!(qnet)
qnet.which != 1 ? error("qnet which not correctly assigned") : nothing
qnet.hybrid[1].k != 2 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
qnet.hybrid[1].typeHyb != 1 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
qnet.hybrid[1].prev.number != -5 ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

eliminateHybridization!(qnet) #error here
size(qnet.hybrid,1) != 0 || qnet.numHybrids != 0 ? error("qnet should not have hybrid nodes anymore") : nothing
qnet.t1 != 0.2-log(1-0.1*(1-exp(-1.1))) ? error("internal edge length not correctly updated") : nothing

#---
node = qnet.hybrid[1];

updateSplit!(qnet)
qnet.split != [1,1,2,2] ? error("qnet.split not correctly assigned") : nothing

updateFormula!(qnet)
qnet.formula != [2,1,2] ? error("qnet.formula not correctly assigned") : nothing

calculateExpCF!(qnet)
qnet.expCF != [1/3*exp(-qnet.t1),1-2/3*exp(-qnet.t1),1/3*exp(-qnet.t1)] ? error("qnet.expCF wrongly calculated") : nothing

