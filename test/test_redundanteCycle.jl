# test of code to remove redundante cycles after extracting quartet
# Claudia May 2015


include("/Users/Clauberry/Documents/phylo/software/CFimplementation/julia/git_laptop/CFnetworks/types.jl")
include("/Users/Clauberry/Documents/phylo/software/CFimplementation/julia/git_laptop/CFnetworks/functions.jl")

tree = "(((((((1,2),3),4),5),(6,7)),(8,9)),10);"


seed = 2738
currT0 = readTopologyUpdate(tree);
Random.seed!(seed)
besttree = deepcopy(currT0);
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
success
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
success
printEdges(besttree)
writeTopologyLevel1(besttree,true)
net = deepcopy(besttree);

q1 = Quartet(1,["1","2","3","4"],[0.5,0.4,0.1]);
extractQuartet!(net,q1);
qnet=deepcopy(q1.qnet);
printEdges(qnet)
printNodes(qnet)

## redundantCycle!(qnet,qnet.node[14])
## redundantCycle!(qnet,qnet.node[8])
## printEdges(qnet)
## printNodes(qnet)
identifyQuartet!(qnet)
eliminateHybridization!(qnet)
updateSplit!(qnet)
updateFormula!(qnet)
calculateExpCF!(qnet)
