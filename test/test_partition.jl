# test to see if update partition works
# Claudia May 2015


include("/Users/Clauberry/Documents/phylo/software/CFimplementation/julia/git_laptop/CFnetworks/types.jl")
include("/Users/Clauberry/Documents/phylo/software/CFimplementation/julia/git_laptop/CFnetworks/functions.jl")

tree = "(((((((1,2),3),4),5),(6,7)),(8,9)),10);"

#seed = 2738
seed = 56326

currT0 = readTopologyUpdate(tree);
srand(seed)
besttree = deepcopy(currT0);
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
success
printEdges(besttree)
printNodes(besttree)
writeTopology(besttree,true)
net = deepcopy(besttree);
length(net.partition)
[n.number for n in net.partition[1].edges]
[n.number for n in net.partition[2].edges]
[n.number for n in net.partition[3].edges]
[n.number for n in net.partition[4].edges]
[n.number for n in net.partition[5].edges]
[n.number for n in net.partition[6].edges]


success = false
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
success
printEdges(besttree)
printNodes(besttree)
writeTopology(besttree,true)
net = deepcopy(besttree);
length(net.partition)
[n.number for n in net.partition[1].edges]
[n.number for n in net.partition[2].edges]
[n.number for n in net.partition[3].edges]
[n.number for n in net.partition[4].edges]
[n.number for n in net.partition[5].edges]
[n.number for n in net.partition[6].edges]
[n.number for n in net.partition[7].edges]
[n.number for n in net.partition[8].edges]
[n.number for n in net.partition[9].edges]
[n.number for n in net.partition[10].edges]

deleteHybridizationUpdate!(net,net.node[21]);
length(net.partition)
length(net.partition) == 6
# 15,11,10,[9,7,5,3,1,2,4,6,8],17,14
[n.number for n in net.partition[1].edges] == [15]
[n.number for n in net.partition[2].edges] == [11]
[n.number for n in net.partition[3].edges] == [10]
[n.number for n in net.partition[4].edges] == [17]
[n.number for n in net.partition[5].edges] == [14]
[n.number for n in net.partition[6].edges] == [3,1,2,21,8,6,4,9,7]
printNodes(net)

deleteHybridizationUpdate!(net,net.node[18]);
length(net.partition) == 0
