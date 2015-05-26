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
[n.number for n in net.partition[1]]
[n.number for n in net.partition[2]]
[n.number for n in net.partition[3]]
[n.number for n in net.partition[4]]
[n.number for n in net.partition[5]]
[n.number for n in net.partition[6]]

success = false
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
success
printEdges(besttree)
printNodes(besttree)
writeTopology(besttree,true)
net = deepcopy(besttree);
length(net.partition)
[n.number for n in net.partition[1]]
[n.number for n in net.partition[2]]
[n.number for n in net.partition[3]]
[n.number for n in net.partition[4]]
[n.number for n in net.partition[5]]
[n.number for n in net.partition[6]]
[n.number for n in net.partition[7]]
[n.number for n in net.partition[8]]
[n.number for n in net.partition[9]]
[n.number for n in net.partition[10]]
