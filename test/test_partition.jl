# test to see if update partition works
# Claudia May 2015

tree = "(((((((1,2),3),4),5),(6,7)),(8,9)),10);"

#seed = 2738
seed = 56326

currT0 = readTopologyLevel1(tree);
Random.seed!(seed);
besttree = deepcopy(currT0);
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
success
#printEdges(besttree)
#printNodes(besttree)
writeTopologyLevel1(besttree,true)
net = deepcopy(besttree);
length(net.partition)
[n.number for n in net.partition[1].edges] == [15] || error("wrong partition")
[n.number for n in net.partition[2].edges] == [11] || error("wrong partition")
[n.number for n in net.partition[3].edges] == [10] || error("wrong partition")
[n.number for n in net.partition[4].edges] == [9,7,5,3,1,2,4,6,8] || error("wrong partition")
[n.number for n in net.partition[5].edges] == [17] || error("wrong partition")
[n.number for n in net.partition[6].edges] == [14] || error("wrong partition")


success = false
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
success
#printEdges(besttree)
#printNodes(besttree)
@test_logs writeTopologyLevel1(besttree,true)
net = deepcopy(besttree);
length(net.partition)
[n.number for n in net.partition[1].edges] == [15] || error("wrong partition")
[n.number for n in net.partition[2].edges] == [11] || error("wrong partition")
[n.number for n in net.partition[3].edges] == [10] || error("wrong partition")
[n.number for n in net.partition[4].edges] == [17] || error("wrong partition")
[n.number for n in net.partition[5].edges] == [14] || error("wrong partition")
[n.number for n in net.partition[6].edges] == [3,1,2] || error("wrong partition")
[n.number for n in net.partition[7].edges] == [21] || error("wrong partition")
[n.number for n in net.partition[8].edges] == [8] || error("wrong partition")
[n.number for n in net.partition[9].edges] == [6] || error("wrong partition")
[n.number for n in net.partition[10].edges] == [4] || error("wrong partition")

deleteHybridizationUpdate!(net,net.node[21]);
length(net.partition)
length(net.partition) == 6 || error("wrong partition")
# 15,11,10,[9,7,5,3,1,2,4,6,8],17,14
[n.number for n in net.partition[1].edges] == [15] || error("wrong partition")
[n.number for n in net.partition[2].edges] == [11] || error("wrong partition")
[n.number for n in net.partition[3].edges] == [10] || error("wrong partition")
[n.number for n in net.partition[4].edges] == [17] || error("wrong partition")
[n.number for n in net.partition[5].edges] == [14] || error("wrong partition")
[n.number for n in net.partition[6].edges] == [3,1,2,21,8,6,4,9,7] || error("wrong partition")
#printNodes(net)

deleteHybridizationUpdate!(net,net.node[18]);
length(net.partition) == 0 || error("wrong partition")
