# test to see if update partition works
# Claudia May 2015

tree = "(((((((1,2),3),4),5),(6,7)),(8,9)),10);"

# earlier seeds: 2738, 56326 up to v0.14.2
seed = 41

currT0 = readTopologyLevel1(tree);
Random.seed!(seed);
besttree = deepcopy(currT0);
successful,_ = PhyloNetworks.addHybridizationUpdate!(besttree);
@test successful
@test_logs PhyloNetworks.writeTopologyLevel1(besttree,true)
net = deepcopy(besttree);
length(net.partition)
@test Set([e.number for e in p.edges] for p in net.partition) ==
    Set([[15], [11], [10], [9,7,5,3,1,2,4,6,8], [17], [14]])

successful = false
successful,_ = PhyloNetworks.addHybridizationUpdate!(besttree);
@test successful
@test_logs PhyloNetworks.writeTopologyLevel1(besttree,true)
net = deepcopy(besttree);
@test length(net.partition) == 9
@test Set([e.number for e in p.edges] for p in net.partition) ==
    Set([[14], [11], [10], [17], [15], [22,8,9], [3,1,2], [4], [6]])

deleteHybridizationUpdate!(net,net.node[21]);
@test length(net.partition) == 6
@test Set([e.number for e in p.edges] for p in net.partition) ==
    Set([[14], [11], [10], [17], [15], [22,8,9,3,1,2,4,21,5]])

deleteHybridizationUpdate!(net,net.node[18]);
@test length(net.partition) == 0
