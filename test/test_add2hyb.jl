# adding more than one hybrid
# claudia may 2015

tree = "(((((((1,2),3),4),5),(6,7)),(8,9)),10);"
currT0 = readTopologyLevel1(tree);
#printEdges(currT0)
besttree = deepcopy(currT0);

Random.seed!(16);
successful,hybrid,flag,nocycle,flag2,flag3 = PhyloNetworks.addHybridizationUpdate!(besttree);
@test all([successful, flag, flag2, flag3, !nocycle, hybrid.hybrid])
@test hybrid.number == 11
@test sum(e.inCycle==11 for e in besttree.edge) == 4
@test hybrid.k == 4
@test !hybrid.isBadDiamondI
@test !hybrid.isBadDiamondII

successful,hybrid,flag,nocycle,flag2,flag3 = PhyloNetworks.addHybridizationUpdate!(besttree); #will add a bad triangle
# seed was chosen such that we tried to add a bad triangle. We should notice.
@test all([!successful, flag, !flag2, flag3, !nocycle])
@test hybrid.k == 3
@test hybrid.isVeryBadTriangle
ed = hybridEdges(hybrid)
@test ed[1].isMajor
@test ed[1].gamma > 0.5
@test ed[1].hybrid

# did not recognize as bad diamond II
tree = "(6,(5,#H7:0.0):9.970714072991349,(3,(((2,1):0.2950382234364404,4):0.036924483697671304)#H7:0.00926495670648208):1.1071489442240392);"
net = readTopologyLevel1(tree);
net.node[10].isBadDiamondII || error("does not recognize as bad diamond II")
