# test to see that deleteHybridizationUpdate undoes all attributes
# prompted by Cecile finding cases when containRoot was not updated
# Claudia December 2015

if !isdefined(:individualtest) individualtest = false; end

if(individualtest)
    include("../src/types.jl")
    include("../src/functions.jl")
    const DEBUG = true
end

if isdefined(:PhyloNetworks)
    PhyloNetworks.CHECKNET || error("need CHECKNET==true in PhyloNetworks to test snaq in test_correctLik.jl")
else
    CHECKNET || error("need CHECKNET==true to test snaq in test_correctLik.jl")
end

tree = "(((((((1,2),3),4),5),(6,7)),(8,9)),10);"

#seed = 2738
seed = 56326

currT0 = readTopologyLevel1(tree);
## printEdges(currT0)
## printNodes(currT0)
## writeTopologyLevel1(currT0)
checkNet(currT0)
srand(seed);
besttree = deepcopy(currT0);

# ===== first hybridization ==========================
println("starting first hybridization")
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
success || error("not able to place first hybridization")
#printEdges(besttree)
#printNodes(besttree)
println(writeTopologyLevel1(besttree))
net = deepcopy(besttree);
# test contain root
!net.edge[15].containRoot || error("edge 15 wrong contain Root")
!net.edge[19].containRoot || error("edge 19 wrong contain Root")
!net.edge[20].containRoot || error("edge 20 wrong contain Root")

# test inCycle
net.edge[12].inCycle == 11 || error("wrong incycle")
net.edge[13].inCycle == 11 || error("wrong incycle")
net.edge[16].inCycle == 11 || error("wrong incycle")
net.edge[18].inCycle == 11 || error("wrong incycle")
net.edge[19].inCycle == 11 || error("wrong incycle")
net.edge[20].inCycle == 11 || error("wrong incycle")

net.node[12].inCycle == 11 || error("wrong incycle")
net.node[13].inCycle == 11 || error("wrong incycle")
net.node[16].inCycle == 11 || error("wrong incycle")
net.node[17].inCycle == 11 || error("wrong incycle")
net.node[19].inCycle == 11 || error("wrong incycle")
net.node[20].inCycle == 11 || error("wrong incycle")

# test partition
length(net.partition) == 6 || error("wrong partition")
[n.number for n in net.partition[1].edges] == [15] || error("wrong partition")
[n.number for n in net.partition[2].edges] == [11] || error("wrong partition")
[n.number for n in net.partition[3].edges] == [10] || error("wrong partition")
[n.number for n in net.partition[4].edges] == [9,7,5,3,1,2,4,6,8] || error("wrong partition")
[n.number for n in net.partition[5].edges] == [17] || error("wrong partition")
[n.number for n in net.partition[6].edges] == [14] || error("wrong partition")

# ===== second hybridization ==========================
println("starting second hybridization")
success = false
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
success || error("could not add second hybridization")
#printEdges(besttree)
#printNodes(besttree)
println(writeTopologyLevel1(besttree,true))
net = deepcopy(besttree);

# test contain root
!net.edge[15].containRoot || error("edge 15 wrong contain Root")
!net.edge[19].containRoot || error("edge 19 wrong contain Root")
!net.edge[20].containRoot || error("edge 20 wrong contain Root")
!net.edge[3].containRoot || error("edge 15 wrong contain Root")
!net.edge[1].containRoot || error("edge 19 wrong contain Root")
!net.edge[2].containRoot || error("edge 20 wrong contain Root")
!net.edge[23].containRoot || error("edge 15 wrong contain Root")
!net.edge[22].containRoot || error("edge 19 wrong contain Root")

# test inCycle
net.edge[12].inCycle == 11 || error("wrong incycle")
net.edge[13].inCycle == 11 || error("wrong incycle")
net.edge[16].inCycle == 11 || error("wrong incycle")
net.edge[18].inCycle == 11 || error("wrong incycle")
net.edge[19].inCycle == 11 || error("wrong incycle")
net.edge[20].inCycle == 11 || error("wrong incycle")

net.node[12].inCycle == 11 || error("wrong incycle")
net.node[13].inCycle == 11 || error("wrong incycle")
net.node[16].inCycle == 11 || error("wrong incycle")
net.node[17].inCycle == 11 || error("wrong incycle")
net.node[19].inCycle == 11 || error("wrong incycle")
net.node[20].inCycle == 11 || error("wrong incycle")

net.edge[9].inCycle == 13 || error("wrong incycle")
net.edge[7].inCycle == 13 || error("wrong incycle")
net.edge[5].inCycle == 13 || error("wrong incycle")
net.edge[22].inCycle == 13 || error("wrong incycle")
net.edge[23].inCycle == 13 || error("wrong incycle")

net.node[5].inCycle == 13 || error("wrong incycle")
net.node[7].inCycle == 13 || error("wrong incycle")
net.node[9].inCycle == 13 || error("wrong incycle")
net.node[21].inCycle == 13 || error("wrong incycle")
net.node[22].inCycle == 13 || error("wrong incycle")


# test partition
length(net.partition) == 10 || error("wrong partition")
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

## # ===== identify containRoot for net.node[21]
## net0=deepcopy(net);
## hybrid = net.node[21];
## hybrid.isBadDiamondI
## nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,hybrid);
## [n.number for n in edgesInCycle] == [23,9,7,5,22] || error("wrong identified cycle")
## [n.number for n in nodesInCycle] == [13,14,-4,-5,-6] || error("wrong identified cycle")
## edgesRoot = identifyContainRoot(net,hybrid);
## [n.number for n in edgesRoot] == [3,1,2] || error("wrong identified contain root")
## [e.containRoot for e in edgesRoot] == [false,false,false] || error("edges root should be false contain root")
## edges = hybridEdges(hybrid);
## [e.number for e in edges] == [22,23,3] || error("wrong identified edges")
## undoGammaz!(hybrid,net);
## othermaj = getOtherNode(edges[1],hybrid);
## othermaj.number == -6 || error("wrong othermaj")
## edgesmaj = hybridEdges(othermaj);
## [e.number for e in edgesmaj] == [22,5,4] || error("wrong identified edges")
## edgesmaj[3].containRoot || error("wrong edgesmaj[3] contain root")
## undoContainRoot!(edgesRoot);
## [e.containRoot for e in edgesRoot] == [true,true,true] || error("edges root should be false contain root")
## deleteHybrid!(hybrid,net,true, false)
## [e.containRoot for e in edgesRoot] == [true,true,true] || error("edges root should be false contain root")
## printEdges(net)

# ================= delete second hybridization =============================
println("starting deletion of second hybridization")
deleteHybridizationUpdate!(net,net.node[21], false,false);
length(net.partition) == 6 || error("wrong partition")
# 15,11,10,[9,7,5,3,1,2,4,6,8],17,14
[n.number for n in net.partition[1].edges] == [15] || error("wrong partition")
[n.number for n in net.partition[2].edges] == [11] || error("wrong partition")
[n.number for n in net.partition[3].edges] == [10] || error("wrong partition")
[n.number for n in net.partition[4].edges] == [17] || error("wrong partition")
[n.number for n in net.partition[5].edges] == [14] || error("wrong partition")
[n.number for n in net.partition[6].edges] == [3,1,2,8,6,4,9,7,5] || error("wrong partition")
#printNodes(net)
#printEdges(net)
#printPartitions(net)
println(writeTopologyLevel1(net))

# test contain root
!net.edge[15].containRoot || error("edge 15 wrong contain Root")
!net.edge[19].containRoot || error("edge 19 wrong contain Root")
!net.edge[20].containRoot || error("edge 20 wrong contain Root")
net.edge[1].containRoot || error("edge 1 wrong contain Root")
net.edge[2].containRoot || error("edge 2 wrong contain Root")
net.edge[3].containRoot || error("edge 3 wrong contain Root")

# test inCycle
net.edge[12].inCycle == 11 || error("wrong incycle")
net.edge[13].inCycle == 11 || error("wrong incycle")
net.edge[16].inCycle == 11 || error("wrong incycle")
net.edge[18].inCycle == 11 || error("wrong incycle")
net.edge[19].inCycle == 11 || error("wrong incycle")
net.edge[20].inCycle == 11 || error("wrong incycle")

net.node[12].inCycle == 11 || error("wrong incycle")
net.node[13].inCycle == 11 || error("wrong incycle")
net.node[16].inCycle == 11 || error("wrong incycle")
net.node[17].inCycle == 11 || error("wrong incycle")
net.node[19].inCycle == 11 || error("wrong incycle")
net.node[20].inCycle == 11 || error("wrong incycle")


# =============== delete first hybridization ===================
println("starting deletion of first hybridization")
deleteHybridizationUpdate!(net,net.node[19]);
checkNet(net)
length(net.partition) == 0 || error("wrong partition")
#printEdges(net)

