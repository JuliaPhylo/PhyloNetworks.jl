# test to start debugging adding more than one hybrid
# claudia may 2015

seed=1234
tree = "(((((((1,2),3),4),5),(6,7)),(8,9)),10);"
currT0 = readTopologyLevel1(tree);
#printEdges(currT0)
besttree = deepcopy(currT0);
Random.seed!(1234);
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
success || error("added hybrid not successful")
flag || error("added hybrid not successful")
flag2 || error("added hybrid not successful")
flag3 || error("added hybrid not successful")
!nocycle || error("added hybrid not successful")
hybrid.hybrid || error("added hybrid not successful")
hybrid.number == 11 || error("added hybrid not successful")
besttree.edge[12].inCycle == 11 || error("incycle not correct")
besttree.edge[13].inCycle == 11 || error("incycle not correct")
besttree.edge[18].inCycle == 11 || error("incycle not correct")
besttree.edge[20].inCycle == 11 || error("incycle not correct")
hybrid.k == 4 || error("not correct hybrid.k")
!hybrid.isBadDiamondI || error("thinks its bad diamond I")
!hybrid.isBadDiamondII || error("thinks its bad diamond II")
#printEdges(besttree)
#printNodes(besttree)
#writeTopology(besttree,true)
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree); #will add a bad triangle
!success || error("added bad triangle and did not notice")
flag || error("added bad triangle and did not notice")
!flag2 || error("added bad triangle and did not notice")
flag3 || error("added bad triangle and did not notice")
!nocycle || error("added bad triangle and did not notice")
hybrid.k == 3 || error("added bad triangle and did not notice")
hybrid.isVeryBadTriangle || error("added bad triangle and did not notice")
ed = hybridEdges(hybrid)
ed[1].isMajor || error("in bad triangle, update major hybrid not working, and it should")
ed[1].gamma > 0.5 || error("in bad triangle, update major hybrid not working, and it should")
ed[1].hybrid || error("in bad triangle, update major hybrid not working, and it should")

# did not recognize as bad diamond II
tree = "(6,(5,#H7:0.0):9.970714072991349,(3,(((2,1):0.2950382234364404,4):0.036924483697671304)#H7:0.00926495670648208):1.1071489442240392);"
net = readTopologyLevel1(tree);
net.node[10].isBadDiamondII || error("does not recognize as bad diamond II")

#print("NO ERRORS!!")

# printEdges(besttree)
# printNodes(besttree)
# writeTopology(besttree,true)

#----add hyb by parts
# include("../src/types.jl")
# include("../src/functions.jl")
# seed=1234
# tree = "(((((((1,2),3),4),5),(6,7)),(8,9)),10);"
# currT0 = readTopologyUpdate(tree);
# printEdges(currT0)
# besttree = deepcopy(currT0);
# Random.seed!(1234);
# success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
# printEdges(besttree)
# ##Random.seed!(5678);
# hybrid = addHybridization!(besttree)
# printEdges(besttree)
# flag=false
# nocycle=true
# flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(besttree,hybrid)
# printEdges(besttree)
# updateMajorHybrid!(besttree,hybrid)
# besttree.root=9
# plot(besttree)
# writeTopology(besttree)
# flag2, edgesGammaz = updateGammaz!(besttree,hybrid,false)
# flag3, edgesRoot = updateContainRoot!(besttree,hybrid)
# plot(besttree)

