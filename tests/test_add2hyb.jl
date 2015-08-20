# test to start debugging adding more than one hybrid
# claudia may 2015


include("../src/types.jl")
include("../src/functions.jl")

seed=1234
tree = "(((((((1,2),3),4),5),(6,7)),(8,9)),10);"
currT0 = readTopologyUpdate(tree);
printEdges(currT0)
besttree = deepcopy(currT0);
srand(1234)
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
printEdges(besttree)
printNodes(besttree)
writeTopology(besttree,true)
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
printEdges(besttree)
printNodes(besttree)
writeTopology(besttree,true)

#----add hyb by parts
include("../src/types.jl")
include("../src/functions.jl")
seed=1234
tree = "(((((((1,2),3),4),5),(6,7)),(8,9)),10);"
currT0 = readTopologyUpdate(tree);
printEdges(currT0)
besttree = deepcopy(currT0);
srand(1234)
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(besttree);
printEdges(besttree)
hybrid = addHybridization!(besttree,false);
printEdges(besttree)
flag=false
nocycle=true
flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(besttree,hybrid);
printEdges(besttree)
updateMajorHybrid!(besttree,hybrid);
flag2, edgesGammaz = updateGammaz!(besttree,hybrid,false);
flag3, edgesRoot = updateContainRoot!(besttree,hybrid);
