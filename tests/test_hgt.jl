# test with the network in HGTinconsistency
# Claudia April 2015

include("../types.jl")
include("../functions.jl")

truenetwork = "((((1,2),((3,4))#H1),(#H1,5)),6);"
treefile = "1.ms"

d = readInputData(treefile, "HGTtableCF.txt"); #will get list of all quartets (allQuartets.txt), and obsCF (HGTtableCF.txt)
length(d.quartet)

# for expCF
net = readTopologyUpdate(truenetwork);
printEdges(net)
extractQuartet!(net,d) #extract qnet from that list and calculate expCF
df2 = writeExpCF(d.quartet)
writetable("HGT_truenet_expCF.csv",df2)

# compare table with expCF and obsCF: very different!

# optTopLevel with expCF and starting tree 1_astral.out, no branches updated
df2 = readtable("HGT_truenet_expCF.csv")
d2 = readTableCF(df2); #expCF

currT0 = readTopologyUpdate("1_astral.out");
printEdges(currT0)
srand(1234)
currT = deepcopy(currT0);
addHybridizationUpdate!(currT); #add hybrid at random (different test would be to start with the tree)
printEdges(currT)

@time optTopLevel!(currT,d2,1,false)

# debug1hgtBad.txt : could not find bug!

net2 = readTopologyUpdate("(5,((3,(2,1):0.9184973864504608):0.0,(4)#H7:0.7994139650661589::1.0):1.2422144462061508,(6,#H7:0.0002521174970950474::0.0):0.20997432883063805);");
printEdges(net2)
qnet = QuartetNetwork(net2);
printNodes(qnet)
printEdges(qnet)

deleteLeaf!(qnet,qnet.node[1])

qnet0 = deepcopy(qnet);
deleteLeaf!(qnet,qnet.node[9])
