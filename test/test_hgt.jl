# test with the network in HGTinconsistency
# Claudia April 2015

include("../src/types.jl")
include("../src/functions.jl")

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

# ----------------- optTopLevel with expCF and starting tree 1_astral.out, no branches updated --------
include("../src/types.jl")
include("../src/functions.jl")

df2 = readtable("HGT_truenet_expCF.csv")
d2 = readTableCF(df2); #expCF

truenetwork = "((((1,2),((3,4))#H1),(#H1,5)),6);"
net = readTopologyUpdate(truenetwork);
printEdges(net)
@time optBL!(net,d2) #loglik~1e.-15 in 0.17secs

currT0 = readTopologyUpdate("1_astral.out");
printEdges(currT0)
Random.seed!(1234) #found right network in 135secs, wrong BL, right gamma (debug3hgtGood.txt)
Random.seed!(4568) # bug with movedownlevel, will leave for later (debug4hgt.txt)
Random.seed!(11233) #found right network in 119.34secs (debug4hgtBad.txt, by mistake)
currT = deepcopy(currT0);
addHybridizationUpdate!(currT); #add hybrid at random (different test would be to start with the tree)
printEdges(currT)

@time optTopLevel!(currT,d2,1)

# ----------------- optTopLevel with expCF and starting tree 1_astral.out, branches updated --------
include("../src/types.jl")
include("../src/functions.jl")

df2 = readtable("HGT_truenet_expCF.csv")
d2 = readTableCF(df2); #expCF

truenetwork = "((((1,2),((3,4))#H1),(#H1,5)),6);"
net = readTopologyUpdate(truenetwork);
printEdges(net)
@time optBL!(net,d2) #loglik~1e.-15 in 0.17secs

currT0 = readTopologyUpdate("1_astral.out");
x = updateBL!(currT0,d2)
printEdges(currT0)
Random.seed!(1234) #right network in 20secs(debug12hgt)
Random.seed!(4568) #movedownlevel: debug13hgtBad
Random.seed!(11233) #very close to right network in 135secs (debug8hgtgood.txt)
currT = deepcopy(currT0);
addHybridizationUpdate!(currT); #add hybrid at random (different test would be to start with the tree)
printEdges(currT)

@time optTopLevel!(currT,d2,1);

# ----------------- optTopLevel with obsCF and starting tree 1_astral.out, no branches updated --------
include("../src/types.jl")
include("../src/functions.jl")

df = readtable("HGTtableCF.txt") #from 1.ms
d = readTableCF(df); #obsCF

truenetwork = "((((1,2),((3,4))#H1),(#H1,5)),6);"
net = readTopologyUpdate(truenetwork);
printEdges(net)
@time optBL!(net,d,true) #loglik ~144.74 for true net

currT0 = readTopologyUpdate("1_astral.out");
printEdges(currT0)
Random.seed!(1234) #local max found debug5hgtBad
Random.seed!(4568) #local max2 found debug6hgtBad
Random.seed!(11233) #local max found debug7hgtBad
currT = deepcopy(currT0);
addHybridizationUpdate!(currT); #add hybrid at random (different test would be to start with the tree)
printEdges(currT)

@time optTopLevel!(currT,d,1);

# ----------------- optTopLevel with obsCF and starting tree 1_astral.out, branches updated --------
include("../src/types.jl")
include("../src/functions.jl")

df = readtable("HGTtableCF.txt") #from 1.ms
d = readTableCF(df); #obsCF

truenetwork = "((((1,2),((3,4))#H1),(#H1,5)),6);"
net = readTopologyUpdate(truenetwork);
printEdges(net)
@time optBL!(net,d,true) #loglik ~144.74 for true net

currT0 = readTopologyUpdate("1_astral.out");
x = updateBL!(currT0,d)
printEdges(currT0)
Random.seed!(1234) #local max found debug9hgt.txt
Random.seed!(4568) #local max2 found debug9hgt.txt
Random.seed!(11233) #local max3 found debug9hgt.txt
currT = deepcopy(currT0);
addHybridizationUpdate!(currT); #add hybrid at random (different test would be to start with the tree)
printEdges(currT)

@time optTopLevel!(currT,d,1);


# -----------------------
# bug in Random.seed!(4568)
net = readTopologyUpdate("((4,#H7:9.99670403892172::0.43454301575229803):1.5467254857425556,((6,(5)#H7:2.512064322645178::0.565456984247702):9.221085796210835,(2,1):0.38003642076628485):0.0,3);");
Random.seed!(4568)
flag = proposedTop!(:nni,net,true,0,100, zeros(Int,18), zeros(Int,6))
df = readtable("HGTtableCF.txt") #from 1.ms
d = readTableCF(df); #obsCF

qnet=QuartetNetwork(net);
extractQuartet!(net,d)
printEdges(qnet)
printNodes(qnet)
deleteLeaf!(qnet,qnet.node[1]) #leaf1
qnet0=deepcopy(qnet);
deleteLeaf!(qnet,qnet.node[2]) #leaf4

