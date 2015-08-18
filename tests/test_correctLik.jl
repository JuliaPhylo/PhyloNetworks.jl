# test to see if the likelihood is correctly calculated
# and if the networks are correctly estimated
# Claudia August 2015

# -------------------5taxon tree------------------

include("../types.jl")
include("../functions.jl")

df = readtable("Tree_output.csv")
d = readTableCF(df)

# starting tree:
tree = "((6,4),(7,8),10);"
currT = readTopologyUpdate(tree);
printEdges(currT)

estTree = snaq(currT,d,hmax=0,returnNet=true);

approxEq(estTree.loglik,0.0) || error("not correct tree estimated")

# ------------------5taxon network 1 hybridization: Case H-----------------
# starting topology: Case G
tree = "((((6:0.1,4:1.5)1:0.2,(7)11#H1)5:0.1,(11#H1,8)),10:0.1);" # Case G
currT = readTopologyUpdate(tree);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output.csv")
d = readTableCF(df)

estNet = snaq(currT,d,hmax=1,returnNet=true)




