# test to see if the likelihood is correctly calculated
# and if the networks are correctly estimated
# Claudia August 2015

# -------------------5taxon tree------------------

include("src/types.jl")
include("src/functions.jl")

df = readtable("Tree_output.csv")
d = readTableCF(df)

# starting tree:
tree = "((6,4),(7,8),10);"
currT = readTopologyUpdate(tree);
printEdges(currT)

extractQuartet!(currT,d)
calculateExpCFAll!(d)
lik = logPseudoLik(d)

approxEq(lik,193.7812623319291) || error("not correct likelihood calculated with tree")

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

extractQuartet!(currT,d)
calculateExpCFAll!(d)
lik = logPseudoLik(d)

approxEq(lik,50.17161079450669) || error("not correct likelihood calculated with tree")

estNet = snaq(currT,d,hmax=1,runs=1,returnNet=true);

approxEq(estNet.loglik,0.0,0.01,1.e18) || estNet.loglik < 0.2 || error("not correct net estimated")


