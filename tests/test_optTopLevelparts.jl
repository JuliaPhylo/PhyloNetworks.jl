# test the parts in optTopLevel
# Claudia February 2015

# ------------------5taxon network 1 hybridization-----------------
# starting topology: Case G
include("../case_g_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output.csv")
d = readDataCF(df)
#df2 = readtable("CaseG_output.csv")

currloglik,currxmin = optBL!(currT,d)
newT = deepcopy(currT);
count = 0
N = 10
#move = whichMove(currT)
move = :CHdir
move = :MVorigin
move = :MVtarget
move = :nni

flag = proposedTop!(move,newT,true,count,N)
flag
printEdges(newT)
printNodes(newT)

newloglik, newxmin = optBL!(newT,d)
newloglik - currloglik

currT = deepcopy(newT);
currloglik = newloglik
currxmin = newxmin

# -------------------5taxon tree------------------
include("../tree_example_read.jl");
printEdges(net)
q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(net,d)

df = writeExpCF(d.quartet)
writetable("CaseG_output.csv",df)

