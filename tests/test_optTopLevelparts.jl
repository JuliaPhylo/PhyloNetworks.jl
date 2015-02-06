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
include("tree_example.jl");
