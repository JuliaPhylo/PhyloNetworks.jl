# test the parts in optTopLevel
# Claudia February 2015

# ------------------5taxon network 1 hybridization-----------------
# starting topology: Case G
include("../examples/case_g_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output.csv")
d = readDataCF(df)
#df2 = readtable("CaseG_output.csv")

optBL!(currT,d)
newT = deepcopy(currT);
count = 0
N = 100
move = whichMove(currT)
move = :CHdir
move = :MVorigin
move = :MVtarget
move = :nni

flag = proposedTop!(move,newT,true,count,N, zeros(Int,18), zeros(Int,6))
printEdges(newT)
printNodes(newT)
count([e.hybrid for e in newT.edge]) == 2 || error("there are not 2 hybrid edges")
newT.hybrid[1].k

optBL!(newT,d)
newloglik - currloglik
changeDirectionUpdate!(newT.node[5],newT,false)

currT = deepcopy(newT);
currloglik = newloglik
currxmin = newxmin

# ------------------5taxon network 1 hybridization-----------------
# starting topology: Case F
include("../case_f_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output.csv")
d = readDataCF(df)

currloglik,currxmin = optBL!(currT,d)
updateParameters!(currT)
updateLik!(currT,currloglik)
newT = deepcopy(currT);
count = 0
N = 100
move = whichMove(currT)
move = :CHdir
move = :MVorigin
move = :MVtarget
move = :nni

flag = proposedTop!(move,newT,true,count,N)
printEdges(newT)
printNodes(newT)
count([e.hybrid for e in newT.edge]) == 2 || error("there are not 2 hybrid edges")
newT.hybrid[1].k

newloglik, newxmin = optBL!(newT,d)
newloglik - currloglik

currT = deepcopy(newT);
currloglik = newloglik
currxmin = newxmin

# -------------------5taxon tree------------------

include("../src/types.jl")
include("../src/functions.jl")

df = readtable("Tree_output.csv")
d = readDataCF(df)

# starting tree:
tree = "((6,4),(7,8),10);"
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
currT = readTopologyUpdate("prueba_tree.txt");
printEdges(currT)


optBL!(currT,d)
newT = deepcopy(currT);
count = 0
N = 10
move = :nni

flag = proposedTop!(move,newT,true,count,N, zeros(Int,18), zeros(Int,6))
flag
printEdges(newT)
printNodes(newT)

newloglik, newxmin = optBL!(newT,d)
newloglik - currloglik
`
currT = deepcopy(newT);
currloglik = newloglik
currxmin = newxmin
