# test afterOptBL
# Claudia March 2015
# use cases in test_optBL and test_optTopLevel

# test_optBL -------------------------------------------------------------
# should not do anything, should leave net unchanged

# CASE G
include("../src/types.jl")
include("../src/functions.jl")

df = readtable("CaseG_output.csv")
d2 = readTableCF(df)

tree = "((((6,4)1,(7)11#H1:::0.8)5,(11#H1:::0.2,8)),10);" # Case G different starting branch lengths
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

net.ht
realht = [0.1,0.2,0.1,1.0]

@time fmin,xmin=optBL!(net,d2)
# with 100*logPseudoLik and t in (0,10): takes longer, but finds it!
#got 0.0 at [0.10085,0.19991,0.1001,0.98557] after 1123 iterations (returned FTOL_REACHED)
#elapsed time: 53.20413795 seconds (107820188 bytes allocated, 0.13% gc time)
#(2.169777681982341e-9,[0.100853,0.199907,0.100098,0.985569])

f=afterOptBL!(net,d2)
all(f) || error("afterOptBL not correctly for optBL in CASE G")

#--------
# Case H
include("../src/types.jl")
include("../src/functions.jl")

df = readtable("CaseH_output.csv")
d2 = readTableCF(df)

# starting ht (gamma,t3,t5,t7)
ht = [0.2,1.,1.,1.]

tree = string("((((6:0.1,4:1.5):",string(ht[2]),",#H1:::",string(ht[1]),"),7:0.2):",string(ht[4]),",(8)#H1:::",string(1-ht[1]),",10:0.1);") # Case H
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

net.ht
realht = [0.1,0.1,1.,0.1]

@time fmin,xmin=optBL!(net,d2)
# with 100*logPseudoLik and t in (0,10):
#got 0.0 at [0.1,0.1,1.0,0.1] after 208 iterations (returned FTOL_REACHED)
#elapsed time: 6.892366587 seconds (20646048 bytes allocated)
#(2.4010737722561867e-13,[0.0999999,0.0999992,1.0,0.1])

f=afterOptBL!(net,d2)
all(f) || error("afterOptBL not correctly for optBL in CASE H")

# -----------
# CASE J
include("../src/types.jl")
include("../src/functions.jl")

df = readtable("CaseJ_output.csv")
d2 = readTableCF(df)

# starting ht (gamma,t3,t5,t7)
ht = [0.2,1.,1.,1.]

tree = string("((((6)#H1:::",string(1-ht[1]),",4:1.5):",string(ht[2]),",7:0.2):",string(ht[3]),",8:0.1,(#H1:::",string(ht[1]),",10));") # Case J
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

net.ht
realht = [0.1,0.2,0.1,1.0]

@time fmin,xmin=optBL!(net,d2)
# with 100*logPseudoLik and t in (0,10):
#got 0.0 at [0.1,0.2,0.1,1.00003] after 249 iterations (returned FTOL_REACHED)
#elapsed time: 6.656956697 seconds (23010204 bytes allocated, 0.56% gc time)
#(5.280560068867562e-12,[0.0999983,0.199999,0.0999999,1.00003])

f=afterOptBL!(net,d2)
all(f) || error("afterOptBL not correctly for optBL in CASE J")

# -----------
# CASE F
include("../src/types.jl")
include("../src/functions.jl")

df = readtable("CaseF_output.csv")
d2 = readTableCF(df)

# starting ht (gamma,t4,t5,t9)
ht = [0.1,1.,1.,1.]

tree = string("(((6:0.1,(4)11#H1:::",string(1-ht[1]),")1:",string(ht[3]),",(11#H1:::",string(ht[1]),",7))5:",string(ht[4])",8:0.1,10:0.1);") # Case F: bad diamond I
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

net.ht
realht = [0.1,0.127,0.0285]

@time fmin,xmin=optBL!(net,d2)
# with 100*logPseudoLik and t in (0,10):
# got 0.0 at [0.1,0.12689,0.02855] after 3356 iterations (returned FTOL_REACHED)
# elapsed time: 90.405727385 seconds (281649264 bytes allocated, 0.24% gc time)
# (8.157692309723136e-12,[0.0999999,0.126889,0.028549])

f=afterOptBL!(net,d2)
all(f) || error("afterOptBL not correctly for optBL in CASE F")

# ---------------
# CASE I
include("../src/types.jl")
include("../src/functions.jl")

df = readtable("CaseI_output.csv")
d2 = readTableCF(df)

# starting ht (gamma,t4,t6,t9,t10)
ht = [0.2,0.0,2.0,2.0,2.0]
ht = [0.2,0.0,0.5,0.5,0.5]

tree = string("((((8,10):",string(ht[2]),")#H1:::",string(1-ht[1]),",7):",string(ht[3]),",6,(4,#H1:",string(ht[4]),"::",string(ht[1]),"):",string(ht[5]),");") # Case I Bad diamond II
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

net.ht
realht = [0.1,2.0,1.0,1.0,1.0]

@time fmin,xmin=optBL!(net,d2)
# with 100*logPseudoLik and t in (0,10):
# ht = [0.2,0.0,2.0,2.0,2.0]
#got 0.0 at [0.1,2.0,1.0,0.99994,1.0] after 753 iterations (returned FTOL_REACHED)
#elapsed time: 47.213517434 seconds (80265788 bytes allocated, 0.16% gc time)
#(1.2568189823957448e-12,[0.0999998,2.0,1.0,0.999941,1.0])

f=afterOptBL!(net,d2)
all(f) || error("afterOptBL not correctly for optBL in CASE I")


# test_optTopLevelparts.jl------------------------------------------------------------------

# starting topology: Case G
include("../examples/case_g_example.jl");
currT = deepcopy(net);
printEdges(currT)
printNodes(currT)

# real network: Case H
df = readtable("CaseH_output.csv")
d = readTableCF(df)

# real network: Case H
df = readtable("CaseH_output.csv")
d = readTableCF(df)

currloglik,currxmin = optBL!(currT,d, false)
isValid(currT)
success,flagh,flagt,flaghz = afterOptBL!(currT,d)
reject = afterOptBLAll!(currT,d)


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

newloglik, newxmin = optBL!(newT,d)
reject = afterOptBLAll!(newT,d)
isValid(newT)
afterOptBL!(newT,d)
currT=deepcopy(newT);


# starting topology: Case F
include("../examples/case_f_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output.csv")
d = readTableCF(df)

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

