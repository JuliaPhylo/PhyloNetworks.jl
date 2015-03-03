# test for the optimization of branch lengths
# with Case g
# Claudia January 2015

## include("../case_g_example.jl");
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("CaseG_output.csv",df)

include("../types.jl")
include("../functions.jl")

df = readtable("CaseG_output.csv")
d2 = readDataCF(df)

tree = "((((6,4)1,(7)11#H1:::0.8)5,(11#H1:::0.2,8)),10);" # Case G different starting branch lengths
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

net.ht
realht = [0.1,0.2,0.1,1.0]

@time fmin,xmin=optBL!(net,d2)
# with usual logPseudoLik:
#got 5.33118 at [0.247,0.18076,0.12026,0.2956] after 81 iterations (returned FTOL_REACHED)
#elapsed time: 8.654081812 seconds (8624832 bytes allocated)
#(5.331178696104555,[0.246999,0.180762,0.120265,0.295604])

# with new logPseudoLik:
#got 0.0 at [0.14088,0.19534,0.10493,0.59544] after 512 iterations (returned FTOL_REACHED)
#elapsed time: 6.552122524 seconds (95796384 bytes allocated, 0.56% gc time)
#(5.522931422394296e-8,[0.140883,0.195339,0.10493,0.595437])

# with 100*logPseudoLik and t in (0,10): takes longer, but finds it!
#got 0.0 at [0.10085,0.19991,0.1001,0.98557] after 1123 iterations (returned FTOL_REACHED)
#elapsed time: 53.20413795 seconds (107820188 bytes allocated, 0.13% gc time)
#(2.169777681982341e-9,[0.100853,0.199907,0.100098,0.985569])

@allocated fmin,xmin=optBL!(net,d2)

# -------- different starting point ------

# (0.12,1.,1.,1.)
tree = "((((6,4)1,(7)11#H1:::0.88)5,(11#H1:::0.12,8)),10);" # Case G different starting branch lengths
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

@time fmin,xmin=optBL!(net,d2)


# (0.1,1.,1.,1.)
tree = "((((6,4)1,(7)11#H1)5,(11#H1,8)),10);" # Case G different starting branch lengths
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

@time fmin,xmin=optBL!(net,d2)


# (0.05,1.,1.,1.)
tree = "((((6,4)1,(7)11#H1:::0.95)5,(11#H1:::0.05,8)),10);" # Case G different starting branch lengths
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

@time fmin,xmin=optBL!(net,d2)


# (0.15,1.,1.,1.)
tree = "((((6,4)1,(7)11#H1:::0.85)5,(11#H1:::0.15,8)),10);" # Case G different starting branch lengths
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

@time fmin,xmin=optBL!(net,d2)


# (0.35,1.,1.,1.)
tree = "((((6,4)1,(7)11#H1:::0.65)5,(11#H1:::0.35,8)),10);" # Case G different starting branch lengths
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

@time fmin,xmin=optBL!(net,d2)


# (0.24,1.,1.,1.)
tree = "((((6,4)1,(7)11#H1:::0.75)5,(11#H1:::0.25,8)),10);" # Case G different starting branch lengths
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)
printNodes(net)

@time fmin,xmin=optBL!(net,d2)

# ==================================================================================================

# test optBL with Case H
# Claudia January 2015


## include("../case_h_example.jl");
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("CaseH_output.csv",df)

include("../types.jl")
include("../functions.jl")

df = readtable("CaseH_output.csv")
d2 = readDataCF(df)

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
# with usual logPseudoLik:
#got 4.34098 at [0.1,0.10001,0.99999,0.1] after 151 iterations (returned FTOL_REACHED)
#elapsed time: 21.076511672 seconds (21355496 bytes allocated, 0.17% gc time)
#(4.340977085797292,[0.0999962,0.100005,0.999989,0.0999984])

# with new logPseudoLik:
#got 0.0 at [0.1,0.1,1.0,0.1] after 192 iterations (returned FTOL_REACHED)
#elapsed time: 0.973165584 seconds (16698980 bytes allocated, 3.44% gc time)
#(4.3376522534621724e-13,[0.1,0.0999979,1.0,0.100001])

# with 100*logPseudoLik and t in (0,10):
#got 0.0 at [0.1,0.1,1.0,0.1] after 208 iterations (returned FTOL_REACHED)
#elapsed time: 6.892366587 seconds (20646048 bytes allocated)
#(2.4010737722561867e-13,[0.0999999,0.0999992,1.0,0.1])

# ==================================================================================================================================

# test optBL with Case J
# Claudia January 2015


## include("../case_j_example.jl");
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("CaseJ_output.csv",df)

include("../types.jl")
include("../functions.jl")

df = readtable("CaseJ_output.csv")
d2 = readDataCF(df)

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
# with usual logPseudoLik:
#got 5.39161 at [0.09997,0.19999,0.1,1.00061] after 246 iterations (returned FTOL_REACHED)
#elapsed time: 59.297791699 seconds (71939204 bytes allocated, 0.11% gc time)
#(5.3916134028946034,[0.0999661,0.199991,0.0999991,1.00061])

# with new logPseudoLik:
#got 0.0 at [0.10002,0.20001,0.1,0.99958] after 260 iterations (returned FTOL_REACHED)
#elapsed time: 1.211645168 seconds (21966196 bytes allocated)
#(9.797186804450541e-12,[0.100022,0.200007,0.100001,0.999579])

# with 100*logPseudoLik and t in (0,10):
#got 0.0 at [0.1,0.2,0.1,1.00003] after 249 iterations (returned FTOL_REACHED)
#elapsed time: 6.656956697 seconds (23010204 bytes allocated, 0.56% gc time)
#(5.280560068867562e-12,[0.0999983,0.199999,0.0999999,1.00003])

# ==================================================================================================================================

# test optBL with Case F Bad Diamond I
# Claudia January 2015


## include("../case_f_example.jl");
## parameters!(net)
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("CaseF_output.csv",df)

include("../types.jl")
include("../functions.jl")

df = readtable("CaseF_output.csv")
d2 = readDataCF(df)

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
# with usual logPseudoLik:
#got 5.40235 at [0.1,0.12689,0.02855] after 116 iterations (returned FTOL_REACHED)
#elapsed time: 15.648447251 seconds (17823396 bytes allocated, 0.18% gc time)
#(5.402353356033268,[0.1,0.126887,0.0285486])

# with new logPseudoLik:
#got 0.0 at [0.1,0.12689,0.02855] after 116 iterations (returned FTOL_REACHED)
#elapsed time: 0.494250509 seconds (9272420 bytes allocated)
#(3.6216219428084243e-12,[0.1,0.126887,0.0285488])

# with new logPseudoLik and new algorithm LN_COBYLA
#got 0.0 at [0.1,0.12689,0.02855] after 3277 iterations (returned FTOL_REACHED)
#elapsed time: 23.001580343 seconds (312597904 bytes allocated, 0.98% gc time)
#(1.920897287842076e-13,[0.0999998,0.126889,0.0285491])

# with 100*logPseudoLik and t in (0,10):
# got 0.0 at [0.1,0.12689,0.02855] after 3356 iterations (returned FTOL_REACHED)
# elapsed time: 90.405727385 seconds (281649264 bytes allocated, 0.24% gc time)
# (8.157692309723136e-12,[0.0999999,0.126889,0.028549])

# ==================================================================================================================================

# test optBL with Case I Bad Diamond II
# Claudia January 2015


## include("../case_i_example.jl");
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("CaseI_output.csv",df)

include("../types.jl")
include("../functions.jl")

df = readtable("CaseI_output.csv")
d2 = readDataCF(df)

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
## with usual logPseudoLik:
# with ht = [0.2,0.0,0.5,0.5,0.5]
#got 3.69564 at [0.09753,1.97129,0.99532,0.46822,1.04214] after 267 iterations (returned FTOL_REACHED)
#elapsed time: 45.70869142 seconds (35023452 bytes allocated, 0.06% gc time)
#(3.6956398763931144,[0.0975269,1.97129,0.995318,0.468224,1.04214])

#ht: [0.2,1.0,2.0,2.0,2.0]
#got 3.69564 at [0.1034,2.03878,1.00653,3.58586,0.94929] after 176 iterations (returned FTOL_REACHED)
#elapsed time: 15.743171805 seconds (17639780 bytes allocated)
#(3.6956415275081724,[0.103403,2.03878,1.00653,3.58586,0.949288])


## with new logPseudoLik:
#ht: [0.2,1.0,2.0,2.0,2.0]
#got 0.0 at [0.10342,2.03906,1.00656,3.60666,0.94899] after 195 iterations (returned FTOL_REACHED)
#elapsed time: 1.11298624 seconds (17529096 bytes allocated, 2.98% gc time)
#(3.3635155002465492e-6,[0.103419,2.03906,1.00656,3.60666,0.948988])

# with ht = [0.2,0.0,0.5,0.5,0.5]
#got 0.0 at [0.09997,1.99969,0.99995,0.99284,1.00042] after 424 iterations (returned FTOL_REACHED)
#elapsed time: 1.733474741 seconds (36994324 bytes allocated, 1.90% gc time)
#(1.8811458391292683e-10,[0.0999737,1.99969,0.99995,0.992835,1.00042])

# with 100*logPseudoLik and t in (0,10):
# ht = [0.2,0.0,2.0,2.0,2.0]
#got 0.0 at [0.1,2.0,1.0,0.99994,1.0] after 753 iterations (returned FTOL_REACHED)
#elapsed time: 47.213517434 seconds (80265788 bytes allocated, 0.16% gc time)
#(1.2568189823957448e-12,[0.0999998,2.0,1.0,0.999941,1.0])

# -------------------5taxon tree------------------
## include("../tree_example_read.jl");
## printEdges(net)
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("Tree_output.csv",df)


include("../types.jl")
include("../functions.jl")

df = readtable("Tree_output.csv")
d2 = readDataCF(df)

# starting tree:
ht = [1.0,1.0]
tree = string("(((6:0.1,4:1.5)1:",string(ht[1]),",7:0.2)5:",string(ht[2]),",8:0.1,10:0.1);") # normal tree
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");


net.ht
realht = [0.2,0.1]

@time fmin,xmin=optBL!(net,d2)
#got 5.34957 at [0.2,0.1] after 28 iterations (returned FTOL_REACHED)
#elapsed time: 5.533522804 seconds (54162200 bytes allocated, 0.41% gc time)
#(5.349567420518451,[0.2,0.1])

#==============================================================================================
#================ Debugging optBL ============================================================


# test optBL with Case I Bad Diamond II
# does not yield correct ht for one starting point

include("../types.jl")
include("../functions.jl")

df = readtable("CaseI_output.csv")
d2 = readDataCF(df)

# starting ht (gamma,t4,t6,t9,t10)
wronght = [0.1,1.0,1.0,3.6,1.0]
ht = wronght
realht = [0.1,2.0,1.0,1.0,1.0]

tree = string("((((8,10):",string(ht[2]),")#H1:::",string(1-ht[1]),",7):",string(ht[3]),",6,(4,#H1:",string(ht[4]),"::",string(ht[1]),"):",string(ht[5]),");") # Case I Bad diamond II
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

net.ht

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(net,d)

wrongdf = writeExpCF(d.quartet)
writetable("CaseI_output_wrong.csv",wrongdf)


# test optBL with Case G
# does not yield correct ht for one starting point
include("../types.jl")
include("../functions.jl")

df = readtable("CaseG_output.csv")
d2 = readDataCF(df)

# starting ht (gamma,t3,t6,t9)
wronght = [0.14,0.2,0.1,0.6]
ht = wronght
realht = [0.1,0.2,0.1,1.0]

tree = string("((((6,4)1:",string(ht[2]),",(7)11#H1:::",string(1-ht[1]),")5:",string(ht[3]),",(11#H1:::",string(ht[1]),",8):",string(ht[4]),"),10);") # Case G different starting branch lengths
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net

net.ht

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(net,d)

wrongdf = writeExpCF(d.quartet)
writetable("CaseG _output_wrong.csv",wrongdf)
