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
#got 5.33118 at [0.247,0.18076,0.12026,0.2956] after 81 iterations (returned FTOL_REACHED)
#elapsed time: 8.654081812 seconds (8624832 bytes allocated)
#(5.331178696104555,[0.246999,0.180762,0.120265,0.295604])

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
#got 4.34098 at [0.1,0.10001,0.99999,0.1] after 151 iterations (returned FTOL_REACHED)
#elapsed time: 21.076511672 seconds (21355496 bytes allocated, 0.17% gc time)
#(4.340977085797292,[0.0999962,0.100005,0.999989,0.0999984])

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
#got 5.39161 at [0.09997,0.19999,0.1,1.00061] after 246 iterations (returned FTOL_REACHED)
#elapsed time: 59.297791699 seconds (71939204 bytes allocated, 0.11% gc time)
#(5.3916134028946034,[0.0999661,0.199991,0.0999991,1.00061])


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
#got 5.40235 at [0.1,0.12689,0.02855] after 116 iterations (returned FTOL_REACHED)
#elapsed time: 15.648447251 seconds (17823396 bytes allocated, 0.18% gc time)
#(5.402353356033268,[0.1,0.126887,0.0285486])

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
# with ht = [0.2,0.0,0.5,0.5,0.5]
#got 3.69564 at [0.09753,1.97129,0.99532,0.46822,1.04214] after 267 iterations (returned FTOL_REACHED)
#elapsed time: 45.70869142 seconds (35023452 bytes allocated, 0.06% gc time)
#(3.6956398763931144,[0.0975269,1.97129,0.995318,0.468224,1.04214])

#ht: [0.2,1.0,2.0,2.0,2.0]
#got 3.69564 at [0.1034,2.03878,1.00653,3.58586,0.94929] after 176 iterations (returned FTOL_REACHED)
#elapsed time: 15.743171805 seconds (17639780 bytes allocated)
#(3.6956415275081724,[0.103403,2.03878,1.00653,3.58586,0.949288])


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
