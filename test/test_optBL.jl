# test for the optimization of branch lengths
# with Case g
# Claudia January 2015

## include("../examples/case_g_example.jl");
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("CaseG_output.csv",df)

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

@time optBL!(net,d2)
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

# with new ftol, xtol values
#got 0.00073 at [0.2474,0.18189,0.11927,0.30014] after 48 iterations (returned XTOL_REACHED)
#elapsed time: 0.023274544 seconds (4256960 bytes allocated)
#

# with old ftol,xtol values
#got 0.0 at [0.13768,0.19568,0.10456,0.61459] after 561 iterations (returned FTOL_REACHED)
#elapsed time: 2.122653096 seconds (118537484 bytes allocated, 4.50% gc time)

@time optBL!(net,d2,false,1e-5,1e-6,1e-3,1e-4)
#got 0.0001 at [0.24708,0.18075,0.12033,0.29593] after 64 iterations (returned SUCCESS)
#elapsed time: 0.028643709 seconds (5502208 bytes allocated)


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

@time optBL!(net,d2)
#got 0.01558 at [0.09049,0.20435,0.09535,0.91173] after 31 iterations (returned XTOL_REACHED)
#elapsed time: 0.014890177 seconds (2772976 bytes allocated)

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


# (0.11,0.22,0.11,1.1)
tree = "((((6,4)1,(7)11#H1:::0.89)5,(11#H1:::0.11,8)),10);" # Case G different starting branch lengths
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


## include("../examples/case_h_example.jl");
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("CaseH_output.csv",df)

include("../src/types.jl")
include("../src/functions.jl")

df = readtable("CaseH_output.csv")
df = readtable("CaseH_output2.csv")
d2 = readTableCF(df)

# starting ht (gamma,t3,t5,t7)
ht = [0.2,1.,1.,1.]
ht = [0.05,0.0001,2.,1.] #strange case un debug18bad.txt

tree = string("((((6:0.1,4:1.5):",string(ht[2]),",#H1:::",string(ht[1]),"):",string(ht[3]),",7:0.2):",string(ht[4]),",(8)#H1:::",string(1-ht[1]),",10:0.1);") # Case H
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

net.ht
realht = [0.1,0.1,1.,0.1]

@time optBL!(net,d2)
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

# with new ftol, xtol (start 0.2,1,1,1)
#got 0.04565 at [0.20363,0.64827,0.40823,0.11374] after 51 iterations (returned XTOL_REACHED)
#elapsed time: 0.057075544 seconds (4462816 bytes allocated, 53.14% gc time)

#got 0.04478 at [0.20178,0.64205,0.41684,0.11366] after 25 iterations (returned SUCCESS)
#elapsed time: 0.016354695 seconds (2226424 bytes allocated)

# with old ftol, xol
#got 0.0 at [0.1,0.1,1.0,0.1] after 176 iterations (returned FTOL_REACHED)
#elapsed time: 0.117537385 seconds (15143800 bytes allocated, 29.95% gc time)

@time optBL!(net,d2,false,1e-5,1e-6,1e-5,1e-6)
#got 1.0e-5 at [0.10043,0.1054,0.99408,0.10017] after 115 iterations (returned FTOL_REACHED)
#elapsed time: 0.055701792 seconds (10398360 bytes allocated)

@time optBL!(net,d2,false,1e-5,1e-6,1e-3,1e-4)
#got 1.0e-5 at [0.10054,0.10558,0.99394,0.10022] after 25 iterations (returned FTOL_REACHED)
#elapsed time: 0.013228257 seconds (2226504 bytes allocated)

# from debug18bad
tree = string("(4,6,(#H2:0.7392085405544356::0.046179825120885414,(7,(10,(8)#H2:0.0::0.9538201748791146):0.9803511144374873):2.212878358589699):0.00000038687);")
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

net.ht
realht = [0.1,0.1,1.,0.1]

@time optBL!(net,d2)


# ==================================================================================================================================

# test optBL with Case J
# Claudia January 2015


## include("../examples/case_j_example.jl");
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("CaseJ_output.csv",df)

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

@time optBL!(net,d2)
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

# with ftol, xtol
#got 0.05209 at [0.26015,0.26161,0.10523,0.24363] after 42 iterations (returned XTOL_REACHED)
#elapsed time: 0.021402455 seconds (3634568 bytes allocated)

# with old ftol, xtol
#got 0.0 at [0.1,0.2,0.1,1.0] after 233 iterations (returned FTOL_REACHED)
#elapsed time: 0.112556852 seconds (19426488 bytes allocated)

@time optBL!(net,d2,false,1e-5,1e-6,1e-3,1e-4)
#got 0.00219 at [0.1303,0.20847,0.10099,0.6283] after 106 iterations (returned FTOL_REACHED)
#elapsed time: 0.060469424 seconds (8887448 bytes allocated)

@time optBL!(net,d2,false,1e-5,1e-6,1e-5,1e-6)
#got 0.00177 at [0.13073,0.20994,0.10121,0.64478] after 26 iterations (returned FTOL_REACHED)
#elapsed time: 0.012715507 seconds (2277912 bytes allocated)

# ==================================================================================================================================

# test optBL with Case F Bad Diamond I
# Claudia January 2015


## include("../examples/case_f_example.jl");
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

@time optBL!(net,d2)
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

# with ftol,xtol
#got 0.25041 at [0.06362,0.15339,0.04918] after 24 iterations (returned XTOL_REACHED)
#elapsed time: 0.045551191 seconds (2403576 bytes allocated)

#with old ftol, xtol
#got 0.0 at [0.1,0.12689,0.02855] after 2917 iterations (returned FTOL_REACHED)
#elapsed time: 1.303746554 seconds (229545768 bytes allocated, 11.98% gc time)

@time optBL!(net,d2,false,1e-5,1e-6,1e-3,1e-4)
#got 5.0e-5 at [0.09974,0.12735,0.02892] after 92 iterations (returned XTOL_REACHED)
#elapsed time: 0.079739994 seconds (7321368 bytes allocated, 45.28% gc time)

# ==================================================================================================================================

# test optBL with Case I Bad Diamond II
# Claudia January 2015


## include("../examples/case_i_example.jl");
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("CaseI_output.csv",df)

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

@time optBL!(net,d2)
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

# with ftol,xtol
# ht = [0.2,0.0,2.0,2.0,2.0]
#got 0.00256 at [0.10406,2.0465,1.00637,3.33775,0.98475] after 58 iterations (returned XTOL_REACHED)
#elapsed time: 0.034984668 seconds (5205128 bytes allocated)

# with old ftol, xtol
#got 0.00034 at [0.1034,2.03863,1.00652,3.5823,0.94937] after 201 iterations (returned FTOL_REACHED)
#elapsed time: 0.094041162 seconds (17632400 bytes allocated)

@time optBL!(net,d2,false,1e-5,1e-6,1e-3,1e-4)
# ht = [0.2,0.0,2.0,2.0,2.0]
#got 0.00043 at [0.10287,2.02652,1.00596,5.93267,0.96293] after 51 iterations (returned SUCCESS)
#elapsed time: 0.027747237 seconds (4558192 bytes allocated)

@time optBL!(net,d2,false,1e-5,1e-6,1e-5,1e-6)
#got 0.00036 at [0.10311,2.03708,1.00602,3.43898,0.95781] after 35 iterations (returned FTOL_REACHED)
#elapsed time: 0.018274282 seconds (3167744 bytes allocated)

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


include("../src/types.jl")
include("../src/functions.jl")

df = readtable("Tree_output.csv")
d2 = readTableCF(df)

# starting tree:
ht = [1.0,1.0]
tree = string("(((6:0.1,4:1.5)1:",string(ht[1]),",7:0.2)5:",string(ht[2]),",8:0.1,10:0.1);") # normal tree
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");


net.ht
realht = [0.2,0.1]

@time optBL!(net,d2)
#got 5.34957 at [0.2,0.1] after 28 iterations (returned FTOL_REACHED)
#elapsed time: 5.533522804 seconds (54162200 bytes allocated, 0.41% gc time)
#(5.349567420518451,[0.2,0.1])

#with ftol, xtol
#got 0.0 at [0.19999,0.09999] after 20 iterations (returned XTOL_REACHED)
#elapsed time: 0.007440742 seconds (1234840 bytes allocated)

#==============================================================================================
#================ Debugging optBL ============================================================
# compare expCF from wrong estimates with real estimates
# ===========================================================================================

# test optBL with Case I Bad Diamond II
# does not yield correct ht for one starting point

include("../src/types.jl")
include("../src/functions.jl")

df = readtable("CaseI_output.csv")
d2 = readTableCF(df)

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
include("../src/types.jl")
include("../src/functions.jl")

df = readtable("CaseG_output.csv")
d2 = readTableCF(df)

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
