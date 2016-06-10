# test for whole optimization on the space of topologies
# for h<= hmax
# Claudia March 2015
# based on test_optTopLevel.jl, but now we use a more identifiable Case G, H to begin with
# and now we use afterOptBLALL

# -------------------5taxon tree------------------

include("../src/types.jl")
include("../src/functions.jl")

df = readtable("Tree_output.csv")
d = readTableCF(df)

# starting tree:
tree = "((6,4),(7,8),10);"
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
currT = readTopologyUpdate("prueba_tree.txt");
printEdges(currT)

epsilon = eps()
N = 100

@time newT = optTopLevel!(currT,epsilon,d,0);
#got 5.34957 at [0.2,0.1] after 28 iterations (returned FTOL_reached)
#loglik_1 = 5.34957
#found minimizer topology at step 1 with -loglik=5.34957 and ht_min=[0.2,0.1]
#elapsed time: 8.910952599 seconds (91065584 bytes allocated, 0.59% gc time)
printEdges(newT)
# forgot to copy, but true tree!



# ------------------5taxon network 1 hybridization: Case H-----------------
# starting topology: Case G

## include("../examples/case_h_example2.jl");
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("CaseH_output2.csv",df)

include("../examples/case_g_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output2.csv")
d = readTableCF(df)

@time optTopLevel!(currT,d,1)


tree = string("(1,2,((7,(8)#H5:5.814544267883624):0.9977876663423212,(#H5:1.0,4):0.0):1.9430580774498776);")
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
newT = readTopologyUpdate("prueba_tree.txt");

printEdges(newT)

newT.ht

q1 = Quartet(1,["1","4","2","7"],[0.5,0.4,0.1]);
q2 = Quartet(2,["1","4","8","7"],[0.5,0.4,0.1]);
q3 = Quartet(3,["8","4","2","7"],[0.5,0.4,0.1]);
q4 = Quartet(4,["1","8","2","7"],[0.5,0.4,0.1]);
q5 = Quartet(5,["1","4","2","8"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(newT,d)

wrongdf = writeExpCF(d.quartet)
writetable("CaseH_output_wrong_optTop_startCaseG.csv",wrongdf)

# ------------------
# starting topology: Case F
include("../examples/case_f_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output2.csv")
d = readTableCF(df)

@time optTopLevel!(currT,d,1)

tree = string("(4,6,(#H2:0.7392085405544356::0.046179825120885414,(7,(10,(8)#H2:0.0::0.9538201748791146):0.9803511144374873):2.212878358589699):0.00000073506);")
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
newT = readTopologyUpdate("prueba_tree.txt");

printEdges(newT)

#include("../examples/case_h_example2.jl");
newT.ht

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

dd = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(newT,dd)

wrongdf = writeExpCF(d.quartet)
writetable("CaseH2_output_wrong_optTop_startCaseF.csv",wrongdf)

#elapsed time: 4.236802192 seconds (274215268 bytes allocated, 4.26% gc time)
#WARNING: newT.loglik 1.4887847814495139 not really close to 0.0, you might need to redo with another starting point


# ------------------5taxon network 1 hybridization: Case F-----------------
# starting topology: Case G
include("../examples/case_g_example.jl");
currT = deepcopy(net);
printEdges(currT)

## include("../examples/case_f_example2.jl");
## parameters!(net)
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("CaseF_output2.csv",df)

# real network: Case F
df = readtable("CaseF_output.csv")
df = readtable("CaseF_output2.csv")
d = readTableCF(df)


@time optTopLevel!(currT,d,1)

printEdges(newT)


# starting topology: Case H
include("../examples/case_h_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case F
df = readtable("CaseF_output.csv")
df = readtable("CaseF_output2.csv")
d = readTableCF(df)


@time optTopLevel!(currT,d,1)

printEdges(newT)
