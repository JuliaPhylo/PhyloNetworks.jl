# test for whole optimization on the space of topologies
# for h<= hmax
# Claudia March 2015
# based on test_optTopLevel.jl, but now we use a more identifiable Case G, H to begin with
# and now we use afterOptBLALL

# -------------------5taxon tree------------------

include("../types.jl")
include("../functions.jl")

df = readtable("Tree_output.csv")
d = readDataCF(df)

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

## include("../case_h_example2.jl");
## q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
## q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
## q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
## q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
## q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

## d = DataCF([q1,q2,q3,q4,q5]);
## extractQuartet!(net,d)

## df = writeExpCF(d.quartet)
## writetable("CaseH_output2.csv",df)

include("../case_g_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output2.csv")
d = readDataCF(df)

epsilon = eps()
N = 100

@time newT = optTopLevel!(currT,epsilon,N,d,1)
printEdges(newT)
# with afterOptBLAll
#elapsed time: 3.727836428 seconds (199955596 bytes allocated, 3.19% gc time)
#WARNING: newT.loglik 1.533449204376918 not really close to 0.0, you might need to redo with another starting point

newT.ht

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(newT,d)

wrongdf = writeExpCF(d.quartet)
writetable("CaseH_output_wrong_optTop_startCaseG.csv",wrongdf)

# ------------------
# starting topology: Case F
include("../case_f_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output2.csv")
d = readDataCF(df)


epsilon = eps()
N = 100

@time newT = optTopLevel!(currT,epsilon,d,1);
printEdges(newT)
#elapsed time: 4.236802192 seconds (274215268 bytes allocated, 4.26% gc time)
#WARNING: newT.loglik 1.4887847814495139 not really close to 0.0, you might need to redo with another starting point


# ------------------5taxon network 1 hybridization: Case F-----------------
# starting topology: Case G
include("../case_g_example.jl");
currT = deepcopy(net);
printEdges(currT)

## include("../case_f_example2.jl");
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
df = readtable("CaseF_output2.csv")
d = readDataCF(df)


epsilon = eps()
N = 100

@time newT = optTopLevel!(currT,epsilon,N,d);
#elapsed time: 11.59461217 seconds (112302760 bytes allocated, 0.26% gc time)
printEdges(newT)


# starting topology: Case H
include("../case_h_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case F
df = readtable("CaseF_output.csv")
d = readDataCF(df)


epsilon = eps()
N = 100

@time newT = optTopLevel!(currT,epsilon,N,d);
#elapsed time: 26.763174225 seconds (142165976 bytes allocated, 0.25% gc time)
printEdges(newT)
