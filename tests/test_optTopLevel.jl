# test for whole optimization on the space of topologies
# with the same number of hybridizations
# Claudia February 2015

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
include("../case_g_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output.csv")
d = readDataCF(df)


epsilon = eps()
N = 100

@time optTopLevel!(currT,epsilon,N,d,1)
printEdges(newT)
# with original optBL
# elapsed time: 347.527745834 seconds (372704328 bytes allocated, 0.06% gc time)
# did not find right network (Case H), came back to starting point (Case G)

# with new added inequality
#elapsed time: 58.711593728 seconds (498732864 bytes allocated, 0.49% gc time)
# did not find the right network (Case H), stopped in bad diamond II case

# with afterOptBLAll
# elapsed time: 3.936120092 seconds (173514144 bytes allocated, 3.14% gc time)

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
df = readtable("CaseH_output.csv")
d = readDataCF(df)


epsilon = eps()
N = 100

@time newT = optTopLevel!(currT,epsilon,N,d);
printEdges(newT)
#elapsed time: 10.644039127 seconds (127552964 bytes allocated, 0.89% gc time)

# ------------------5taxon network 1 hybridization: Case F-----------------
# starting topology: Case G
include("../case_g_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case F
df = readtable("CaseF_output.csv")
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
