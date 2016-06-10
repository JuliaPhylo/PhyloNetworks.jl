# test for whole optimization on the space of topologies
# with the same number of hybridizations
# Claudia February 2015

# -------------------5taxon tree------------------

include("../src/types.jl")
include("../src/functions.jl")

df = readtable("Tree_output.csv")
d = readTableCF(df)

# starting tree:
tree = "((6,4),(7,8),10);"
tree = "((((6,4),7),8),10);" #true tree
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
currT = readTopologyUpdate("prueba_tree.txt");
printEdges(currT)


@time optTopLevel!(currT,d,0)
net = snaq(currT,d,hmax=0);

#old:
#@time optTopLevel!(currT,M,N,d,0);
#got 5.34957 at [0.2,0.1] after 28 iterations (returned FTOL_reached)
#loglik_1 = 5.34957
#found minimizer topology at step 1 with -loglik=5.34957 and ht_min=[0.2,0.1]
#elapsed time: 8.910952599 seconds (91065584 bytes allocated, 0.59% gc time)
printEdges(net)
# forgot to copy, but true tree!



# ------------------5taxon network 1 hybridization: Case H-----------------
# starting topology: Case G
include("../examples/case_g_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output.csv")
d = readTableCF(df)


@time newT = optTopLevel!(currT,d,1)
# with original optBL
# elapsed time: 347.527745834 seconds (372704328 bytes allocated, 0.06% gc time)
# did not find right network (Case H), came back to starting point (Case G)

# with new added inequality
#elapsed time: 58.711593728 seconds (498732864 bytes allocated, 0.49% gc time)
# did not find the right network (Case H), stopped in bad diamond II case

# with afterOptBLAll
# elapsed time: 6.321085645 seconds (340867164 bytes allocated, 3.52% gc time)
# found correct network!!

# with ftol,xtol in optBL
#got 0.11429 at [0.67712,0.92599,0.13123,0.23364] after 219 iterations (returned FTOL_REACHED)
#WARNING: newT.loglik 0.11428632002947303 not really close to 0.0, you might need to redo with another starting point
#END optTopLevel: found minimizer topology at step 100 with -loglik=0.18708 and ht_min=[0.1995,0.62389,0.44779,0.07295]

# with old ftol, xtol
#got 73.2177 at [0.42152,1.0016,0.64937,0.0] after 91 iterations (returned FTOL_REACHED)
#before comparing: newT.loglik 73.21769826979121, currT.loglik 0.09699788946913745
#ends while for 100 with delta 0.7803748966954892
#WARNING: newT.loglik 0.09699788946913745 not really close to 0.0, you might need to redo with another starting point
#END optTopLevel: found minimizer topology at step 100 with -loglik=0.097 and ht_min=[0.96515,0.06804,0.0846]

tree = "(1,2,((4,#H-6::0.29974323628759736):0.0,((8)#H-6:0.2532215227803048::0.7002567637124026,7):0.1265313400455008):0.9651457506039909);"
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
newT = readTopologyUpdate("prueba_tree.txt");
printEdges(newT)


newT.ht

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

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
df = readtable("CaseH_output.csv")
d = readTableCF(df)

@time optTopLevel!(currT,d,1)

printEdges(newT)


# ------------------5taxon network 1 hybridization: Case F-----------------
# starting topology: Case G
include("../examples/case_g_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case F
df = readtable("CaseF_output.csv")
d = readTableCF(df)


@time optTopLevel!(currT,d,1)

printEdges(newT)


# starting topology: Case H
include("../examples/case_h_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case F
df = readtable("CaseF_output.csv")
d = readTableCF(df)


@time optTopLevel!(currT,d,1)

printEdges(newT)
