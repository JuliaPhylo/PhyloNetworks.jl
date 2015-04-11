# test for moveDownLevel
# Claudia April 2015

# starting topology: Case G
include("../examples/case_g_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output.csv")
d = readTableCF(df)

optBL!(currT,d)


include("../examples/case_f_example.jl");
currT = deepcopy(net);
printEdges(currT)

# real network: Case H
df = readtable("CaseH_output.csv")
d = readTableCF(df)

optBL!(currT,d)
