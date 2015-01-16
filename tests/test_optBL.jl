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

fmin,xmin=optBL(net,d2)
