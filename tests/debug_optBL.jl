# code to try to figure out problem with gammaz inequality in bad diamond I
# Claudia March 2015

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

include("../types.jl")
include("../functions.jl")

df = readtable("CaseF_output.csv")
df2 = readtable("CaseF_output2.csv") #longer branches
d = readDataCF(df)
d2 = readDataCF(df2)

# starting ht (gamma,t4,t5,t9)
ht = [0.1,1.,1.,1.]
ht = [0.3,10.,10.,1.] # crashes!

tree = string("(((6:0.1,(4)11#H1:::",string(1-ht[1]),")1:",string(ht[3]),",(11#H1:::",string(ht[1]),",7))5:",string(ht[4])",8:0.1,10:0.1);") # Case F: bad diamond I
f = open("prueba_tree.txt","w")
write(f,tree)
close(f)
net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)

net.ht
realht1 = [0.1,0.127,0.0285]
realht2 = [1.1,0.49,0.2]

@time optBL!(net,d2,true)

# ----------------- NLopt example website -----------------------
# to do here: put website example and change starting point to close to the boundary
# also, do a small example of optimization with inequality x_n+x_n-1<=1 to send to the author of nlopt

