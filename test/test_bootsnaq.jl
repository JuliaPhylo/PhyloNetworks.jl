# test the function to do bootstrap on snaq
# Claudia October 2015

include("../src/types.jl")
include("../src/functions.jl")

d=readTableCF("../examples/tableCF.txt");
T=readStartTop("../examples/startTree.txt",d);
writeTopology(T)
net1=snaq(T,d,filename="net1_snaq");

df=readtable("../examples/tableCFCI.txt",separator=';')
bootnet = bootsnaq(T,df,nrep=4,bestNet=net1);
bootnet = bootsnaq(T,df,nrep=4);
