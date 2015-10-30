# functions in readme file
# Claudia October 2015

include("../src/types.jl")
include("../src/functions.jl")

# read data
d=readTrees2CF("../examples/treefile.txt");
d2=readTrees2CF("../examples/treefile.txt",whichQ="rand",numQ=10);
d=readTableCF("../examples/tableCF.txt");
T=readStartTop("../examples/startTree.txt",d);
T2=readStartTop("../examples/startTree.txt",d2);
writeTopology(T)
# snaq
net1=snaq!(T2,d2,filename="net1_snaq")
net2=snaq!(T,d,hmax=2, filename="net2_snaq");
writeTopology(net1)
writeTopology(net1,di=true)
# topologyMaxQ
net1topo = readTopologyLevel1("(2,(4,(3,(5,(6,#H1)))),(1)#H1);");
topologyMaxQPseudolik!(net1topo,d)
writeTopology(net1topo)
# topologyQ
net1withBL = readTopologyLevel1("(2,(4,(3,(5,(6,#H6:1.0::0.288):5.006):0.518):0.491):1.533,(1)#H6:1.0::0.712);");
topologyQPseudolik!(net1withBL,d)
#plot obsCF vs expCF
df = dfObsExpCF(d)
using Gadfly
p = plot(df,layer(x="obsCF1",y="expCF1",Geom.point,Theme(default_color=colorant"orange")),layer(x="obsCF2",y="expCF2",Geom.point,Theme(default_color=colorant"purple")),layer(x="obsCF3",y="expCF3",Geom.point,Theme(default_color=color("blue"))),layer(x=0:1,y=0:1),Geom.line,Theme(default_color=color("black")))
