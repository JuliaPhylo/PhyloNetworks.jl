# functions in readme file
# Claudia October 2015

include("../src/types.jl")
include("../src/functions.jl")

# read data
try
    d=readTrees2CF("../examples/treefile.txt", CFfile = "ex.txt")
    d2=readTrees2CF("../examples/treefile.txt",whichQ="rand",numQ=10);
    d=readTableCF("../examples/tableCF.txt");
    T=readTopology("../examples/startTree.txt")
    T=readTopologyLevel1("../examples/startTree.txt");
    T2=readTopology("../examples/startTree.txt");
    T3=readTopology("(2,3,(4,(5,(1,6))));")
    writeTopologyLevel1(T)
    # snaq
    net1=snaq!(T,d,filename="net1_snaq")
    net2=snaq!(T,d,hmax=2, filename="net2_snaq");
    writeTopologyLevel1(net1)
    writeTopologyLevel1(net1,di=true)
    # topologyMaxQ
    net1topo = readTopologyLevel1("(2,(4,(3,(5,(6,#H1)))),(1)#H1);");
    topologyMaxQPseudolik!(net1topo,d)
    writeTopologyLevel1(net1topo)
    # topologyQ
    net1withBL = readTopologyLevel1("(2,(4,(3,(5,(6,#H6:1.0::0.288):5.006):0.518):0.491):1.533,(1)#H6:1.0::0.712);");
    topologyQPseudolik!(net1withBL,d)
catch
    error("error in README")
end

#plot obsCF vs expCF
## df = dfObsExpCF(d)
## using Gadfly
## p = plot(df,layer(x="obsCF1",y="expCF1",Geom.point,Theme(default_color=colorant"orange")),layer(x="obsCF2",y="expCF2",Geom.point,Theme(default_color=colorant"purple")),layer(x="obsCF3",y="expCF3",Geom.point,Theme(default_color=color("blue"))),layer(x=0:1,y=0:1),Geom.line,Theme(default_color=color("black")))
