# functions in readme file
# Claudia October 2015

##using PhyloNetworks
using DataFrames

PhyloNetworks.CHECKNET || error("need CHECKNET==true in PhyloNetworks to test snaq in test_correctLik.jl")

try
    d=readTrees2CF("../examples/treefile.txt", CFfile = "none", writeTab=false, writeSummary=false)
    T=readTopologyLevel1("../examples/startTree.txt");
    writeTopologyLevel1(T)
    # snaq
    net1=snaq!(T,d,filename="", runs=2)
    otherNet = undirectedOtherNetworks(net1)
    net2=snaq!(T,d,hmax=2, filename="", runs=2)
    d=readTableCF("../examples/tableCF.txt", summaryfile="");
    # snaq
    net1=snaq!(T,d,filename="", runs=2)
    net2=snaq!(T,d,hmax=2, filename="", runs=2)
    rootatnode!(net1, "4")
    writeTopologyLevel1(net1)
    writeTopologyLevel1(net1,di=true)
    # topologyMaxQ
    net1topo = readTopology("(2,(4,(3,(5,(6,#H1)))),(1)#H1);")
    net1par = topologyMaxQPseudolik!(net1topo,d)
    # topologyQ
    net1withBL = readTopology("(2,(4,(3,(5,(6,#H6:1.0::0.3):5.006):0.518):0.491):1.533,(1)#H6:1.0::0.7);");
    topologyQPseudolik!(net1withBL,d)
catch
    error("error in README with input trees")
end

try
    T=readTopologyLevel1("../examples/startTree.txt");
    d=readTableCF("../examples/tableCF.txt", summaryfile="")
    # snaq
    net1=snaq!(T,d,filename="", runs=2)
    net2=snaq!(T,d,hmax=2, filename="", runs=2)
    rootatnode!(net1, "4")
    writeTopologyLevel1(net1)
    writeTopologyLevel1(net1,di=true)
    # topologyMaxQ
    net1topo = readTopology("(2,(4,(3,(5,(6,#H1)))),(1)#H1);")
    net1par = topologyMaxQPseudolik!(net1topo,d)
    # topologyQ
    net1withBL = readTopology("(2,(4,(3,(5,(6,#H6:1.0::0.3):5.006):0.518):0.491):1.533,(1)#H6:1.0::0.7);");
    topologyQPseudolik!(net1withBL,d)
catch
    error("error in README with CF table")
end

println("NO ERRORS!")

if(false)
try
    df = readtable("../examples/tableCFCI.csv")
    bootnet = bootsnaq(T, df, hmax=1, nrep=3, runs=2) # need to make sure does not produce files
    BStable, tree1 = treeEdgesBootstrap(bootnet,net1)
    BSn, BSe, BSc, BSgam, BSedgenum = hybridBootstrapSupport(bootnet, net1)
catch
    error("error in README for bootstrap")
end
end

#plot obsCF vs expCF
## df = fittedQuartetCF(d)
## using Gadfly
## p = plot(df,layer(x="obsCF1",y="expCF1",Geom.point,Theme(default_color=colorant"orange")),layer(x="obsCF2",y="expCF2",Geom.point,Theme(default_color=colorant"purple")),layer(x="obsCF3",y="expCF3",Geom.point,Theme(default_color=color("blue"))),layer(x=0:1,y=0:1),Geom.line,Theme(default_color=color("black")))
