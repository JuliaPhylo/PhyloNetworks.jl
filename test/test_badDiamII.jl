# updateGammaz does not recognize this as bad diamond II
# Claudia April 2015

if !isdefined(:individualtest) individualtest = false; end

if(individualtest)
    include("../src/types.jl")
    include("../src/functions.jl")
end

if isdefined(:PhyloNetworks)
    PhyloNetworks.CHECKNET || error("need CHECKNET==true in PhyloNetworks to test snaq in test_correctLik.jl")
else
    CHECKNET || error("need CHECKNET==true to test snaq in test_correctLik.jl")
end

tree = "(6,(5,#H7:0.0):9.970714072991349,(3,(((2,1):0.2950382234364404,4):0.036924483697671304)#H7:0.00926495670648208):1.1071489442240392);"
net = readTopologyLevel1(tree);
checkNet(net)
printNodes(net)
printEdges(net)
net.node[10].number == 3 || error("wrong hybrid")
net.node[10].hybrid || error("does not know it is hybrid")
net.node[10].isBadDiamondII || error("does not know it is bad diamond II")
##plot(net,showEdgeNumber=true)
[e.inCycle for e in net.edge] == [-1,-1,3,3,-1,-1,-1,-1,-1,-1,3,3] || error("error in incycle")
[e.containRoot for e in net.edge] == [true,true,false,true,true,false,false,false,false,false,false,true] || error("error in contain root")
[e.istIdentifiable for e in net.edge] == [false,false,true,true,false,false,false,true,false,false,true,true] || error("istIdentifiable not correct")
(net.edge[3].hybrid && net.edge[11].hybrid) || error("hybrid edges wrong")

df=DataFrame(t1=["1","1","2","2","1","2","2","2","2","2","1","2","2","3","2"],
             t2=["3","3","3","3","3","1","5","1","1","1","5","1","1","5","3"],
             t3=["5","5","5","5","6","5","6","3","6","5","6","3","3","6","6"],
             t4=["6","4","4","6","4","6","4","5","4","4","4","6","4","4","4"],
             CF1234=[0.565,0.0005,0.0005,0.565,0.00003,0.99986,0.0410167,1,0.99987,1,0.040167,0.998667,1,0.073167,0.00003],
             CF1324=[0.0903,0.8599,0.8599,0.0903,0.8885,0.00006,0.263,0,0.00006,0,0.2630,0.00006,0,0.0424667,0.8885])
df[:CF1423] = 1-df[:CF1234]-df[:CF1324]
d = readTableCF(df)

net2 = topologyMaxQPseudolik!(net,d)
340 < net2.loglik < 340.5 || error("wrong loglik")
net2.edge[3].istIdentifiable || error("wrong hybrid is t identifiable")
9.95 < net2.edge[3].length < 9.99 || error("wrong bl estimated")
net2.edge[11].istIdentifiable || error("wrong hybrid is t identifiable")
net2.edge[11].length < 0.01 || error("wrong bl estimated")
net2.edge[10].length == 0 || error("tree edge in bad diamond II not 0")
printEdges(net2)


## testing readTopology----------------------------------------------------------------
## will remove this from the test because readTopology should not have to worry about istIdentifiable
## we move into updateGammaz
## tree = "(6,(5,#H7:0.0):9.970714072991349,(3,(((2,1):0.2950382234364404,4):0.036924483697671304)#H7:0.00926495670648208):1.1071489442240392);"
## net2 = readTopology(tree)
## printEdges(net2)

## for e in net2.edge
##     if(e.node[1].leaf || e.node[2].leaf)
##         !e.istIdentifiable || error("ext edge should not identifiable")
##     else
##         e.istIdentifiable || error("int edge should not identifiable")
##     end
## end

## ## this case works fine
## tree = "((((8,10))#H1,7),6,(4,#H1));" # Case I Bad diamond I
## net = readTopologyLevel1(tree)
## checkNet(net)
## net2 = readTopology(tree)
## printEdges(net2)

## for e in net2.edge
##     if(e.node[1].leaf || e.node[2].leaf)
##         !e.istIdentifiable || error("ext edge should not identifiable")
##     else
##         e.istIdentifiable || error("int edge should not identifiable")
##     end
## end
