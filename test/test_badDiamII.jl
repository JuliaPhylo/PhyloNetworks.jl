# updateGammaz does not recognize this as bad diamond II
# Claudia April 2015

PhyloNetworks.CHECKNET || error("need CHECKNET==true in PhyloNetworks to test snaq in test_correctLik.jl")

@testset "test: bad diamond, max pseudo lik" begin
global tree, net, df, d

tree = "(6,(5,#H7:0.0):9.970714072991349,(3,(((2,1):0.2950382234364404,4):0.036924483697671304)#H7:0.00926495670648208):1.1071489442240392);"
net = readTopologyLevel1(tree);
checkNet(net)
#printNodes(net)
#printEdges(net)
@test net.node[10].number == 3 # or: wrong hybrid
@test net.node[10].hybrid # or: does not know it is hybrid
@test net.node[10].isBadDiamondII # or: does not know it is bad diamond II
##plot(net,showEdgeNumber=true)
@test [e.inCycle for e in net.edge] == [-1,-1,3,3,-1,-1,-1,-1,-1,-1,3,3] # or: error in incycle
@test [e.containRoot for e in net.edge] == [true,true,false,true,true,false,false,false,false,false,false,true] # or: error in contain root
@test [e.istIdentifiable for e in net.edge] == [false,false,true,true,false,false,false,true,false,false,true,true] # or: istIdentifiable not correct
@test (net.edge[3].hybrid && net.edge[11].hybrid) # or: hybrid edges wrong")

df=DataFrame(t1=["1","1","2","2","1","2","2","2","2","2","1","2","2","3","2"],
             t2=["3","3","3","3","3","1","5","1","1","1","5","1","1","5","3"],
             t3=["5","5","5","5","6","5","6","3","6","5","6","3","3","6","6"],
             t4=["6","4","4","6","4","6","4","5","4","4","4","6","4","4","4"],
             CF1234=[0.565,0.0005,0.0005,0.565,0.00003,0.99986,0.0410167,1,0.99987,1,0.040167,0.998667,1,0.073167,0.00003],
             CF1324=[0.0903,0.8599,0.8599,0.0903,0.8885,0.00006,0.263,0,0.00006,0,0.2630,0.00006,0,0.0424667,0.8885])
df[!,:CF1423] = 1.0 .- df[!,:CF1234] .- df[!,:CF1324]
d = readTableCF(df)

Random.seed!(345);
net2 = topologyMaxQPseudolik!(net,d, ftolRel=1e-5,ftolAbs=1e-6,xtolRel=1e-3,xtolAbs=1e-4)
@test net2.loglik < 340.5 # loglik ~ 339.95
@test net2.edge[3].istIdentifiable # or: wrong hybrid is t identifiable
# @test 9.95 < net2.edge[3].length < 9.99
# why it's not a good idea to test the branch length estimate:
# BL~4.92 when loglik~339.95
# BL>9.95 when loglik>340
# BL~0.0  when loglik~339.88 (using tolerances of 1e-8)
@test net2.edge[11].istIdentifiable # or: wrong hybrid is t identifiable
@test net2.edge[11].length < 0.01 # or: wrong bl estimated
@test net2.edge[10].length == 0.0 # or: tree edge in bad diamond II not 0
#printEdges(net2)

@test_logs show(devnull, net2)
@test_logs [show(devnull, net2.node[i]) for i in [1,3,10]];
@test_logs [show(devnull, net2.edge[i]) for i in [1,3,11]];
@test_logs show(devnull, d)
@test_logs show(devnull, d.quartet[1])
@test_logs show(devnull, d.quartet[1].qnet)
@test tipLabels(d.quartet) == ["1","2","3","4","5","6"]
a = (@test_logs fittedQuartetCF(d, :wide));
@test size(a) == (15,10)
a = (@test_logs fittedQuartetCF(d, :long));
@test size(a) == (45,7)

end

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
