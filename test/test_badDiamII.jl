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


## testing readTopology----------------------------------------------------------------
tree = "(6,(5,#H7:0.0):9.970714072991349,(3,(((2,1):0.2950382234364404,4):0.036924483697671304)#H7:0.00926495670648208):1.1071489442240392);"
net2 = readTopology(tree)
printEdges(net2)

for(e in net2.edge)
    if(e.node[1].leaf || e.node[2].leaf)
        !e.istIdentifiable || error("ext edge should not identifiable")
    else
        e.istIdentifiable || error("int edge should not identifiable")
    end
end

## this case works fine
tree = "((((8,10))#H1,7),6,(4,#H1));" # Case I Bad diamond I
net = readTopologyLevel1(tree)
checkNet(net)
net2 = readTopology(tree)
printEdges(net2)

for(e in net2.edge)
    if(e.node[1].leaf || e.node[2].leaf)
        !e.istIdentifiable || error("ext edge should not identifiable")
    else
        e.istIdentifiable || error("int edge should not identifiable")
    end
end
