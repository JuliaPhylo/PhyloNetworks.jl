# functions to try to debug cui3.out
# converted to test function

PhyloNetworks.CHECKNET || error("need CHECKNET==true in PhyloNetworks to test snaq in test_correctLik.jl")

@testset "test: level-1 partitions" begin

global net, flag, flag2, flag3, nocycle
text = "(Xmayae,((Xhellerii,(((Xclemenciae_F2,Xmonticolus):1.458,(((((Xmontezumae,(Xnezahuacoyotl)#H26:0.247::0.804):0.375,((Xbirchmanni_GARC,Xmalinche_CHIC2):0.997,Xcortezi):0.455):0.63,(#H26:0.0::0.196,((Xcontinens,Xpygmaeus):1.932,(Xnigrensis,Xmultilineatus):1.401):0.042):2.439):2.0)#H7:0.787::0.835,(Xmaculatus,(Xandersi,(Xmilleri,((Xxiphidium,#H7:9.563::0.165):1.409,(Xevelynae,(Xvariatus,(Xcouchianus,(Xgordoni,Xmeyeri):0.263):3.532):0.642):0.411):0.295):0.468):0.654):1.022):0.788):1.917)#H27:0.149::0.572):0.668,Xalvarezi):0.257,(Xsignum,#H27:1.381::0.428):4.669);"

#net = readTopology("cui3.out")
net = readTopology(text)
#printPartitions(net)
#printEdges(net)
cleanBL!(net)
cleanAfterRead!(net,false)
length(net.hybrid)
@test [n.number for n in net.hybrid] == [7,15,25] # or: wrong hybrid nodes"

# ----------------- update everything but partition ------------------
@testset "udpate hybrid 1" begin
n=net.hybrid[1];
flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,n);
@test flag
@test [n.number for n in nodesInCycle] == [7,-16,-10,-11,-12]
@test [n.number for n in edgesInCycle] == [16,24,15,9,8]
flag2, edgesGammaz = updateGammaz!(net,n,false);
@test flag2
@test isempty([n.number for n in edgesGammaz])
@test !n.isBadDiamondI
@test !n.isBadDiamondII
flag3, edgesRoot = updateContainRoot!(net,n);
@test flag3
@test [n.number for n in edgesRoot] == [7]
end

@testset "udpate hybrid 2" begin
n=net.hybrid[2];
flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,n);
@test flag
@test [n.number for n in nodesInCycle] == [15,-24,-23,-22,-21,-20,-8]
@test [n.number for n in edgesInCycle] == [31,32,42,43,44,45,26]
flag2, edgesGammaz = updateGammaz!(net,n,false);
@test flag2
@test isempty([n.number for n in edgesGammaz])
@test !n.isBadDiamondI
@test !n.isBadDiamondII
flag3, edgesRoot = updateContainRoot!(net,n);
@test flag3
@test [n.number for n in edgesRoot] == [25,15,9,6,14,12,10,11,13,24,23,19,17,18,22,20,21]
end

@testset "udpate hybrid 3" begin
n=net.hybrid[3];
flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,n);
@test flag
@test [n.number for n in nodesInCycle] == [25,-29,-2,-3,-4]
@test [n.number for n in edgesInCycle] == [53,54,51,49,48]
flag2, edgesGammaz = updateGammaz!(net,n,false);
@test flag2
@test isempty([n.number for n in edgesGammaz])
@test !n.isBadDiamondI
@test !n.isBadDiamondII
flag3, edgesRoot = updateContainRoot!(net,n);
@test flag3
@test [n.number for n in edgesRoot] == [47,5,3,4,46,45,27,44,28,43,29,42,32,30,41,33,40,34,39,35,38,36,37]

@test isempty(net.partition)
end

# -------------- update partition ---------------------
@testset "partition hybrid 1" begin
n=net.hybrid[1];
nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,n);
@test [n.number for n in edgesInCycle] == [16,24,15,9,8]
@test [n.number for n in nodesInCycle] == [7,-16,-10,-11,-12]
updatePartition!(net,nodesInCycle)
@test length(net.partition) == 5
@test [n.number for n in net.partition[1].edges] == [7]
@test [n.number for n in net.partition[2].edges] == [23,19,17,18,22,20,21]
@test [n.number for n in net.partition[3].edges] == [25]
@test [n.number for n in net.partition[4].edges] == [14,12,10,11,13]
@test [n.number for n in net.partition[5].edges] == [6]
end

@testset "partition hybrid 2" begin
n=net.hybrid[2];
nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,n);
@test [n.number for n in nodesInCycle] == [15,-24,-23,-22,-21,-20,-8]
@test [n.number for n in edgesInCycle] == [31,32,42,43,44,45,26]
updatePartition!(net,nodesInCycle)
@test length(net.partition) == 11
@test [n.number for n in net.partition[6].edges] == [30]
@test [n.number for n in net.partition[7].edges] == [41,33,40,34,39,35,38,36,37]
@test [n.number for n in net.partition[8].edges] == [29]
@test [n.number for n in net.partition[9].edges] == [28]
@test [n.number for n in net.partition[10].edges] == [27]
@test [n.number for n in net.partition[11].edges] == [46,5,3,4,47]
end

@testset "partition hybrid 3" begin
n=net.hybrid[3];
nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,n);
@test [n.number for n in nodesInCycle] == [25,-29,-2,-3,-4]
@test [n.number for n in edgesInCycle] == [53,54,51,49,48]
updatePartition!(net,nodesInCycle)
@test length(net.partition) == 15
@test [n.number for n in net.partition[12].edges] == [52]
@test [n.number for n in net.partition[13].edges] == [1]
@test [n.number for n in net.partition[14].edges] == [50]
@test [n.number for n in net.partition[15].edges] == [2]
end
#printPartitions(net)
end
