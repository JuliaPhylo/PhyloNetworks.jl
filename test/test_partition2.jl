# functions to try to debug cui3.out
# converted to test function

if !isdefined(:individualtest) individualtest = false; end

if(individualtest)
    include("../src/types.jl")
    include("../src/functions.jl")
    const DEBUG = true
end

if isdefined(:PhyloNetworks)
    PhyloNetworks.CHECKNET || error("need CHECKNET==true in PhyloNetworks to test snaq in test_correctLik.jl")
else
    CHECKNET || error("need CHECKNET==true to test snaq in test_correctLik.jl")
end

text = "(Xmayae,((Xhellerii,(((Xclemenciae_F2,Xmonticolus):1.458,(((((Xmontezumae,(Xnezahuacoyotl)#H26:0.247::0.804):0.375,((Xbirchmanni_GARC,Xmalinche_CHIC2):0.997,Xcortezi):0.455):0.63,(#H26:0.0::0.196,((Xcontinens,Xpygmaeus):1.932,(Xnigrensis,Xmultilineatus):1.401):0.042):2.439):2.0)#H7:0.787::0.835,(Xmaculatus,(Xandersi,(Xmilleri,((Xxiphidium,#H7:9.563::0.165):1.409,(Xevelynae,(Xvariatus,(Xcouchianus,(Xgordoni,Xmeyeri):0.263):3.532):0.642):0.411):0.295):0.468):0.654):1.022):0.788):1.917)#H27:0.149::0.572):0.668,Xalvarezi):0.257,(Xsignum,#H27:1.381::0.428):4.669);"


#net = readTopology("cui3.out")
net = readTopology(text)
printPartitions(net)
printEdges(net)
cleanBL!(net)
cleanAfterRead!(net,false)
length(net.hybrid)
[n.number for n in net.hybrid] == [7,15,25] || error("wrong hybrid nodes")

# ----------------- update everything but partition ------------------
n=net.hybrid[1];
flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,n);
flag || error("wrong update in hybrid 1")
[n.number for n in nodesInCycle] == [7,-16,-10,-11,-12] || error("wrong update in hybrid 1")
[n.number for n in edgesInCycle] == [16,24,15,9,8] || error("wrong update in hybrid 1")
flag2, edgesGammaz = updateGammaz!(net,n,false);
flag2 || error("wrong update in hybrid 1")
isempty([n.number for n in edgesGammaz]) || error("wrong update in hybrid 1")
!n.isBadDiamondI || error("wrong update in hybrid 1")
!n.isBadDiamondII || error("wrong update in hybrid 1")
flag3, edgesRoot = updateContainRoot!(net,n);
flag3 || error("wrong update in hybrid 1")
[n.number for n in edgesRoot] == [7] || error("wrong update in hybrid 1")

n=net.hybrid[2];
flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,n);
flag || error("wrong update in hybrid 2")
[n.number for n in nodesInCycle] == [15,-24,-23,-22,-21,-20,-8] || error("wrong update in hybrid 2")
[n.number for n in edgesInCycle] == [31,32,42,43,44,45,26] || error("wrong update in hybrid 2")
flag2, edgesGammaz = updateGammaz!(net,n,false);
flag2 || error("wrong update in hybrid 2")
isempty([n.number for n in edgesGammaz]) || error("wrong update in hybrid 2")
!n.isBadDiamondI || error("wrong update in hybrid 2")
!n.isBadDiamondII || error("wrong update in hybrid 2")
flag3, edgesRoot = updateContainRoot!(net,n);
flag3 || error("wrong update in hybrid 2")
[n.number for n in edgesRoot] == [25,15,9,6,14,12,10,11,13,24,23,19,17,18,22,20,21] || error("wrong update in hybrid 2")


n=net.hybrid[3];
flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,n);
flag || error("wrong update in hybrid 3")
[n.number for n in nodesInCycle] == [25,-29,-2,-3,-4] || error("wrong update in hybrid 3")
[n.number for n in edgesInCycle] == [53,54,51,49,48] || error("wrong update in hybrid 3")
flag2, edgesGammaz = updateGammaz!(net,n,false);
flag2 || error("wrong update in hybrid 3")
isempty([n.number for n in edgesGammaz]) || error("wrong update in hybrid 3")
!n.isBadDiamondI || error("wrong update in hybrid 3")
!n.isBadDiamondII || error("wrong update in hybrid 3")
flag3, edgesRoot = updateContainRoot!(net,n);
flag3 || error("wrong update in hybrid 3")
[n.number for n in edgesRoot] == [47,5,3,4,46,45,27,44,28,43,29,42,32,30,41,33,40,34,39,35,38,36,37] || error("wrong update in hybrid 3")

isempty(net.partition) || error("wrong partition after initial update")

# -------------- update partition ---------------------
n=net.hybrid[1];
nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,n);
[n.number for n in edgesInCycle] == [16,24,15,9,8] || error("wrong partition hybrid 1")
[n.number for n in nodesInCycle] == [7,-16,-10,-11,-12] || error("wrong partition hybrid 1")
updatePartition!(net,nodesInCycle)
length(net.partition) == 5 || error("wrong partition hybrid 1")
[n.number for n in net.partition[1].edges] == [7] || error("wrong partition hybrid 1")
[n.number for n in net.partition[2].edges] == [23,19,17,18,22,20,21] || error("wrong partition hybrid 1")
[n.number for n in net.partition[3].edges] == [25] || error("wrong partition hybrid 1")
[n.number for n in net.partition[4].edges] == [14,12,10,11,13] || error("wrong partition hybrid 1")
[n.number for n in net.partition[5].edges] == [6] || error("wrong partition hybrid 1")

n=net.hybrid[2];
nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,n);
[n.number for n in nodesInCycle] == [15,-24,-23,-22,-21,-20,-8] || error("wrong partition in hybrid 2")
[n.number for n in edgesInCycle] == [31,32,42,43,44,45,26] || error("wrong partition in hybrid 2")
updatePartition!(net,nodesInCycle)
length(net.partition) == 11 || error("wrong partition hybrid 2")
[n.number for n in net.partition[6].edges] == [30] || error("wrong partition hybrid 2")
[n.number for n in net.partition[7].edges] == [41,33,40,34,39,35,38,36,37] || error("wrong partition hybrid 2")
[n.number for n in net.partition[8].edges] == [29] || error("wrong partition hybrid 2")
[n.number for n in net.partition[9].edges] == [28] || error("wrong partition hybrid 2")
[n.number for n in net.partition[10].edges] == [27] || error("wrong partition hybrid 2")
[n.number for n in net.partition[11].edges] == [46,5,3,4,47] || error("wrong partition hybrid 2")

n=net.hybrid[3];
nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,n);
[n.number for n in nodesInCycle] == [25,-29,-2,-3,-4] || error("wrong partition in hybrid 3")
[n.number for n in edgesInCycle] == [53,54,51,49,48] || error("wrong partition in hybrid 3")
updatePartition!(net,nodesInCycle)
length(net.partition) == 15 || error("wrong partition hybrid 2")
[n.number for n in net.partition[12].edges] == [52] || error("wrong partition hybrid 3")
[n.number for n in net.partition[13].edges] == [1] || error("wrong partition hybrid 3")
[n.number for n in net.partition[14].edges] == [50] || error("wrong partition hybrid 3")
[n.number for n in net.partition[15].edges] == [2] || error("wrong partition hybrid 3")


printPartitions(net)
