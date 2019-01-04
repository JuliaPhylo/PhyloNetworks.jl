## julia script to test undirectedOtherNetworks
## function that will give all the networks obtained by
## moving the hybrid nodes inside each cycle
## Claudia April 2016

@testset "test: move hybrid around cycle" begin
global net
net = readTopology("(((B)#H2,((C,#H2)S1,(A)S2)S3)S4);")
@test_logs hybridatnode!(net, 7)
@test_throws ErrorException hybridatnode!(net, 8)

## very important semicolon at the end because show changes containRoot
net = readTopologyLevel1("(6,(((1,2):1.8702617089780738,(3,(4)#H7:0.042493238243074646::0.9015570666393798):10.0):10.0,(5,#H7:1.1830905092006163::0.09844293336062017):10.0):0.5);");
out = "6"
##plot(net, showNodeNumber=true, showEdgeNumber=true)
checkNet(net)
otherNet1 = undirectedOtherNetworks(net, insideSnaq=true);
length(otherNet1) == 4 || error("wrong number of other networks")
otherNet1[1].hybrid[1].number == -8 || error("wrong new hybrid")
otherNet1[2].hybrid[1].number == -3 || error("wrong new hybrid")
otherNet1[3].hybrid[1].number == -4 || error("wrong new hybrid")
otherNet1[4].hybrid[1].number == -6 || error("wrong new hybrid")
checkNet(otherNet1[1], false,checkPartition=false)
checkNet(otherNet1[2], false,checkPartition=false)
checkNet(otherNet1[3], false,checkPartition=false)
checkNet(otherNet1[4], false,checkPartition=false)
##i=4
##plot(otherNet1[i], showEdgeNumber=true)
otherNet2 = undirectedOtherNetworks(net, outgroup=out, insideSnaq=true);
length(otherNet2) == 3 || error("wrong number of other networks")
otherNet2[1].hybrid[1].number == -8 || error("wrong new hybrid")
otherNet2[2].hybrid[1].number == -4 || error("wrong new hybrid")
otherNet2[3].hybrid[1].number == -6 || error("wrong new hybrid")
checkNet(otherNet2[1], false,checkPartition=false)
checkNet(otherNet2[2], false,checkPartition=false)
checkNet(otherNet2[3], false,checkPartition=false)
##i=3
##plot(otherNet2[i])



net = readTopologyLevel1("(4:1.0,(3:1.0,((1:1.0,2:1.0):1.87)#H7:10.0::0.51):0.042,(5:1.0,(#H7:10.0::0.49,6:1.0):10.0):1.183);");
out = "4"
checkNet(net)
##plot(net, showNodeNumber=true, showEdgeNumber=true)
otherNet1 = undirectedOtherNetworks(net, insideSnaq=true);
length(otherNet1) == 4 || error("wrong number of other networks")
otherNet1[1].hybrid[1].number == -7 || error("wrong new hybrid")
otherNet1[2].hybrid[1].number == -6 || error("wrong new hybrid")
otherNet1[3].hybrid[1].number == -2 || error("wrong new hybrid")
otherNet1[4].hybrid[1].number == -3 || error("wrong new hybrid")
checkNet(otherNet1[1], false,checkPartition=false)
checkNet(otherNet1[2], false,checkPartition=false)
checkNet(otherNet1[3], false,checkPartition=false)
checkNet(otherNet1[4], false,checkPartition=false)
##i=4
##plot(otherNet1[i], showEdgeNumber=true)
otherNet2 = undirectedOtherNetworks(net, outgroup=out, insideSnaq=true);
length(otherNet2) == 3 || error("wrong number of other networks")
otherNet2[1].hybrid[1].number == -7 || error("wrong new hybrid")
otherNet2[2].hybrid[1].number == -6 || error("wrong new hybrid")
otherNet2[3].hybrid[1].number == -3 || error("wrong new hybrid")
checkNet(otherNet2[1], false,checkPartition=false)
checkNet(otherNet2[2], false,checkPartition=false)
checkNet(otherNet2[3], false,checkPartition=false)
##i=3
##plot(otherNet2[i])


net = readTopologyLevel1("(6,((5,#H7:0.0::0.29999):10.0,(((2,#H9:0.0::0.15902):0.02332,(1,(3)#H9:0.9417::0.84098):0.12997):1.2018,(4)#H7:0.01722::0.70001):9.99428):0.24593);");
out = "6"
checkNet(net)
##plot(net, showNodeNumber=true, showEdgeNumber=true)
otherNet1 = undirectedOtherNetworks(net, insideSnaq=true);
length(otherNet1) == 6 || error("wrong number of other networks")
checkNet(otherNet1[1], false,checkPartition=false)
checkNet(otherNet1[2], false,checkPartition=false)
checkNet(otherNet1[3], false,checkPartition=false)
checkNet(otherNet1[4], false,checkPartition=false)
checkNet(otherNet1[5], false,checkPartition=false)
checkNet(otherNet1[6], false,checkPartition=false)
otherNet1[1].hybrid[2].number == -7 && otherNet1[1].hybrid[1].number == 3 || error("wrong new hybrid")
otherNet1[2].hybrid[2].number == -6 && otherNet1[2].hybrid[1].number == 3 || error("wrong new hybrid")
otherNet1[3].hybrid[2].number == -8 && otherNet1[3].hybrid[1].number == 3 || error("wrong new hybrid")
otherNet1[4].hybrid[1].number == 5 && otherNet1[4].hybrid[2].number == -4 || error("wrong new hybrid")
otherNet1[5].hybrid[1].number == 5 && otherNet1[5].hybrid[2].number == -3 || error("wrong new hybrid")
otherNet1[6].hybrid[1].number == 5 && otherNet1[6].hybrid[2].number == -5 || error("wrong new hybrid")
##i=6
##plot(otherNet1[i], showNodeNumber=true)
otherNet2 = undirectedOtherNetworks(net, outgroup=out, insideSnaq=true);
length(otherNet2) == 4 || error("wrong number of other networks")
checkNet(otherNet2[1], false,checkPartition=false)
checkNet(otherNet2[2], false,checkPartition=false)
checkNet(otherNet2[3], false,checkPartition=false)
checkNet(otherNet2[4], false,checkPartition=false)
otherNet2[1].hybrid[2].number == -7 && otherNet2[1].hybrid[1].number == 3 || error("wrong new hybrid")
otherNet2[2].hybrid[1].number == 3 && otherNet2[2].hybrid[2].number == -8 || error("wrong new hybrid")
otherNet2[3].hybrid[1].number == 5 && otherNet2[3].hybrid[2].number == -4 || error("wrong new hybrid")
otherNet2[4].hybrid[1].number == 5 && otherNet2[4].hybrid[2].number == -5 || error("wrong new hybrid")
##i=3
##plot(otherNet2[i])

## net = readTopologyLevel1("(Xgordoni,Xmeyeri,(Xcouchianus,(Xvariatus,(Xevelynae,((Xxiphidium,#H25:9.992::0.167):1.383,(Xmilleri,(Xandersi,(Xmaculatus,((((Xhellerii,(Xalvarezi,Xmayae):0.327):0.259,Xsignum):1.866,(Xclemenciae_F2,Xmonticolus):1.461):0.786,((((Xmontezumae,(Xnezahuacoyotl)#H26:0.247::0.807):0.372,((Xbirchmanni_GARC,Xmalinche_CHIC2):1.003,Xcortezi):0.454):0.63,((Xcontinens,Xpygmaeus):1.927,((Xnigrensis,Xmultilineatus):1.304,#H26:0.0::0.193):0.059):2.492):2.034)#H25:0.707::0.833):1.029):0.654):0.469):0.295):0.41):0.646):3.509):0.263);")

end
