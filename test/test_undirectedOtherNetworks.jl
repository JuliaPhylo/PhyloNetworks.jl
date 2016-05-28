## julia script to test undirectedOtherNetworks
## function that will give all the networks obtained by
## moving the hybrid nodes inside each cycle
## Claudia April 2016

include("../src/types.jl")
include("../src/functions.jl")

net = readTopologyLevel1("(6,(((1,2):1.8702617089780738,(3,(4)#H7:0.042493238243074646::0.9015570666393798):10.0):10.0,(5,#H7:1.1830905092006163::0.09844293336062017):10.0):0.5);")
out = "6"

net = readTopologyLevel1("(4:1.0,(3:1.0,((1:1.0,2:1.0):1.87)#H7:10.0::0.51):0.042,(5:1.0,(#H7:10.0::0.49,6:1.0):10.0):1.183);")
## error found when new node is -1: inCycle is meaningless

net = readTopologyLevel1("(6,((5,#H7:0.0::0.29999):10.0,(((2,#H9:0.0::0.15902):0.02332,(1,(3)#H9:0.9417::0.84098):0.12997):1.2018,(4)#H7:0.01722::0.70001):9.99428):0.24593);")
net = readTopologyLevel1("(Xgordoni,Xmeyeri,(Xcouchianus,(Xvariatus,(Xevelynae,((Xxiphidium,#H25:9.992::0.167):1.383,(Xmilleri,(Xandersi,(Xmaculatus,((((Xhellerii,(Xalvarezi,Xmayae):0.327):0.259,Xsignum):1.866,(Xclemenciae_F2,Xmonticolus):1.461):0.786,((((Xmontezumae,(Xnezahuacoyotl)#H26:0.247::0.807):0.372,((Xbirchmanni_GARC,Xmalinche_CHIC2):1.003,Xcortezi):0.454):0.63,((Xcontinens,Xpygmaeus):1.927,((Xnigrensis,Xmultilineatus):1.304,#H26:0.0::0.193):0.059):2.492):2.034)#H25:0.707::0.833):1.029):0.654):0.469):0.295):0.41):0.646):3.509):0.263);")



plot(net, showNodeNumber=true, showEdgeNumber=true)
otherNet1 = undirectedOtherNetworks(net)
length(otherNet1)
i=4
plot(otherNet1[i])

otherNet2 = undirectedOtherNetworks(net, outgroup=out)
length(otherNet2)
i=3
plot(otherNet2[i])

