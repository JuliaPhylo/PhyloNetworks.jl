# test for hybridatnode
# still not automatic test function
# Claudia April 2016


include("../src/types.jl")
include("../src/functions.jl")

net = readTopologyLevel1("(Gma,(Gch,((Gmo_C7,((Gmo_N2,Gmo_I3):0.037,Gmo_I7):0.181):2.318,#H8:::0.03):2.892):2.347,(Bsa)#H8:::0.97);")
plot(net)
plot(net, showNodeNumber=true)
# want to change hybrid to -3

net2 = hybridatnode(net,-3)
plot(net2)
plot(net2, showEdgeNumber=true)
rootonedge!(net2,13) #Bsa outgroup
plot(net2)


hybridatnode!(net,-3)
plot(net)


net = readTopology("(A:1.0,(B,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,D:0.8):1.3):0.7):0.1;");
plot(net, showNodeNumber=true)
hybridatnode!(net, -7)
plot(net)
