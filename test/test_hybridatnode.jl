# test for hybridatnode
# still not automatic test function
# Claudia April 2016


include("../src/types.jl")
include("../src/functions.jl")

net = readTopologyLevel1("(Gma,(Gch,((Gmo_C7,((Gmo_N2,Gmo_I3):0.037,Gmo_I7):0.181):2.318,#H8:::0.03):2.892):2.347,(Bsa)#H8:::0.97);")
#plot(net)
#plot(net, showNodeNumber=true)
# want to change hybrid to -3

net2 = hybridatnode(net,-3)
plot(net2)
plot(net2, showEdgeNumber=true)
rootonedge!(net2,13) #Bsa outgroup
plot(net2)


