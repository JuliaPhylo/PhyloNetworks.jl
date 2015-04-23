# updateGammaz does not recognize this as bad diamond II
# Claudia April 2015

include("../types.jl")
include("../functions.jl")

tree = "(6,(5,#H7:0.0):9.970714072991349,(3,(((2,1):0.2950382234364404,4):0.036924483697671304)#H7:0.00926495670648208):1.1071489442240392);"
net = readTopologyUpdate(tree);
