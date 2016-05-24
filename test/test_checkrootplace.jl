## script to test function checkRootPlace
## Claudia May 2016

include("../src/types.jl")
include("../src/functions.jl")

net = readTopologyLevel1("(6,((5,#H7:0.0::0.29999):10.0,(((2,#H9:0.0::0.15902):0.02332,(1,(3)#H9:0.9417::0.84098):0.12997):1.2018,(4)#H7:0.01722::0.70001):9.99428):0.24593);")
plot(net)
checkRootPlace!(net,outgroup="6")
plot(net)
checkRootPlace!(net,outgroup="2")
plot(net)
checkRootPlace!(net,outgroup="1")
plot(net)
checkRootPlace!(net,outgroup="3") ## can't do it, does not change net
plot(net)
