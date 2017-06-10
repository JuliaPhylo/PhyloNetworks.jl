# development of parsimony algorithm: by William Sparks

## test of parsimony on networks

importall PhyloNetworks
include("RCall_network.jl")

println("Testing rcall-based plotting function on a network")

function testrcall_network(net)

    net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")

    RPlotNetworks(net)
end
