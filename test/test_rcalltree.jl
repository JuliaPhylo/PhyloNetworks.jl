# development of parsimony algorithm: by William Sparks

## test of rcall-based plotting function on trees

importall PhyloNetworks
include("RCall_tree.jl")

println("Testing RCall-based plotting functions on a tree")

function testrcall_tree(net, tips)

    net = readTopology("(A,(B,(C,D)));")

    RPlot(net)
end

