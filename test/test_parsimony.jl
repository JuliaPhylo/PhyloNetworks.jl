# development of parsimony algorithm: by William Sparks

## test of parsimony on trees

importall PhyloNetworks
include("parsimony.jl")

println("Testing parsimony() on a tree")

function testparsimony_tree(net, tips)

    net = readTopology("(A,(B,(C,D)));")

    tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1)
end

## test of parsimony on networks

println("Testing parsimony() on a network")

function testparsimony_network(net, tips)

    net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")

    tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1)

    parsimony(net, tips)
end

