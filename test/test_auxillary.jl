#= # for local testing, need this:
using Test
using PhyloNetworks
using PhyloPlots
using CSV
=#

@testset "setlengths and setgammas" begin
str_level1_s = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));" # indviduals S1A S1B S1C go on leaf 1
net = readTopology(str_level1_s)

PhyloNetworks.setlengths!([net.edge[1]], [1.1])
@test net.edge[1].length == 1.1
PhyloNetworks.setlengths!([net.edge[3], net.edge[4]], [3.3, 4.4])
@test net.edge[3].length == 3.3
@test net.edge[4].length == 4.4

PhyloNetworks.setgammas!([net.edge[18]], [0.25])
@test net.edge[18].gamma == 0.25
@test net.edge[16].gamma == 0.75

@test PhyloNetworks.getlengths([net.edge[1]]) == [net.edge[1].length]
@test PhyloNetworks.getlengths([net.edge[1], net.edge[5]]) == [net.edge[1].length, net.edge[5].length]

end
