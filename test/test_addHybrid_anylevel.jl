#= # for local testing, need this:
using Test
using PhyloNetworks
using PhyloPlots
using CSV
=#

@testset "addHybrid" begin

str_level1 = "(((S8,S9),(((((S1,S2,S3),S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
net_level1 = readTopology(str_level1);

# check that the function adds a hybrid
net, newhybridnode = PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[3], net_level1.edge[9])
@test net == net_level1
@test writeTopology(net) == "(((S8,S9),(((((S1,S2,S3),S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
@test newhybridnode.hybrid
@test getParents(newhybridnode)[1].gamma < 1

@test !isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[15], net_level1.edge[3]))
@test !isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[2], net_level1.edge[17]))
@test !isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[20], net_level1.edge[16]))


# check that wont add hybrid if edge 1 is a directed descendant of edge 2
@test isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[6], net_level1.edge[9]))

# check that wont add hybrid if it would create a 3 cycle
@test isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[3], net_level1.edge[17]))
@test isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[16], net_level1.edge[17]))
@test isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[16], net_level1.edge[3]))
@test isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[16], net_level1.edge[18]))

end
