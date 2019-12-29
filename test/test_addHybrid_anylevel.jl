#= # for local testing, need this:
using Test
using PhyloNetworks
using PhyloPlots
using CSV
=#

@testset "addhybridedge" begin
str_level1 = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));" # indviduals S1A S1B S1C go on leaf 1
net_level1 = readTopology(str_level1)

# ALLOWED MOVES
net, newhybridnode = PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[3], net_level1.edge[9], true)
@test newhybridnode.hybrid
@test PhyloNetworks.getMajorParentEdge(newhybridnode).gamma >= 0.5
@test PhyloNetworks.getMajorParentEdge(newhybridnode).hybrid
@test PhyloNetworks.getMinorParentEdge(newhybridnode).gamma <= 0.5
@test PhyloNetworks.getMinorParentEdge(newhybridnode).hybrid

net_level1 = readTopology(str_level1);
@test !isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[15], net_level1.edge[3], true))
@test net_level1.numHybrids == 3
@test length(net_level1.edge) == 23

net_level1 = readTopology(str_level1);
@test !isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[2], net_level1.edge[17], true))
@test net_level1.numHybrids == 3
@test length(net_level1.edge) == 23

net_level1 = readTopology(str_level1);
@test !isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[20], net_level1.edge[16], true))
@test net_level1.numHybrids == 3
@test length(net_level1.edge) == 23

# case 3 (good hybrid edge choice leads to a DAG)
net_level1 = readTopology(str_level1);
@test !isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[6], net_level1.edge[20], false)) 

# DIRECTIONAL ERRORS
    # throws error if edge 1 is a directed descendant of edge 2
net_level1 = readTopology(str_level1);
# case 6
@test_throws ErrorException("directional conflict: edge 6 is a directional descendant of edge 15.") PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[6], net_level1.edge[15], true) 
# case 2
@test_throws ErrorException("directional conflict: edge 6 is a directional descendant of edge 18.") PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[6], net_level1.edge[18], true)
# case 3 (bad hybrid edge choice leads to a nonDAG)
@test_throws ErrorException("directional conflict: edge 6 is a directional descendant of edge 20.") PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[6], net_level1.edge[20], true)

# 3-CYCLE ERRORS
    # throws error if addhybridedge would create a 3 cycle
net_level1 = readTopology(str_level1);
@test_throws ErrorException("hybrid between edge 6 and edge 9 would create 3-cycle.") PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[6], net_level1.edge[9], true)
@test_throws ErrorException("hybrid between edge 3 and edge 17 would create 3-cycle.") PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[3], net_level1.edge[17], true)
@test_throws ErrorException("hybrid between edge 16 and edge 17 would create 3-cycle.") PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[16], net_level1.edge[17], true)
@test_throws ErrorException("hybrid between edge 16 and edge 3 would create 3-cycle.") PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[16], net_level1.edge[3], true)
@test_throws ErrorException("hybrid between edge 16 and edge 18 would create 3-cycle.") PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[16], net_level1.edge[18], true)

# NEW HYBRID TO AN EXISTING HYBRID EDGE #? do we want to allow this?
net_level1 = readTopology(str_level1);
@test !isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[9], net_level1.edge[10], false)) 
# TODO this currently creates a stack overflow error in isdirectionaldescendant!
end
