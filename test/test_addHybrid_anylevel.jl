#= # for local testing, need this:
using Test
using PhyloNetworks
using PhyloPlots
using CSV
=#

@testset "addhybridedge! top function" begin
str_tree = "(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);";
tree = readTopology(str_tree)
@test !isnothing(PhyloNetworks.addhybridedge!(tree, true)) # should be able to add a hybrid
@test tree.numHybrids == 1
@test !isnothing(PhyloNetworks.addhybridedge!(tree, true)) # should be able to add a hybrid
@test tree.numHybrids == 2
@test !any([n.hybrid for n in PhyloNetworks.getParents(tree.hybrid[2])]) # tests if network is treechild

tree = readTopology(str_tree)
@test !isnothing(PhyloNetworks.addhybridedge!(tree, false)) # should be able to add a hybrid
@test tree.numHybrids == 1

str_level1 = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));" 
net_level1 = readTopology(str_level1)
@test !isnothing(PhyloNetworks.addhybridedge!(net_level1, true)) # should be able to add a hybrid
@test net_level1.numHybrids == 3
@test !any([n.hybrid for n in PhyloNetworks.getParents(net_level1.hybrid[3])]) # tests if network is treechild

net_level1 = readTopology(str_level1)
@test !isnothing(PhyloNetworks.addhybridedge!(net_level1, false)) # should be able to add a hybrid
@test net_level1.numHybrids == 3
end

@testset "addhybridedge! helper function" begin
str_level1 = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));" 
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

net_level1 = readTopology(str_level1); # case 3 (good hybrid edge choice leads to a DAG)
@test !isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[6], net_level1.edge[20], false))

# NEW HYBRID INTO AN EXISTING HYBRID EDGE 
net_level1 = readTopology(str_level1);
@test !isnothing(PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[9], net_level1.edge[10], false))
end

@testset "edge checking functions" begin
str_level1 = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));" 
net_level1 = readTopology(str_level1)

# 3-CYCLE ERRORS
    # throws error if addhybridedge would create a 3 cycle
net_level1 = readTopology(str_level1);
@test PhyloNetworks.hybrid3cycle(net_level1, net_level1.edge[6], net_level1.edge[9])
@test PhyloNetworks.hybrid3cycle(net_level1, net_level1.edge[3], net_level1.edge[17])
@test PhyloNetworks.hybrid3cycle(net_level1, net_level1.edge[16], net_level1.edge[17])
@test PhyloNetworks.hybrid3cycle(net_level1, net_level1.edge[16], net_level1.edge[3])
@test PhyloNetworks.hybrid3cycle(net_level1, net_level1.edge[16], net_level1.edge[18])

# DIRECTIONAL ERRORS
    # throws error if edge 1 is a directed descendant of edge 2
net_level1 = readTopology(str_level1);
# case 6
@test PhyloNetworks.directionalconflict(net_level1, net_level1.edge[6], net_level1.edge[15], true) 
# case 2
@test PhyloNetworks.directionalconflict(net_level1, net_level1.edge[6], net_level1.edge[18], true)
# case 3 (bad hybrid edge choice leads to a nonDAG)
@test PhyloNetworks.directionalconflict(net_level1, net_level1.edge[6], net_level1.edge[20], true) 
end

@testset "delete a hybridization" begin
str_level1 = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));" 

net = readTopology(str_level1)
net, newhybridnode = PhyloNetworks.addhybridedge!(net, net.edge[3], net.edge[9], true)
# deletes hybrid, adds to blacklist $(treeedge1.number)
PhyloNetworks.removehybridedge!(net, newhybridnode, true, true)
@test net.numHybrids == 2
@test 9 in net.blacklist # edge 9 should be in blacklist

net = readTopology(str_level1)
net, newhybridnode = PhyloNetworks.addhybridedge!(net, net.edge[3], net.edge[9], true)
# deletes hybrid, does not add to blacklist
PhyloNetworks.removehybridedge!(net, newhybridnode, true, false)
@test net.numHybrids == 2
@test !(9 in net.blacklist) # edge 9 should not be in blacklist
end


