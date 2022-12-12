#= # for local testing, need this:
using Test
using PhyloNetworks
using Random
=#

@testset "addhybridedge! top function" begin
# caution: this test has randomness in its choice of edges; may error sometimes and not others
str_tree = "(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);";
tree = readTopology(str_tree)
Random.seed!(5432);
@test !isnothing(PhyloNetworks.addhybridedge!(tree, true, true))
@test tree.numHybrids == 1
@test !isnothing(PhyloNetworks.addhybridedge!(tree, true, true)) # should be able to add a hybrid
@test tree.numHybrids == 2
@test !any([n.hybrid for n in PhyloNetworks.getparents(tree.hybrid[2])]) # tests if network is treechild

str_level1 = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
netl1 = readTopology(str_level1)
@test !isnothing(PhyloNetworks.addhybridedge!(netl1, true, true))
@test netl1.numHybrids == 3
@test !any([n.hybrid for n in PhyloNetworks.getparents(netl1.hybrid[3])]) # tests if network has no hybrid ladder

netl1 = readTopology(str_level1)
newhybridnode, newhybridedge = PhyloNetworks.addhybridedge!(netl1, false, true)
@test !isnothing(newhybridnode)
@test netl1.numHybrids == 3
PhyloNetworks.deletehybridedge!(netl1, PhyloNetworks.getMinorParentEdge(newhybridnode))
@test hardwiredClusterDistance(netl1, readTopology(str_level1), true) == 0
end # of addhybridedge! top function

@testset "addhybridedge! helper function" begin
str_level1 = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
# allowed moves
netl1 = readTopology(str_level1)
newhybridnode, newhybridedge = PhyloNetworks.addhybridedge!(netl1, netl1.edge[3], netl1.edge[9], true, 0.0, 0.2)
@test newhybridnode.hybrid
@test PhyloNetworks.getparentedge(newhybridnode).gamma == 0.8
@test PhyloNetworks.getMinorParentEdge(newhybridnode).gamma == 0.2
netl1 = readTopology(str_level1);
@test !isnothing(PhyloNetworks.addhybridedge!(netl1, netl1.edge[15], netl1.edge[3], true))
@test writeTopology(netl1) == "(((((((S1,S4),(S5)#H1),(#H1,(S6,S7))),#H3))#H2,((S8,S9))#H3),(#H2,S10));"
netl1 = readTopology(str_level1);
@test !isnothing(PhyloNetworks.addhybridedge!(netl1, netl1.edge[2], netl1.edge[17], true))
@test writeTopology(netl1) == "((#H2,S10),(((S8,(S9,#H3)),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2))#H3);"
netl1 = readTopology(str_level1);
@test !isnothing(PhyloNetworks.addhybridedge!(netl1, netl1.edge[20], netl1.edge[16], true))
@test writeTopology(netl1) == "(((S8,S9),(((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2)#H3),((#H2,S10),#H3));"
netl1 = readTopology(str_level1); # good hybrid edge choice leads to a DAG when reverting the direction of edge2
@test !isnothing(PhyloNetworks.addhybridedge!(netl1, netl1.edge[6], netl1.edge[20], false))
@test netl1.root < 19 # the root must have been changed due to changing some edges' directions
@test writeTopology(netl1) == "(#H2,S10,((((S8,S9),((((S5)#H1,((S1,S4),#H3)),(#H1,(S6,S7))))#H2)))#H3);"
# new hybrid into an existing hybrid edge
netl1 = readTopology(str_level1);
@test !isnothing(PhyloNetworks.addhybridedge!(netl1, netl1.edge[9], netl1.edge[10], true))
@test writeTopology(netl1) == "(((S8,S9),((((S6,S7),(#H1)#H3),(((S1,S4),(S5)#H1),#H3)))#H2),(#H2,S10));"
end # of addhybridedge! helper function

@testset "edge checking functions" begin
str_level1 = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
# 3-cycles: throws error if adding a hybrid edge would create a 3-cycle
netl1 = readTopology(str_level1);
@test PhyloNetworks.hybrid3cycle(netl1.edge[6], netl1.edge[9])
@test PhyloNetworks.hybrid3cycle(netl1.edge[3], netl1.edge[17])
@test PhyloNetworks.hybrid3cycle(netl1.edge[16], netl1.edge[17])
@test PhyloNetworks.hybrid3cycle(netl1.edge[16], netl1.edge[3])
@test PhyloNetworks.hybrid3cycle(netl1.edge[16], netl1.edge[18])
@test !PhyloNetworks.hybrid3cycle(netl1.edge[16], netl1.edge[20])

# directional: throws error if the new network would not be a DAG, e.g. if edge 1 is a directed descendant of edge 2
# case 6
nodeS145 = PhyloNetworks.getparent(netl1.edge[6])
@test PhyloNetworks.directionalconflict(netl1, nodeS145, netl1.edge[15], true)
# case 2
@test PhyloNetworks.directionalconflict(netl1, nodeS145, netl1.edge[18], true)
# case 3 (bad hybrid edge choice leads to a nonDAG)
@test PhyloNetworks.directionalconflict(netl1, nodeS145, netl1.edge[20], true)
@test PhyloNetworks.directionalconflict(netl1, nodeS145, netl1.edge[4], false)
@test !PhyloNetworks.directionalconflict(netl1, nodeS145, netl1.edge[4], true)
end # of edge checking functions
