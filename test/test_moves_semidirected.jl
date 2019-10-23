#= # for local testing, need this:
using Test
using PhyloNetworks
using PhyloPlots
=#

@testset "unconstrained NNI moves" begin

str_level1 = "(((8,9),(((((1,2,3),4),(5)#H1),(#H1,(6,7))))#H2),(#H2,10));"
net_level1 = readTopology(str_level1);
str_hybridbelowroot = "((8,9),(((((1,2,3),4),(5)#H1),(#H1,(6,7))))#H2,(#H2,10));"
net_hybridbelowroot = readTopology(str_hybridbelowroot)
# same topology as: rootatnode!(net_level1, -3). edges 1:22
str_nontreechild = "((((Ag,E))#H3,(#H1:7.159::0.056,((M:0.0)#H2:::0.996,(Ak,(#H3:0.08,#H2:0.0::0.004):0.023):0.078):2.49):2.214):0.026,((Az:2.13,As:2.027):1.697)#H1:0.0::0.944,Ap);"
net_nontreechild = readTopology(str_nontreechild);
# problem: the plot has an extra vertical segment, for a clade that's not in the major tree
# --> fix that in PhyloPlots (fixit)
str_hybridladder = "(((B)#H1)#H2,((C,#H2:::0.8),(#H1,(A1,A2))),O);"
net_hybridladder = readTopology(str_hybridladder);

#=
plot(net_level1, :R, showNodeNumber=true, showEdgeNumber=true)
plot(net_hybridbelowroot, :R, showNodeNumber=true, showEdgeNumber=true)
plot(net_nontreechild, :R, showNodeNumber=true, showEdgeNumber=true)
plot(net_hybridladder, :R, showNodeNumber=true, showEdgeNumber=true)
=#

@test isnothing(PhyloNetworks.nni!(net_level1, net_level1.edge[1], 0x01)) # external edge

@testset "edge 3: BB undirected, move $move" for move in 0x01:0x08
    undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[3], move);
    PhyloNetworks.nni!(undoinfo...);
    @test writeTopology(net_level1) == str_level1
end
@testset "edge 17: BB directed, move $move" for move in 0x01:0x02
    undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[16], move);
    PhyloNetworks.nni!(undoinfo...);
    @test writeTopology(net_level1) == str_level1
end
@test_throws Exception PhyloNetworks.nni!(net_level1, net_level1.edge[16], 0x03);
@testset "edge 13: BR directed, move $move" for move in 0x01:0x03
    undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[13], move);
    PhyloNetworks.nni!(undoinfo...);
    @test writeTopology(net_level1) == str_level1
end
end # of testset on unconstrained NNIs

#=
myconstraint = PhyloNetworks.TopologyConstraint(0x03,Set("8"),net_level1)
undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[1], myconstraint)
=#
