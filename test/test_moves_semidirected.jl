#= # for local testing, need this:
using Test
using PhyloNetworks
using PhyloPlots
=#

@testset "unconstrained NNI moves" begin

str_level1 = "(((8,9),(((((1,2,3),4),(5)#H1),(#H1,(6,7))))#H2),(#H2,10));"
net_level1 = readTopology(str_level1);
# same topology as: rootatnode!(net_level1, -3). edges 1:22
str_nontreechild = "((((Ag,E))#H3,(#H1:7.159::0.056,((M:0.0)#H2:::0.996,(Ak,(#H3:0.08,#H2:0.0::0.004):0.023):0.078):2.49):2.214):0.026,((Az:2.13,As:2.027):1.697)#H1:0.0::0.944,Ap);"
net_nontreechild = readTopology(str_nontreechild);
# problem: the plot has an extra vertical segment, for a clade that's not in the major tree
# --> fix that in PhyloPlots (fixit)
str_hybridladder = "(#H2:::0.2,((C,((B)#H1)#H2:::0.8),(#H1,(A1,A2))),O);"
net_hybridladder = readTopology(str_hybridladder);

#=
plot(net_level1, :R, showNodeNumber=true, showEdgeNumber=true)
plot(net_nontreechild, :R, showNodeNumber=true, showEdgeNumber=true)
plot(net_hybridladder, :R, showNodeNumber=true, showEdgeNumber=true)
=#

@test isnothing(PhyloNetworks.nni!(net_level1, net_level1.edge[1], 0x01)) # external edge

#TODO test isconnected function

@testset "level1 edge 3: BB undirected move $move" for move in 0x01:0x08
    undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[3], move);
    #location of v node (number -4)
    if move in [1, 2, 6, 8] #check that edge alpha connected to v
        nodes = []
        for n in net_level1.edge[20].node
            push!(nodes, n.number)
        end
        @test -4 in nodes
    elseif move in [3, 4, 5, 7] #check that edge beta connected to v
        nodes = []
        for n in net_level1.edge[19].node
            push!(nodes, n.number)
        end
        @test -4 in nodes
    end
    #location of u node (number -3)
    if move in [2, 4, 5, 6] #check that edge gamma connected to u
        nodes = []
        for n in net_level1.edge[1].node
            push!(nodes, n.number)
        end
        @test -3 in nodes
    elseif move in [1, 3, 7, 8]  #check that edge delta connected to u
        nodes = []
        for n in net_level1.edge[2].node
            push!(nodes, n.number)
        end
        @test -3 in nodes
    end
    #check directionality
    if move in [3, 4, 5, 7] 
        #keeping α, β or flipping uv keeps node -4 as child of edge 3
        @test PhyloNetworks.getChild(net_level1.edge[3]).number == -4
    else 
        #switching α, β AND flipping uv or doing neither makes node -3 child of edge 3
        @test PhyloNetworks.getChild(net_level1.edge[3]).number == -3
    end
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    #confirm we're back to original topology 
    @test writeTopology(net_level1) == str_level1
end #of level1 edge 3: BB undirected

@testset "level1 edge 16: BB directed move $move" for move in 0x01:0x02
    # e not hybrid, tree parent:  BB case, 2 NNIs if directed, 8 if undirected
    undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[16], move);
    nodes = [n.number for n in net_level1.edge[13].node] #beta's connections
    @test -11 in nodes #check beta connected to u node (number -11)
    nodes = [n.number for n in net_level1.edge[14].node] #gamma's connections
    if move == 1 #check that edge gamma connected to v
        @test -12 in nodes
    else #check that edge gamma connected to u
        @test -11 in nodes
    end
    #check directionality
    #node -11 child of edge 16 in both cases
    @test PhyloNetworks.getChild(net_level1.edge[16]).number == -11
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    #confirm we're back to original topology 
    @test writeTopology(net_level1) == str_level1
end #of level1 edge 16: BB directed

@test_throws Exception PhyloNetworks.nni!(net_level1, net_level1.edge[16], 0x03);

#TODO add BR undirected 

#TODO check this looks like it labels 13 as βu because βu = u.edge[ci[1]] = 13
#This is a four-cycle so several moves will not be allowed
@testset "level1 DAG test edge 13: BR directed move $move" for move in 0x01:0x03
    # e.hybrid and tree parent:  BR case, 3 because e cannot contain the root
    # tests for DAG problems too
    #if there is a path from beta to gamma, resulting graph not a DAG
    #there isnt, but there is a path from alpha to delta
    #moves that switch alpha and beta AND alpha and delta should cause problems #TODO or all?
    undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[13], move);
    #test that move was made
    if move == 1 || (move == 3 && PhyloNetworks.getChild(net_level1.edge[17]).number == -11) #checks that alpha -> u
        #check that edge alpha connected to v
        nodes = []
        for n in net_level1.edge[17].node
            push!(nodes, n.number)
        end
        @test 8 in nodes
    elseif move == 2 || (move == 3 && PhyloNetworks.getChild(net_level1.edge[16]).number == -11) #checks beta -> u
        #check that edge beta connected to v
        nodes = []
        for n in net_level1.edge[16].node
            push!(nodes, n.number)
        end
        @test 8 in nodes
    end
    #location of u node (number -11)
    if move in [3] #check that edge gamma connected to u
        nodes = []
        for n in net_level1.edge[11].node
            push!(nodes, n.number)
        end
        @test -11 in nodes
    elseif move in [1, 2]  #check that edge delta connected to u
        nodes = []
        for n in net_level1.edge[10].node
            push!(nodes, n.number)
        end
        @test -11 in nodes
    end
    #check directionality
    if move == 0x01
        @test PhyloNetworks.getChild(net_level1.edge[13]).number == -11
    else
        @test PhyloNetworks.getChild(net_level1.edge[13]).number == 8
    end
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    #confirm we're back to original topology 
    @test writeTopology(net_level1) == str_level1
end #of level1 edge 13: BR directed

@testset "level1 edge 18: RR (directed) move $move" for move in 0x01:0x04
    # RB case, 4 moves. uv edge cannot contain the root (always directed)
    undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[18], move);
    #test that move was made
    nodes = [n.number for n in net_level1.edge[19].node] #alpha's connections
    if move in [0x01, 0x03]
        @test -6 in nodes #check that edge alpha connected to v
    elseif move in [0x02, 0x04]
        @test !(-6 in nodes)
    end
    nodes = [n.number for n in net_level1.edge[12].node] #gamma's connections
    if move in [0x03, 0x04] #check that edge gamma connected to u (11)
        @test 11 in nodes
    elseif move in [0x01, 0x02]  #check that edge delta connected to u, not gamma
        @test !(11 in nodes)
    end
    #check directionality (should point toward u, node 11)
    @test PhyloNetworks.getChild(net_level1.edge[18]).number == 11
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    #confirm we're back to original topology 
    @test writeTopology(net_level1) == str_level1
end #of level1 edge 18: RR (directed)

@testset "non tree child net edge 3: RB (directed) move $move" for move in 0x01:0x04
    # RB case, 4 moves. uv edge cannot contain the root
    undoinfo = PhyloNetworks.nni!(net_nontreechild, net_nontreechild.edge[3], move);
    #test that move was made
    #location of v node (node -5)
    nodes = [n.number for n in net_nontreechild.edge[4].node] #alpha's connections
    if move in [0x01, 0x03]
        #check that edge alpha connected to v
        @test -5 in nodes
    else
        #check that edge beta connected to v, not u
        @test !(-5 in nodes)
    end
    nodes = [n.number for n in net_nontreechild.edge[2].node] #delta's connections
    if move in [0x01, 0x02] #check that edge delta connected to u
        @test 3 in nodes
    else
        @test !(3 in nodes)
    end
    #check directionality (should point toward u, node 3)
    @test PhyloNetworks.getChild(net_nontreechild.edge[3]).number == 3
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    #confirm we're back to original topology 
    @test writeTopology(net_nontreechild) == str_nontreechild
end #of non tree child net edge 5: RB (directed)

@testset "hybrid ladder net edge 4: RR (directed) move $move" for move in 0x01:0x02
    # RR case, 2 moves. uv edge cannot contain the root (always directed)
    undoinfo = PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[4], move);
    #test that move was made
    nodes = [n.number for n in net_hybridladder.edge[5].node] #alpha's connections
    if move == 0x01 #check that edge alpha connected to v
        @test 4 in nodes
    else move == 0x02 #check that edge beta connected to v, not alpha
        @test !(4 in nodes)
    end
    #check that edge delta connected to u, not gamma
    nodes = [n.number for n in net_hybridladder.edge[3].node]
    @test 1 in nodes
    #check directionality (should point toward u, node 1)
    @test PhyloNetworks.getChild(net_hybridladder.edge[4]).number == 1
    #undo move
    PhyloNetworks.nni!(undoinfo...); 
    #confirm we're back to original topology 
    @test writeTopology(net_hybridladder) == str_hybridladder
end #of hybrid ladder net edge 4: RR (directed)

#TODO add test of edge 1
@testset "hybrid ladder net edge 5: BR undirected move $move" for move in 0x01:0x06
    # BR case, 6 moves. uv edge can contain the root
    undoinfo = PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[5], move);
    nodes = [n.number for n in net_hybridladder.edge[6].node] #alpha
    if move in [0x01, 0x05] || (move == 0x03 && PhyloNetworks.getChild(net_hybridladder.edge[6]).number == -6) || 
        (move == 0x06 && PhyloNetworks.getChild(net_hybridladder.edge[4]).number == -6)
        @test 1 in nodes #check that edge alpha connected to v (node 1)
    else
        @test !(1 in nodes) #check that edge alpha not connected to v (node 1)
    end
    nodes = [n.number for n in net_hybridladder.edge[4].node] #delta's connections
    if move in [0x01, 0x02, 0x04, 0x05]
        @test -4 in nodes #delta connected to u
    else
        @test !(-4 in nodes)
    end
    #check directionality (edge should point toward u, node -4)
    @test PhyloNetworks.getChild(net_hybridladder.edge[5]).number == -4
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    #confirm we're back to original topology 
    @test writeTopology(net_hybridladder) == str_hybridladder
end #of hybrid ladder net edge 5: BR undirected

@testset "hybrid ladder net edge 12: BB undirected (edge below root) move $move" for move in 0x01:0x08
    # BB case, 8 moves. uv edge can contain the root
    undoinfo = PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[12], move);
    nodes = [n.number for n in net_hybridladder.edge[1].node] #alpha
    if move in [0x01, 0x02, 0x06, 0x08]
        @test -3 in nodes #check that edge alpha connected to v
    else
        @test !(-3 in nodes) #check that edge alpha not connected to v
    end
    nodes = [n.number for n in net_hybridladder.edge[11].node] #delta's connections
    if move in [0x01, 0x03, 0x07, 0x08]
        @test -2 in nodes #delta connected to u
    else
        @test !(-2 in nodes)
    end
    #check directionality
    if move in [0x01, 0x02, 0x03, 0x04]
        @test PhyloNetworks.getChild(net_hybridladder.edge[12]).number == -2 #TODO check which direction should this point?
    else
        @test PhyloNetworks.getChild(net_hybridladder.edge[12]).number == -3
    end
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    #confirm we're back to original topology 
    @test writeTopology(net_hybridladder) == str_hybridladder
end #of hybrid ladder net edge 12: BB undirected (edge below root)

@testset "hybrid ladder net edge 1: BR undirected move $move" for move in 0x01:0x03
    # BR case, 6 moves. uv edge can contain the root
    #if there is a path from beta to gamma, resulting graph not a DAG
    #there isnt, but there is a path from alpha to gamma and delta
    #moves that switch alpha and beta (4,5,6) should cause problems #TODO check
    undoinfo = PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[1], move);
    nodes = [n.number for n in net_hybridladder.edge[12].node] #alpha
    if move == 0x01 || (move == 0x03 && PhyloNetworks.getChild(net_hybridladder.edge[12]).number == -2)
        @test 1 in nodes #check that edge alpha connected to v
    else
        @test !(1 in nodes) #check that edge alpha not connected to v
    end
    nodes = [n.number for n in net_hybridladder.edge[4].node] #delta's connections
    if move in [0x01, 0x02]
        @test -2 in nodes #delta connected to u
    else
        @test !(-2 in nodes)
    end
    #check directionality (edge should point toward u, node -2)
    @test PhyloNetworks.getChild(net_hybridladder.edge[1]).number == -2
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    #confirm we're back to original topology 
    @test writeTopology(net_hybridladder) == str_hybridladder
end #of hybrid ladder net edge 1: BR undirected

@testset "hybrid ladder net edge 1: DAG errors BR undirected move $move" for move in 0x04:0x06
    # BR case, 6 moves. uv edge can contain the root
    #if there is a path from beta to gamma, resulting graph not a DAG
    #there isnt, but there is a path from alpha to gamma and delta
    #moves that switch alpha and beta (4,5,6) should cause problems #TODO fix
    @test isnothing(PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[1], move))
end

#for edge 17, u is connected in a 4 cycle with gamma
#if gamma connected to u, would create a 3 cycle
#TODO our checks don't catch this type of case yet
@testset "net_level1 no3cyle check edge 17: BB directed move 1" begin
    @test isnothing(PhyloNetworks.nni!(net_level1, net_level1.edge[17], 0x01))
end #of net_level1 no3cyle check edge 17: BB directed

@testset "net_level1 no3cyle check edge 13: BR directed moves 1" begin
    @test isnothing(PhyloNetworks.nni!(net_level1, net_level1.edge[13], 0x01))
end #of net_level1 no3cyle check edge 17: BB directed

@testset "test isdescendant and isconnected functions" begin
    @test isdescendant(-9, -6) == true
    @test isdescendant(-9, -4) == false
    @test isconnected(-7, -6) == true
    @test isconnect(-7, -3) == false
end

end # of testset on unconstrained NNIs


myconstraint = PhyloNetworks.TopologyConstraint(0x03,Set("8"),net_level1)
undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[1], myconstraint)