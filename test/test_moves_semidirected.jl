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

@testset "level1 edge 3: BB undirected move $move" for move in 0x01:0x08
    undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[3], move);
    #location of v node (number -4)
    nodes = [n.number for n in net_level1.edge[20].node] #α's connections
    if move in [1, 2, 6, 8] #check that edge α connected to v
        @test -4 in nodes
    else
        @test !(-4 in nodes)
    end
    #location of u node (number -3)
    nodes = [n.number for n in net_level1.edge[1].node]
    if move in [2, 4, 5, 6] #check that edge γ connected to u
        @test -3 in nodes
    else
        @test !(-3 in nodes)
    end
    # check directionality
    if move in [3, 4, 5, 7] 
        # keeping α, β or flipping uv keeps node -4 as child of edge 3
        @test PhyloNetworks.getChild(net_level1.edge[3]).number == -4
    else 
        # switching α, β AND flipping uv or doing neither makes node -3 child of edge 3
        @test PhyloNetworks.getChild(net_level1.edge[3]).number == -3
    end
    # undo move
    PhyloNetworks.nni!(undoinfo...);
    # confirm we're back to original topology 
    @test writeTopology(net_level1) == str_level1
end # of level1 edge 3: BB undirected

@testset "level1 test edge 13: BR directed move $move" for move in 0x01:0x03
    # e.hybrid and tree parent:  BR case, 3 moves because e cannot contain the root
    # 3cycle test: α connected to γ so move 1 will create 3 cycles.
    if move == 0x01
        @test isnothing(PhyloNetworks.nni!(net_level1, net_level1.edge[13], move)) # would create a 3cycle
    # move 1 would create a 3 cycle, but should work if we don't forbid 3cycles
        undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[13], move, false); #no3cycle=false
        nodes = [n.number for n in net_level1.edge[17].node]
        @test 8 in nodes # check that edge α connected to v
        nodes = [n.number for n in net_level1.edge[11].node]
        @test !(-11 in nodes)
        @test PhyloNetworks.getChild(net_level1.edge[13]).number == -11
        PhyloNetworks.nni!(undoinfo...); # undo move 
        @test writeTopology(net_level1) == str_level1 # confirm we're back to original topology
    else
        undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[13], move);
        # test that move was made
        nodes = [n.number for n in net_level1.edge[17].node]
        if move == 0x03
            @test 8 in nodes # check that edge α connected to v
        else
            @test !(8 in nodes)
        end
        #location of u node (number -11)
        nodes = [n.number for n in net_level1.edge[11].node]
        if move == 0x03 # check that edge γ connected to u
            @test -11 in nodes
        else # check that edge δ connected to u
            @test !(-11 in nodes)
        end
        # check directionality
        @test PhyloNetworks.getChild(net_level1.edge[13]).number == 8
        PhyloNetworks.nni!(undoinfo...); # undo move
        @test writeTopology(net_level1) == str_level1 # confirm we're back to original topology
    end
end # of level1 edge 13: BR directed

@testset "level1 edge 16: BB directed move $move" for move in 0x01:0x02
    # e not hybrid, tree parent:  BB case, 2 NNIs if directed, 8 if undirected
    undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[16], move);
    nodes = [n.number for n in net_level1.edge[13].node] #β's connections
    @test -11 in nodes #check β connected to u node (number -11)
    nodes = [n.number for n in net_level1.edge[14].node] #γ's connections
    if move == 1 #check that edge γ connected to v
        @test -12 in nodes
    else #check that edge γ connected to u
        @test -11 in nodes
    end
    # check directionality node -11 child of edge 16 in both cases
    @test PhyloNetworks.getChild(net_level1.edge[16]).number == -11
    PhyloNetworks.nni!(undoinfo...); # undo move
    @test writeTopology(net_level1) == str_level1 #confirm we're back to original topology 
end #of level1 edge 16: BB directed

# BB directed has two allowed moves, so 0x03 should throw an exception
@test_throws Exception PhyloNetworks.nni!(net_level1, net_level1.edge[16], 0x03);

@testset "level1 edge 17: BB directed with 4cycle move $move" for move in 0x01:0x02
    # no3cycle test
        # move 2 should fail because of 3cycle problems
        # β is connected in a 4 cycle with γ so a move that makes γ and β a pair
        # would create a 3 cycle 
    if move == 0x01
        undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[17], move);
        nodes = [n.number for n in net_level1.edge[12].node] #β's connections
        @test -6 in nodes #check β connected to u node (number -6)
        nodes = [n.number for n in net_level1.edge[13].node] #γ's connections
        #check that edge γ connected to v
        @test -11 in nodes
        # check directionality node -6 child of edge 17 in both cases
        @test PhyloNetworks.getChild(net_level1.edge[17]).number == -6
        PhyloNetworks.nni!(undoinfo...); # undo move
        @test writeTopology(net_level1) == str_level1 #confirm we're back to original topology 
    else
        @test isnothing(PhyloNetworks.nni!(net_level1, net_level1.edge[17], move)) 
        undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[17], move, false) #should work if we dont check for 3cycles 
        nodes = [n.number for n in net_level1.edge[12].node] #β's connections
        @test -6 in nodes #check β connected to u node (number -6)
        nodes = [n.number for n in net_level1.edge[13].node] #γ's connections
        #check that edge γ connected to u
        @test -6 in nodes
        # check directionality node -11 child of edge 17 in both cases
        @test PhyloNetworks.getChild(net_level1.edge[17]).number == -6
        PhyloNetworks.nni!(undoinfo...); # undo move
        @test writeTopology(net_level1) == str_level1 #confirm we're back to original topology 
    end
end # of level1 edge 17: BB directed

@testset "level1 edge 18: RR (directed) move $move" for move in 0x01:0x04
    # RB case, 4 moves. uv edge cannot contain the root (always directed)
    undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[18], move);
    #test that move was made
    nodes = [n.number for n in net_level1.edge[19].node] #α's connections
    if move in [0x01, 0x03]
        @test -6 in nodes #check that edge α connected to v
    elseif move in [0x02, 0x04]
        @test !(-6 in nodes)
    end
    nodes = [n.number for n in net_level1.edge[12].node] #γ's connections
    if move in [0x03, 0x04] #check that edge γ connected to u (11)
        @test 11 in nodes
    elseif move in [0x01, 0x02]  #check that edge δ connected to u, not γ
        @test !(11 in nodes)
    end
    #check directionality (should point toward u, node 11)
    @test PhyloNetworks.getChild(net_level1.edge[18]).number == 11
    PhyloNetworks.nni!(undoinfo...);
    @test writeTopology(net_level1) == str_level1
end # of level1 edge 18: RR (directed)

@testset "non tree child net edge 3: RB (directed) move $move" for move in 0x01:0x04
    # RB case, 4 moves. uv edge cannot contain the root
    undoinfo = PhyloNetworks.nni!(net_nontreechild, net_nontreechild.edge[3], move);
    #test that move was made
    #location of v node (node -5)
    nodes = [n.number for n in net_nontreechild.edge[4].node] #α's connections
    if move in [0x01, 0x03]
        #check that edge α connected to v
        @test -5 in nodes
    else
        #check that edge β connected to v, not u
        @test !(-5 in nodes)
    end
    nodes = [n.number for n in net_nontreechild.edge[2].node] #δ's connections
    if move in [0x01, 0x02] #check that edge δ connected to u
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

@testset "hybrid ladder net edge 1: BR undirected (u at root, potential nonDAG, 3cycle) move $move" for move in 0x01:0x06
    # BR case, 6 moves. uv edge can contain the root. u at root
    # DAG check: u is root, so if α -> γ, moves 1 and 3 will create a nonDAG
        # α -> γ so moves 1 and 3 will fail. 
        # moves 1 and 3 in the notes correspond to 1, 1' and 3, 3' (1, 4 and 3, 6)
    # 3cycle check: could create 3 cycle (α connected to γ) so moves 1 and 5 forbidden
    if move in [0x01, 0x03, 0x05, 0x06] # DAG check
        @test isnothing(PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[1], move, false))
    end
    if move in [0x01, 0x05] # 3cycle check
        @test isnothing(PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[1], move))
    end
    if move == 0x04
        undoinfo = PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[1], move, false);
        nodes = [n.number for n in net_hybridladder.edge[12].node] # α
        @test !(1 in nodes) # check that edge α is connected to v
        nodes = [n.number for n in net_hybridladder.edge[4].node] # δ's connections
        @test -2 in nodes # check that edge δ is connected to u
        # check directionality (edge should point toward u, node -2)
        @test PhyloNetworks.getChild(net_hybridladder.edge[1]).number == 1
        PhyloNetworks.nni!(undoinfo...);
        @test writeTopology(net_hybridladder) == str_hybridladder
    elseif move == 0x02
        undoinfo = PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[1], move);
        nodes = [n.number for n in net_hybridladder.edge[12].node] # α
        @test !(1 in nodes) # check that edge α not connected to v
        nodes = [n.number for n in net_hybridladder.edge[4].node] # δ's connections
        @test -2 in nodes # δ connected to u
        # check directionality (edge should point toward u, node -2)
        @test PhyloNetworks.getChild(net_hybridladder.edge[1]).number == 1
        PhyloNetworks.nni!(undoinfo...);
        @test writeTopology(net_hybridladder) == str_hybridladder
    end
end # of hybrid ladder net edge 1: BR undirected

@testset "hybrid ladder net edge 4: RR (directed) move $move" for move in 0x01:0x02
    # RR case, 2 moves. uv edge cannot contain the root (always directed)
    undoinfo = PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[4], move);
    #test that move was made
    nodes = [n.number for n in net_hybridladder.edge[5].node] #α's connections
    if move == 0x01 #check that edge α connected to v
        @test 4 in nodes
    else move == 0x02 #check that edge β connected to v, not α
        @test !(4 in nodes)
    end
    #check that edge δ connected to u, not γ
    nodes = [n.number for n in net_hybridladder.edge[3].node]
    @test 1 in nodes
    #check directionality (should point toward u, node 1)
    @test PhyloNetworks.getChild(net_hybridladder.edge[4]).number == 1
    #undo move
    PhyloNetworks.nni!(undoinfo...); 
    #confirm we're back to original topology 
    @test writeTopology(net_hybridladder) == str_hybridladder
end #of hybrid ladder net edge 4: RR (directed)

@testset "hybrid ladder net edge 5: BR undirected move $move" for move in 0x01:0x06
    # BR case, 6 moves. uv edge can contain the root
    # 3 cycle test: α connected to γ, α -> u
    #   moves 1, 5 forbidden 
    # DAG test:
    #   no path from α -> γ or β -> γ so all moves should work
    if move in [0x01, 0x05]
        @test isnothing(PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[5], move)) #3cycle forbidden
        undoinfo = PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[5], move, false) #3cycles allowed
    else
        undoinfo = PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[5], move);
    end
    #tests for all moves:
    nodes = [n.number for n in net_hybridladder.edge[6].node] #α
    if move in [0x01, 0x03, 0x05, 0x06]
        @test 1 in nodes #check that edge α connected to v
    else
        @test !(1 in nodes) #check that edge α not connected to v
    end
    nodes = [n.number for n in net_hybridladder.edge[4].node] #δ's connections
    if move in [0x01, 0x02, 0x04, 0x05]
        @test -4 in nodes #δ connected to u
    else
        @test !(-4 in nodes)
    end
    #check directionality 
    if move in [0x01, 0x05]
        #(edge should point toward u, node -4)
        @test PhyloNetworks.getChild(net_hybridladder.edge[5]).number == -4
    else
        @test PhyloNetworks.getChild(net_hybridladder.edge[5]).number == 1
    end
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    #confirm we're back to original topology 
    @test writeTopology(net_hybridladder) == str_hybridladder
end #of hybrid ladder net edge 5: BR undirected

@testset "hybrid ladder net edge 12: BB undirected (edge below root) move $move" for move in 0x01:0x08
    # BB case, 8 moves. uv edge can contain the root. no flip.
    # no3cycle: moves 1, 4, 5, 8 would create a 3 cycle because α is connected to γ
    if move in [0x01, 0x04, 0x05, 0x08]
        @test isnothing(PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[12], move))
        undoinfo = PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[12], move, false);
    else
        undoinfo = PhyloNetworks.nni!(net_hybridladder, net_hybridladder.edge[12], move);
    end
    nodes = [n.number for n in net_hybridladder.edge[1].node] #α
    if move in [0x01, 0x02, 0x06, 0x08]
        @test -3 in nodes #check that edge α connected to v
    else
        @test !(-3 in nodes) #check that edge α not connected to v
    end
    nodes = [n.number for n in net_hybridladder.edge[11].node] #δ's connections
    if move in [0x01, 0x03, 0x07, 0x08]
        @test -2 in nodes #δ connected to u
    else
        @test !(-2 in nodes)
    end
    #check directionality
    @test PhyloNetworks.getChild(net_hybridladder.edge[12]).number == -3
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    # confirm we're back to original topology 
    @test writeTopology(net_hybridladder) == str_hybridladder
end # of hybrid ladder net edge 12: BB undirected (edge below root)

@testset "test isdescendant and isconnected functions" begin
    net_level1 = readTopology(str_level1);
    @test PhyloNetworks.isdescendant(net_level1.node[7], net_level1.node[17]) == true # nodes -9, -6
    @test PhyloNetworks.isdescendant(net_level1.node[7], net_level1.node[3]) == false # nodes -9, -4
    @test PhyloNetworks.isdescendant(net_level1.node[15], net_level1.node[17]) == true # node -12, -6
    @test PhyloNetworks.isconnected(net_level1.node[12], net_level1.node[17]) == true #node -7, -6
    @test PhyloNetworks.isconnected(net_level1.node[12], net_level1.node[19]) == false # node -7, -3 
end

#TODO add these after adding contraint functionality
# myconstraint = PhyloNetworks.TopologyConstraint(0x03,Set("8"),net_level1)
# undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[1], myconstraint)

end # of testset on unconstrained NNIs

