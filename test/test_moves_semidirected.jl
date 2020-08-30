#= # for local testing, need this:
using Test
using PhyloNetworks
using Random
=#

@testset "unconstrained NNI moves" begin

str_level1 = "(((S8,S9),(((((S1,S2,S3),S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
net_level1 = readTopology(str_level1); # this network has a polytomy at node -9
# same topology as: rootatnode!(net_level1, -3). edges 1:22
str_nontreechild = "((((Ag,E))#H3,(#H1:7.159::0.056,((M:0.0)#H2:::0.996,(Ak,(#H3:0.08,#H2:0.0::0.004):0.023):0.078):2.49):2.214):0.026,((Az:2.13,As:2.027):1.697)#H1:0.0::0.944,Ap);"
net_nontreechild = readTopology(str_nontreechild);
# problem: the plot has an extra vertical segment, for a clade that's not in the major tree
# --> fix that in PhyloPlots (fixit)
str_hybridladder = "(#H2:::0.2,((C,((B)#H1)#H2:::0.8),(#H1,(A1,A2))),O);"
net_hybridladder = readTopology(str_hybridladder);

#=
using PhyloPlots
plot(net_level1, :R, showNodeNumber=true, showEdgeNumber=true)
plot(net_nontreechild, :R, showNodeNumber=true, showEdgeNumber=true)
plot(net_hybridladder, :R, showNodeNumber=true, showEdgeNumber=true)
=#

@test isnothing(nni!(net_level1, net_level1.edge[1], 0x01, true, true)) # external edge

@testset "level1 edge 3: BB undirected move $move" for move in 0x01:0x08
    undoinfo = nni!(net_level1, net_level1.edge[3], move, true, true);
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
    nni!(undoinfo...);
    # confirm we're back to original topology
    @test writeTopology(net_level1) == str_level1
end # of level1 edge 3: BB undirected

@testset "level1 test edge 13: BR directed move $move" for move in 0x01:0x03
    # e.hybrid and tree parent:  BR case, 3 moves because e cannot contain the root
    # 3cycle test: α connected to γ so move 1 will create 3 cycles.
    if move == 0x01
        @test isnothing(nni!(net_level1, net_level1.edge[13], move, true, true)) # would create a 3cycle
    # move 1 would create a 3 cycle, but should work if we don't forbid 3cycles
        undoinfo = nni!(net_level1, net_level1.edge[13], move, true, false); #no3cycle=false
        nodes = [n.number for n in net_level1.edge[17].node]
        @test 8 in nodes # check that edge α connected to v
        nodes = [n.number for n in net_level1.edge[11].node]
        @test !(-11 in nodes)
        @test PhyloNetworks.getChild(net_level1.edge[13]).number == -11
        nni!(undoinfo...); # undo move
        @test writeTopology(net_level1) == "(((S8,S9),(((((S1,S2,S3),S4),(S5:0.0)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
        # confirm we're back to original topology (but now with 0.0 branch lengths)
    else
        undoinfo = nni!(net_level1, net_level1.edge[13], move, true, false);
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
        nni!(undoinfo...); # undo move
        @test writeTopology(net_level1) == "(((S8,S9),(((((S1,S2,S3),S4),(S5:0.0)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"# confirm we're back to original topology
    end
end # of level1 edge 13: BR directed

@testset "level1 edge 16: BB directed move $move" for move in 0x01:0x02
    # e not hybrid, tree parent:  BB case, 2 NNIs if directed, 8 if undirected
    undoinfo = nni!(net_level1, net_level1.edge[16], move, true, true);
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
    nni!(undoinfo...); # undo move
    @test writeTopology(net_level1) == "(((S8,S9),(((((S1,S2,S3),S4),(S5:0.0)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
        # confirm we're back to original topology
end #of level1 edge 16: BB directed

# BB directed has two allowed moves, so 0x03 should throw an exception
@test_throws Exception nni!(net_level1, net_level1.edge[16], 0x03);

@testset "level1 edge 17: BB directed with 4cycle move $move" for move in 0x01:0x02
    # no3cycle test
        # move 2 should fail because of 3cycle problems
        # β is connected in a 4 cycle with γ so a move that makes γ and β a pair
        # would create a 3 cycle
    if move == 0x01
        undoinfo = nni!(net_level1, net_level1.edge[17], move, true, true);
        nodes = [n.number for n in net_level1.edge[12].node] #β's connections
        @test -6 in nodes #check β connected to u node (number -6)
        nodes = [n.number for n in net_level1.edge[13].node] #γ's connections
        #check that edge γ connected to v
        @test -11 in nodes
        # check directionality node -6 child of edge 17 in both cases
        @test PhyloNetworks.getChild(net_level1.edge[17]).number == -6
        nni!(undoinfo...); # undo move
        @test writeTopology(net_level1) == "(((S8,S9),(((((S1,S2,S3),S4),(S5:0.0)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
    else
        @test isnothing(nni!(net_level1, net_level1.edge[17], move, true, true))
        undoinfo = nni!(net_level1, net_level1.edge[17], move, true, false) #should work if we dont check for 3cycles
        nodes = [n.number for n in net_level1.edge[12].node] #β's connections
        @test -6 in nodes #check β connected to u node (number -6)
        nodes = [n.number for n in net_level1.edge[13].node] #γ's connections
        #check that edge γ connected to u
        @test -6 in nodes
        # check directionality node -11 child of edge 17 in both cases
        @test PhyloNetworks.getChild(net_level1.edge[17]).number == -6
        nni!(undoinfo...); # undo move
        @test writeTopology(net_level1) == "(((S8,S9),(((((S1,S2,S3),S4),(S5:0.0)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
    end
end # of level1 edge 17: BB directed

@testset "level1 edge 18: RR (directed) move $move" for move in 0x01:0x04
    # RB case, 4 moves. uv edge cannot contain the root (always directed)
    undoinfo = nni!(net_level1, net_level1.edge[18], move, true, true);
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
    nni!(undoinfo...);
        # edge constrained at 0.0 now
    @test writeTopology(net_level1) == "(((S8,S9),(((((S1,S2,S3),S4),(S5:0.0)#H1),(#H1,(S6,S7))):0.0)#H2),(#H2,S10));"
end # of level1 edge 18: RR (directed)

@testset "non tree child net edge 3: RB (directed) move $move" for move in 0x01:0x04
    # RB case, 4 moves. uv edge cannot contain the root
    undoinfo = nni!(net_nontreechild, net_nontreechild.edge[3], move, true, true);
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
    nni!(undoinfo...);
        # keep constrained edges at 0.0, but otherwise topology completely restored
    @test writeTopology(net_nontreechild) == "((((Ag,E):0.0)#H3,(#H1:7.159::0.056,((M:0.0)#H2:::0.996,(Ak,(#H3:0.08,#H2:0.0::0.004):0.023):0.078):2.49):2.214):0.026,((Az:2.13,As:2.027):1.697)#H1:0.0::0.944,Ap);"
end #of non tree child net edge 5: RB (directed)

@testset "hybrid ladder net edge 1: BR undirected (u at root, potential nonDAG, 3cycle) move $move" for move in 0x01:0x06
    # BR case, 6 moves. uv edge can contain the root. u at root
    # DAG check: u is root, so if α -> γ, moves 1 and 3 will create a nonDAG
        # α -> γ so moves 1 and 3 will fail.
        # moves 1 and 3 in the notes correspond to 1, 1' and 3, 3' (1, 4 and 3, 6)
    # 3cycle check: could create 3 cycle (α connected to γ) so moves 1 and 5 forbidden
    if move in [0x01, 0x03, 0x05, 0x06] # DAG check
        @test isnothing(nni!(net_hybridladder, net_hybridladder.edge[1], move, false, true))
    end
    if move in [0x01, 0x05] # 3cycle check
        @test isnothing(nni!(net_hybridladder, net_hybridladder.edge[1], move, true, true))
    end
    if move == 0x04
        undoinfo = nni!(net_hybridladder, net_hybridladder.edge[1], move, false, true);
        nodes = [n.number for n in net_hybridladder.edge[12].node] # α
        @test !(1 in nodes) # check that edge α is connected to v
        nodes = [n.number for n in net_hybridladder.edge[4].node] # δ's connections
        @test -2 in nodes # check that edge δ is connected to u
        # check directionality (edge should point toward u, node -2)
        @test PhyloNetworks.getChild(net_hybridladder.edge[1]).number == 1
        nni!(undoinfo...);
        @test writeTopology(net_hybridladder)== "(#H2:::0.2,((C,((B)#H1:0.0)#H2:::0.8),(#H1,(A1,A2))),O);" # restored but edge below hybrid node constrained at 0.0
    elseif move == 0x02
        undoinfo = nni!(net_hybridladder, net_hybridladder.edge[1], move, true, true);
        nodes = [n.number for n in net_hybridladder.edge[12].node] # α
        @test !(1 in nodes) # check that edge α not connected to v
        nodes = [n.number for n in net_hybridladder.edge[4].node] # δ's connections
        @test -2 in nodes # δ connected to u
        # check directionality (edge should point toward u, node -2)
        @test PhyloNetworks.getChild(net_hybridladder.edge[1]).number == 1
        nni!(undoinfo...);
        @test writeTopology(net_hybridladder) == "(#H2:::0.2,((C,((B)#H1:0.0)#H2:::0.8),(#H1,(A1,A2))),O);" # restored but edge below hybrid node constrained at 0.0
    end
end # of hybrid ladder net edge 1: BR undirected

@testset "hybrid ladder net edge 4: RR (directed) move $move" for move in 0x01:0x02
    # RR case, 2 moves. uv edge cannot contain the root (always directed)
    undoinfo = nni!(net_hybridladder, net_hybridladder.edge[4], move, false, true);
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
    nni!(undoinfo...);
    @test writeTopology(net_hybridladder) == "(#H2:::0.2,((C,((B)#H1:0.0)#H2:::0.8),(#H1,(A1,A2))),O);" # restored but edge below hybrid node constrained at 0.0
end #of hybrid ladder net edge 4: RR (directed)

@testset "hybrid ladder net edge 5: BR undirected move $move" for move in 0x01:0x06
    # BR case, 6 moves. uv edge can contain the root
    # 3 cycle test: α connected to γ, α -> u
    #   moves 1, 5 forbidden 
    # DAG test:
    #   no path from α -> γ or β -> γ so all moves should work
    if move in [0x01, 0x05]
        @test isnothing(nni!(net_hybridladder, net_hybridladder.edge[5], move, false, true)) # 3-cycles forbidden
        undoinfo = nni!(net_hybridladder, net_hybridladder.edge[5], move, false, false) # 3-cycles allowed
    else
        undoinfo = nni!(net_hybridladder, net_hybridladder.edge[5], move, false, true);
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
    nni!(undoinfo...);
    @test writeTopology(net_hybridladder) == "(#H2:::0.2,((C,((B)#H1:0.0)#H2:::0.8),(#H1,(A1,A2))),O);" # restored but edge below hybrid node constrained at 0.0
end #of hybrid ladder net edge 5: BR undirected

@testset "hybrid ladder net edge 12: BB undirected (edge below root) move $move" for move in 0x01:0x08
    # BB case, 8 moves. uv edge can contain the root. no flip.
    # no3cycle: moves 1, 4, 5, 8 would create a 3 cycle because α is connected to γ
    if move in [0x01, 0x04, 0x05, 0x08]
        @test isnothing(nni!(net_hybridladder, net_hybridladder.edge[12], move, true, true))
        undoinfo = nni!(net_hybridladder, net_hybridladder.edge[12], move, true, false);
    else
        undoinfo = nni!(net_hybridladder, net_hybridladder.edge[12], move, true, true);
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
    nni!(undoinfo...);
    @test writeTopology(net_hybridladder) == "(#H2:::0.2,((C,((B)#H1:0.0)#H2:::0.8),(#H1,(A1,A2))),O);" # restored but edge below hybrid node constrained at 0.0
end # of hybrid ladder net edge 12: BB undirected (edge below root)

@testset "test isdescendant and isconnected functions" begin
    net_level1 = readTopology(str_level1);
    @test  PhyloNetworks.isdescendant(net_level1.node[7], net_level1.node[17])  # nodes -9, -6
    @test !PhyloNetworks.isdescendant(net_level1.node[7], net_level1.node[3])   # nodes -9, -4
    @test  PhyloNetworks.isdescendant(net_level1.node[15], net_level1.node[17]) # nodes -12, -6
    @test !PhyloNetworks.isdescendant(net_level1.node[12], net_level1.node[12])
    @test  PhyloNetworks.isconnected(net_level1.node[12], net_level1.node[17])  # nodes -7, -6
    @test !PhyloNetworks.isconnected(net_level1.node[12], net_level1.node[19])  # nodes -7, -3
    # mess up the direction of some tree edges, then check descendence relationships with isdescendant_undirected
    for i in [4,5,6,7,9,10,12,17,3,20] net_level1.edge[i].isChild1 = !net_level1.edge[i].isChild1; end
    @test  PhyloNetworks.isdescendant_undirected(net_level1.node[7], net_level1.node[17], net_level1.edge[18])
    @test !PhyloNetworks.isdescendant_undirected(net_level1.node[7], net_level1.node[3], net_level1.edge[3])
    @test  PhyloNetworks.isdescendant_undirected(net_level1.node[15], net_level1.node[17], net_level1.edge[18])
    @test !PhyloNetworks.isdescendant_undirected(net_level1.node[12], net_level1.node[12], net_level1.edge[12])
end

end # of testset on unconstrained NNIs

@testset "constrained NNI moves" begin
# subsets for: species constraints; move root (species & clade constraints);
#              clade constraints (not yet: TODO)

@testset "species constraints" begin # multiple individuals from each species
str_level1_s = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));" # indviduals S1A S1B S1C go on leaf 1
net_level1_s = readTopology(str_level1_s)
#=
plot(net_level1_s, :R, showNodeNumber=true, showEdgeNumber=true)
=#

# test breakedge! function
newnode, newedge = PhyloNetworks.breakedge!(net_level1_s.edge[4], net_level1_s);
@test length(net_level1_s.node) == 20
@test length(net_level1_s.edge) == 21
@test newnode.edge[1].number == 4
@test newnode.edge[2].number == 21
@test PhyloNetworks.getParent(net_level1_s.edge[4]) === newnode
@test PhyloNetworks.getChild(newedge) === newnode

# test addleaf! function
net_level1_s = readTopology(str_level1_s)
PhyloNetworks.addleaf!(net_level1_s, net_level1_s.node[4], "S1A");
@test !net_level1_s.node[findfirst([n.number == 3 for n in net_level1_s.node])].leaf
PhyloNetworks.addleaf!(net_level1_s, net_level1_s.node[4], "S1B");
PhyloNetworks.addleaf!(net_level1_s, net_level1_s.node[4], "S1C");
@test net_level1_s.edge[21].containRoot == false # check containRoot on edge 4 and exterior edges
@test net_level1_s.edge[22].containRoot == false
@test PhyloNetworks.getChild(net_level1_s.edge[21]).name == "S1A"
@test PhyloNetworks.getChild(net_level1_s.edge[22]).name == "S1B"
# test addleaf! on edge
net_level1_s = readTopology(str_level1_s)
PhyloNetworks.addleaf!(net_level1_s, net_level1_s.edge[4], "S1A");
@test length(net_level1_s.node) == 21
@test net_level1_s.node[21].leaf

# test addindividuals! function
net_level1_s = readTopology(str_level1_s)
PhyloNetworks.addindividuals!(net_level1_s, "S1", ["S1A", "S1B", "S1C"])
@test !net_level1_s.node[findfirst([n.number == 3 for n in net_level1_s.node])].leaf
@test length(net_level1_s.node[findfirst([n.number == 3 for n in net_level1_s.node])].edge) == 4
# spaces in name
net_level1_s = readTopology(str_level1_s)
@test_logs (:warn, r"^species S 1 not") PhyloNetworks.addindividuals!(net_level1_s, "S 1", ["S1A", "S1B", "S1C"])
@test writeTopology(net_level1_s) == str_level1_s # network unchanged
@test_logs (:warn, r"^Spaces in \"S1 A\" may cause errors") PhyloNetworks.addindividuals!(net_level1_s, "S1", ["S1 A", "S1B", "S1C"])
@test writeTopology(net_level1_s) == "(((S8,S9),(((((S1_A,S1B,S1C)S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
# # test mapindividuals! function
net_level1_s = readTopology(str_level1_s)
# in net env
filename = joinpath(@__DIR__, "..","examples","mappingIndividuals.csv")
# filename = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "mappingIndividuals.csv"))
net_level1_i, c_species = PhyloNetworks.mapindividuals(net_level1_s, filename)
@test string(c_species[1]) == "Species constraint, on tips: S1A, S1B, S1C\n stem edge number 4\n crown node number 3"
@test c_species[1].taxonnames == ["S1A","S1B","S1C"]
@test c_species[1].taxonnums == Set([11,12,13])
@test writeTopology(net_level1_i) == "(((S8,S9),(((((S1A,S1B,S1C)S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"

# test clade constraint contructor
c_clade = PhyloNetworks.TopologyConstraint(0x02, ["S1A","S1B","S1C","S4"], net_level1_i)
@test string(c_clade) == "Clade constraint, on tips: S1A, S1B, S1C, S4\n stem edge number 6\n crown node number -8"

# test errors in species constructor
@test_throws ErrorException PhyloNetworks.TopologyConstraint(0x02, ["S1A"], net_level1_i) # only 1 tip
@test_throws ErrorException PhyloNetworks.TopologyConstraint(0x02, ["S1A", "TypoTaxa"], net_level1_i) # typo
@test_throws ErrorException PhyloNetworks.TopologyConstraint(0x02, ["S1A", "S8"], net_level1_i) # not a clade

# NNIs under species constraints
Random.seed!(1234);
# no nni on stem edge for species example
@test isnothing(nni!(net_level1_i , net_level1_i.edge[4], true, true, c_species))
@testset "NNI, 1 species constraint, net level 1, edge $ei" for ei in [8,3,9]
    # 8: BR directed, 3: BB undirected, 9: BB directed, 15: RB directed
    # note: there are no cases of RR directed in net_level1_i
    undoinfo = nni!(net_level1_i , net_level1_i.edge[ei], true, true, c_species);
    @test undoinfo !== nothing
    nni!(undoinfo...);
        # restored to original network, except that edges below hybrid nodes will now have length 0.0
    @test writeTopology(net_level1_i) == "(((S8,S9),(((((S1A,S1B,S1C)S1,S4),(S5:0.0)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
end
# TODO: BR case edge 8: nni move 3 causes problems. hybrid node 6 has 0 or 2+ major hybrid parents
end # of species constraints

@testset "test move root & constraint checking under species & clade constraints" begin
# "(((S8,S9),(((((S1A,S1B,S1C)S1,S4),#H1),((S5)#H1,(S6,S7))))#H2),(#H2,S10));"
netl1_i = readTopology("(((((S1A,S1B,S1C)S1,S4),#H1),((S5)#H1,(S6,S7))));")
con = [PhyloNetworks.TopologyConstraint(0x01, ["S1A","S1B","S1C"], netl1_i),
       PhyloNetworks.TopologyConstraint(0x02, ["S5","S6","S7"], netl1_i)]
Random.seed!(765);
@test PhyloNetworks.moveroot!(netl1_i, con) # only 2 options
writeTopology(netl1_i) == "(((S1A,S1B,S1C)S1,S4),#H1,(((S5)#H1,(S6,S7))));" # now unrooted
@test PhyloNetworks.moveroot!(netl1_i, con) # only 1 option
writeTopology(netl1_i) == "((S1A,S1B,S1C)S1,S4,(#H1,(((S5)#H1,(S6,S7)))));"
netl1_i.root = 14 # back to original rooted network. This node is still of degree 2
@test !PhyloNetworks.checkspeciesnetwork!(netl1_i, con) # false: root *at* clade crown
@test netl1_i.root == 13 # now unrooted (via removedegree2nodes!), root was moved, con[2] stem edge was deleted too...
netl1_i.root = 7; directEdges!(netl1_i) # move root strictly above clade crown
con[2] = PhyloNetworks.TopologyConstraint(0x02, ["S5","S6","S7"], netl1_i)
@test PhyloNetworks.checkspeciesnetwork!(netl1_i, con) # now fine: root *above* clade crown
undoinfo = nni!(netl1_i,netl1_i.edge[8],0x01,false,false);
@test !PhyloNetworks.checkspeciesnetwork!(netl1_i, con)
nni!(undoinfo...);
@test PhyloNetworks.checkspeciesnetwork!(netl1_i, con)
undoinfo = nni!(netl1_i,netl1_i.edge[8],0x03,false,false) # creates a 2-cycle
@test netl1_i.numEdges == 13
PhyloNetworks.deletehybridedge!(netl1_i, netl1_i.edge[10])
@test netl1_i.numEdges == 10 # 2-cycle removed
netl1_i = readTopology("(((S1A,S1B,S1C),S4),#H1,((S5)#H1,(S6,S7)));")
undoinfo = nni!(netl1_i,netl1_i.edge[12],0x03,false,false) # 4-cycle now
@test nni!(netl1_i,netl1_i.edge[12],0x02,true,true) === nothing # would create a 3-cycle
end

#=
str_level1 = "(((S8,S9),(((((S1,S2,S3),S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
net_level1 = readTopology(str_level1); #polytomy at node -9 for leaves 3, 4, 5

str_nontreechild = "((((Ag,E))#H3,(#H1:7.159::0.056,((M:0.0)#H2:::0.996,(Ak,(#H3:0.08,#H2:0.0::0.004):0.023):0.078):2.49):2.214):0.026,((Az:2.13,As:2.027):1.697)#H1:0.0::0.944,Ap);"
net_nontreechild = readTopology(str_nontreechild);

str_polytomy_species = "(((S8,S9),(((((S1,S2,S3),S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
net_species = readTopology(str_polytomy_species);
=#

#=
@testset "clade contraints" begin
c_nontree_cladebelowhybrid = PhyloNetworks.TopologyConstraint(0x02, ["Az", "As"], net_nontreechild)
c_level1_clade = PhyloNetworks.TopologyConstraint(0x02, ["S1", "S2", "S3"], net_level1)
# this test below returns nothing but shouldn't
nni!(net_nontreechild, net_nontreechild.edge[20], true, [c_nontree_cladebelowhybrid])
# this test errors
@test isnothing(nni!(net_nontreechild, net_nontreechild.edge[18], true, [c_nontree_cladebelowhybrid]))

@test PhyloNetworks.checkspeciesnetwork!(net_nontreechild, [c_nontree_cladebelowhybrid])

PhyloNetworks.addindividuals!(net_level1_s, "S1", ["S1A", "S1B", "S1C"])
@test_throws ErrorException PhyloNetworks.checkspeciesnetwork!(net_level1_i, [c_species])
@test PhyloNetworks.cladesviolated(net_level1_i, c_species)
@test_throws ErrorException PhyloNetworks.checkspeciesnetwork!(net_level1_s, [c_level1_species])
@test PhyloNetworks.checknetwork(net_level1_i, c_species) #TODO will want to remove if we remove this function
end # of testset on checknetwork functions for species constraints

=#

end # of constrained NNI moves
