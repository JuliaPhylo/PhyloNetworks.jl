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
    #assign labels
    # uv = net_level1.edge[3]
    # u = PhyloNetworks.getParent(uv)
    # v = PhyloNetworks.getChild(uv)
    # labs = [PhyloNetworks.edgerelation(e, u, uv) for e in u.edge]
    # pi = findfirst(isequal(:parent), labs)
    # αu = u.edge[pi]
    # ci = findall(isequal(:child), labs)  # vector. there should be 1 or 2
    # βu = u.edge[ci[1]] #edge 19
    # labs = [PhyloNetworks.edgerelation(e, v, uv) for e in v.edge]
    # ci = findall(isequal(:child), labs)
    # vγ = v.edge[ci[1]] # γ = getChild(vγ)
    # vδ = v.edge[ci[2]]
    #do move
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
end

@testset "edge 16: BB directed, move $move" for move in 0x01:0x02
    # e not hybrid, tree parent:  BB case, 2 NNIs if directed, 8 if undirected
    undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[16], move);
    #check beta connected to u node (number -11)
    nodes = []
    for n in net_level1.edge[13].node
        push!(nodes, n.number)
    end
    @test -11 in nodes
    nodes = []
    for n in net_level1.edge[14].node
        push!(nodes, n.number)
    end
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
end
@test_throws Exception PhyloNetworks.nni!(net_level1, net_level1.edge[16], 0x03);

#TODO check this looks like it labels 13 as βu because βu = u.edge[ci[1]] = 13
@testset "edge 13: BR directed, move $move" for move in 0x01:0x03
    # e.hybrid and tree parent:   BR case, 3 because e cannot contain the root
    undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[13], move);
    #test that move was made
    #location of v node
    if move == 1 || (move == 3 && PhyloNetworks.getChild(net_level1.edge[17]) == -11)
        #check that edge alpha connected to v
        nodes = []
        for n in net_level1.edge[17].node
            push!(nodes, n.number)
        end
        @test 8 in nodes
    elseif move == 2 || (move == 3 && PhyloNetworks.getChild(net_level1.edge[16]) == -11) 
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
    @test PhyloNetworks.getChild(net_level1.edge[13]).number == 8
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    #confirm we're back to original topology 
    @test writeTopology(net_level1) == str_level1
end

#hybrid at root
@testset "net_hybridbelowroot edge 20: BR undirected, move $move" for move in 0x01:0x06
    #in this case, α->u
    undoinfo = PhyloNetworks.nni!(net_hybridbelowroot, net_hybridbelowroot.edge[20], move);
    #test that move was made
    #connections to gamma: 
    nodes = []
    for n in net_hybridbelowroot.edge[19].node
        push!(nodes, n.number)
    end
    if move in [1, 2, 4, 5] #should be connected to v
        @test 11 in nodes
    elseif move in [3, 6] #should be connected to u 
        #both α->u and β->u cases
        @test -12 in nodes
    end
    #connections to alpha (edge 22)
    nodes = []
    for n in net_hybridbelowroot.edge[22].node
        push!(nodes, n.number)
    end
    if move in [1, 3, 5] #should be connected to v 
        @test 11 in nodes
    else #should be connected to u (for 6 too since α->u)
        @test -12 in nodes
    end
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    #confirm we're back to original topology 
    @test writeTopology(net_hybridbelowroot) == str_hybridbelowroot
end

#hybrid at root
@testset "net_hybridbelowroot edge 22: BB undirected, move $move" for move in 0x01:0x06
    #in this case, α->u
    undoinfo = PhyloNetworks.nni!(net_hybridbelowroot, net_hybridbelowroot.edge[22], move);
    #test that move was made
    #connections to gamma
    nodes = []
    for n in net_hybridbelowroot.edge[20].node
        push!(nodes, n.number)
    end
    if move in [1, 3, 7, 8] #gamma connected to v
        @test -12 in nodes
    else #gamma connected to u
        @test -2 in nodes
    end
    #connections to alpha
    nodes = []
    for n in net_hybridbelowroot.edge[3].node
        push!(nodes, n.number)
    end
    if move in [1, 2, 6, 8] #alpha connected to v
        @test -12 in nodes
    else #alpha connected to u
        @test -2 in nodes
    end
    #undo move
    PhyloNetworks.nni!(undoinfo...);
    #confirm we're back to original topology 
    @test writeTopology(net_hybridbelowroot) == str_hybridbelowroot
end


#TODO test no3cycle option by adding a 4 cycle network and trying moves
#net_level1 has a four-cycle at hybrid edge 17
@testset "net_level1 no3cyle check edge 13: BB directed, move 1" for move in 0x01:0x02
    @test_throws Exception PhyloNetworks.nni!(net_level1, net_level1.edge[17], 0x01);
end #of net_level1 has a four-cycle at hybrid edge 17

end # of testset on unconstrained NNIs


myconstraint = PhyloNetworks.TopologyConstraint(0x03,Set("8"),net_level1)
undoinfo = PhyloNetworks.nni!(net_level1, net_level1.edge[1], myconstraint)