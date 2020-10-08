# using PhyloNetworks
# using PhyloPlots
# using Test
# using Debugger

@testset "tree component" begin

    treestr = "(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);"
    tree = readTopology(treestr)
    rcomp = updateroot!(tree, tree_edge_components(tree))
    @test length(rcomp) == 7
    @test all([e.containRoot for e = tree.edge])

    netstr = "(#H1:::0.1,#H2:::0.2,(((b)#H1)#H2,a));"
    net = readTopology(netstr)
    rcomp = updateroot!(net, tree_edge_components(net))
    # plot(net, :R, showNodeNumber=true, showEdgeNumber=true);
    @test Set(n.number for n in rcomp) == Set([-2 -3 4])
    @test Set(e.number for e in net.edge if e.containRoot) == Set([7 6 5 2 1])

    # test semidirected cycle case
    net.edge[1].isChild1 = !net.edge[1].isChild1
    net.edge[7].isChild1 = !net.edge[7].isChild1
    net.edge[7].hybrid = true
    net.edge[7].gamma  = net.edge[4].gamma
    net.edge[4].gamma  = 1.0
    net.edge[4].hybrid = false
    net.root = 5
    net.node[6].hybrid = true
    net.node[6].name = "H3"
    net.node[2].hybrid = false
    net.hybrid[1] = net.node[6]
    mem = tree_edge_components(net)
    @test_throws PhyloNetworks.RootMismatch updateroot!(net, mem)
    try updateroot!(net, mem)
    catch e
        @test occursin("Semidirected cycle", sprint(showerror, e))
    end

    # test multiple entry points case
    str_level1 = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
    netl1 = readTopology(str_level1)
    # TODO rename P
    P = PhyloNetworks
    root = netl1.node[P.getIndexNode(-2, netl1)]
    e1 = netl1.edge[P.getIndexEdge(20, netl1)]
    e2 = netl1.edge[P.getIndexEdge(17, netl1)]
    P.deleteEdge!(netl1, e1)
    P.deleteEdge!(netl1, e2)
    P.deleteNode!(netl1, root)
    P.removeEdge!(netl1.node[P.getIndexNode(-12, netl1)], e1)
    P.removeEdge!(netl1.node[P.getIndexNode(-3, netl1)], e2)
    mem = tree_edge_components(netl1)
    @test_throws PhyloNetworks.RootMismatch updateroot!(netl1, mem)
    try updateroot!(netl1, mem)
    catch e
        @test occursin("common ancestor", sprint(showerror, e))
    end

    # test undirected cycle case
    netl1 = readTopology(str_level1)
    n1 = netl1.node[P.getIndexNode(-6, netl1)]
    n2 = netl1.node[P.getIndexNode(-8, netl1)]
    e = P.Edge(length(netl1.edge)+1,1.0,false,1.0)
    P.setNode!(e, n1)
    P.setNode!(e, n2)
    push!(n1.edge, e)
    push!(n2.edge, e)
    push!(net.edge, e)
    @test_throws PhyloNetworks.RootMismatch tree_edge_components(netl1)
    try tree_edge_components(netl1)
    catch e
        @test occursin("Undirected cycle", sprint(showerror, e))
    end

end
