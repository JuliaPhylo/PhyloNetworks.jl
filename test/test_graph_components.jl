# using Test, PhyloNetworks
# using PhyloPlots
# using Debugger

@testset "tree component" begin

    treestr = "(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);"
    tree = readTopology(treestr)
    for e in tree.edge e.containRoot=false; end # wrong, on purpose
    @test collect(values(treeedgecomponents(tree))) == repeat([1], inner=7)
    rcomp = checkroot!(tree)
    @test rcomp == 1
    @test all([e.containRoot for e = tree.edge])

    netstr = "(#H1:::0.1,#H2:::0.2,(((b)#H1)#H2,a));"
    net = readTopology(netstr)
    for e in net.edge e.containRoot=false; end # wrong, on purpose
    node2comp = treeedgecomponents(net) # e.g. [1,1,1,2,2,3] or [2,2,2,1,1,3]
    compsize = [count(isequal(i), values(node2comp)) for i in 1:3]
    @test sort(compsize) == [1,2,3]
    rcompID = checkroot!(net, node2comp)
    @test compsize[rcompID] == 3
    # plot(net, :R, showNodeNumber=true, showEdgeNumber=true);
    rcomp = keys(filter(p -> p.second == rcompID, node2comp))
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
    mem = treeedgecomponents(net)
    nnum2comp = Dict(n.number => uc for (n,uc) in mem)
    @test nnum2comp[4] == nnum2comp[-3]
    @test nnum2comp[1] == nnum2comp[2] == nnum2comp[3]
    @test length(unique(nnum2comp[nn] for nn in [4,1,-2])) == 3
    @test_throws PhyloNetworks.RootMismatch checkroot!(net, mem)
    try checkroot!(net, mem)
    catch e
        @test occursin("Semidirected cycle", sprint(showerror, e))
    end

    # test multiple entry points case
    str_level1 = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
    netl1 = readTopology(str_level1)
    P = PhyloNetworks # binding P: local to test set
    root = netl1.node[19] # 19 = findfirst(n -> n.number == -2, netl1.node)
    e1 = netl1.edge[20]   # 20 = findfirst(e -> e.number == 20, netl1.edge)
    e2 = netl1.edge[17]   # 17 = findfirst(e -> e.number == 17, netl1.edge)
    P.deleteEdge!(netl1, e1)
    P.deleteEdge!(netl1, e2)
    P.deleteNode!(netl1, root)
    P.removeEdge!(netl1.node[18], e1) # 18 = P.getIndexNode(-12, netl1)
    P.removeEdge!(netl1.node[16], e2) # 16 = P.getIndexNode(-3, netl1)
    # nodenumber2UC = Dict(n.number => uc for (n,uc) in treeedgecomponents(netl1))
    mem = treeedgecomponents(netl1)
    @test_throws PhyloNetworks.RootMismatch checkroot!(netl1, mem)
    try checkroot!(netl1, mem)
    catch e
        @test occursin("no common ancestor", sprint(showerror, e))
    end

    # test undirected cycle case
    netl1 = readTopology(str_level1)
    n1 = netl1.node[14] # 14 = P.getIndexNode(-6, netl1)
    n2 = netl1.node[6]  #  6 = P.getIndexNode(-8, netl1)
    e = P.Edge(21,1.0,false,1.0) # 21 = length(netl1.edge) + 1
    P.setNode!(e, n1)
    P.setNode!(e, n2)
    push!(n1.edge, e)
    push!(n2.edge, e)
    push!(net.edge, e)
    @test_throws PhyloNetworks.RootMismatch treeedgecomponents(netl1)
    try treeedgecomponents(netl1)
    catch e
        @test occursin("Undirected cycle", sprint(showerror, e))
    end

end
