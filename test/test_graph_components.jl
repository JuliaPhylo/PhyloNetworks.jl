# using Test, PhyloNetworks
# using PhyloPlots
# using Debugger

@testset "Tarjan's biconnected components" begin

    net = readnewick("(A,(B,(C,D)));");
    a = biconnectedcomponents(net);
    @test [[e.number for e in b] for b in a] == [[1],[2],[3],[4],[5],[6],]
    net = readnewick("(((A,(((C,(D)#H2),(E,#H2)))#H1),(B,#H1)),F);");
    a = biconnectedcomponents(net);
    @test [[e.number for e in b] for b in a] == [[1],[2],[3],[6],
        [8, 7, 4, 5],[9],[12],[14, 13, 10, 11],[15],[16]]
    net = readnewick("(((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2));");
    a = biconnectedcomponents(net);
    @test [[e.number for e in b] for b in a] == [[1],
        [2],[5],[6],[12],[10, 14, 13, 7, 8, 9, 3, 4, 11]]
    net = readnewick("((((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2)),(((F)#H3,G),(H,#H3)));");
    a = biconnectedcomponents(net);
    @test [[e.number for e in b] for b in a] == [[1],[2],[5],[6],[12],
        [10, 14, 13, 7, 8, 9, 3, 4, 11],[15],[16],[20],[18],
        [22, 21, 17, 19],[23]]
    a = biconnectedcomponents(net, true);
    @test [[e.number for e in b] for b in a] == [[10, 14, 13, 7, 8, 9, 3, 4, 11],
        [22, 21, 17, 19]]
    # net with hybrid ladder; 3 degree-2 nodes; root edge above LSA
    net = readnewick("((((((((((((Ae_caudata))#H1,#H2),(Ae_umbellulata,#H1)),Ae_comosa),((((Ae_searsii)#H2,#H3),#H4)))),(((Ae_speltoides_Tr223,Ae_speltoides_Tr251))#H3,(Ae_mutica)#H4))),S_vavilovii)));")
    a = biconnectedcomponents(net, false);
    a = [sort!([e.number for e in b]) for b in a]
    @test length(a) == 14
    @test a[14] == [31] # root edge
    @test a[10] == [3,4,5, 7,8,9, 11, 13,14,15,16,17,18,19,20, 24, 26,27]

    net = readnewick("((((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2)),(((F)#H3,G),(H,#H3)));");
    r,major,minor = PhyloNetworks.blobinfo(net, false);
    @test [n.number for n in r] == [-5, 3, -8, 6, -10, -3, -2, 9, -14, -12, -11, -2]
    r,major,minor = PhyloNetworks.blobinfo(net);
    @test [n.number for n in r] == [-3,-11,-2]
    @test [[e.number for e in h] for h in major] == [[7, 3],[17],[]]
    @test [[e.number for e in h] for h in minor] == [[13,9],[21],[]]
    forest, blobs = blobdecomposition(net);
    @test length(blobs)==3
    @test writenewick(forest) == "(dummy -3,dummy -11);"
    s = IOBuffer()
    writesubtree!(s, blobs[1], nothing, false, true)
    @test String(take!(s)) == "(((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2));"
    writesubtree!(s, blobs[2], nothing, false, true)
    @test String(take!(s)) == "(((F)#H3,G),(H,#H3));"
    writesubtree!(s, blobs[3], nothing, false, true)
    @test String(take!(s)) == "(dummy -3,dummy -11);"

    # h=2, 2 non-trivial blobs above the LSA, LSA = the single tip
    net = readnewick("((((t1)#H22:::0.8,#H22))#H10:::0.7,#H10);")
    a = biconnectedcomponents(net,true);
    @test [[e.number for e in b] for b in a] == [[3,2], [6,5]]
    lsa, lsaind = PhyloNetworks.leaststableancestor(net)
    @test (lsa.number,lsaind) == (2,4)

    # h=3, one level-1 blob above the LSA, one level-2 blob below including a 2-cycle
    net = readnewick("((((((t2,#H25:::0.3))#H22:::0.8,#H22),(t1)#H25:::0.7))#H10:::0.6,#H10);")
    a = biconnectedcomponents(net,true);
    [[e.number for e in b] for b in a] == [ [5,8,2,3,4,6], [11,10]]

    aentry = PhyloNetworks.biconnectedcomponent_entrynodes(net, a)
    @test [n.number for n in aentry] == [-4,-2]
    aexit = PhyloNetworks.biconnectedcomponent_exitnodes(net, a)
    @test [[n.number for n in ae] for ae in aexit] == [[2,-7],[5]]

    # balanced tree + extra root edge
    net = readnewick("((((t1,t2),(t3,t4))));")
    _, lsaindex = PhyloNetworks.leaststableancestor(net)
    @test net.vec_node[lsaindex].number == -4
    # LSA = root & entry to non-trivial blob
    lsa, _ = PhyloNetworks.leaststableancestor(readnewick("(#H2:::0.2,((b)#H2,a));"))
    @test lsa.number == -2

    # level, process_biconnectedcomponents!
    net = readnewick("(((((#H25)#H22:::0.8,#H22),((t2:0.1,t1))#H25:::0.7)));")
    @test_throws "no biconnected components stored" PhyloNetworks.leaststableancestor(net, false, false)
    PhyloNetworks.process_biconnectedcomponents!(net)
    @test length(net.partition) == 5
    checkpart(i) = (net.partition[i].cycle, sort!([e.number for e in net.partition[i].edges]))
    @test checkpart(1) == ([1,2], [9])
    @test checkpart(2) == ([2,5], [1,2,3,4,8])
    @test checkpart(3) == ([5,6], [7])
    @test checkpart(4) == ([6], [5])
    @test checkpart(5) == ([6], [6])
    @testset for i in net.numedges
        @test net.edge[i] ∈ net.partition[net.edge[i].inte1].edges
    end
    @test PhyloNetworks.entrynode_preindex.(net.partition) == [1,2,5,6,6]
    @test collect.(PhyloNetworks.exitnodes_preindex.(net.partition)) == [[2],[5],[6],[],[]]
    @test length(PhyloNetworks.exitnodes_preindex(net.partition[5])) == 0
    @test getlevel(net,false,false) == 2
    @test PhyloNetworks.istrivial.(net.partition[[1,2,4]]) == [true, false,true]
    @test PhyloNetworks.ispendent.(net.partition[[1,2,4]]) == [false,false,true]
    lsa, _ = PhyloNetworks.leaststableancestor(net, false, false)
    @test lsa.number == -8
    PhyloNetworks.empty!.(net.partition)
end

@testset "tree component" begin

    treestr = "(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);"
    tree = readnewick(treestr)
    for e in tree.edge e.containroot=false; end # wrong, on purpose
    @test collect(values(treeedgecomponents(tree))) == repeat([1], inner=7)
    rcomp = checkroot!(tree)
    @test rcomp == 1
    @test all([e.containroot for e = tree.edge])

    netstr = "(#H1:::0.1,#H2:::0.2,(((b)#H1)#H2,a));"
    net = readnewick(netstr)
    for e in net.edge e.containroot=false; end # wrong, on purpose
    node2comp = treeedgecomponents(net) # e.g. [1,1,1,2,2,3] or [2,2,2,1,1,3]
    compsize = [count(isequal(i), values(node2comp)) for i in 1:3]
    @test sort(compsize) == [1,2,3]
    rcompID = checkroot!(net, node2comp)
    @test compsize[rcompID] == 3
    rcomp = keys(filter(p -> p.second == rcompID, node2comp))
    @test Set(n.number for n in rcomp) == Set([-2 -3 4])
    @test Set(e.number for e in net.edge if e.containroot) == Set([7 6 5 2 1])
    for e in net.edge e.containroot=true; end # wrong, on purpose
    net.rooti = 2 # node number 1, also H1: not in the root component on purpose
    checkroot!(net)
    @test [e.number for e in net.edge if e.containroot] == [1,2,5,6,7]
    @test net.rooti == 5 # i.e. node number -3. only 1 other option: 6, i.e. nn -2

    # test semidirected cycle case
    net.edge[1].ischild1 = !net.edge[1].ischild1
    net.edge[7].ischild1 = !net.edge[7].ischild1
    net.edge[7].hybrid = true
    net.edge[7].gamma  = net.edge[4].gamma
    net.edge[4].gamma  = 1.0
    net.edge[4].hybrid = false
    net.rooti = 5
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
    netl1 = readnewick(str_level1)
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
    netl1 = readnewick(str_level1)
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

@testset "ToB" begin
# binary, with 2-cycle -> degree-2 node
net = readnewick("((S1:0.1,(((S2,(S3)#H1),(#H1,S4:0.4)):0.2)#H2:0.2),(#H2,(#H3,(S5)#H3)));");
@test writenewick(treeofblobs(net)) == "((S5),S1:0.1,(S4:0.4,S2,S3):0.2);"
# rooted, chain of 2 degree-2 blobs to 1 taxon
net = readnewick("((((b)#H1,#H1))#H2,#H2);")
@test writenewick(treeofblobs(net)) == "((b));"
# extra root edge; non-binary: 1 block's entry = side of other
net = readnewick("((#H3,((a)#H3,(b1,#H1),(b2)#H1)));")
@test writenewick(treeofblobs(net)) == "(b2,b1,a);"
# unrooted, 2 blocks at root
net = readnewick("(((a1)#H2,(((#H2),a2))#H1),#H1,#H3,((b1)#H3,b2));")
@test writenewick(treeofblobs(net)) == "(b2,b1,a2,a1);"
end
end
