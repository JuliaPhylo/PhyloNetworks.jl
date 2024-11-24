#= # for local testing, need this:
using Test
using PhyloNetworks
using PhyloPlots
using CSV
=#

@testset "auxiliary" begin

@testset "show, setlengths, setgammas" begin
originalstdout = stdout
redirect_stdout(devnull) # requires julia v1.6
@test_nowarn PhyloNetworks.citation()
redirect_stdout(originalstdout)

str_level1_s = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));" # indviduals S1A S1B S1C go on leaf 1
net = readnewick(str_level1_s)
net0 = readnewick(str_level1_s)

s = IOBuffer()
show(s, net0)
@test String(take!(s)) == """
HybridNetwork, Rooted Network
20 edges
19 nodes: 8 tips, 2 hybrid nodes, 9 internal tree nodes.
tip labels: S8, S9, S1, S4, ...
(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));
"""
show(s,net0.node[8])
@test String(take!(s)) == """
PhyloNetworks.Node:
 number:6
 name:H1
 hybrid node
 attached to 3 edges, numbered: 7 8 10
"""
show(s,net0.edge[10])
@test String(take!(s)) == """
PhyloNetworks.EdgeT{PhyloNetworks.Node}:
 number:10
 length:-1.0
 minor hybrid edge with gamma=-1.0
 attached to 2 node(s) (parent first): -10 6
"""

@test PhyloNetworks.isEqual(net, net0)
@test getchild(net0.hybrid[1]).name == "S5"
@test_throws "0 children instead of 1" getchild(net0.leaf[1])
@test_throws "did not find a partner" getpartneredge(net0.edge[1])
@test PhyloNetworks.getIndex(net0.edge[10], net) == 10 # note: edge in net0, searched in net
@test PhyloNetworks.getIndexNode(net0.edge[10], net0.node[8]) == 1
@test PhyloNetworks.getIndexHybrid(net0.node[8], net) == 1
@test PhyloNetworks.getconnectingedge(net0.node[8], net0.node[13]) === net0.edge[10]
PhyloNetworks.deleteIntNode!(net0, getroot(net0))
@test net0.numedges == 19
@test net0.numnodes == 18
e1,e2,e3 = PhyloNetworks.hybridEdges(net0.node[13])
# broken: bc utilities to update `intn1` (cycle number for level-1 nets) are in SNaQ.
# fix later with blob decomposition for general networks
# @test (e1.number, e2.number, e3.number) == (10,14,13)
# edge 14 is part of the cycle from hybrid 6 (H1), yet net.edge[14].inte1 is -1
net0.hybrid[1].name = ""
PhyloNetworks.assignhybridnames!(net0)
@test net0.hybrid[1].name == "H6"
ee = PhyloNetworks.adjacentedges(net0.edge[10])
@test [e.number for e in ee] == [7,8,10,13,14]

setlengths!([net.edge[1]], [1.1])
@test net.edge[1].length == 1.1
setlengths!([net.edge[3], net.edge[1]], [3.3, missing])
@test net.edge[3].length == 3.3
@test net.edge[1].length == -1
setlength!(net.edge[3], -1)
@test net.edge[3].length == -1
@test_throws "non negative" setlength!(net.edge[3], -0.5)

PhyloNetworks.setmultiplegammas!([net.edge[18]], [0.25])
@test net.edge[18].gamma == 0.25
@test net.edge[16].gamma == 0.75

@test PhyloNetworks.getlengths([net.edge[1]]) == [net.edge[1].length]
@test PhyloNetworks.getlengths([net.edge[1], net.edge[5]]) == [net.edge[1].length, net.edge[5].length]
end

@testset "hashybridladder, istreechild" begin
tree = readnewick("(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);");
@test !PhyloNetworks.hashybridladder(tree)
PhyloNetworks.addhybridedge!(tree, tree.edge[5], tree.edge[1], true)
PhyloNetworks.addhybridedge!(tree, tree.edge[2], tree.edge[1], true)
@test PhyloNetworks.hashybridladder(tree)

# 2 biconnected components at the root,
# degree-2 node above hybrid in directed part: not tree-child.
# after suppression: weakly tree-child but not rooted tree-child,
net = readnewick("(((a1)#H2,(((#H2),a2))#H1),#H1,#H3,((b1)#H3,b2));")
@test getlevel(net) == 2
@test istreechild(net) == (false, false, false)
removedegree2nodes!(net)
@test istreechild(net) == (false, true, false)
rootatnode!(net, -3)
@test istreechild(net) == (true, true, false)
PhyloNetworks.deletehybridedge!(net, net.edge[8])
@test istreechild(net) == (true, true, true)
end # of testing hashybridladder

@testset "shrink edges and cycles" begin
# shrink 1 edge, and hassinglechild
net = readnewick("((A:2.0,(((B1,B2):1.0)0.01)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5);")
@test !PhyloNetworks.shrinkedge!(net, net.edge[10])
@test_throws Exception PhyloNetworks.shrinkedge!(net, net.edge[6]) # hybrid edge
@test_throws Exception PhyloNetworks.shrinkedge!(net, net.edge[3]) # external edge
@test hassinglechild(net.hybrid[1])
@test hassinglechild(net.node[5]) && !net.node[5].hybrid # degree-2 tree node
@test !PhyloNetworks.shrinkedge!(net, net.edge[5])
@test  PhyloNetworks.shrinkedge!(net, net.edge[4])
@test !hassinglechild(net.hybrid[1])
# shrink cycles
net = readnewick("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
@test !shrink2cycles!(net)
@test !shrink3cycles!(net)
PhyloNetworks.addhybridedge!(net, net.edge[7], net.edge[4], true, 0.1, 0.2)
@test !shrink2cycles!(net)
@test shrink3cycles!(net) # tree case
@test writenewick(net) == "(((A:2.0,(B:1.0)#H1:0.1::0.9):1.37,(C:0.6,#H1:1.0::0.1):0.9):0.6,D:2.0);"
PhyloNetworks.addhybridedge!(net, net.edge[3], net.edge[6], true, 0.3, 0.4)
@test shrink3cycles!(net) # hybrid case
@test writenewick(net, round=true, digits=5) == "(((A:2.0,(B:1.0)#H1:0.13191::0.94):1.37,(C:0.6,#H1:1.0::0.06):0.9):0.6,D:2.0);"
net0 = readnewick("((((((a:1)#H1:1::.9)#H2:1::.8)#H3:1::.7,#H3:0.5):1,#H2:1):1,(#H1:1,b:1):1,c:1);")
net = deepcopy(net0) # new 2/3 cycles appear when some are shrunk
@test shrink2cycles!(net)
@test writenewick(net, round=true) == "((#H1:1.0::0.1,b:1.0):1.0,c:1.0,(a:1.0)#H1:4.48::0.9);"
@test shrink3cycles!(net0)
writenewick(net0, round=true) == "(c:1.1,a:5.132,b:1.9);"
# non-tree-child network: w shape
net0 = readnewick("((a:1,#H1:.1::.1):1,(((b:.5)#H3:1)#H1:1,(#H3:0.8::.4)#H2:1):1,(#H2:.2::.0,c:1):1);")
PhyloNetworks.deletehybridedge!(net0, net0.edge[2], false,true,false,false)
net = deepcopy(net0) # delete edge 10 (with γ=0) then shrink 2-cycle
PhyloNetworks.deletehybridedge!(net, net.edge[7], false,true,false,false)
@test shrink3cycles!(net)
@test writenewick(net) == "(a:2.0,c:2.0,b:3.42);"
# shrink 3-cycle, which deletes edge 10 (with γ=0) : same result in the end
@test shrink3cycles!(net0)
@test writenewick(net) == "(a:2.0,c:2.0,b:3.42);"
end

@testset "readwrite.jl" begin
exdir = joinpath(@__DIR__,"..","examples")
# exdir = joinpath(dirname(pathof(PhyloNetworks)), "..","examples")
ns = readmultinewick(joinpath(exdir,"net1.networks"), false) # not fast
@test length(ns) == 5
ns = (@test_logs (:warn, r"won't erase with") (:warn, r"^skipped phylogeny on line 15") readnexus_treeblock(
    joinpath(exdir,"test_reticulatetreeblock.nex")))
@test length(ns) == 3
@test tiplabels(ns[2]) == ["tax4","tax3","tax2","tax1"]
# net with a hybrid leaf
net = readnewick("((#H1:::0.1,b),c,#H1:::0.9);")
PhyloNetworks.addChild!(net, net.node[1])
PhyloNetworks.addChild!(net, net.node[1])
@test writenewick(net) == "((#H1:::0.1,b),c,(H1_4:0.0,H1_5:0.0)#H1:::0.9);"
PhyloNetworks.expandChild!(net, net.node[1])
@test writenewick(net) == "((#H1:::0.1,b),c,((H1_5:0.0,H1_4:0.0)H1_6:0.0)#H1:::0.9);"
end

end # of set of auxiliary test sets
