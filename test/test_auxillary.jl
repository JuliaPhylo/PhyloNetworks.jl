#= # for local testing, need this:
using Test
using PhyloNetworks
using PhyloPlots
using CSV
=#

@testset "auxiliary" begin

@testset "setlengths and setgammas" begin
originalstdout = stdout
redirect_stdout(open("/dev/null", "w")) # not portable to Windows
@test_nowarn PhyloNetworks.citation()
redirect_stdout(originalstdout)

str_level1_s = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));" # indviduals S1A S1B S1C go on leaf 1
net = readTopology(str_level1_s)

PhyloNetworks.setlengths!([net.edge[1]], [1.1])
@test net.edge[1].length == 1.1
PhyloNetworks.setlengths!([net.edge[3], net.edge[4]], [3.3, 4.4])
@test net.edge[3].length == 3.3
@test net.edge[4].length == 4.4

PhyloNetworks.setmultiplegammas!([net.edge[18]], [0.25])
@test net.edge[18].gamma == 0.25
@test net.edge[16].gamma == 0.75

@test PhyloNetworks.getlengths([net.edge[1]]) == [net.edge[1].length]
@test PhyloNetworks.getlengths([net.edge[1], net.edge[5]]) == [net.edge[1].length, net.edge[5].length]
end

@testset "hashybridladder" begin
tree = readTopology("(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);");
@test !PhyloNetworks.hashybridladder(tree)
PhyloNetworks.addhybridedge!(tree, tree.edge[5], tree.edge[1], true)
PhyloNetworks.addhybridedge!(tree, tree.edge[2], tree.edge[1], true)
@test PhyloNetworks.hashybridladder(tree)
end # of testing hashybridladder

@testset "shrink 2/3 cycles" begin
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
@test !PhyloNetworks.shrink2cycles!(net)
@test !PhyloNetworks.shrink3cycles!(net)
PhyloNetworks.addhybridedge!(net, net.edge[7], net.edge[4], true, 0.1, 0.2)
@test !PhyloNetworks.shrink2cycles!(net)
@test PhyloNetworks.shrink3cycles!(net) # tree case
@test writeTopology(net) == "(((A:2.0,(B:1.0)#H1:0.1::0.9):1.37,(C:0.6,#H1:1.0::0.1):0.9):0.6,D:2.0);"
PhyloNetworks.addhybridedge!(net, net.edge[3], net.edge[6], true, 0.3, 0.4)
@test PhyloNetworks.shrink3cycles!(net) # hybrid case
@test writeTopology(net, round=true, digits=5) == "(((A:2.0,(B:1.0)#H1:0.13191::0.94):1.37,(C:0.6,#H1:1.0::0.06):0.9):0.6,D:2.0);"
net0 = readTopology("((((((a:1)#H1:1::.9)#H2:1::.8)#H3:1::.7,#H3:0.5):1,#H2:1):1,(#H1:1,b:1):1,c:1);")
net = deepcopy(net0) # new 2/3 cycles appear when some are shrunk
@test PhyloNetworks.shrink2cycles!(net)
@test writeTopology(net, round=true) == "((#H1:1.0::0.1,b:1.0):1.0,c:1.0,(a:1.0)#H1:4.48::0.9);"
@test PhyloNetworks.shrink3cycles!(net0)
writeTopology(net0, round=true) == "(c:1.1,a:5.132,b:1.9);"
# non-tree-child network: w shape
net0 = readTopology("((a:1,#H1:.1::.1):1,(((b:.5)#H3:1)#H1:1,(#H3:0.8::.4)#H2:1):1,(#H2:.2::.0,c:1):1);")
PhyloNetworks.deletehybridedge!(net0, net0.edge[2], false,true,false,false)
net = deepcopy(net0) # delete edge 10 (with γ=0) then shrink 2-cycle
PhyloNetworks.deletehybridedge!(net, net.edge[7], false,true,false,false)
@test PhyloNetworks.shrink3cycles!(net)
@test writeTopology(net) == "(a:2.0,c:2.0,b:3.42);"
# shrink 3-cycle, which deletes edge 10 (with γ=0) : same result in the end
@test PhyloNetworks.shrink3cycles!(net0)
@test writeTopology(net) == "(a:2.0,c:2.0,b:3.42);"
end

end # of set of auxiliary test sets
