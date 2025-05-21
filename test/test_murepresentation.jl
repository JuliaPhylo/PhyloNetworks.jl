@testset "μ distances" begin
@testset "on μ-representations" begin

net1 = readnewick("((C,(B)#H1),(#H1,A));")
μ1 = (@test_logs (:warn, r"^leaf C") PN.node_murepresentation(net1, ["D","A","B"]))
@test string(μ1) == """PhyloNetworks.NodeMuRepresentation
3 taxa in μ-vectors: ["D", "A", "B"]
4 nodes, map node number => μ0 to hybrids and μ-vector to taxa:
  -5 => μ0=1 μ=[0, 1, 1]
  -3 => μ0=1 μ=[0, 0, 1]
  -2 => μ0=2 μ=[0, 1, 2]
  3 => μ0=1 μ=[0, 0, 1]"""

# net1 as subnetwork: net2 = (net1, (D)), and degree-2 node in net2
net2 = readnewick("(((A,(B)#H1),(#H1,C)),(D));")
labels = ["D","A","B","C"] # different order than from both net1 & net2
μ1 = (@test_logs PN.node_murepresentation(net1, labels; preorder=false))
μ2 = (@test_logs PN.node_murepresentation(net2, labels))
@test μ1 != μ2
@test PN.mudistance(μ1, μ2) == 2 # root and degree-2 node above D
@test ((@test_logs (:warn,r"leaf") (:warn,r"leaf") (:warn,r"leaf") mudistance_semidirected(
    net1, net2; preorder=false, labels=["A","B"]))) == 2

μ2o = PN.edge_murepresentation(net2, labels; preorder=false)
@test string(μ2o) == """PhyloNetworks.EdgeMuRepresentation
4 taxa in μ-vectors: ["D", "A", "B", "C"]
root μ-vector: μ0=2 μ=[1, 1, 2, 1]
maps edge number => μ-entry:
3 tree edges in the root component
  4 => μ0=1 μ=[0, 1, 1, 0]; μ0=1 μ=[1, 0, 1, 1]
  7 => μ0=1 μ=[0, 0, 1, 1]; μ0=1 μ=[1, 1, 1, 0]
  8 => μ0=2 μ=[0, 1, 2, 1]; μ0=0 μ=[1, 0, 0, 0]
2 edges in the directed part:
  5 => incident; μ0=1 μ=[0, 0, 1, 0]
  3 => incident; μ0=1 μ=[0, 0, 1, 0]"""

deleteleaf!(net2, "D")
μ2 = (@test_logs PN.node_murepresentation(net2, labels))
@test μ1 == μ2
@test PN.mudistance(μ1, μ2) == 0

μ1 = PN.edge_murepresentation(net1, labels; preorder=false)
μ2 = PN.edge_murepresentation(net2, labels; preorder=false)
@test μ2 == μ1
@test μ2 != μ2o
end

@testset "on more complex networks" begin

net1 = readnewick("((((b,a))#H1,#H1)h2,c,d)R;")
net2 = readnewick("(((((a,b))#H1,c)h2,#H1),d)R;")
@test mudistance_rooted([net1, net2]) == [0 3; 3 0]
@test mudistance_semidirected([net1, net2]; preorder=false, userootμ=true) == [0 2; 2 0]

labs = ["a","b","c","d"]
μ1  = PN.edge_murepresentation(net1, labs; preorder=false)
@test string(μ1) == """
PhyloNetworks.EdgeMuRepresentation
4 taxa in μ-vectors: ["a", "b", "c", "d"]
root μ-vector: μ0=2 μ=[2, 2, 1, 1]
maps edge number => μ-entry:
1 tree edges in the root component
  6 => μ0=2 μ=[2, 2, 0, 0]; μ0=0 μ=[0, 0, 1, 1]
3 edges in the directed part:
  5 => incident; μ0=1 μ=[1, 1, 0, 0]
  4 => incident; μ0=1 μ=[1, 1, 0, 0]
  3 => tree; μ0=0 μ=[1, 1, 0, 0]"""
net1r = deepcopy(net1); rootatnode!(net1r, 4);
μ1r = PN.edge_murepresentation(net1r, labs)
@test string(μ1r.mumap_rootcomp[6]) == "(μ0=0 μ=[0, 0, 1, 1], μ0=2 μ=[2, 2, 0, 0])"
@test mudistance_semidirected([net1, net1r]; userootμ=true) == [0 0; 0 0]

# networks with differen leaf sets
deleteleaf!(net1r,"c")
@test mudistance_rooted(net1, net1r) == 3
@test mudistance_semidirected(net1, net1r) == 1
@test mudistance_semidirected(net1, net1r; preorder=false, userootμ=true) == 2

# non-orchard network, with all tag types
nwkstr1 = "((b1:1,((a1:1)#H12:3,#H1:0)r1:2):1,((x1:.5,((#H21:1,#H3:.5)w1:1)#H1:1):1,((((c:1)#H3:.5,#H12:1)w2:1)#H2:1,x2:.5):1):1,((#H2:0,(a2:1)#H21:3)r2:2,b2:1):1);"
net1 = readnewick(nwkstr1)
# plot(net1, useedgelength=true, showedgenumber=true, shownodelabel=true);
μ1 = PN.edge_murepresentation(net1) # default labels
@test tiplabels(μ1) == ["b1","a1","x1","c","x2","a2","b2"]
@test μ1.mu_root == PN.MuVector(14, [1,3,1,4,1,3,1])
rv = μ1.muvec_rootcomp
@test length(rv) == 7
@test [(m[1].mu_hybs, m[2].mu_hybs) for m in rv] == [(3,11),(3,11),(4,10),(4,10),(4,10),(4,10),(6,8)]
@test rv[1][1].mu_tips == [0,0,1,1,0,1,0]
@test rv[1][2].mu_tips == [1,3,0,3,1,2,1]
@test rv[7][1].mu_tips == [0,1,1,2,1,1,0]
@test rv[7][2].mu_tips == [1,2,0,2,0,2,1]
dv = μ1.muvec_directed
@test length(dv) == 12
@test [m[1] for m in dv] == vcat(repeat([PN.mutag_i],6),repeat([PN.mutag_t],2),
    repeat([PN.mutag_h],4))

# swap b1 and b2: different topology, yet same hardwired clusters etc
nwkstr2 = replace(replace(replace(nwkstr1,"b1"=>"tmp"), "b2"=>"b1"), "tmp"=>"b2")
net2 = readnewick(nwkstr2)
# @test hardwiredclusterdistance(net1, net2, true) == 0
@test mudistance_rooted(net1, net2) == 0
@test mudistance_semidirected(net1, net2) == 0
# test networks with different leaf sets
net2 = readnewick("(((A,(B)#H1),(#H1,C)),(D));")
μ1 = (@test_logs PN.node_murepresentation(net1, tiplabels(net1)))
μ2 = (@test_logs PN.node_murepresentation(net2, tiplabels(net2)))
@test μ1 != μ2
@test µ2 != μ1
end
end
