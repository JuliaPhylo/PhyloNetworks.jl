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

# net1 subnetwork: net2 = (net1, (D)), and degree-2 node in net2
net2 = readnewick("(((A,(B)#H1),(#H1,C)),(D));")
labels = ["D","A","B","C"] # different order than from both net1 & net2
μ1 = (@test_logs PN.node_murepresentation(net1, labels, false)) # net1 already preordered
μ2 = (@test_logs PN.node_murepresentation(net2, labels))
@test μ1 != μ2
@test PN.mudistance(μ1, μ2) == 2 # root and degree-2 node above D

μ2o = PN.edge_murepresentation(net2, labels, false)
@test string(μ2o) == """PhyloNetworks.EdgeMuRepresentation
4 taxa in μ-vectors: ["D", "A", "B", "C"]
root μ-vector: μ0=2 μ=[1, 1, 2, 1]
maps edge number => μ-entry:
4 tree edges in the root component
  4 => μ0=1 μ=[0, 1, 1, 0]; μ0=1 μ=[1, 0, 1, 1]
  7 => μ0=1 μ=[0, 0, 1, 1]; μ0=1 μ=[1, 1, 1, 0]
  10 => μ0=0 μ=[1, 0, 0, 0]; μ0=2 μ=[0, 1, 2, 1]
  8 => μ0=2 μ=[0, 1, 2, 1]; μ0=0 μ=[1, 0, 0, 0]
2 edges in the directed part:
  5 => incident; μ0=1 μ=[0, 0, 1, 0]
  3 => incident; μ0=1 μ=[0, 0, 1, 0]"""

deleteleaf!(net2, "D")
μ2 = (@test_logs PN.node_murepresentation(net2, labels))
@test μ1 == μ2
@test PN.mudistance(μ1, μ2) == 0

μ1 = PN.edge_murepresentation(net1, labels, false)
μ2 = PN.edge_murepresentation(net2, labels, false)
@test μ2 == μ1
@test μ2 != μ2o
end

@testset "on more complex networks" begin

net1 = readnewick("((((b,a))#H1,#H1)h2,c,d)R;")
net2 = readnewick("(((((a,b))#H1,c)h2,#H1),d)R;")
@test mudistance_rooted(net1, net2) == 3

labs = ["a","b","c","d"]
μ1  = PN.edge_murepresentation(net1 , labs, false)
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
@test mudistance_semidirected(net1, net1r; preorder=false) == 0

# non-orchard network, with all tag types
end
end
