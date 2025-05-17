@testset "μ distances" begin
@testset "on the μ-representations" begin

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
  4 => μ0=1 μ=[0, 1, 1, 0]; μ0=1 μ=[0, 1, 1, 0]
  7 => μ0=1 μ=[0, 0, 1, 1]; μ0=1 μ=[0, 0, 1, 1]
  10 => μ0=0 μ=[1, 0, 0, 0]; μ0=0 μ=[1, 0, 0, 0]
  8 => μ0=2 μ=[0, 1, 2, 1]; μ0=2 μ=[0, 1, 2, 1]
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


@testset "on networks" begin
    net1 = readnewick("(((a,b)#H1,#H1)h2,c,d)R;")
    labels = ["a", "b", "c", "d"]

    num2name = Dict{Int, String}()
    for i in net1.node
        num2name[i.number] = i.name
    end
    μ1 = PN.node_murepresentation(net1, labels)
    name_map = Dict{String, Tuple{Vararg{Int}}}()
    
    μ_map = μ1.mu_map
    for (i, value) in μ_map
        name_map[num2name[i]] = value
    end
    println("net1",name_map)
    @test name_map == Dict{String, Tuple{Vararg{Int64}}}( "c" => (0, 0, 1, 0),"H1" => (1, 1, 0, 0), "h2" => (2, 2, 0, 0), "b" => (0, 1, 0, 0), "R" => (2, 2, 1, 1), "a" => (1, 0, 0, 0), "d" => (0, 0, 0, 1))

    net2 = readnewick("(((a,b)#H1,c)h2,#H1,d)R;")
    labels = ["a", "b", "c", "d"]

    num2name = Dict{Int, String}()
    for i in net2.node
        num2name[i.number] = i.name
    end
    μ2 = PN.node_murepresentation(net2, labels)
    name_map = Dict{String, Tuple{Vararg{Int}}}()
    μ_map = μ2.mu_map
    for (i, value) in μ_map
        name_map[num2name[i]] = value
    end
    println("net2",name_map)
    @test name_map == Dict{String, Tuple{Vararg{Int64}}}( "c" => (0, 0, 1, 0),"H1" => (1, 1, 0, 0), "h2" => (1, 1, 1, 0), "b" => (0, 1, 0, 0), "R" => (2, 2, 1, 1), "a" => (1, 0, 0, 0), "d" => (0, 0, 0, 1))
    
    @test mudistance_rooted(net1, net2) == 2
end
end
