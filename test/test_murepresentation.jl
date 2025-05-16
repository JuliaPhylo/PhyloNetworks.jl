
@testset "mu representation equality test" begin
net1 = readnewick("((C,(B)#H1),(#H1,A));")
net2 = readnewick("((A,(B)#H1),(#H1,C));")

μ1 = (@test_logs (:warn, r"^leaf C") PN.node_murepresentation(net1, ["D","A","B"]))
@test string(μ1) == """PhyloNetworks.NodeMuRepresentation
3 taxa in μ-vectors: ["D", "A", "B"]
4 nodes, map node number => μ0 to hybrids and μ-vector to taxa:
  -5 => μ0=1 μ=[0, 1, 1]
  -3 => μ0=1 μ=[0, 0, 1]
  -2 => μ0=2 μ=[0, 1, 2]
  3 => μ0=1 μ=[0, 0, 1]"""

labels = ["D","A","B","C"]
μ1 = PN.node_murepresentation(net1, labels)
μ2 = PN.node_murepresentation(net2, labels)

node_distance = mudistance_rooted(net1, net2)

@test μ1 == μ2
@test node_distance == 0

μ1 = PN.edge_murepresentation(net1, labels)
μ2 = PN.edge_murepresentation(net2, labels)

edge_distance = mudistance_rooted(net1, net2)

@test μ1 == μ2
@test edge_distance == 0
end


@testset "non zero distance test" begin
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
