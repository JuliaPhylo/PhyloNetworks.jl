# @testset "node_mu_rep invariant checks for various networks" begin

#     # Helper function to check μ-vector invariants
#     function check_mu_invariants(newick_str::String)
#         net = readnewick(newick_str)
#         labels = tipLabels(net)
#         mu = node_mu(net, labels).mu_map

#         # Map each leaf node to its index in labels
#         label_map = Dict{Node, Int}(leaf => i for (i, leaf) in enumerate(net.leaf))

#         for node in net.node
#             if node.leaf
#                 @test mu[node.number][label_map[node]] == 1
#                 @test sum(mu[node.number]) == 1
#             else
#                 leaf_count = sum(mu[node.number])
#                 if node == net.node[1]  # root node after preorder!
#                     @test leaf_count == length(net.leaf)
#                 else
#                     @test leaf_count ≥ 1
#                 end
#             end
#         end
#     end

#     # All networks tested under one parent set
#     newicks = [
#         "((a,b)mu1,c,d)r;",
#         "((C:0.9,(B:0.2)#H1:0.7::0.6):0.6,(#H1:0.6,A:1):0.5);"
#     ]

#     for newick in newicks
#         check_mu_invariants(newick)
#     end
# end


@testset "node mu equality test" begin
    net1 = readnewick("((C,(B)#H1),(#H1,A));")
    net2 = readnewick("((B,(C)#H1),(#H1,A));")

    labels = ["A", "B", "C"]  # must match all tips used in both networks

    μ1 = PN.node_mu(net1, labels)
    μ2 = PN.node_mu(net2, labels)
    println(μ1== μ2)
    @test μ1 == μ2  

    newick1 = "((C,(B)#H1),(#H1,A));"
    newick2 = "((B,(C)#H1),(#H1,A));"

    labels = ["A", "B", "C"]  # must match all tips used in both networks

    μ1 = node_mu(net1, labels)
    μ2 = node_mu(net2, labels)
    @test μ1 == μ2 

end
