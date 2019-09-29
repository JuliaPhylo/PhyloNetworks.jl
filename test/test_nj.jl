@testset "Neighbour joining implementation" begin

    D = CSV.read(joinpath(@__DIR__, "..", "examples", "caudata_dist.txt"))
    tree = nj(D)
    # read results from ape implementation of nj
    apetree = readTopology(joinpath(@__DIR__, "..", "examples", "caudata_dist_nj.txt"))
    @test hardwiredClusterDistance(tree, apetree, false) == 0
    # also check branch lengths (more or less)
    @test sort!([e.length for e in tree.edge]) ≈ sort!([e.length for e in apetree.edge])

    # example where Ints are converted to Floats, and there's a < 0 edge length
    df = DataFrame(s1=[0, 5,9,9,5], s2=[5,0,10,10,6], s3=[9,10,0,8,4],
                   s4=[9,10,8,0,0], s5=[5,6, 4, 0,0])
    tree = (@test_logs (:info, r"have negative lengths") nj(df))
    @test writeTopology(tree) == "(((s1:2.0,s2:3.0):3.0,s3:4.0):2.0,s4:2.0,s5:-2.0);"
    # example with no names argument, force_nonnegative_edges
    D = [0 5 9 9 6.5; 5 0 10 10 7.5; 9 10 0 8 5.5; 9 10 8 0 1.5; 6.5 7.5 5.5 1.5 0]
    tree = (@test_logs (:info, r"reset to 0") PhyloNetworks.nj!(D, force_nonnegative_edges=true))
    @test writeTopology(tree) == "(((1:2.0,2:3.0):3.0,3:4.0):2.0,4:2.0,5:0.0);"
end
