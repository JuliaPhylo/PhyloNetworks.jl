@testset "Neighbour joining implementation" begin

    D = CSV.read(joinpath(@__DIR__, "..", "examples", "caudata_dist.txt");
              types=[Float64 for i in 1:197])
    tree = nj(D)

    # read results from ape implementation of nj
    apetree = readTopology(joinpath(@__DIR__, "..", "examples", "caudata_dist_nj.txt"))

    @test hardwiredClusterDistance(tree, apetree, false) == 0

end
