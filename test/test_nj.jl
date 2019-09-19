@testset "Neighbour joining implementation" begin

    D = CSV.read(joinpath(@__DIR__, "..", "examples", "caudata_dist.txt");
              types=[Float64 for i in 1:197])
    D = convert(Matrix, D)
    tree = nj!(copy(D))

    # read results from ape implementation of nj
    apetree = readTopology(joinpath(@__DIR__, "..", "examples", "caudata_dist_nj.txt"))
    for i in 1:apetree.numTaxa
        apetree.leaf[i].name = ""
    end

    @test hardwiredClusterDistance(tree, apetree, false) == 0

end
