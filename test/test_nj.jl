@testset "neighbour joining" begin

@testset "neighbour joining tree" begin
    D = CSV.read(joinpath(@__DIR__,"..","examples","caudata_dist.txt"),  DataFrame)
    tree = nj(D)
    # read results from ape implementation of nj
    apetree = readnewick(joinpath(@__DIR__, "..", "examples", "caudata_dist_nj.txt"))
    @test hardwiredclusterdistance(tree, apetree, false) == 0
    # also check branch lengths (more or less)
    @test sort!([e.length for e in tree.edge]) ≈ sort!([e.length for e in apetree.edge])

    # example where Ints are converted to Floats, and there's a < 0 edge length
    df = DataFrame(s1=[0, 5,9,9,5], s2=[5,0,10,10,6], s3=[9,10,0,8,4],
                   s4=[9,10,8,0,0], s5=[5,6, 4, 0,0])
    tree = (@test_logs (:info, r"negative length") nj(df))
    @test writenewick(tree) == "(((s1:2.0,s2:3.0):3.0,s3:4.0):2.0,s5:-2.0,s4:2.0);"
    # example with no names argument, force_nonnegative_edges
    D = [0 5 9 9 6.5; 5 0 10 10 7.5; 9 10 0 8 5.5; 9 10 8 0 1.5; 6.5 7.5 5.5 1.5 0]
    tree = (@test_logs (:info, r"reset to 0") PhyloNetworks.nj!(D, force_nonnegative_edges=true))
    @test writenewick(tree) == "(((1:2.0,2:3.0):3.0,3:4.0):2.0,5:0.0,4:2.0);"
end

@testset "neighbour joining tree of blobs" begin
    D = CSV.read(joinpath(@__DIR__,"..","examples","caudata_dist.txt"),  DataFrame)
    tree = PhyloNetworks.nj_treeofblobs(D)
    # read results from ape implementation of nj
    apetree = readnewick(joinpath(@__DIR__, "..", "examples", "caudata_dist_nj.txt"))
    @test hardwiredclusterdistance(tree, apetree, false) == 0

    # example with two blobs
    onet = readnewick("(((a:6,(b:2,(h:1.5)#H1:0::0.65):1):2,(#H1:0,c:3):3):1,((d:1,(e:2.5,(hh:3)#H2:0::0.75):1):1,(#H2:0,f:2):2):0.5);");
    oD = pairwisetaxondistancematrix(onet)
    njtob = PhyloNetworks.nj_treeofblobs!(copy(oD), tiplabels(onet))
    otob = treeofblobs(onet)
    @test hardwiredclusterdistance(otob, njtob, false) == 0

    # true blob vs no blob
    @test PN.cutblobpvalue(oD, BitVector([1,1,1,1,0,0,0,0]), 1) == 1
    @test PN.cutblobpvalue(oD, BitVector([1,1,1,0,1,0,0,0]), 1) ≈ 0 atol=1e-10
end

# @testset "Mixture Models" begin
#     mm = PN.mixtureChi2Diagonal(10)
#     @test cdf(mm, [5.4, 12.3, 23.4]) ≈ [0.0166, 0.3522, 0.9173] atol=1e-4
# end


end