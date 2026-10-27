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

@testset "nj tob from distances" begin
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

    # with 2-cycle (invisible from the distance)
    onet = readnewick("((S1:0.1,(((S2:1,(S3:2)#H1:0.6::0.65):0.2,(#H1:0.6,S4:0.4):0.2):0.2)#H2:0.4):0.3,(#H2:0.2::0.75,(#H3:0.1::0.35,(S5:0.2)#H3:0.1):0.4):0.5);");
    oD = pairwisetaxondistancematrix(onet)
    njtob = PhyloNetworks.nj_treeofblobs!(copy(oD), tiplabels(onet))
    otob = treeofblobs(onet)
    otob = removedegree2nodes!(otob) # remove degree 2 node leading to S5
    @test hardwiredclusterdistance(otob, njtob, false) == 0

    # Level 2, galled network
    onet = readnewick("(((x1:0.1,#H1:0.7::0.7):0.6,((((x2:0.3)#H1:0.3,x3:0.9):0.4,#H2:0.3::0.45):0.3,((x4:0.4)#H2:0.3,x5:0.3):0.3):0.3):2.0,((x6:1,(x7:2.5,(x8:3)#H3:0::0.75):1):1,(#H3:0,x9:2):2):0.5);")
    oD = pairwisetaxondistancematrix(onet)
    njtob = PhyloNetworks.nj_treeofblobs!(copy(oD), tiplabels(onet))
    otob = treeofblobs(onet)
    @test hardwiredclusterdistance(otob, njtob, false) == 0
end

# @testset "Mixture Models" begin
#     mm = PN.mixtureChi2Diagonal(10)
#     @test cdf(mm, [5.4, 12.3, 23.4]) ≈ [0.0166, 0.3522, 0.9173] atol=1e-4
# end

@testset "nj tob from f2" begin
    D = [0.0        0.0705476  0.134613   0.29065    0.283362;
        0.0705476  0.0        0.0970743  0.283194   0.275375;
        0.134613   0.0970743  0.0        0.295239   0.282258;
        0.29065    0.283194   0.295239   0.0        0.0557102;
        0.283362   0.275375   0.282258   0.0557102  0.0]

    Qf2 = [  
        1.28002e-5   7.14817e-6  -1.82795e-6   -2.79157e-6    9.44815e-7  -1.14182e-6  -1.35661e-6  -5.3093e-7   -8.51401e-7   -3.43044e-7;
        7.14817e-6   4.11061e-5   2.16668e-6    1.87709e-6    1.62578e-5  -2.71216e-6  -3.13008e-6   2.35249e-6  -2.53742e-7    1.16926e-6;
        -1.82795e-6   2.16668e-6   0.000107486   8.75079e-5   -8.77585e-7   7.28136e-5   5.93078e-5   5.14487e-5   4.27387e-5   -3.95933e-6;
        -2.79157e-6   1.87709e-6   8.75079e-5    0.000103984  -7.8349e-7    5.97909e-5   7.3265e-5    4.29741e-5   5.38306e-5   -2.6987e-6;
        9.44815e-7   1.62578e-5  -8.77585e-7   -7.8349e-7     2.29756e-5  -1.84632e-6  -1.62229e-6   2.45254e-7  -1.58019e-6    2.58864e-7;
        -1.14182e-6  -2.71216e-6   7.28136e-5    5.97909e-5   -1.84632e-6   9.82021e-5   8.16645e-5   6.24272e-5   5.17748e-5   -3.82119e-6;
        -1.35661e-6  -3.13008e-6   5.93078e-5    7.3265e-5    -1.62229e-6   8.16645e-5   9.80959e-5   5.11546e-5   6.4475e-5    -3.69453e-6;
        -5.3093e-7    2.35249e-6   5.14487e-5    4.29741e-5    2.45254e-7   6.24272e-5   5.11546e-5   0.00010627   8.53185e-5   -1.62183e-6;
        -8.51401e-7  -2.53742e-7   4.27387e-5    5.38306e-5   -1.58019e-6   5.17748e-5   6.4475e-5    8.53185e-5   0.000103266  -3.96531e-6;
        -3.43044e-7   1.16926e-6  -3.95933e-6   -2.6987e-6     2.58864e-7  -3.82119e-6  -3.69453e-6  -1.62183e-6  -3.96531e-6    1.29104e-5
        ]

    njtob = PhyloNetworks.nj_treeofblobs!(copy(D), nloci = 10000)
    njtob2 = PhyloNetworks.nj_treeofblobs!(copy(D), copy(Qf2), nloci = 10)
    @test hardwiredclusterdistance(njtob, njtob2, false) == 0
end

end