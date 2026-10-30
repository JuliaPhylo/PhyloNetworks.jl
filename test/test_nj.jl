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
    D = [
        0.0        0.0733858  0.129473   0.290621   0.281591;
        0.0733858  0.0        0.0931066  0.281698   0.273439;
        0.129473   0.0931066  0.0        0.286105   0.2779;
        0.290621   0.281698   0.286105   0.0        0.0556367;
        0.281591   0.273439   0.2779     0.0556367  0.0
    ]

    Qf2 = [  
        1.50478e-6   8.09507e-7  -7.36444e-8  -1.30472e-7   1.07442e-7  -1.82893e-7  -2.03363e-7  -9.49698e-8  -1.54482e-7  -1.83797e-8
        8.09507e-7   3.67832e-6   8.96178e-8   2.11283e-8   1.34691e-6  -3.55667e-7  -4.55563e-7   2.49094e-8  -8.8153e-9    6.3331e-9
        -7.36444e-8   8.96178e-8   1.03447e-5   8.15767e-6  -8.78664e-8   7.04807e-6   5.43957e-6   5.16223e-6   3.88522e-6  -2.74406e-7
        -1.30472e-7   2.11283e-8   8.15767e-6   9.86822e-6  -1.64891e-7   5.45756e-6   6.81228e-6   3.87791e-6   5.01919e-6  -3.02898e-7
        1.07442e-7   1.34691e-6  -8.78664e-8  -1.64891e-7   2.11966e-6  -6.30336e-8  -1.33292e-7  -8.63752e-8  -8.02527e-8   9.62365e-9
        -1.82893e-7  -3.55667e-7   7.04807e-6   5.45756e-6  -6.30336e-8   9.86575e-6   7.81126e-6   6.193e-6     4.74428e-6  -2.78907e-7
        -2.03363e-7  -4.55563e-7   5.43957e-6   6.81228e-6  -1.33292e-7   7.81126e-6   9.45751e-6   4.69395e-6   5.97689e-6  -3.00158e-7
        -9.49698e-8   2.49094e-8   5.16223e-6   3.87791e-6  -8.63752e-8   6.193e-6     4.69395e-6   1.0102e-5    7.87003e-6  -2.58748e-7
        -1.54482e-7  -8.8153e-9    3.88522e-6   5.01919e-6  -8.02527e-8   4.74428e-6   5.97689e-6   7.87003e-6   9.63801e-6  -2.21982e-7
        -1.83797e-8   6.3331e-9   -2.74406e-7  -3.02898e-7   9.62365e-9  -2.78907e-7  -3.00158e-7  -2.58748e-7  -2.21982e-7   1.19921e-6
        ]

    njtob = PhyloNetworks.nj_treeofblobs!(copy(D), nloci = 100000)
    njtob2 = PhyloNetworks.nj_treeofblobs!(copy(D), copy(Qf2))
    @test hardwiredclusterdistance(njtob, njtob2, false) == 0

    njtob3 = PhyloNetworks.nj_treeofblobs!(copy(D), copy(Qf2); conditionalCutTest=false)
    @test PN.isEqual(njtob2,njtob3)

end

end