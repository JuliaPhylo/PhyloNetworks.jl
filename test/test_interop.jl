@testset "testing interoperability, matrix-based net" begin

# tree, some edge lengths missing
tree1 = readTopology("(A,(B:1.0,(C:1.0,D:1.0):1.0):1.0);");
@test_logs PhyloNetworks.resetNodeNumbers!(tree1);
tree1.edge[3].number = 50
@test_logs (:warn, r"^resetting edge numbers") PhyloNetworks.resetEdgeNumbers!(tree1);
@test tree1.edge[3].number == 3
@test PhyloNetworks.majoredgematrix(tree1) == [5 1; 5 6; 6 2; 6 7; 7 3; 7 4]
@test all(PhyloNetworks.majoredgelength(tree1) .===  [missing, 1.0,1.0,1.0,1.0,1.0])
@test PhyloNetworks.minorreticulationmatrix(tree1) == Array{Int64,2}(undef, 0,2)
@test size(PhyloNetworks.minorreticulationlength(tree1)) == (0,)
@test size(PhyloNetworks.minorreticulationgamma(tree1)) == (0,)

# network, h=1, some missing gamma values
net1 = (@test_logs (:warn, r"^third colon : without gamma value") readTopology("(((A:4.0,(B:1.0)#H1:1.1::):0.5,(C:0.6,#H1:1.0):1.0):3.0,D:5.0);"));
@test_logs PhyloNetworks.resetNodeNumbers!(net1);
@test PhyloNetworks.majoredgematrix(net1) == [5 6; 5 4; 6 8; 6 7; 7 3; 8 1; 8 9; 9 2]
@test PhyloNetworks.majoredgelength(net1) ==  [3.,5.,.5,1.,.6,4.,1.1,1.]
@test PhyloNetworks.minorreticulationmatrix(net1) == [7 9]
@test PhyloNetworks.minorreticulationlength(net1) == [1.]
@test all(PhyloNetworks.minorreticulationgamma(net1) .=== [missing])

# network, h=2 hybridizations
s = "(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,
      (((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);";
net2 = readTopology(s);
@test_logs PhyloNetworks.resetNodeNumbers!(net2);
@test PhyloNetworks.majoredgematrix(net2) == [13 15; 13 14; 13 12; 14 10; 14 11; 15 18;
  15 16; 16 17; 16 9; 17 24; 17 8; 18 1; 18 19; 19 20; 20 21; 20 23; 21 2; 21 22; 22 3;
  23 4; 24 25; 25 26; 25 7; 26 5; 26 6]
@test all(PhyloNetworks.majoredgelength(net2) .===
  [5.943,1.863,missing,missing,missing,.026,.723,.187,missing,0,missing,missing,
  2.214,2.490,.078,missing,missing,.023,.08,0,1.697,2.11,2.027,.002,.023])

#a.isnull[[3,4,5,9,11,12,16,17]] == [true for i in 1:8]
#@test convert(Array, a[[1,2,6,7,8,10,13,14,15]]) == [5.943,1.863, .026,.723,.187, 0, 2.214,2.490,.078]
#@test convert(Array, a[18:end]) == [.023,.08,0,1.697,2.11,2.027,.002,.023]
@test PhyloNetworks.minorreticulationmatrix(net2) == [19 24; 22 23]
@test PhyloNetworks.minorreticulationlength(net2) == [7.159,0.]
@test PhyloNetworks.minorreticulationgamma(net2) == [0.056,0.004]

end # of testset for interop
