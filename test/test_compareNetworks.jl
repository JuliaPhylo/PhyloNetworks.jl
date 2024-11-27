if !(@isdefined doalltests) doalltests = false; end

@testset "test sets: compareNetworks" begin
global net, tree
#----------------------------------------------------------#
#   testing functions to delete edges and nodes            #
#----------------------------------------------------------#
@testset "deletehybridedge!" begin

PN = PhyloNetworks # local to testset

# with nofuse=true
netstr = "(((A:4.0,(B:1.0)#H1:1.1::0.9):0.5,(C:0.6,#H1:1.0::0.1):1.0):3.0,D:5.0);"
net = readnewick(netstr)
@test_logs PN.deletehybridedge!(net, net.edge[6], true);
@test writenewick(net) == "(((A:4.0,(B:1.0)H1:1.1):0.5,(C:0.6):1.0):3.0,D:5.0);"
@test net.edge[3].gamma == 0.9
@test net.node[3].name == "H1"
net = readnewick(netstr)
@test_logs PN.deletehybridedge!(net, net.edge[3], true);
@test writenewick(net) == "(((A:4.0):0.5,(C:0.6,(B:1.0)H1:1.0):1.0):3.0,D:5.0);"
@test  net.edge[5].gamma == 0.1
@test !net.edge[5].hybrid
@test  net.edge[5].ismajor
@test  net.node[3].name == "H1"
@test !net.node[3].hybrid

# example of network with one hybrid edge connected to the root
#  3 edges below the root (1 of them hybrid):
netstr = "((Adif:1.0,(Aech:0.122,#H6:10.0::0.047):10.0):1.614,Aten:1.0,((Asub:1.0,Agem:1.0):0.0)#H6:5.062::0.953);";
net = readnewick(netstr);
@test_logs PN.deletehybridedge!(net, net.edge[10]); # will be rooted
@test writenewick(net) == "((Adif:1.0,(Aech:0.122,(Asub:1.0,Agem:1.0):10.0):10.0):1.614,Aten:1.0);"
net = readnewick(netstr);
@test_logs PN.deletehybridedge!(net, net.edge[10], false, true); # unrooted
@test writenewick(net) == "(Adif:1.0,(Aech:0.122,(Asub:1.0,Agem:1.0):10.0):10.0,Aten:2.614);"
net = readnewick(netstr);
@test_logs PN.deletehybridedge!(net, net.edge[3]);
@test writenewick(net) == "((Adif:1.0,Aech:10.122):1.614,Aten:1.0,(Asub:1.0,Agem:1.0):5.062);"
# 2 edges below the root (1 of them hybrid):
netstr = "((Adif:1.0,(Aech:0.122,#H6:10.0::0.047):10.0):1.614,((Asub:1.0,Agem:1.0):0.0)#H6:5.062::0.953);";
net = readnewick(netstr);
@test_logs PN.deletehybridedge!(net, net.edge[9]);
@test writenewick(net) == "(Adif:1.0,(Aech:0.122,(Asub:1.0,Agem:1.0):10.0):10.0);"
net = readnewick(netstr);
@test_logs PN.deletehybridedge!(net, net.edge[9], true);
@test writenewick(net) == "(Adif:1.0,(Aech:0.122,((Asub:1.0,Agem:1.0):0.0)H6:10.0):10.0);"
net = readnewick(netstr);
@test_logs PN.deletehybridedge!(net, net.edge[9], true, false, false,
    true, true); # last: keeporiginalroot=true to keep the root of degree 1
@test writenewick(net) == "((Adif:1.0,(Aech:0.122,((Asub:1.0,Agem:1.0):0.0)H6:10.0):10.0):1.614);"

if doalltests
net=readnewick("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.1361111):10.0);");
PN.deletehybridedge!(net, net.edge[11]);
writenewick(net) == "(4,((1,2):2.069,(6,5):3.423):0.265,3);" ||
 error("deletehybridedge! didn't work on 11th edge")
net=readnewick("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.1361111):10.0);");
PN.deletehybridedge!(net, net.edge[4]);
writenewick(net) == "(4,((6,5):3.423,1):0.265,(3,2):10.0);" ||
 error("deletehybridedge! didn't work on 4th edge")
end

# example with wrong attributed inChild1
net=readnewick("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.1361111):10.0);");
net.edge[5].ischild1 = false;
@test_logs PN.deletehybridedge!(net, net.edge[4]);
@test writenewick(net) == "(4,((6,5):3.423,1):0.265,(3,2):10.0);"

if doalltests
net = readnewick("((Adif:1.0,(Aech:0.122,#H6:10.0::0.047):10.0):1.614,Aten:1.0,((Asub:1.0,Agem:1.0):0.0)#H6:5.062::0.953);");
net.edge[5].ischild1 = false # edge incident to the root, now pointing towards the root
# tests that deletehybridedge! and writenewick don't use ischild1,
# now in conflict with net.rooti: node -2 instead of -3 as would ischild1 suggest
PN.deletehybridedge!(net, net.edge[10]);
@test writenewick(net) == "((Adif:1.0,(Aech:0.122,(Asub:1.0,Agem:1.0):10.0):10.0):1.614,Aten:1.0);"
end

# example with simplify=false
net0 = readnewick("((((((a:1)#H1:1::.9)#H2:1::.8)#H3:1::.7,#H3:0.5):1,#H2:1):1,(#H1:1,b:1):1,c:1);")
net = deepcopy(net0)
@test writenewick(PN.deletehybridedge!(net, net.edge[5]), round=true) ==
  "((#H1:1.0::0.1,b:1.0):1.0,c:1.0,(a:1.0)#H1:3.0::0.9);"
@test writenewick(PN.deletehybridedge!(net0, net0.edge[5],false,true,false,false), round=true) ==
  "((#H2:1.0::0.2,((a:1.0)#H1:1.0::0.9)#H2:3.0::0.8):1.0,(#H1:1.0::0.1,b:1.0):1.0,c:1.0);"

# level-2 degree-2 blob simplifying to parallel edges, + polytomy below
net0 = readnewick("(africa_east:0.003,((#H3:0::0.003,(non_africa_west:0.2,non_africa_east:0.2)#H1:0.3::1):0.3,(#H1:0::0)#H3:0::0.997)H2:0);");
PN.deletehybridedge!(net0, net0.edge[2], false, false, false, true, false)
@test all(!n.hybrid for n in net0.node)
@test all(e.containroot for e in net0.edge)
@test writenewick(net0) == "(africa_east:0.003,(non_africa_west:0.2,non_africa_east:0.2)H1:0.0);"

end # of testing deletehybridedge!

@testset "deleteleaf! and hardwiredclusterdistance" begin

cui2str = "(Xgordoni,Xmeyeri,(Xcouchianus,(Xvariatus,(Xevelynae,((Xxiphidium,#H25:9.992::0.167):1.383,(Xmilleri,(Xandersi,(Xmaculatus,((((Xhellerii,(Xalvarezi,Xmayae):0.327):0.259,Xsignum):1.866,(Xclemenciae_F2,Xmonticolus):1.461):0.786,((((Xmontezumae,(Xnezahuacoyotl)#H26:0.247::0.807):0.372,((Xbirchmanni_GARC,Xmalinche_CHIC2):1.003,Xcortezi):0.454):0.63,((Xcontinens,Xpygmaeus):1.927,((Xnigrensis,Xmultilineatus):1.304,#H26:0.0::0.193):0.059):2.492):2.034)#H25:0.707::0.833):1.029):0.654):0.469):0.295):0.41):0.646):3.509):0.263);"
cui3str = "(Xmayae,((Xhellerii,(((Xclemenciae_F2,Xmonticolus):1.458,(((((Xmontezumae,(Xnezahuacoyotl)#H26:0.247::0.804):0.375,((Xbirchmanni_GARC,Xmalinche_CHIC2):0.997,Xcortezi):0.455):0.63,(#H26:0.0::0.196,((Xcontinens,Xpygmaeus):1.932,(Xnigrensis,Xmultilineatus):1.401):0.042):2.439):2.0)#H7:0.787::0.835,(Xmaculatus,(Xandersi,(Xmilleri,((Xxiphidium,#H7:9.563::0.165):1.409,(Xevelynae,(Xvariatus,(Xcouchianus,(Xgordoni,Xmeyeri):0.263):3.532):0.642):0.411):0.295):0.468):0.654):1.022):0.788):1.917)#H27:0.149::0.572):0.668,Xalvarezi):0.257,(Xsignum,#H27:1.381::0.428):4.669);"

if doalltests
net3  = readnewick(cui3str);
net2  = readnewick(cui2str);
# major tree, root with outgroup then delete 2 leaves:
tree2 = majortree(net2);        tree3 = majortree(net3);
rootatnode!(tree2,"Xmayae");    rootatnode!(tree3,"Xmayae");
deleteleaf!(tree2,"Xhellerii"); deleteleaf!(tree3,"Xhellerii");
deleteleaf!(tree2,"Xsignum");   deleteleaf!(tree3,"Xsignum");
hardwiredclusterdistance(tree2, tree3, false) == 0 ||
  error("HWDist not 0, major tree, root then prune");
# major tree, delete 2 leaves then root with outgroup:
tree2 = majortree(net2);        tree3 = majortree(net3);
deleteleaf!(tree2,"Xhellerii"); deleteleaf!(tree3,"Xhellerii");
deleteleaf!(tree2,"Xsignum");   deleteleaf!(tree3,"Xsignum");
hardwiredclusterdistance(tree2, tree3, false) == 0 || error("HWD not 0, major tree - 2 taxa");
@test hardwiredclusterdistance(tree2, tree3, true) == 21
rootatnode!(tree3,"Xmaculatus");
hardwiredclusterdistance(tree2, tree3, true) == 15 || error("rooted RF dist not 15");
rootatnode!(tree2,"Xgordoni");
hardwiredclusterdistance(tree2, tree3, true) == 16 || error("rooted RF dist not 16");
hardwiredclusterdistance(tree2, tree3, false) == 0 || error("HWD not 0, major tree - 2 taxa");
end

# network: delete 2 leaves
net2  = readnewick(cui2str);
@test_logs deleteleaf!(net2,"Xhellerii");
@test_logs deleteleaf!(net2,"Xsignum");
@test_logs rootatnode!(net2,"Xmayae");
net3  = readnewick(cui3str);
@test_logs deleteleaf!(net3,"Xhellerii");
@test_logs deleteleaf!(net3,"Xsignum");
@test hardwiredclusterdistance(net2, net3, true) == 3
@test_logs rootatnode!(net3,"Xmayae");
@test hardwiredclusterdistance(net2, net3, true) == 4
@test hardwiredclusterdistance(net2, net3, false) == 4
@test_logs deleteleaf!(net3,"Xmayae"; unroot=true);    #plot(net3);
@test net3.numhybrids == 2
# using simplify=false in deleteleaf!
net3  = readnewick(cui3str);
deleteleaf!(net3,"Xhellerii"); deleteleaf!(net3,"Xsignum");
deleteleaf!(net3,"Xmayae", simplify=false);
@test net3.numhybrids==3
@test net3.numnodes == 47
net3  = readnewick(cui3str);
deleteleaf!(net3,"Xhellerii"); deleteleaf!(net3,"Xsignum");
deleteleaf!(net3,"Xmayae", simplify=false, unroot=true);
@test net3.numhybrids==3
@test net3.numnodes == 46

end # of testset for deleteleaf! and hardwiredclusterdistance

#----------------------------------------------------------#
#   testing functions to display trees / subnetworks       #
#----------------------------------------------------------#

@testset "deletehybridthreshold!" begin

if doalltests
net21 = readnewick("(A,((B,#H1:::0.5),(C,(D)#H1)));");
deletehybridthreshold!(net21,0.2);
writenewick(net21) == "(A,((B,#H1:::0.5),(C,(D)#H1:::0.5)));" ||
 error("deletehybridthreshold! didn't work on net21, gamma=0.2")
deletehybridthreshold!(net21,0.5);
writenewick(net21) == "(A,((C,D),B));" ||
 error("deletehybridthreshold! didn't work on net21, gamma=0.5")

net22 = readnewick("(A,((B,#H1:::0.2),(C,(D)#H1:::0.8)));");
deletehybridthreshold!(net22,0.3);
writenewick(net22) == "(A,((C,D),B));" ||
 error("deletehybridthreshold! didn't work on net22, gamma=0.3")

net31 = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(C:0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
net32 = deepcopy(net31); for e in net32.edge e.length=-1.0; end
# plot(net31)
deletehybridthreshold!(net31,0.3);
writenewick(net31) == "(A:1.0,((C:0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("deletehybridthreshold! didn't work on net31, gamma=0.3")
deletehybridthreshold!(net32,0.3)
writenewick(net32) == "(A,((C,D),B));" ||
 error("deletehybridthreshold! didn't work on net32, gamma=0.3")

net42 = readnewick("(A:1.0,((B:1.1,#H1:0.2):1.2,(C:0.9,(D:0.8)#H1:0.3):1.3):0.7):0.1;");
deletehybridthreshold!(net42,0.5);
writenewick(net42) == "(A:1.0,((C:0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("deletehybridthreshold! didn't work on net42, gamma=0.5")
end

net5 = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
# plot(net5)
@test_logs deletehybridthreshold!(net5,0.5);  # both H1 and H2 eliminated
@test writenewick(net5) == "(A:1.0,((((C:0.52,E:0.52):0.6,F:1.5):0.9,D:1.1):1.3,B:2.3):0.7);"
net5 = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
@test_logs deletehybridthreshold!(net5,0.3);  # H2 remains
@test writenewick(net5) == "(A:1.0,((((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,D:1.1):1.3,B:2.3):0.7);"

end # of testset, deletehybridthreshold!

@testset "displayednetworks! and displayedtrees" begin

net3 = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(C:0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
net31 = PhyloNetworks.displayednetworks!(net3, net3.node[6]); #H1 = 6th node
@test writenewick(net31) == "(A:1.0,((B:1.1,D:1.0):1.2,C:2.2):0.7);"
@test writenewick(net3)  == "(A:1.0,((C:0.9,D:1.1):1.3,B:2.3):0.7);"

if doalltests
net3 = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(C:0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
a = displayedtrees(net3, 0.2);
@test length(a) == 2
@test writenewick(a[1]) == "(A:1.0,((C:0.9,D:1.1):1.3,B:2.3):0.7);"
@test writenewick(a[2]) == "(A:1.0,((B:1.1,D:1.0):1.2,C:2.2):0.7);"
end

net5 = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
a = displayedtrees(net5, 0.5);
@test length(a) == 1
@test writenewick(a[1]) == "(A:1.0,((((C:0.52,E:0.52):0.6,F:1.5):0.9,D:1.1):1.3,B:2.3):0.7);"
a = displayedtrees(net5, 0.1);
@test length(a) == 4
@test writenewick(a[1]) == "(A:1.0,((((C:0.52,E:0.52):0.6,F:1.5):0.9,D:1.1):1.3,B:2.3):0.7);"
@test writenewick(a[2]) == "(A:1.0,((B:1.1,D:1.0):1.2,((C:0.52,E:0.52):0.6,F:1.5):2.2):0.7);"
@test writenewick(a[3]) == "(A:1.0,((((F:0.7,E:0.51):0.8,C:1.12):0.9,D:1.1):1.3,B:2.3):0.7);"
@test writenewick(a[4]) == "(A:1.0,((B:1.1,D:1.0):1.2,((F:0.7,E:0.51):0.8,C:1.12):2.2):0.7);"

net = readnewick("(((A:4.0,(B:1.0)#H1:1.1::0.9):0.5,(C:0.6,#H1:1.0::0.1):1.0):3.0,D:5.0);")
trees = (@test_logs displayedtrees(net,0.0; nofuse=true));
@test writenewick(trees[1])=="(((A:4.0,(B:1.0)H1:1.1):0.5,(C:0.6):1.0):3.0,D:5.0);"
@test writenewick(trees[2])=="(((A:4.0):0.5,(C:0.6,(B:1.0)H1:1.0):1.0):3.0,D:5.0);"
@test PhyloNetworks.inheritanceweight.(trees) ≈ [log(0.9), log(0.1)]

end # of testset, displayednetworks! & displayedtrees

@testset "majortree and displayednetworkat!" begin

net5 = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
@test writenewick(majortree(net5)) == "(A:1.0,((((C:0.52,E:0.52):0.6,F:1.5):0.9,D:1.1):1.3,B:2.3):0.7);"
@test_logs displayednetworkat!(net5, net5.hybrid[1]);
@test writenewick(net5) == "(A:1.0,((((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,D:1.1):1.3,B:2.3):0.7);"
net = readnewick("((((B)#H1)#H2,((D,C,#H2)S1,(#H1,A)S2)S3)S4);") # missing γ's, level 2
@test writenewick(majortree(net)) == "(((D,C)S1,A)S3,B)S4;"
@test writenewick(majortree(net; nofuse=true)) == "(((B)H1)H2,((D,C)S1,(A)S2)S3)S4;"
setgamma!(net.edge[8], 0.8)
@test writenewick(majortree(net)) == "((D,C)S1,(A,B)S2)S3;"
@test writenewick(majortree(net; nofuse=true)) == "((D,C)S1,((B)H1,A)S2)S3;"
@test writenewick(majortree(net; nofuse=true, keeporiginalroot=true)) == "(((D,C)S1,((B)H1,A)S2)S3)S4;"
@test writenewick(majortree(net; keeporiginalroot=true)) == "(((D,C)S1,(A,B)S2)S3)S4;"

# net6 below: hybrid ladder H2 -> H1; and H2 child of root
# using multgammas=true, to test multiplygammas and how it is used
net6 = readnewick("(#H2:::0.2,((C,((B)#H1:::0.6)#H2:::0.8),(#H1,(A1,A2))),O);")
tre6 = displayedtrees(net6, 0.5; nofuse=false, multgammas=true)[1]
@test tre6.edge[2].gamma ≈ 0.48
displayednetworkat!(net6, net6.node[4], false, false, true) # nofuse, unroot=false, multgammas=true
@test net6.edge[3].gamma ≈ 0.6
net6 = readnewick("(#H2:::0.2,((C,((B)#H1:::0.6)#H2:::0.8),(#H1,(A1,A2))),O);")
displayednetworkat!(net6, net6.node[3], false, false, true) # nofuse, unroot=false, multgammas=true
@test net6.edge[5].gamma ≈ 0.4
@test net6.edge[3].gamma ≈ 0.48

end # of testset, majortree & displayednetworkat!

#----------------------------------------------------------#
#   testing functions to compare trees                     #
#----------------------------------------------------------#

@testset "tree2Matrix" begin

if doalltests
net5 = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
tree = displayedtrees(net5, 0.0);
taxa = tiplabels(net5)
M1 = PhyloNetworks.tree2Matrix(tree[1], taxa, rooted=false);
M2 = PhyloNetworks.tree2Matrix(tree[2], taxa, rooted=false);
M1 ==
[15 0 0 1 1 1 1;
 12 0 0 1 1 1 0;
  8 0 0 1 1 0 0] || error("bad M1, from net5")
M2 ==
[ 4 0 1 0 0 0 1;
 12 0 0 1 1 1 0;
  8 0 0 1 1 0 0] || error("bad M2, from net5")
hardwiredclusterdistance(tree[1], tree[2], true) == 2 || error("wrong dist between rooted trees 1 and 2");
hardwiredclusterdistance(tree[1], tree[2], false)== 2 || error("wrong dist between unrooted trees 1 and 2");
hardwiredclusters(net5, taxa) ==
[16 0 1 1 1 1 1 10;
  4 0 1 0 0 0 1 10;
  3 0 0 0 0 0 1 11;
 15 0 0 1 1 1 1 10;
 12 0 0 1 1 1 0 10;
  8 0 0 1 1 0 0 10;
  7 0 0 0 1 0 0 11;
 11 0 0 0 1 1 0 10] || error("wrong matrix of hardwired clusters for net5");
hardwiredclusterdistance(net5, tree[2], true) == 4 || error("wrong dist between net5 and tree 2");
end

## check things with R
# library(ape); set.seed(15); phy <- rmtree(10,8); write.tree(phy,"tenRtrees.tre")
# library(phangorn); RF.dist(phy,rooted=F)
#     1  2  3  4  5  6  7  8  9
# 2   8
# 3  10 10
# 4  10 10 10
# 5  10  8  8 10
# 6  10 10  8 10  8
# 7  10 10 10 10 10 10
# 8  10  8  8 10 10  8 10
# 9  10 10  8  8 10 10 10 10
# 10  8  6 10  8  8 10 10 10  8
# i=10; j=1; RF.dist(phy[[i]],phy[[j]],rooted=T) # 10
# i=10; j=2; RF.dist(phy[[i]],phy[[j]],rooted=T) #  8
## now checking if we get the same results in julia
# phy = readmultinewick("tenRtrees.tre")
# or
phy1 = readnewick("((t8:0.8623136566,(((t6:0.141187073,t2:0.7767125128):0.9646669542,t4:0.8037273993):0.447443719,t5:0.7933459524):0.8417851452):0.7066285675,(t1:0.0580010619,(t7:0.6590069213,t3:0.1069735419):0.5657461432):0.3575631182);");
phy2 = readnewick("((t3:0.9152618761,t4:0.4574306419):0.7603277895,(((t1:0.4291725352,t8:0.3302786439):0.3437780738,(t5:0.8438980761,(t6:0.6667000714,t2:0.7141199473):0.01087239943):0.752832541):0.2591188031,t7:0.7685037958):0.9210739341);");
phy3 = readnewick("(((t7:0.3309174306,t6:0.8330178803):0.7741786113,(((t2:0.4048132468,t8:0.6809111023):0.6810255498,(t4:0.6540613638,t5:0.2610215396):0.8490990005):0.6802781771,t3:0.2325445588):0.911911567):0.94644987,t1:0.09404937108);");
phy10= readnewick("((t4:0.1083955287,((t1:0.8376079942,t8:0.1745392387):0.6178579947,((t6:0.3196466176,t2:0.9228881211):0.3112748025,t7:0.05162345758):0.7137957355):0.5162231021):0.06693460606,(t5:0.005652675638,t3:0.2584615161):0.7333540542);");

@test hardwiredclusterdistance(phy1, phy10, false)== 8
@test hardwiredclusterdistance(phy2, phy10, false)== 6
@test hardwiredclusterdistance(phy3, phy10, false)==10
@test hardwiredclusterdistance(phy1, phy2,  false)== 8
@test hardwiredclusterdistance(phy1, phy3,  false)==10
@test hardwiredclusterdistance(phy2, phy3,  false)==10
@test hardwiredclusterdistance(phy1, phy10, true) ==10
@test hardwiredclusterdistance(phy2, phy10, true) == 8
# or: wrong RF distance between some of the trees

end # of testset, tree2Matrix

@testset "ladderpartition" begin
tree = readnewick("(O,A,((B1,B2),(E,(C,D))));")
PhyloNetworks.resetnodenumbers!(tree; checkpreorder=true, type=:postorder)
below, above = PhyloNetworks.ladderpartition(tree)
@test below == [[[1]],[[2]],[[3]],[[4]],[[5]],[[6]],[[7]],
  [[3],[4]], [[6],[7]], [[5],[6,7]], [[3,4],[5,6,7]], [[1],[2],[3,4,5,6,7]]]
@test above == [[[2],[3,4,5,6,7]], [[1],[3, 4, 5, 6, 7]], [[4], [5,6,7],[1],[2]],
  [[3],[5,6,7],[1],[2]], [[6,7],[3,4],[1],[2]], [[7],[5],[3,4],[1],[2]],
  [[6],[5],[3,4],[1],[2]], [[5,6,7],[1],[2]], [[5],[3,4],[1],[2]], [[3,4],[1],[2]], [[1],[2]], []]
end

#----------------------------------------------------------#
#   testing function to compare networks                   #
#   with hardwired clusters                                #
#   used for detection of given hybridization event        #
#----------------------------------------------------------#

@testset "displayedtrees, hardwiredclusters, hardwiredclusterdistance, displayednetworkat!" begin

estnet = readnewick("(6,((5,#H7:0.0::0.402):8.735,((1,2):6.107,((3,4):1.069)#H7:9.509::0.598):6.029):0.752);")
# originally from "../msSNaQ/simulations/estimatedNetworks/baseline/nloci10/1_julia.out"
trunet = readnewick("((((1,2),((3,4))#H1),(#H1,5)),6);");
@test hardwiredclusterdistance(majortree(trunet), majortree(estnet),false) == 0 # false: unrooted
truminor = minortreeat(trunet, 1); # (1:1.0,2:1.0,((5:1.0,(3:1.0,4:1.0):2.0):1.0,6:1.0):2.0);
estminor = minortreeat(estnet, 1); # (5:1.0,(3:1.0,4:1.0):1.069,(6:1.0,(1:1.0,2:1.0):10.0):8.735);
@test writenewick(truminor) == "(((5,(3,4)),(1,2)),6);"
@test writenewick(estminor) == "(6,((5,(3,4):1.069):8.735,(1,2):12.136):0.752);"
@test hardwiredclusterdistance(truminor, estminor, false) == 0 # false: unrooted
# so the hybrid edge was estimated correctly!!
rootatnode!(trunet, -8)
@test hardwiredclusterdistance(estnet, trunet, true) == 3
# next: testing hardwiredclusterdistance_unrooted, via the option rooted=false
@test hardwiredclusterdistance(estnet, trunet, false) == 0
h0est = readnewick("(((2:0.01,1:0.01):0.033,(3:0.0154,4:0.0149):0.0186):0.0113,6:0.0742,5:0.0465);")
truenet = readnewick("((((1,2),((3,4))#H1),(#H1,5)),6);")
h1est = readnewick("(5:0.0,6:0.0,(((2:0.0)#H1:0.0::0.95,1:0.0):0.0,((4:0.0,3:0.0):0.0,#H1:0.0::0.05):0.0):0.0);")
@test hardwiredclusterdistance(h0est, truenet, false) == 2
@test hardwiredclusterdistance(truenet, h0est, false) == 2
@test hardwiredclusterdistance(h1est, truenet, false) == 4
@test hardwiredclusterdistance(truenet, h1est, false) == 4

net5 = readnewick("(A,((B,#H1:::0.2),(((C,(E)#H2:::0.7),(#H2:::0.3,F)),(D)#H1:::0.8)));");
tree = displayedtrees(net5, 0.0);
taxa = tiplabels(net5)
@test hardwiredclusters(tree[1], taxa) ==
[16 0 1 1 1 1 1 10;
 15 0 0 1 1 1 1 10;
 12 0 0 1 1 1 0 10;
  8 0 0 1 1 0 0 10]
@test hardwiredclusters(tree[2], taxa) ==
[16 0 1 1 1 1 1 10;
  4 0 1 0 0 0 1 10;
 12 0 0 1 1 1 0 10;
  8 0 0 1 1 0 0 10]
@test hardwiredclusters(net5, taxa) ==
[16 0 1 1 1 1 1 10;
  4 0 1 0 0 0 1 10;
  3 0 0 0 0 0 1 11;
 15 0 0 1 1 1 1 10;
 12 0 0 1 1 1 0 10;
  8 0 0 1 1 0 0 10;
  7 0 0 0 1 0 0 11;
 11 0 0 0 1 1 0 10]

if doalltests
trunet = readnewick("(((1,2),((3,4))#H1),(#H1,5),6);"); # unrooted
taxa = tiplabels(trunet)
hardwiredclusters(trunet, taxa) ==
[8 1 1 1 1 0 0 10
 3 1 1 0 0 0 0 10
 7 0 0 1 1 0 0 11
 6 0 0 1 1 0 0 10
11 0 0 1 1 1 0 10] || error("wrong hardwired cluster matrix for unrooted trunet");
trunet = readnewick("((((1,2),((3,4))#H1),(#H1,5)),6);"); # rooted: good!
hardwiredclusters(trunet, taxa) ==
[12 1 1 1 1 1 0 10;
  8 1 1 1 1 0 0 10;
  3 1 1 0 0 0 0 10;
  7 0 0 1 1 0 0 11;
  6 0 0 1 1 0 0 10;
 11 0 0 1 1 1 0 10] || error("wrong hardwired cluster matrix for trunet");
hardwiredclusters(estnet, taxa) ==
[13 1 1 1 1 1 0 10
  4 0 0 1 1 1 0 10
  3 0 0 1 1 0 0 11
 10 0 0 1 1 0 0 10
 12 1 1 1 1 0 0 10
  7 1 1 0 0 0 0 10] || error("wrong hardwired cluster matrix for estnet")
hardwiredclusterdistance(trunet,estnet,true) == 0 ||
 error("trunet and estnet should be found to be at HWDist 0");

net51 = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;")
directedges!(net51); # doing this avoids a warning on the next line:
displayednetworkat!(net51, net51.hybrid[1]) # H2
# "WARNING: node -3 being the root is contradicted by ischild1 of its edges."
rootatnode!(net51, "A");
writenewick(net51) == "(A:1.0,((((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,D:1.1):1.3,B:2.3):0.7);" ||
  error("wrong net51 after displayednetworkat!");
net52 = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
displayednetworkat!(net52, net52.hybrid[2])
writenewick(net52) == "(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,E:0.52):0.6,F:1.5):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7);" ||
  error("wrong net52 after displayednetworkat!");
taxa = tiplabels(net52) # order: A B C E F D
hardwiredclusters(net51, taxa) ==
[16 0 1 1 1 1 1 10;
 15 0 0 1 1 1 1 10;
 12 0 0 1 1 1 0 10;
  8 0 0 1 1 0 0 10;
  7 0 0 0 1 0 0 11;
 11 0 0 0 1 1 0 10] || error("wrong hardwired clusters for net51");
hardwiredclusters(net52, taxa) ==
[16 0 1 1 1 1 1 10;
  4 0 1 0 0 0 1 10;
  3 0 0 0 0 0 1 11;
 15 0 0 1 1 1 1 10;
 12 0 0 1 1 1 0 10;
  8 0 0 1 1 0 0 10] || error("wrong hardwired clusters for net52");
hardwiredclusterdistance(net51,net52,true) == 4 ||
 error("wrong HWDist between net51 and net52");
end

end # of testset: displayedtrees, hardwiredclusters, hardwiredclusterdistance, displayednetworkat!

@testset "hardwiredcluster! on single nodes" begin

net5 = "(A,((B,#H1),(((C,(E)#H2),(#H2,F)),(D)#H1)));" |> readnewick |> directedges! ;
taxa = net5 |> tiplabels # ABC EF D
m = hcat([true,false,false,false,false,false],
[false,true,false,false,false,false],
[false,false,false,false,false,true],
[false,true,false,false,false,true],
[false,false,true,false,false,false],
[false,false,false,true,false,false],
[false,false,false,true,false,false],
[false,false,true,true,false,false],
[false,false,false,true,false,false],
[false,false,false,false,true,false],
[false,false,false,true,true,false],
[false,false,true,true,true,false],
[false,false,false,false,false,true],
[false,false,false,false,false,true],
[false,false,true,true,true,true],
[false,true,true,true,true,true])
for i = 1:16
  @test hardwiredcluster(net5.edge[i], taxa) == m[:,i]
end
end # of testset, hardwiredcluster! on single nodes

end
