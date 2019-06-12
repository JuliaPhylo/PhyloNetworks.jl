# test of deleteHybridEdge!, functions to extract displayed trees/subnetworks,
#      used to compare networks with the hardwired cluster distance.
# Cecile March 2016

if !(@isdefined doalltests) doalltests = false; end
@testset "test sets: compareNetworks" begin
global net, tree
#----------------------------------------------------------#
#   testing functions to delete edges and nodes            #
#----------------------------------------------------------#
@testset "testing deleteHybridEdge!" begin

# with keepNodes=true
netstr = "(((A:4.0,(B:1.0)#H1:1.1::0.9):0.5,(C:0.6,#H1:1.0::0.1):1.0):3.0,D:5.0);"
net = readTopology(netstr)
@test_logs PhyloNetworks.deleteHybridEdge!(net, net.edge[6], true);
@test writeTopology(net) == "(((A:4.0,(B:1.0)H1:1.1):0.5,(C:0.6):1.0):3.0,D:5.0);"
@test net.edge[3].gamma == 0.9
@test net.node[3].name == "H1"
net = readTopology(netstr)
@test_logs PhyloNetworks.deleteHybridEdge!(net, net.edge[3], true);
@test writeTopology(net) == "(((A:4.0):0.5,(C:0.6,(B:1.0)H1:1.0):1.0):3.0,D:5.0);"
@test  net.edge[5].gamma == 0.1
@test !net.edge[5].hybrid
@test  net.edge[5].isMajor
@test  net.node[3].name == "H1"
@test !net.node[3].hybrid

# example of network with one hybrid edge connected to the root
#  3 edges below the root (1 of them hybrid):
netstr = "((Adif:1.0,(Aech:0.122,#H6:10.0::0.047):10.0):1.614,Aten:1.0,((Asub:1.0,Agem:1.0):0.0)#H6:5.062::0.953);";
net = readTopology(netstr);
@test_logs PhyloNetworks.deleteHybridEdge!(net, net.edge[10]);
@test writeTopology(net) == "(Adif:1.0,(Aech:0.122,(Asub:1.0,Agem:1.0):10.0):10.0,Aten:2.614);"
net = readTopology(netstr);
@test_logs PhyloNetworks.deleteHybridEdge!(net, net.edge[3]);
@test writeTopology(net) == "((Adif:1.0,Aech:10.122):1.614,Aten:1.0,(Asub:1.0,Agem:1.0):5.062);"
# 2 edges below the root (1 of them hybrid):
netstr = "((Adif:1.0,(Aech:0.122,#H6:10.0::0.047):10.0):1.614,((Asub:1.0,Agem:1.0):0.0)#H6:5.062::0.953);";
net = readTopology(netstr);
@test_logs PhyloNetworks.deleteHybridEdge!(net, net.edge[9]);
@test writeTopology(net) == "(Adif:1.0,(Aech:0.122,(Asub:1.0,Agem:1.0):10.0):10.0);"
net = readTopology(netstr);
@test_logs PhyloNetworks.deleteHybridEdge!(net, net.edge[9], true);
@test writeTopology(net) == "((Adif:1.0,(Aech:0.122,((Asub:1.0,Agem:1.0):0.0)H6:10.0):10.0):1.614);"

if doalltests
net=readTopology("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.1361111):10.0);");
# plot(net, showEdgeNumber=true, showNodeNumber=true)
deleteHybridEdge!(net, net.edge[11]);
writeTopologyLevel1(net) == "(4,((1,2):2.069,(6,5):3.423):0.265,3);" ||
 error("deleteHybridEdge! didn't work on 11th edge")
net=readTopology("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.1361111):10.0);");
deleteHybridEdge!(net, net.edge[4]);
writeTopologyLevel1(net) == "(4,((6,5):3.423,1):0.265,(3,2):10.0);" ||
 error("deleteHybridEdge! didn't work on 4th edge")
end

# example with wrong attributed inChild1
net=readTopology("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.1361111):10.0);");
net.edge[5].isChild1 = false;
@test_logs deleteHybridEdge!(net, net.edge[4]);
@test writeTopologyLevel1(net) == "(4,((6,5):3.423,1):0.265,(3,2):10.0);"
# or: deleteHybridEdge! didn't work on 4th edge when isChild1 was outdated

if doalltests
net = readTopology("((Adif:1.0,(Aech:0.122,#H6:10.0::0.047):10.0):1.614,Aten:1.0,((Asub:1.0,Agem:1.0):0.0)#H6:5.062::0.953);");
net.edge[5].isChild1 = false # edge 5 from -1 to -2
deleteHybridEdge!(net, net.edge[10]);
println("a warning is expected: \"node -1 being the root is contradicted by isChild1 of its edges.\"")
writeTopologyLevel1(net) == "(Adif:1.0,(Aech:0.122,(Asub:1.0,Agem:1.0):10.0):10.0,Aten:2.614);" ||
 error("deleteHybridEdge! didn't work on 10th edge after isChild1 was changed")
# plot(net, showEdgeNumber=true, showNodeNumber=true)
end

end # of testing deleteHybridEdge!

@testset "testing deleteleaf! and hardwiredClusterDistance" begin

cui2str = "(Xgordoni,Xmeyeri,(Xcouchianus,(Xvariatus,(Xevelynae,((Xxiphidium,#H25:9.992::0.167):1.383,(Xmilleri,(Xandersi,(Xmaculatus,((((Xhellerii,(Xalvarezi,Xmayae):0.327):0.259,Xsignum):1.866,(Xclemenciae_F2,Xmonticolus):1.461):0.786,((((Xmontezumae,(Xnezahuacoyotl)#H26:0.247::0.807):0.372,((Xbirchmanni_GARC,Xmalinche_CHIC2):1.003,Xcortezi):0.454):0.63,((Xcontinens,Xpygmaeus):1.927,((Xnigrensis,Xmultilineatus):1.304,#H26:0.0::0.193):0.059):2.492):2.034)#H25:0.707::0.833):1.029):0.654):0.469):0.295):0.41):0.646):3.509):0.263);"
cui3str = "(Xmayae,((Xhellerii,(((Xclemenciae_F2,Xmonticolus):1.458,(((((Xmontezumae,(Xnezahuacoyotl)#H26:0.247::0.804):0.375,((Xbirchmanni_GARC,Xmalinche_CHIC2):0.997,Xcortezi):0.455):0.63,(#H26:0.0::0.196,((Xcontinens,Xpygmaeus):1.932,(Xnigrensis,Xmultilineatus):1.401):0.042):2.439):2.0)#H7:0.787::0.835,(Xmaculatus,(Xandersi,(Xmilleri,((Xxiphidium,#H7:9.563::0.165):1.409,(Xevelynae,(Xvariatus,(Xcouchianus,(Xgordoni,Xmeyeri):0.263):3.532):0.642):0.411):0.295):0.468):0.654):1.022):0.788):1.917)#H27:0.149::0.572):0.668,Xalvarezi):0.257,(Xsignum,#H27:1.381::0.428):4.669);"

net3  = readTopology(cui3str);
net2  = readTopology(cui2str);

if doalltests
# major tree, root with outgroup then delete 2 leaves:
tree2 = majorTree(net2);        tree3 = majorTree(net3);
rootatnode!(tree2,"Xmayae");    rootatnode!(tree3,"Xmayae");
deleteleaf!(tree2,"Xhellerii"); deleteleaf!(tree3,"Xhellerii");
deleteleaf!(tree2,"Xsignum");   deleteleaf!(tree3,"Xsignum");
hardwiredClusterDistance(tree2, tree3, false) == 0 ||
  error("HWDist not 0, major tree, root then prune");
# major tree, delete 2 leaves then root with outgroup:
tree2 = majorTree(net2);        tree3 = majorTree(net3);
deleteleaf!(tree2,"Xhellerii"); deleteleaf!(tree3,"Xhellerii");
deleteleaf!(tree2,"Xsignum");   deleteleaf!(tree3,"Xsignum");
hardwiredClusterDistance(tree2, tree3, false) == 0 || error("HWD not 0, major tree - 2 taxa");
hardwiredClusterDistance(tree2, tree3, true) == 20 || error("rooted RF dist not 20");
rootatnode!(tree3,"Xmaculatus");
hardwiredClusterDistance(tree2, tree3, true) == 15 || error("rooted RF dist not 15");
rootatnode!(tree2,"Xgordoni");
hardwiredClusterDistance(tree2, tree3, true) == 16 || error("rooted RF dist not 16");
hardwiredClusterDistance(tree2, tree3, false) == 0 || error("HWD not 0, major tree - 2 taxa");
end

# network: delete 2 leaves
net2  = readTopology(cui2str);
@test_logs deleteleaf!(net2,"Xhellerii");
@test_logs deleteleaf!(net2,"Xsignum");
# earlier warning: """node 13 is a leaf. Will create a new node if needed, to set taxon "Xmayae" as outgroup."""
@test_logs rootatnode!(net2,"Xmayae");
net3  = readTopology(cui3str);
@test_logs deleteleaf!(net3,"Xhellerii");
@test_logs deleteleaf!(net3,"Xsignum");
# earlier warning: """node 13 is a leaf. Will create a new node if needed, to set taxon "Xmayae" as outgroup."""
@test_logs rootatnode!(net2,"Xmayae");
@test hardwiredClusterDistance(net2, net3, true) == 3
@test_logs deleteleaf!(net3,"Xmayae");    #plot(net3);
@test net3.numHybrids == 2
# using simplify=false in deleteleaf!
net3  = readTopology(cui3str);
deleteleaf!(net3,"Xhellerii"); deleteleaf!(net3,"Xsignum");
deleteleaf!(net3,"Xmayae", simplify=false);
@test net3.numHybrids==3 # or: deleteleaf wrong on mayae with simplify=false
# plot(net3); # looks weird though: k=2 cycle at the root. 3 root edges:
# one to a leaf, 1 major & 1 minor hybrid edge to the same child.

end # of testset for deleteleaf! and hardwiredClusterDistance

#----------------------------------------------------------#
#   testing functions to display trees / subnetworks       #
#----------------------------------------------------------#

@testset "testing deleteHybridThreshold!" begin

if doalltests
net21 = readTopology("(A,((B,#H1),(C,(D)#H1)));");
# manual bug fix to get good gamma's (0.5 not 1.0) and major/minor
net21.edge[3].gamma = 0.5;
net21.edge[7].gamma = 0.5;
deleteHybridThreshold!(net21,0.2);
writeTopologyLevel1(net21) == "(A,((B,#H1:::0.5),(C,(D)#H1:::0.5)));" ||
 error("deleteHybridThreshold! didn't work on net21, gamma=0.2")
deleteHybridThreshold!(net21,0.5);
writeTopologyLevel1(net21) == "(A,((C,D),B));" ||
 error("deleteHybridThreshold! didn't work on net21, gamma=0.5")

net22 = readTopology("(A,((B,#H1:::0.2),(C,(D)#H1:::0.8)));");
deleteHybridThreshold!(net22,0.3);
writeTopologyLevel1(net22) == "(A,((C,D),B));" ||
 error("deleteHybridThreshold! didn't work on net22, gamma=0.3")

net31 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(C:0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
net32 = deepcopy(net31); for e in net32.edge e.length=-1.0; end
# plot(net31)
deleteHybridThreshold!(net31,0.3);
writeTopologyLevel1(net31) == "(A:1.0,((C:0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("deleteHybridThreshold! didn't work on net31, gamma=0.3")
deleteHybridThreshold!(net32,0.3)
writeTopologyLevel1(net32) == "(A,((C,D),B));" ||
 error("deleteHybridThreshold! didn't work on net32, gamma=0.3")

net42 = readTopology("(A:1.0,((B:1.1,#H1:0.2):1.2,(C:0.9,(D:0.8)#H1:0.3):1.3):0.7):0.1;");
deleteHybridThreshold!(net42,0.5);
writeTopologyLevel1(net42) == "(A:1.0,((C:0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("deleteHybridThreshold! didn't work on net42, gamma=0.5")
end

net5 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
# plot(net5)
@test_logs deleteHybridThreshold!(net5,0.5);  # both H1 and H2 eliminated
@test writeTopologyLevel1(net5) == "(A:1.0,((((C:0.52,E:0.52):0.6,F:1.5):0.9,D:1.1):1.3,B:2.3):0.7);"
# or: deleteHybridThreshold! didn't work on net5, gamma=0.5
net5 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
@test_logs deleteHybridThreshold!(net5,0.3);  # H2 remains
@test writeTopologyLevel1(net5) == "(A:1.0,((((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,D:1.1):1.3,B:2.3):0.7);"
# or: deleteHybridThreshold! didn't work on net5, gamma=0.3

end # of testset, deleteHybridThreshold!

@testset "testing displayedNetworks! and displayedTrees" begin

net3 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(C:0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
net31 = displayedNetworks!(net3, net3.node[6]); #H1 = 6th node
@test writeTopologyLevel1(net31) == "(A:1.0,((B:1.1,D:1.0):1.2,C:2.2):0.7);"
# or: displayedNetworks! didn't work on net3, minor at 6th node
@test writeTopologyLevel1(net3)  == "(A:1.0,((C:0.9,D:1.1):1.3,B:2.3):0.7);"
# or: displayedNetworks! didn't work on net3, major at 6th node

if doalltests
net3 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(C:0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
a = displayedTrees(net3, 0.2);
length(a) == 2 ||
 error("displayedTrees didn't work on net3, gamma=0.2: output not of length 2")
writeTopologyLevel1(a[1]) == "(A:1.0,((C:0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("displayedTrees didn't work on net3, gamma=0.2: 1st tree wrong")
writeTopologyLevel1(a[2]) == "(A:1.0,((B:1.1,D:1.0):1.2,C:2.2):0.7);" ||
 error("displayedTrees didn't work on net3, gamma=0.2: 2nd tree wrong")
end

net5 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
a = displayedTrees(net5, 0.5);
@test length(a) == 1 # or: displayedTrees didn't work on net5, gamma=0.5
@test writeTopologyLevel1(a[1]) == "(A:1.0,((((C:0.52,E:0.52):0.6,F:1.5):0.9,D:1.1):1.3,B:2.3):0.7);"
# or: displayedTrees didn't work on net5, gamma=0.5
a = displayedTrees(net5, 0.1);
@test length(a) == 4 # or: displayedTrees didn't work on net5, gamma=0.1
@test writeTopologyLevel1(a[1]) == "(A:1.0,((((C:0.52,E:0.52):0.6,F:1.5):0.9,D:1.1):1.3,B:2.3):0.7);"
@test writeTopologyLevel1(a[2]) == "(A:1.0,((B:1.1,D:1.0):1.2,((C:0.52,E:0.52):0.6,F:1.5):2.2):0.7);"
@test writeTopologyLevel1(a[3]) == "(A:1.0,((((F:0.7,E:0.51):0.8,C:1.12):0.9,D:1.1):1.3,B:2.3):0.7);"
@test writeTopologyLevel1(a[4]) == "(A:1.0,((B:1.1,D:1.0):1.2,((F:0.7,E:0.51):0.8,C:1.12):2.2):0.7);"

net = readTopology("(((A:4.0,(B:1.0)#H1:1.1::0.9):0.5,(C:0.6,#H1:1.0::0.1):1.0):3.0,D:5.0);")
trees = (@test_logs displayedTrees(net,0.0; keepNodes=true));
@test writeTopology(trees[1])=="(((A:4.0,(B:1.0)H1:1.1):0.5,(C:0.6):1.0):3.0,D:5.0);"
@test writeTopology(trees[2])=="(((A:4.0):0.5,(C:0.6,(B:1.0)H1:1.0):1.0):3.0,D:5.0);"
@test PhyloNetworks.inheritanceWeight.(trees) ≈ [log(0.9), log(0.1)]

end # of testset, displayedNetworks! & displayedTrees

@testset "testing majorTree and displayedNetworkAt!" begin

net5 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
@test writeTopology(majorTree(net5)) == "(A:1.0,((((C:0.52,E:0.52):0.6,F:1.5):0.9,D:1.1):1.3,B:2.3):0.7);"
@test_logs displayedNetworkAt!(net5, net5.hybrid[1]);
@test writeTopology(net5) == "(A:1.0,((((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,D:1.1):1.3,B:2.3):0.7);"
net = readTopology("((((B)#H1)#H2,((D,C,#H2)S1,(#H1,A)S2)S3)S4);") # missing γ's, level 2
@test writeTopology(majorTree(net)) == "(((D,C)S1,A)S3,B)S4;"
setGamma!(net.edge[8], 0.8)
@test writeTopology(majorTree(net)) == "((D,C)S1,(A,B)S2)S3;"

end # of testset, majorTree & displayedNetworkAt!

#----------------------------------------------------------#
#   testing functions to compare trees                     #
#----------------------------------------------------------#

@testset "testing tree2Matrix" begin

if doalltests
net5 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
tree = displayedTrees(net5, 0.0);
taxa = tipLabels(net5);
M1 = tree2Matrix(tree[1], taxa, rooted=false);
M2 = tree2Matrix(tree[2], taxa, rooted=false);
M1 ==
[15 0 0 1 1 1 1;
 12 0 0 1 1 1 0;
  8 0 0 1 1 0 0] || error("bad M1, from net5")
M2 ==
[ 4 0 1 0 0 0 1;
 12 0 0 1 1 1 0;
  8 0 0 1 1 0 0] || error("bad M2, from net5")
hardwiredClusterDistance(tree[1], tree[2], true) == 2 || error("wrong dist between rooted trees 1 and 2");
hardwiredClusterDistance(tree[1], tree[2], false)== 2 || error("wrong dist between unrooted trees 1 and 2");
hardwiredClusters(net5, taxa) ==
[16 0 1 1 1 1 1 10;
  4 0 1 0 0 0 1 10;
  3 0 0 0 0 0 1 11;
 15 0 0 1 1 1 1 10;
 12 0 0 1 1 1 0 10;
  8 0 0 1 1 0 0 10;
  7 0 0 0 1 0 0 11;
 11 0 0 0 1 1 0 10] || error("wrong matrix of hardwired clusters for net5");
hardwiredClusterDistance(net5, tree[2], true) == 4 || error("wrong dist between net5 and tree 2");
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
# phy = readInputTrees("tenRtrees.tre")
# or
phy1 = readTopology("((t8:0.8623136566,(((t6:0.141187073,t2:0.7767125128):0.9646669542,t4:0.8037273993):0.447443719,t5:0.7933459524):0.8417851452):0.7066285675,(t1:0.0580010619,(t7:0.6590069213,t3:0.1069735419):0.5657461432):0.3575631182);");
phy2 = readTopology("((t3:0.9152618761,t4:0.4574306419):0.7603277895,(((t1:0.4291725352,t8:0.3302786439):0.3437780738,(t5:0.8438980761,(t6:0.6667000714,t2:0.7141199473):0.01087239943):0.752832541):0.2591188031,t7:0.7685037958):0.9210739341);");
phy3 = readTopology("(((t7:0.3309174306,t6:0.8330178803):0.7741786113,(((t2:0.4048132468,t8:0.6809111023):0.6810255498,(t4:0.6540613638,t5:0.2610215396):0.8490990005):0.6802781771,t3:0.2325445588):0.911911567):0.94644987,t1:0.09404937108);");
phy10= readTopology("((t4:0.1083955287,((t1:0.8376079942,t8:0.1745392387):0.6178579947,((t6:0.3196466176,t2:0.9228881211):0.3112748025,t7:0.05162345758):0.7137957355):0.5162231021):0.06693460606,(t5:0.005652675638,t3:0.2584615161):0.7333540542);");

@test hardwiredClusterDistance(phy1, phy10, false)== 8
@test hardwiredClusterDistance(phy2, phy10, false)== 6
@test hardwiredClusterDistance(phy3, phy10, false)==10
@test hardwiredClusterDistance(phy1, phy2,  false)== 8
@test hardwiredClusterDistance(phy1, phy3,  false)==10
@test hardwiredClusterDistance(phy2, phy3,  false)==10
@test hardwiredClusterDistance(phy1, phy10, true) ==10
@test hardwiredClusterDistance(phy2, phy10, true) == 8
# or: wrong RF distance between some of the trees

end # of testset, tree2Matrix

#----------------------------------------------------------#
#   testing function to compare networks                   #
#   with hardwired clusters                                #
#   used for detection of given hybridization event        #
#----------------------------------------------------------#

@testset "test displayedTrees, hardwiredClusters, hardwiredClusterDistance, displayedNetworkAt!" begin

estnet = readTopology("(6,((5,#H7:0.0::0.402):8.735,((1,2):6.107,((3,4):1.069)#H7:9.509::0.598):6.029):0.752);")
# originally from "../msSNaQ/simulations/estimatedNetworks/baseline/nloci10/1_julia.out"
trunet = readTopology("((((1,2),((3,4))#H1),(#H1,5)),6);");
@test hardwiredClusterDistance(majorTree(trunet), majorTree(estnet),false) == 0 # false: unrooted
truminor = minorTreeAt(trunet, 1); # (1:1.0,2:1.0,((5:1.0,(3:1.0,4:1.0):2.0):1.0,6:1.0):2.0);
estminor = minorTreeAt(estnet, 1); # (5:1.0,(3:1.0,4:1.0):1.069,(6:1.0,(1:1.0,2:1.0):10.0):8.735);
@test writeTopology(truminor) == "((((1,2),(3,4)),5),6);"
@test writeTopology(estminor) == "(6,(((1,2):6.107,(3,4):10.578):6.029,5):0.752);"
@test hardwiredClusterDistance(truminor, estminor, false) == 0 # false: unrooted
# so the hybrid edge was estimated correctly!!

net5 = readTopology("(A,((B,#H1:::0.2),(((C,(E)#H2:::0.7),(#H2:::0.3,F)),(D)#H1:::0.8)));");
tree = displayedTrees(net5, 0.0);
taxa = tipLabels(net5);
@test hardwiredClusters(tree[1], taxa) ==
[16 0 1 1 1 1 1 10;
 15 0 0 1 1 1 1 10;
 12 0 0 1 1 1 0 10;
  8 0 0 1 1 0 0 10]
@test hardwiredClusters(tree[2], taxa) ==
[16 0 1 1 1 1 1 10;
  4 0 1 0 0 0 1 10;
 12 0 0 1 1 1 0 10;
  8 0 0 1 1 0 0 10]
@test hardwiredClusters(net5, taxa) ==
[16 0 1 1 1 1 1 10;
  4 0 1 0 0 0 1 10;
  3 0 0 0 0 0 1 11;
 15 0 0 1 1 1 1 10;
 12 0 0 1 1 1 0 10;
  8 0 0 1 1 0 0 10;
  7 0 0 0 1 0 0 11;
 11 0 0 0 1 1 0 10]

if doalltests
trunet = readTopology("(((1,2),((3,4))#H1),(#H1,5),6);"); # unrooted
taxa = tipLabels(trunet);
hardwiredClusters(trunet, taxa) ==
[8 1 1 1 1 0 0 10
 3 1 1 0 0 0 0 10
 7 0 0 1 1 0 0 11
 6 0 0 1 1 0 0 10
11 0 0 1 1 1 0 10] || error("wrong hardwired cluster matrix for unrooted trunet");
trunet = readTopology("((((1,2),((3,4))#H1),(#H1,5)),6);"); # rooted: good!
hardwiredClusters(trunet, taxa) ==
[12 1 1 1 1 1 0 10;
  8 1 1 1 1 0 0 10;
  3 1 1 0 0 0 0 10;
  7 0 0 1 1 0 0 11;
  6 0 0 1 1 0 0 10;
 11 0 0 1 1 1 0 10] || error("wrong hardwired cluster matrix for trunet");
hardwiredClusters(estnet, taxa) ==
[13 1 1 1 1 1 0 10
  4 0 0 1 1 1 0 10
  3 0 0 1 1 0 0 11
 10 0 0 1 1 0 0 10
 12 1 1 1 1 0 0 10
  7 1 1 0 0 0 0 10] || error("wrong hardwired cluster matrix for estnet")
hardwiredClusterDistance(trunet,estnet,true) == 0 ||
 error("trunet and estnet should be found to be at HWDist 0");

net51 = readTopologyLevel1("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;")
directEdges!(net51); # doing this avoids a warning on the next line:
displayedNetworkAt!(net51, net51.hybrid[1]) # H2
# "WARNING: node -3 being the root is contradicted by isChild1 of its edges."
rootatnode!(net51, "A");
writeTopologyLevel1(net51) == "(A:0.5,((((C:0.52,(E:0.5)#H2:0.02):0.6,(#H2:0.01,F:0.7):0.8):0.9,D:1.1):1.3,B:2.3):0.5);" ||
  error("wrong net51 after displayedNetworkAt!");
net52 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
displayedNetworkAt!(net52, net52.hybrid[2])
writeTopologyLevel1(net52) == "(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,E:0.52):0.6,F:1.5):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7);" ||
  error("wrong net52 after displayedNetworkAt!");
taxa = tipLabels(net52); # order: A B C E F D
hardwiredClusters(net51, taxa) ==
[16 0 1 1 1 1 1 10;
 15 0 0 1 1 1 1 10;
 12 0 0 1 1 1 0 10;
  8 0 0 1 1 0 0 10;
  7 0 0 0 1 0 0 11;
 11 0 0 0 1 1 0 10] || error("wrong hardwired clusters for net51");
hardwiredClusters(net52, taxa) ==
[16 0 1 1 1 1 1 10;
  4 0 1 0 0 0 1 10;
  3 0 0 0 0 0 1 11;
 15 0 0 1 1 1 1 10;
 12 0 0 1 1 1 0 10;
  8 0 0 1 1 0 0 10] || error("wrong hardwired clusters for net52");
hardwiredClusterDistance(net51,net52,true) == 4 ||
 error("wrong HWDist between net51 and net52");
end

end # of testset: displayedTrees, hardwiredClusters, hardwiredClusterDistance, displayedNetworkAt!

@testset "testing hardwiredCluster! on single nodes" begin

net5 = "(A,((B,#H1),(((C,(E)#H2),(#H2,F)),(D)#H1)));" |> readTopology |> directEdges! ;
taxa = net5 |> tipLabels # ABC EF D
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
  @test hardwiredCluster(net5.edge[i], taxa) == m[:,i]
end
end # of testset, hardwiredCluster! on single nodes

end
