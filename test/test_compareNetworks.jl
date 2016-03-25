# test of deleteHybridEdge!, functions to extract displayed trees/subnetworks,
#      used to compare networks with the hardwired cluster distance.
# Cecile March 2016

using PhyloNetworks

#---- testing deleteHybridEdge! ---------------#
println("\n\nTesting deleteHybridEdge!")

# example of network with one hybrid edge connected to the root:
net = readTopology("((Adif:1.0,(Aech:0.122,#H6:10.0::0.047):10.0):1.614,Aten:1.0,((Asub:1.0,Agem:1.0):0.0)#H6:5.062::0.953);");
# plot(net, showEdgeNumber=true, showNodeNumber=true)
PhyloNetworks.deleteHybridEdge!(net, net.edge[10]);
writeTopology(net) == "(Adif:1.0,(Aech:0.122,(Asub:1.0,Agem:1.0):10.0):10.0,Aten:2.614);" ||
 error("deleteHybridEdge! didn't work on 10th edge")
net = readTopology("((Adif:1.0,(Aech:0.122,#H6:10.0::0.047):10.0):1.614,Aten:1.0,((Asub:1.0,Agem:1.0):0.0)#H6:5.062::0.953);");
PhyloNetworks.deleteHybridEdge!(net, net.edge[3]);
writeTopology(net) == "((Adif:1.0,Aech:10.122):1.614,Aten:1.0,(Asub:1.0,Agem:1.0):5.062);" ||
 error("deleteHybridEdge! didn't work on 3rd edge")
# plot(net, showEdgeNumber=true, showNodeNumber=true)

net=readTopology("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.1361111):10.0);");
# plot(net, showEdgeNumber=true, showNodeNumber=true)
PhyloNetworks.deleteHybridEdge!(net, net.edge[11]);
writeTopology(net) == "(4,((1,2):2.069,(6,5):3.423):0.265,3);" ||
 error("deleteHybridEdge! didn't work on 11th edge")
net=readTopology("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.1361111):10.0);");
PhyloNetworks.deleteHybridEdge!(net, net.edge[4]);
writeTopology(net) == "(4,((6,5):3.423,1):0.265,(3,2):10.0);" ||
 error("deleteHybridEdge! didn't work on 4th edge")

# example with wrong attributed inChild1
net=readTopology("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.1361111):10.0);");
net.edge[5].isChild1 = false;
PhyloNetworks.deleteHybridEdge!(net, net.edge[4]);
writeTopology(net) == "(4,((6,5):3.423,1):0.265,(3,2):10.0);" ||
 error("deleteHybridEdge! didn't work on 4th edge when isChild1 was outdated")

net = readTopology("((Adif:1.0,(Aech:0.122,#H6:10.0::0.047):10.0):1.614,Aten:1.0,((Asub:1.0,Agem:1.0):0.0)#H6:5.062::0.953);");
net.edge[5].isChild1 = false # edge 5 from -1 to -2
PhyloNetworks.deleteHybridEdge!(net, net.edge[10]);
# WARNING: node -1 being the root is contradicted by isChild1 of its edges.
writeTopology(net) == "(Adif:1.0,(Aech:0.122,(Asub:1.0,Agem:1.0):10.0):10.0,Aten:2.614);" ||
 error("deleteHybridEdge! didn't work on 10th edge after isChild1 was changed")

# plot(net, showEdgeNumber=true, showNodeNumber=true)

#----------------------------------------------------------#
#   testing functions to display trees / subnetworks       #
#----------------------------------------------------------#

println("\n\nTesting deleteHybridThreshold!")

net21 = readTopology("(A,((B,#H1),(C,(D)#H1)));");
# manual bug fix to get good gamma's (0.5 not 1.0) and major/minor
net21.edge[3].gamma = 0.5;
net21.edge[7].gamma = 0.5;
deleteHybridThreshold!(net21,0.2);
writeTopology(net21) == "(A,((B,#H1:::0.5),(C,(D)#H1:::0.5)));" ||
 error("deleteHybridThreshold! didn't work on net21, gamma=0.2")
deleteHybridThreshold!(net21,0.5);
writeTopology(net21) == "(A,((C,D),B));" ||
 error("deleteHybridThreshold! didn't work on net21, gamma=0.5")

net22 = readTopology("(A,((B,#H1:::0.2),(C,(D)#H1:::0.8)));");
deleteHybridThreshold!(net22,0.3);
writeTopology(net22) == "(A,((C,D),B));" ||
 error("deleteHybridThreshold! didn't work on net22, gamma=0.3")

net31 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(C:0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
net32 = deepcopy(net31); for (e in net32.edge) e.length=-1.0; end
# plot(net31)
deleteHybridThreshold!(net31,0.3);
writeTopology(net31) == "(A:1.0,((C:0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("deleteHybridThreshold! didn't work on net31, gamma=0.3")
deleteHybridThreshold!(net32,0.3)
writeTopology(net32) == "(A,((C,D),B));" ||
 error("deleteHybridThreshold! didn't work on net32, gamma=0.3")

net42 = readTopology("(A:1.0,((B:1.1,#H1:0.2):1.2,(C:0.9,(D:0.8)#H1:0.3):1.3):0.7):0.1;");
deleteHybridThreshold!(net42,0.5);
writeTopology(net42) == "(A:1.0,((C:0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("deleteHybridThreshold! didn't work on net42, gamma=0.5")

net5 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
# plot(net5)
deleteHybridThreshold!(net5,0.5);  # both H1 and H2 eliminated
writeTopology(net5) == "(A:1.0,((((C:0.52,E:0.52):0.6,F:1.5):0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("deleteHybridThreshold! didn't work on net5, gamma=0.5")
net5 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
deleteHybridThreshold!(net5,0.3);  # H2 remains
writeTopology(net5) == "(A:1.0,((((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("deleteHybridThreshold! didn't work on net5, gamma=0.3")

println("\n\nTesting displayedNetworks! and displayedTrees")

net3 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(C:0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
net31 = PhyloNetworks.displayedNetworks!(net3, net3.node[6]); #H1 = 6th node
writeTopology(net31) == "(A:1.0,((B:1.1,D:1.0):1.2,C:2.2):0.7);" ||
 error("displayedNetworks! didn't work on net3, minor at 6th node")
writeTopology(net3)  == "(A:1.0,((C:0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("displayedNetworks! didn't work on net3, major at 6th node")

net3 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(C:0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
a = displayedTrees(net3, 0.2);
length(a) == 2 ||
 error("displayedTrees didn't work on net3, gamma=0.2: output not of length 2")
writeTopology(a[1]) == "(A:1.0,((C:0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("displayedTrees didn't work on net3, gamma=0.2: 1st tree wrong")
writeTopology(a[2]) == "(A:1.0,((B:1.1,D:1.0):1.2,C:2.2):0.7);" ||
 error("displayedTrees didn't work on net3, gamma=0.2: 2nd tree wrong")

net5 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
a = displayedTrees(net5, 0.5);
length(a) == 1 ||
 error("displayedTrees didn't work on net5, gamma=0.5: output not of length 1")
writeTopology(a[1]) == "(A:1.0,((((C:0.52,E:0.52):0.6,F:1.5):0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("displayedTrees didn't work on net5, gamma=0.5: wrong tree")
a = displayedTrees(net5, 0.1);
length(a) == 4 ||
 error("displayedTrees didn't work on net5, gamma=0.1: output not of length 4")
(writeTopology(a[1]) == "(A:1.0,((((C:0.52,E:0.52):0.6,F:1.5):0.9,D:1.1):1.3,B:2.3):0.7);" &&
 writeTopology(a[2]) == "(A:1.0,((B:1.1,D:1.0):1.2,((C:0.52,E:0.52):0.6,F:1.5):2.2):0.7);" &&
 writeTopology(a[3]) == "(A:1.0,((((F:0.7,E:0.51):0.8,C:1.12):0.9,D:1.1):1.3,B:2.3):0.7);" &&
 writeTopology(a[4]) == "(A:1.0,((B:1.1,D:1.0):1.2,((F:0.7,E:0.51):0.8,C:1.12):2.2):0.7);") ||
 error("displayedTrees didn't work on net5, gamma=0.1: one or more tree(s) is/are wrong")

println("\n\nTesting majorTree and displayedNetworkAt!")

writeTopology(majorTree(net5)) == "(A:1.0,((((C:0.52,E:0.52):0.6,F:1.5):0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("majorTree didn't work on net5")

net5 = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
displayedNetworkAt!(net5, net5.hybrid[1]);
writeTopology(net5) == "(A:1.0,((((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,D:1.1):1.3,B:2.3):0.7);" ||
 error("displayedNetworkAt! didn't work on net5, 1st hybrid")
