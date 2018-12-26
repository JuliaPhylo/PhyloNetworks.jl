# changed node numbering on 5/28 when readSubTree will start internal
# node numbers at -2 to avoid error in undirectedOtherNetworks

if !(@isdefined doalltests) doalltests = false; end

@testset "test: auxiliary" begin
global net
net = readTopology("((((B:102.3456789)#H1)#H2,((D:0.00123456789,C,#H2:::0.123456789)S1,(#H1,A_coolname)S2)S3)S4);")
# using PhyloPlots; plot(net, :R, showEdgeNumber=true, showNodeNumber=true);
s = IOBuffer()
@test_logs printEdges(s, net)
@test String(take!(s)) == """
edge parent child  length  hybrid isMajor gamma   containRoot inCycle istIdentitiable
1    2      1      102.346 false  true    1       true        -1      false
2    3      2              true   true            false       -1      true 
3    10     3              true   true    0.8765  false       -1      true 
4    6      4      0.001   false  true    1       true        -1      false
5    6      5              false  true    1       true        -1      false
6    6      3              true   false   0.1235  false       -1      true 
7    9      6              false  true    1       true        -1      true 
8    8      2              true   false           false       -1      true 
9    8      7              false  true    1       true        -1      false
10   9      8              false  true    1       true        -1      true 
11   10     9              false  true    1       true        -1      true 
"""
close(s); s = IOBuffer()
@test_logs printNodes(s, net)
@test String(take!(s)) == """
node leaf  hybrid hasHybEdge name       inCycle edges'numbers
1    true  false  false      B          -1      1   
2    false true   true       #H1        -1      1    2    8   
3    false true   true       #H2        -1      2    3    6   
4    true  false  false      D          -1      4   
5    true  false  false      C          -1      5   
6    false false  true       S1         -1      4    5    6    7   
7    true  false  false      A_coolname -1      9   
8    false false  true       S2         -1      8    9    10  
9    false false  false      S3         -1      7    10   11  
10   false false  true       S4         -1      3    11  
"""
@test_throws ErrorException PhyloNetworks.getMinorParent(net.node[1])
@test PhyloNetworks.getMajorParent(net.node[1]).number == 2
@test PhyloNetworks.getMinorParent(net.node[3]).number == 6
@test PhyloNetworks.getMajorParent(net.node[3]).number == 10
@test PhyloNetworks.getMajorParentEdge(net.node[6]).number == 7
@test_throws ErrorException PhyloNetworks.getMinorParentEdge(net.node[6])
@test PhyloNetworks.getMajorParentEdge(net.node[2]).number == 2
@test PhyloNetworks.getMinorParentEdge(net.node[2]).number == 8
@test [n.number for n in PhyloNetworks.getChildren(net.node[4])] == [] # leaf
@test [n.number for n in PhyloNetworks.getChildren(net.node[2])] == [1] # hybrid node
@test [n.number for n in PhyloNetworks.getChildren(net.node[9])] == [6,8] # tree node
@test [n.number for n in PhyloNetworks.getChildren(net.node[10])] == [3,9] # at root
@test [n.number for n in PhyloNetworks.getChildren(net.node[6])] == [4,5,3] # polytomy
@test PhyloNetworks.getParent(net.edge[8]).number == 8
@test [n.number for n in PhyloNetworks.getParents(net.node[3])] == [10, 6]
@test [n.number for n in PhyloNetworks.getParents(net.node[6])] == [9]
@test_throws ErrorException deleteleaf!(net, net.node[9])
n = deepcopy(net)
@test_logs deleteleaf!(n, n.node[7])
@test n.numNodes == 8; @test n.numEdges == 9;
@test_logs deleteleaf!(net, net.node[7], simplify=false)
deleteleaf!(net, 4, simplify=false); deleteleaf!(net, 5, simplify=false)
@test net.numNodes == 5; @test net.numEdges == 6;
end

@testset "testing directEdges! and re-rootings" begin
# on a tree, then on a network with h=2"

if doalltests
tre = readTopology("(((((((1,2),3),4),5),(6,7)),(8,9)),10);");
tre.edge[1].isChild1=false; tre.edge[17].isChild1=false
directEdges!(tre)
tre.edge[1].isChild1  || error("directEdges! didn't correct the direction of 1st edge")
tre.edge[17].isChild1 || error("directEdges! didn't correct the direction of 17th edge")
for i=1:18 tre.edge[i].containRoot=false; end;
tre.root = 9;
directEdges!(tre);
!tre.edge[9].isChild1 || error("directEdges! didn't correct the direction of 9th edge")
for i=1:18
 tre.edge[i].containRoot || error("directEdges! didn't correct containRoot of $(i)th edge.")
end
tre = readTopology("(((((((1,2),3),4),5),(6,7)),(8,9)),10);");
rootatnode!(tre, -9); ## clau: previously -8
end

global net
net = readTopology("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);");
# 5th node = node number -7 (clau: previously -6).
net.root = 5
@test_logs directEdges!(net);
@test !net.edge[12].isChild1 # or: directEdges! didn't correct the direction of 12th edge
@test !net.edge[23].isChild1 # or: directEdges! didn't correct the direction of 23th edge"
@test [!net.edge[i].containRoot for i in [8;collect(13:17)]] == [true for i in 1:6]
# or: "directEdges! didn't correct containRoot below a hyb node, $(i)th edge."
@test [net.edge[i].containRoot for i in [9,5,18,2]] == [true for i in 1:4]
# or: "directEdges! didn't correct containRoot of hyb edges."
# plot(net, showNodeNumber=true, showEdgeLength=false, showEdgeNumber=true)
@test_logs rootatnode!(net, -10); # or: rootatnode! complained, node -10
@test_throws PhyloNetworks.RootMismatch rootatnode!(net, "M"; verbose=false);
# println("the rootmismatch about node 5 is good and expected.")
@test_logs rootonedge!(net, 9); # or: rootonedge! complained, edge 9
@test_logs PhyloNetworks.fuseedgesat!(27, net);
# earlier warning: """node 1 is a leaf. Will create a new node if needed, to set taxon "Ag" as outgroup."""
@test_logs rootatnode!(net, "Ag"); # need for new node
# earlier: """node 1 is a leaf. Will create a new node if needed, to set taxon "Ag" as outgroup."""
@test_logs rootatnode!(net, "Ag"); # no need this time
@test length(net.node) == 27 # or: wrong # of nodes after rootatnode! twice on same outgroup
# earlier: """node 10 is a leaf. Will create a new node if needed, to set taxon "Ap" as outgroup."""
@test_logs rootatnode!(net, "Ap");
@test length(net.node) == 27 # or: wrong # of nodes, after 3rd rooting with outgroup
@test_logs rootonedge!(net, 5);

# example with one hybridization below another

if doalltests
net = readTopology("((((((((1,2),3),4),(5)#H1),(#H1,(6,7))))#H2,(8,9)),(#H2,10));");
# sum([!e.containRoot for e in net.edge]) # only 4.
directEdges!(net); # or error("directEdges! says that the root position is incompatible with hybrids")
sum([!e.containRoot for e in net.edge]) == 16 ||
 error("directEdges! wrong on net with 2 stacked hybrids");
plot(net, showEdgeNumber=true, showEdgeLength=false, showNodeNumber=true);
net = readTopology("((((((((1,2),3),4),(5)#H1),(#H1,(6,7))))#H2,(8,9)),(#H2,10));");
net.root=19; # node number -13 (clau: previously -12)
directEdges!(net); # or error("directEdges! says that the root position is incompatible with hybrids");
end

net = readTopology("((((((((1,2),3),4),(5)#H1),(#H1,(6,7))))#H2,(8,9)),(#H2,10));");
net.root=15; # node number -5 (clau: previously -4)
@test_throws PhyloNetworks.RootMismatch directEdges!(net);
# occursin(r"non-leaf node 9 had 0 children",e.msg))
@test_logs rootatnode!(net, -13); # or: rootatnode complained...
@test_throws PhyloNetworks.RootMismatch rootatnode!(net, -5); # verbose = true this time
# occursin(r"non-leaf node 9 had 0 children", e.msg))
@test_throws PhyloNetworks.RootMismatch rootatnode!(net,"#H2"; verbose=false); #try rethrow();
# occursin(r"hybrid edge 17 conflicts", e.msg))
# earlier: """node 12 is a leaf. Will create a new node if needed, to set taxon "10" as outgroup."""
@test_logs rootatnode!(net,"10");

end # of testset for directEdges! and re-rootings

@testset "testing preorder!" begin
# on a tree, then on a network with h=2
global net
tre = readTopology("(((((((1,2),3),4),5),(6,7)),(8,9)),10);");
net = readTopology("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);");

if doalltests
preorder!(tre)
## clau: previously [-1,10,-2,-9,9,8,-3,-8,7,6,-4,5,-5,4,-6,3,-7,2,1];
nodeN = [-2,10,-3,-10,9,8,-4,-9,7,6,-5,5,-6,4,-7,3,-8,2,1];
for i=1:length(tre.node)
  tre.nodes_changed[i].number==nodeN[i] ||
    error("node pre-ordered $i is node number $(tre.nodes_changed[i].number) instead of $(nodeN[i])")
end
end

@test_logs preorder!(net)
## clau previously: [-1,14,-14,13,12,-2,-9,11,-10,10,-3,-4,-5,-6,-7,5,6,4,3,2,-12,9,-13,8,7,1];
nodeN = [-2,14,-15,13,12,-3,-10,11,-11,10,-4,-5,-6,-7,-8,5,6,4,3,2,-13,9,-14,8,7,1];
@test [n.number for n in net.nodes_changed] == nodeN

cui3str = "(Xmayae,((Xhellerii,(((Xclemenciae_F2,Xmonticolus):1.458,(((((Xmontezumae,(Xnezahuacoyotl)#H26:0.247::0.804):0.375,((Xbirchmanni_GARC,Xmalinche_CHIC2):0.997,Xcortezi):0.455):0.63,(#H26:0.0::0.196,((Xcontinens,Xpygmaeus):1.932,(Xnigrensis,Xmultilineatus):1.401):0.042):2.439):2.0)#H7:0.787::0.835,(Xmaculatus,(Xandersi,(Xmilleri,((Xxiphidium,#H7:9.563::0.165):1.409,(Xevelynae,(Xvariatus,(Xcouchianus,(Xgordoni,Xmeyeri):0.263):3.532):0.642):0.411):0.295):0.468):0.654):1.022):0.788):1.917)#H27:0.149::0.572):0.668,Xalvarezi):0.257,(Xsignum,#H27:1.381::0.428):4.669);"
net3  = readTopology(cui3str);
deleteleaf!(net3,"Xhellerii"); deleteleaf!(net3,"Xsignum");
deleteleaf!(net3,"Xmayae", simplify=false);
# now: net3 has a 2-cycle
directEdges!(net3)
preorder!(net3)
@test [n.number for n in net3.nodes_changed] == [-3,25,-6,-8,-20,-21,-22,-23,-25,-26,-27,-28,24,23,22,21,20,-24,15,-10,-16,-17,-19,14,13,-18,12,11,-11,-14,10,-15,9,8,-12,7,6,5,19,18,17,16,-7,4,3,26]

end # of testset for preorder!

@testset "testing cladewiseorder!" begin # on a tree then network with h=2
# cladewiseorder! is used for plotting: to avoid crossing edges in main tree

if doalltests
cladewiseorder!(tre)
nodeN = collect(19:-1:1);
for i=1:length(tre.node)
  tre.cladewiseorder_nodeIndex[i]==nodeN[i] ||
    error("node clade-wise ordered $i is $(tre.cladewiseorder_nodeIndex[i])th node instead of $(nodeN[i])th")
end
end

@test_logs cladewiseorder!(net)
nodeN = collect(26:-1:1);
@test net.cladewiseorder_nodeIndex == nodeN
end # of testset for cladewiseorder!

@testset "testing rotate!" begin
# to change the order of children edges at a given node

if doalltests
net = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
rotate!(net, -5) ## clau: previously -4
[e.number for e in net.node[13].edge] == [14,12,15] || error("rotate didn't work at node -5"); ## clau: previously -4
plot(net); # just to check no error.
end

global net
net=readTopology("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.136):10.0);");
@test_logs rotate!(net, -2, orderedEdgeNum=[1,12,9])
@test [e.number for e in net.node[12].edge] == [1,12,9] # or: rotate didn't work at node -2

end # of testset for rotate
