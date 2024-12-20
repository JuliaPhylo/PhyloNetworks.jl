if !(@isdefined doalltests) doalltests = false; end

@testset "manipulateNet" begin

@testset "test: auxiliary" begin
global net
net = readnewick("((((B:102.3456789)#H1)#H2,((D:0.00123456789,C,#H2:::0.123456789)S1,(#H1,A_coolname)S2)S3)S4);")
s = IOBuffer()
@test_logs printedges(s, net)
@test String(take!(s)) == """
edge parent child  length  hybrid ismajor gamma   containroot i_cycle
1    2      1      102.346 false  true    1       false       -1     
2    3      2              true   true            false       -1     
3    10     3              true   true    0.8765  true        -1     
4    6      4      0.001   false  true    1       true        -1     
5    6      5              false  true    1       true        -1     
6    6      3              true   false   0.1235  true        -1     
7    9      6              false  true    1       true        -1     
8    8      2              true   false           true        -1     
9    8      7              false  true    1       true        -1     
10   9      8              false  true    1       true        -1     
11   10     9              false  true    1       true        -1     
"""
close(s); s = IOBuffer()
@test_logs printnodes(s, net)
@test String(take!(s)) == """
node leaf  hybrid name       i_cycle edges'numbers
1    true  false  B          -1      1   
2    false true   H1         -1      1    2    8   
3    false true   H2         -1      2    3    6   
4    true  false  D          -1      4   
5    true  false  C          -1      5   
6    false false  S1         -1      4    5    6    7   
7    true  false  A_coolname -1      9   
8    false false  S2         -1      8    9    10  
9    false false  S3         -1      7    10   11  
10   false false  S4         -1      3    11  
"""
close(s);
originalstdout = stdout
redirect_stdout(devnull)
printedges(net) # method without io argument
printnodes(net)
redirect_stdout(originalstdout)

@test_throws ErrorException getparent(net.node[net.rooti])
@test_throws ErrorException getparentedge(net.node[net.rooti])
@test_throws ErrorException getparentedgeminor(net.node[net.rooti])
@test_throws ErrorException getparentminor(net.node[1])
@test getparent(net.node[1]).number == 2
@test getparentminor(net.node[3]).number == 6
@test getparent(net.node[3]).number == 10
@test getparentedge(net.node[6]).number == 7
@test_throws ErrorException getparentedgeminor(net.node[6])
@test getparentedge(net.node[2]).number == 2
@test getparentedgeminor(net.node[2]).number == 8
@test [n.number for n in getchildren(net.node[4])] == [] # leaf
@test [n.number for n in getchildren(net.node[2])] == [1] # hybrid node
@test [n.number for n in getchildren(net.node[9])] == [6,8] # tree node
@test [n.number for n in getchildren(net.node[10])] == [3,9] # at root
@test [n.number for n in getchildren(net.node[6])] == [4,5,3] # polytomy
@test getparent(net.edge[8]).number == 8
@test [n.number for n in getparents(net.node[3])] == [10, 6]
@test [n.number for n in getparents(net.node[6])] == [9]
@test_throws ErrorException deleteleaf!(net, net.node[9])
n = deepcopy(net)
@test_logs deleteleaf!(n, n.node[7])
@test n.numnodes == 8; @test n.numedges == 9;
@test_logs deleteleaf!(net, net.node[7], simplify=false)
deleteleaf!(net, 4, simplify=false); deleteleaf!(net, 5, simplify=false)
@test net.numnodes == 5; @test net.numedges == 6;

## test nodeheights functions
net = readnewick("(((C:1,(A:1)#H1:1.5::0.7):1,(#H1:0.3::0.3,E:2.0):2.2):1.0,O:5.2);")
nh0 = [0.0,5.2,1.0,3.2,5.2,2.0,3.5,4.5,3.0]
@test getnodeheights(net) == nh0
@test getnodeheights_majortree(net) == nh0
net.edge[5].length = -1 # remove edge length: make it missing
nh = (@test_logs (:warn,"some hybrid edge length is missing") getnodeheights(net))
@test nh == nh0
@test net.edge[5].length == -1 # Make sure we don't mutate the broken edge length
@test getnodeheights!(net, false) == nh0 # no warning, though not tested
(x->x.length=-1).(view(net.edge, [3,5])) # make both hybrid edges missing
nh0 = [0,5.2,1,3.2,5.2,2,3.2,4.2,3]
@test (@test_logs (:warn, r"missing$") getnodeheights_average(net, false)) == nh0
getnodeheights!(net) == nh0
@test net.edge[3].length ≈ 1.2
@test net.edge[5].length ≈ 0.0 atol=1e-12
@test_logs getnodeheights_average(net, false) == nh0 # no more warning
@test istimeconsistent(net, false)
net.edge[3].length = -1 # make major edge missing
net.edge[5].length = 0.5
nh_major = [0.0,5.2,1.0,3.2,5.2,2.0,3.7,4.7,3.0]
@test (@test_logs (:warn, r"major hybrid edge missing a length") match_mode=:any getnodeheights_majortree(net))==nh_major
net.edge[3].length = 1.2 # make major known again for testing inconsistency
net.edge[5].length = 7 # time-*in*consistent
@test_throws "not time consistent" getnodeheights!(net, false)
@test !istimeconsistent(net, false)
nh = (@test_logs (:warn, r"not time consistent$") getnodeheights_average(net, false))
@test all(nh .≈ [0,5.2,1,3.2,5.2,2, 5.3, 6.3,3])
nh = (@test_logs getnodeheights_majortree(net, false, warn=false))
@test all(nh .≈ [0,5.2,1,3.2,5.2,2, 3.2, 4.2,3])
# re-assign lengths to switch candidate heights between major and minor
net.edge[3].length = 8.2
net.edge[5].length = 0.0
nh = (@test_logs getnodeheights_majortree(net, false, warn=false))
@test all(nh .≈ [0,5.2,1,3.2,5.2,2,10.2,11.2,3])

# below: 3 taxa, h=2, hybrid ladder but no 2-cycle. pruning t9 removes both hybrids.
nwkstring = "((t7:0.23,#H19:0.29::0.47):0.15,(((#H23:0.02::0.34)#H19:0.15::0.53,(t9:0.06)#H23:0.02::0.66):0.09,t6:0.17):0.21);"
net = readnewick(nwkstring)
@test_throws Exception deleteleaf!(net, "t1") # no leaf named t1
net.node[1].name = "t9"
@test_throws Exception deleteleaf!(net, "t9") # 2+ leaves named t1
@test_throws Exception deleteleaf!(net, 10, index=true) # <10 nodes
net.node[1].name = "t7" # back to original name
deleteleaf!(net, "t9", nofuse=true)
@test writenewick(net) == "((t7:0.23):0.15,(t6:0.17):0.21);"
deleteleaf!(net, "t7")
@test writenewick(net) == "(t6:0.17);"
@test_logs (:warn, r"^Only 1 node") deleteleaf!(net, "t6")
@test isempty(net.edge)
@test isempty(net.node)
end

@testset "testing directedges! and re-rootings" begin
# on a tree, then on a network with h=2"

if doalltests
tre = readnewick("(((((((1,2),3),4),5),(6,7)),(8,9)),10);");
tre.edge[1].ischild1=false; tre.edge[17].ischild1=false
directedges!(tre)
@test tre.edge[1].ischild1
@test tre.edge[17].ischild1
for i=1:18 tre.edge[i].containroot=false; end;
tre.rooti = 9;
directedges!(tre);
@test !tre.edge[9].ischild1
for i=1:18
 tre.edge[i].containroot || error("directedges! didn't correct containroot of $(i)th edge.")
end
tre = readnewick("(((((((1,2),3),4),5),(6,7)),(8,9)),10);");
rootatnode!(tre, -9); ## clau: previously -8
end

global net
net = readnewick("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);");
# 5th node = node number -7 (clau: previously -6).
net.rooti = 5
@test_logs directedges!(net)
@test !net.edge[12].ischild1
@test !net.edge[23].ischild1
@test [!net.edge[i].containroot for i in [8;collect(13:17)]] == [true for i in 1:6]
# or: "directedges! didn't correct containroot below a hyb node, $(i)th edge."
@test [net.edge[i].containroot for i in [9,5,18,2]] == [true for i in 1:4]
# or: "directedges! didn't correct containroot of hyb edges."
@test_logs rootatnode!(net, -10); # or: rootatnode! complained, node -10
@test_throws PhyloNetworks.RootMismatch rootatnode!(net, "M");
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
# net with 2 stacked hybrids
net = readnewick("((((((((1,2),3),4),(5)#H1),(#H1,(6,7))))#H2,(8,9)),(#H2,10));");
# sum([!e.containroot for e in net.edge]) # only 4.
directedges!(net); # or error("directedges! says that the root position is incompatible with hybrids")
@test sum([!e.containroot for e in net.edge]) == 16
plot(net, showedgenumber=true, showedgelength=false, shownodenumber=true);
net = readnewick("((((((((1,2),3),4),(5)#H1),(#H1,(6,7))))#H2,(8,9)),(#H2,10));");
net.rooti=19; # node number -13 (clau: previously -12)
directedges!(net)
end

net = readnewick("((((((((1,2),3),4),(5)#H1),(#H1,(6,7))))#H2,(8,9)),(#H2,10));");
net.rooti=15; # node number -5 (clau: previously -4)
@test_throws PhyloNetworks.RootMismatch directedges!(net);
# occursin(r"non-leaf node 9 had 0 children",e.msg))
@test_logs rootatnode!(net, -13); # or: rootatnode complained...
@test_throws PhyloNetworks.RootMismatch rootatnode!(net, -5);
# occursin(r"non-leaf node 9 had 0 children", e.msg))
@test_throws PhyloNetworks.RootMismatch rootatnode!(net,"H2"); #try rethrow();
# occursin(r"hybrid edge 17 conflicts", e.msg))
# earlier: """node 12 is a leaf. Will create a new node if needed, to set taxon "10" as outgroup."""
@test_logs rootatnode!(net,"10");

end # of testset for directedges! and re-rootings

@testset "testing preorder!" begin
# on a tree, then on a network with h=2
global net
tre = readnewick("(((((((1,2),3),4),5),(6,7)),(8,9)),10);");
net = readnewick("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);");

if doalltests
preorder!(tre)
## clau: previously [-1,10,-2,-9,9,8,-3,-8,7,6,-4,5,-5,4,-6,3,-7,2,1];
nodeN = [-2,10,-3,-10,9,8,-4,-9,7,6,-5,5,-6,4,-7,3,-8,2,1];
for i=1:length(tre.node)
  tre.vec_node[i].number==nodeN[i] ||
    error("node pre-ordered $i is node number $(tre.vec_node[i].number) instead of $(nodeN[i])")
end
end

@test_logs preorder!(net)
## clau previously: [-1,14,-14,13,12,-2,-9,11,-10,10,-3,-4,-5,-6,-7,5,6,4,3,2,-12,9,-13,8,7,1];
nodeN = [-2,14,-15,13,12,-3,-10,11,-11,10,-4,-5,-6,-7,-8,5,6,4,3,2,-13,9,-14,8,7,1];
@test [n.number for n in net.vec_node] == nodeN

cui3str = "(Xmayae,((Xhellerii,(((Xclemenciae_F2,Xmonticolus):1.458,(((((Xmontezumae,(Xnezahuacoyotl)#H26:0.247::0.804):0.375,((Xbirchmanni_GARC,Xmalinche_CHIC2):0.997,Xcortezi):0.455):0.63,(#H26:0.0::0.196,((Xcontinens,Xpygmaeus):1.932,(Xnigrensis,Xmultilineatus):1.401):0.042):2.439):2.0)#H7:0.787::0.835,(Xmaculatus,(Xandersi,(Xmilleri,((Xxiphidium,#H7:9.563::0.165):1.409,(Xevelynae,(Xvariatus,(Xcouchianus,(Xgordoni,Xmeyeri):0.263):3.532):0.642):0.411):0.295):0.468):0.654):1.022):0.788):1.917)#H27:0.149::0.572):0.668,Xalvarezi):0.257,(Xsignum,#H27:1.381::0.428):4.669);"
net3  = readnewick(cui3str);
deleteleaf!(net3,"Xhellerii"); deleteleaf!(net3,"Xsignum");
deleteleaf!(net3,"Xmayae", simplify=false, unroot=true);
# now: net3 has a 2-cycle
directedges!(net3)
preorder!(net3)
@test [n.number for n in net3.vec_node] == [-3,25,-6,-8,-20,-21,-22,-23,-25,-26,-27,-28,24,23,22,21,20,-24,15,-10,-16,-17,-19,14,13,-18,12,11,-11,-14,10,-15,9,8,-12,7,6,5,19,18,17,16,-7,4,3,26]

end # of testset for preorder!

@testset "testing cladewiseorder!" begin # on a tree then network with h=2
# cladewiseorder! is used for plotting: to avoid crossing edges in main tree

if doalltests
cladewiseorder!(tre)
nodeN = collect(19:-1:1);
for i=1:length(tre.node)
  tre.vec_int1[i]==nodeN[i] ||
    error("node clade-wise ordered $i is $(tre.vec_int1[i])th node instead of $(nodeN[i])th")
end
end

@test_logs cladewiseorder!(net)
nodeN = collect(26:-1:1);
@test net.vec_int1 == nodeN
end # of testset for cladewiseorder!

@testset "testing rotate!" begin
# to change the order of children edges at a given node

if doalltests
net = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
rotate!(net, -5) ## clau: previously -4
[e.number for e in net.node[13].edge] == [14,12,15] || error("rotate didn't work at node -5"); ## clau: previously -4
plot(net); # just to check no error.
end

global net
net=readnewick("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.136):10.0);");
@test_logs rotate!(net, -2, orderedEdgeNum=[1,12,9])
@test [e.number for e in net.node[12].edge] == [1,12,9] # or: rotate didn't work at node -2

end # of testset for rotate

@testset "other in manipulateNet" begin

net0 = readnewick("((((C:0.9)I1:0.1)I3:0.1,((A:1.0)I2:0.4)I3:0.6):1.4,(((B:0.2)H1:0.6)I2:0.5)I3:2.1);");
removedegree2nodes!(net0, true) # true: to keep the root of degree-2
@test writenewick(net0, round=true) == "((C:1.1,A:2.0):1.4,B:3.4);"

#--- delete above least stable ancestor ---#
# 3 blobs above LSA: cut edge + level-2 blob + cut edge.
net = readnewick("(((((#H25)#H22:::0.8,#H22),((t2:0.1,t1))#H25:::0.7)));")
PhyloNetworks.deleteaboveLSA!(net)
@test writenewick(net) == "(t2:0.1,t1);"
# 2 separate non-trivial clades + root edge
net = readnewick("((((t1,t2):0.1,(t3:0.3,t4)):0.4));")
PhyloNetworks.deleteaboveLSA!(net)
@test writenewick(net) == "((t1,t2):0.1,(t3:0.3,t4));"
# 1 tip + 2-cycle above it, no extra root edge above, hybrid LSA
net = readnewick("((t1)#H22:::0.8,#H22);")
PhyloNetworks.deleteaboveLSA!(net)
@test writenewick(net) == "(t1)H22;"
@test isempty(net.hybrid)
@test [n.name for n in net.vec_node] == ["H22","t1"]

# deleteleaf! and removedegree2nodes! when the root starts a 2-cycle
net = readnewick("((a,(((b)#H1,#H1))#H2),(#H2));")
deleteleaf!(net, "a")
@test writenewick(net) == "((#H2),(((b)#H1,#H1))#H2);"
removedegree2nodes!(net)
@test writenewick(net) == "((((b)#H1,#H1))#H2,#H2);"
net = readnewick("((a,(((b)#H1,#H1))#H2),#H2);") # no degree-2 node adjacent to root this time
deleteleaf!(net, "a", unroot=true)
@test writenewick(net) == "(b)H1;"

# unzip_canonical! and rezip_canonical!
netstr = "(#H2:0.1::0.2,((C:0.2,((B:0.3)#H1:0.4)#H2:0.5::0.8):0.6,(#H1:0.7,((A1:0.8)#H3:0.01,(A2:0.9,#H3:0.02):0.03):1.0):1.1):1.2,O:1.3);"
net = readnewick(netstr)
# 2 reticulation in a hybrid ladder, and another isolated reticulation
undoinfo = PhyloNetworks.unzip_canonical!(net)
@test all(getchildedge(h).length == 0.0 for h in net.hybrid) # unzipped
@test writenewick(net, round=true) == "(#H2:0.8::0.2,((C:0.2,((B:0.0)#H1:0.0)#H2:1.2::0.8):0.6,(#H1:1.0,((A1:0.0)#H3:0.81,(A2:0.9,#H3:0.82):0.03):1.0):1.1):1.2,O:1.3);"
PhyloNetworks.rezip_canonical!(undoinfo...)
@test writenewick(net, round=true) == netstr

end # of testset for other functions in manipulateNet

end # of overall testset for this file
