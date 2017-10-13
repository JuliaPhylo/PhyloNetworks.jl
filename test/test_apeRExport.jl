using RCall # could be moved to runtests.jl

# checking for absence of errors, for the most part
@testset "Testing `apeRExport` and `sexp` functions" begin
R"library(ape)";
R"cat(‘loaded ape library’)"

# on a tree with some edge lengths missing
s = "(A,(B:1.0,(C:1.0,D:1.0):1.0):1.0);";
net = readTopology(s);
tree1 = apeRExport(net);
# check for correct unrooted topology. with the 'score' method,
# we also check for identical lengths of (unrooted) internal edges.
R"tree2 = read.tree(text = $s)";
@test convert(Bool, R"dist.topo($tree1, tree2, method='score') == 0")
# testing sexp: through the use of $net directly, later
net = readTopology(s);
@test convert(Bool, R"dist.topo($net, tree2, method='score') == 0")

# network, h=1, some missing gamma values
#s = "(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0):1.0):1.0,D:1.0);";
#net = readTopology(s);
#tree1 = apeRExport(net);
## commented out tests using h>0 networks, until new version of ape on CRAN
#R"tree2 = read.evonet(text = $s)"
#@test convert(Bool, R"dist.topo($tree1, tree2, method='score') == 0")
#net = readTopology(s) # next: use sexp
#@test convert(Bool, R"dist.topo($net, tree2, method='score') == 0")

# network, h=1, mainTree=true; minor hybrid edge length missing
#s = "(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:::0.1):1.0):1.0,D:1.0);";
#s2 = "(((A:1.0,B:2.0):1.0,C:2.0):1.0,D:1.0);"; # main tree
net = readTopology(s);
net1 = apeRExport(net; mainTree=true);
R"tree2 = read.tree(text = $s2)";
@test convert(Bool, R"dist.topo($net1, tree2, method='score') == 0")
net = readTopology(s);
net = majorTree(net);
R"tree2 = read.tree(text = $s2)" ## fixit: uncomment when CRAN has evonet functions
@test convert(Bool, R"dist.topo($net, tree2, method='score') == 0")

# on a network, h=1, with useEdgeLength=false
#s = "(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);";
#net = readTopology(s);
#net1 = apeRExport(net);
#R"tree2 = read.evonet(text = $s)" ## fixit: uncomment when CRAN has evonet functions
#@test convert(Bool, R"dist.topo($net1, tree2, method='score') == 0")
#net = readTopology(s)
#@test convert(Bool, R"dist.topo($net, tree2, method='score') == 0")

# on a network with two hybridizations
#s = "(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,
#      (((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);";
#net = readTopology(s);
#net1 = apeRExport(net);
#R"tree2 = read.evonet(text = $s)" ## fixit: uncomment when CRAN has evonet functions
#@test convert(Bool, R"dist.topo($net1, tree2, method='score') == 0")
#net = readTopology(s)
#@test convert(Bool, R"dist.topo($net, tree2, method='score') == 0")
end
