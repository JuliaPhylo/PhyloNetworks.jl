R"library(ape)";
# checking for absence of errors, for the most part

#info("Testing `apeRExport` and `sexp` functions")

@testset "Testing `apeRExport` and `sexp` functions" begin

# on a tree with some edge lengths missing
s = "(A,(B:1.0,(C:1.0,D:1.0):1.0):1.0);"
net = readTopology("(A,(B:1.0,(C:1.0,D:1.0):1.0):1.0);")
tree1 = apeRExport(net)
# check for correct topology:
R"tree2 = read.tree(text = $s)"
@test convert(Bool, R"dist.topo($tree1, tree2) == 0")

s = "(A,(B:1.0,(C:1.0,D:1.0):1.0):1.0);"
net = readTopology("(A,(B:1.0,(C:1.0,D:1.0):1.0):1.0);")
R"tree2 = read.tree(text = $s)"
@test convert(Bool, R"dist.topo($net, tree2) == 0")

end

# network, h=1, some missing gamma values
#s = "(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0):1.0):1.0,D:1.0);"
#net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0):1.0):1.0,D:1.0);") 
#tree1 = apeRExport(net)
#R"tree2 = read.evonet(text = $s)"
#@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")

#s = "(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0):1.0):1.0,D:1.0);"
#net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0):1.0):1.0,D:1.0);") 
#tree1 = sexp(net)
#R"tree2 = read.evonet(text = $s)"
#@test convert(Bool, R"dist.topo($net, tree2) == 0")
#@test convert(Bool, R"isTRUE(all.equal(tree2,$net))")

# network, h=1, mainTree=true; minor hybrid edge length missing
#s = "(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:::0.1):1.0):1.0,D:1.0);"
#net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:::0.1):1.0):1.0,D:1.0);") 
#tree1 = apeRExport(net; mainTree=true)
#R"tree2 = read.evonet(text = $s)"
#@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")

#s = "(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:::0.1):1.0):1.0,D:1.0);"
#net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:::0.1):1.0):1.0,D:1.0);") 
#tree1 = sexp(net; mainTree=true)
#R"tree2 = read.evonet(text = $s)"
#@test convert(Bool, R"dist.topo($net, tree2) == 0")

#on a network, h=1, with useEdgeLength=false 
#net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);")
#tree1 = apeRExport(net; useEdgeLength=false)
#R"tree2 = read.evonet(text = $s)"
#@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")

#net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);")
#tree1 = sexp(net; useEdgeLength=false)
#R"tree2 = read.evonet(text = $s)"
#@test convert(Bool, R"dist.topo($net, tree2) == 0")

#on a network with two hybridizations
#net = readTopology("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);")
#s = "(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);"
#tree1 = apeRExport(net)
#R"tree2 = read.evonet(text = $s)"
#@test convert(Bool, R"isTRUE(all.equal(tree2,$tree1))")

#net = readTopology("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);")
#s = "(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);"
#tree1 = sexp(net)
#R"tree2 = read.evonet(text = $s)"
#@test convert(Bool, R"dist.topo($net, tree2) == 0")
