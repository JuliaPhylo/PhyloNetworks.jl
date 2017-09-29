@testset "Testing `apeRExport` function" begin

R"library(ape)"

#on a tree
net = readTopology("(A,(B,(C,D)));")
phy = apeRExport(net)
@test R"plot.evonet($phy)"


    
#on a network
net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
phy = apeRExport(net)
@test R"plot.evonet($phy)"

end

@testset "Testing `sexp` function on `HybridNetwork` objects" begin

R"library(ape)"

#on a tree
net = readTopology("(A,(B,(C,D)));")
phy = sexp(net)
@test R"plot.evonet($phy)"

    
#on a network
net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
phy = sexp(net)
@test R"plot.evonet($phy)"
end
