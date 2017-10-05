R"library(ape)"

@testset "Testing `apeRExport` function" begin
    #checking for absence of errors, not testing for correctness

#on a tree
net = readTopology("(A,(B,(C,D)));")
phy = apeRExport(net)
R"plot.evonet($phy)"    
            
net = readTopology("(A:1.0,(B:1.0,(C:1.0,D:1.0):1.0):1.0);")
phy = apeRExport(net; useEdgeLength=true)
R"plot.evonet($phy)"
        
#on a network
net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
phy = apeRExport(net)
R"plot.evonet($phy)"    
        
net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);")
phy = apeRExport(net; useEdgeLength=true)
R"plot.evonet($phy)"

net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);")
phy = apeRExport(net; mainTree=true)
R"plot.evonet($phy)"

net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);")
phy = apeRExport(net; mainTree=true, useEdgeLength=true)
R"plot.evonet($phy)"

end

@testset "Testing `sexp` function on `HybridNetwork` objects" begin
    #checking for absence of errors, not testing for correctness

#on a tree
net = readTopology("(A,(B,(C,D)));")
phy = sexp(net)
R"plot.evonet($phy)"
        
net = readTopology("(A:1.0,(B:1.0,(C:1.0,D:1.0):1.0):1.0);")
phy = sexp(net)
R"plot.evonet($phy)"

#on a network
net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
phy = sexp(net)
R"plot.evonet($phy)"
        
net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);")
phy = sexp(net)
R"plot.evonet($phy)"
        
end
