R"library(ape)"

info("Testing `apeRExport` function")
    #checking for absence of errors, not testing for correctness

#on a tree
net = readTopology("(A:1.0,(B:1.0,(C:1.0,D:1.0):1.0):1.0);")
phy = apeRExport(net; useEdgeLength=true)
R"plot.evonet($phy)"
        
#on a network
net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);")
phy = apeRExport(net; useEdgeLength=true)
R"plot.evonet($phy)"

net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);")
phy = apeRExport(net; mainTree=true)
R"plot.evonet($phy)"

net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);")
phy = apeRExport(net; mainTree=true, useEdgeLength=true)
R"plot.evonet($phy)"

info("Testing `apeRExport` function")
    #checking for absence of errors, not testing for correctness

#on a tree
net = readTopology("(A:1.0,(B:1.0,(C:1.0,D:1.0):1.0):1.0);")
phy = sexp(net)
R"plot.evonet($phy)"

#on a network
net = readTopology("(((A:1.0,(B:1.0)#H1:1.0::0.9):1.0,(C:1.0,#H1:1.0::0.1):1.0):1.0,D:1.0);")
phy = sexp(net)
R"plot.evonet($phy)"
