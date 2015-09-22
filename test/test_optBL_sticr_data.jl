# test of optBL with the example in STICR panmixia test
# Claudia January 2015

include("types.jl")
include("functions.jl")

df = readtable("/Users/Clauberry/Documents/phylo/software/CFimplementation/sticr/quartetCF.csv")
d2 = readTableCF(df)

net = readTopologyUpdate("/Users/Clauberry/Documents/phylo/software/CFimplementation/sticr/tree.tre")
# fixit: cannot read _ inside the taxon name, need to change readsubtree
printEdges(net)


net.ht
realht = [0.1,0.2,0.1,1.0]

@time fmin,xmin=optBL(net,d2)


