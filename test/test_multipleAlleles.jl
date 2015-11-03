# functions to test the functions for multiple alleles per species
# Claudia October 2015

include("../src/types.jl")
include("../src/functions.jl")

treefile = "(6,(5,(7,(3,4))));"
tree = readTopologyUpdate(treefile);
printEdges(tree)
printNodes(tree)
repSpecies=["7"]
expandLeaves!(repSpecies,tree)
writeTopology(tree)

df = readtable("CFtable1.csv")
alleleDF=DataFrame(allele=["1","2"],species=["7","7"])
newdf = mapAllelesCFtable!(alleleDF,df,false,"",tree)
writeTopology(tree)
