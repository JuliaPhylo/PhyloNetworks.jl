# functions to test the functions for multiple alleles per species
# Claudia October 2015

include("../src/types.jl")
include("../src/functions.jl")

df = readtable("CFtable1.csv")

alleleDF=DataFrame(allele=["1","2"],species=["7","7"])
newdf = mapAllelesCFtable(alleleDF,df,false,"")

