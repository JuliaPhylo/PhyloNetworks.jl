# test all possibilities of readInputData
# Claudia April 2015

include("../types.jl")
include("../functions.jl")

d=readInputData("1.ms"); #tableCF0.txt
d=readInputData("1.ms",:rand,10); #tableCF3.txt
d=readInputData("1.ms",[1,2,3,4,5]); #tableCF4.txt

d=readInputData("1.ms","allQuartets.txt"); #tableCF1.txt
d=readInputData("1.ms","allQuartets.txt",:rand,10); #tableCF2.txt
d=readInputData("1.ms","allQuartets.txt",true,"try4.txt"); #try4.txt

descData(d)

