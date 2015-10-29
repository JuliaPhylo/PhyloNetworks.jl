# test for tree2Matrix
# Claudia October 2015

include("../src/types.jl")
include("../src/functions.jl")

T = readTopology("(A,(B,(C,D)));");
S = ["A","B","C","D"]
M = tree2Matrix(T,S)
