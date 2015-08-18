module PhyloNetworks

using Base.Collections # for updateInCycle with priority queue
using DataFrames # for rep function and read/write csv tables
using NLopt # for branch lengths optimization
using GraphViz #for visualization

const DEBUG = true
const DEBUGC = false #more detailed prints
const CHECKNET = true #for debugging only

# export part

include("types.jl")
include("functions.jl")

end #module
