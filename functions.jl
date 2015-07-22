# functions written in classes.jl and moved here after tested
# for pseudolikelihood implementation (Stage2)
# Claudia August 2014
#
# in julia: include("functions.jl")

# tests of functions in examples_classes.jl outside git_laptop

# needed modules:
#using DataStructures # for updateInCycle with queue
using Base.Collections # for updateInCycle with priority queue
using DataFrames # for rep function and read/write csv tables
using NLopt # for branch lengths optimization
using GraphViz #for visualization

const DEBUG = true
const DEBUGC = false #more detailed prints
const CHECKNET = true #for debugging only

include("src/auxiliary.jl")

include("src/update.jl")

include("src/undo.jl")

include("src/addHybrid.jl")

include("src/deleteHybrid.jl")

include("src/moves.jl")

include("src/readwrite.jl")

include("src/readData.jl")

include("src/optimization.jl")

include("src/pseudolik.jl")

#include("src/visualization.jl")