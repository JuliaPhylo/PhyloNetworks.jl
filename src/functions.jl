# functions written in classes.jl and moved here after tested
# for pseudolikelihood implementation (Stage2)
# Claudia August 2014
#
# in julia: include("functions.jl")

# for development and tests only. 

# needed modules:

using Base.Collections # for updateInCycle with priority queue
using DataFrames # for rep function and read/write csv tables
using GLM # for the lm function
using NLopt # for branch lengths optimization
#using GraphViz #for visualization
using Gadfly # for plotsGadfly, for visualization

import Base.show
import Gadfly.plot

const DEBUG = false
const DEBUGC = false #more detailed prints
const CHECKNET = false #for debugging only
const REDIRECT = false # for debugging to a file later


include("auxiliary.jl")
include("update.jl")
include("undo.jl")
include("addHybrid.jl")
include("deleteHybrid.jl")
include("moves.jl")
include("readwrite.jl")
include("readData.jl")
include("optimization.jl")
include("pseudolik.jl")
include("descriptive.jl")
include("manipulateNet.jl")
include("bootstrap.jl")
include("multipleAlleles.jl")
include("plotsGraphViz.jl")
include("plotsGadfly.jl")
include("compareNetworks.jl")
include("traits.jl")
