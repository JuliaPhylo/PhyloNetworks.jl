VERSION >= v"0.4" || warn("you need to update Julia to current version: 0.4 or higher")
__precompile__()

module PhyloNetworks

using Base.Collections # for updateInCycle with priority queue
using DataFrames # for rep function and read/write csv tables
using NLopt # for branch lengths optimization
using GraphViz #for visualization

import Base.show

const DEBUG = false #for debugging only
const DEBUGC = false #more detailed prints
const CHECKNET = false #for debugging only
const REDIRECT = false # changed for debugging to a file

export
HybridNetwork,
DataCF,
Quartet,
readTopology,
readTopologyLevel1,
tipLabels,
writeTopology,
deleteLeaf!,
printEdges,
printNodes,
readTrees2CF,
readTableCF,
readTableCF!,
readInputTrees,
summarizeDataCF,
readStartTop,
snaq!,
readSnaqNetwork,
snaqDebug,
topologyMaxQPseudolik!,
topologyQPseudolik!,
root!,
dfObsExpCF,
plotPhylonet,
generalExport #for graphviz test

# export part

include("types.jl")
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
#include("visualization.jl")
include("descriptive.jl")
include("bootstrap.jl")
include("genExport.jl")
include("drawTraversal.jl")
include("plotPhylonet.jl")


end #module
