VERSION >= v"0.4" || warn("you need to update Julia to current version: 0.4 or higher")
__precompile__()

module PhyloNetworks

using Base.Collections # for updateInCycle with priority queue
using DataFrames # for rep function and read/write csv tables
using GLM # for the lm function
using NLopt # for branch lengths optimization
# using GraphViz #for visualization
using Gadfly # for plots
using ColorTypes # used by Gadfly already. To resolve data type names (Colorant)

import Base.show
import Gadfly.plot

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
readInputTrees,
summarizeDataCF,
snaq!,
readSnaqNetwork,
snaqDebug,
topologyMaxQPseudolik!,
topologyQPseudolik!,
root!,
preorder!,
cladewiseorder!,
dfObsExpCF,
plotNetGraphViz,
plot,
setLength!,
setGamma!,
mapAllelesCFtable,
deleteHybridThreshold!,
displayedTrees,
majorTree,
minorTreeAt,
displayedNetworkAt!,
hardwiredClusters,
hardwiredClusterDistance,
treeEdgesBootstrap,
hybridDetection,
summarizeHFdf,
bootsnaq

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
include("descriptive.jl")
include("manipulateNet.jl")
include("bootstrap.jl")
include("multipleAlleles.jl")
include("plotsGraphViz.jl")
include("plotsGadfly.jl")
include("compareNetworks.jl")
include("traits.jl")


end #module
