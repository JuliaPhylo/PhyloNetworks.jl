VERSION >= v"0.5" || warn("you need to update Julia to current version: 0.5 or higher")
__precompile__()

module PhyloNetworks

using DataStructures # for updateInCycle with priority queue
using DataFrames # for functions to read/write tables, and names()
using GLM # for the lm function
using NLopt # for branch lengths optimization
# using GraphViz #for visualization
using Gadfly # for plots
using ColorTypes # used by Gadfly already. To resolve data type names (Colorant)
using StatsBase: sample
using Combinatorics.combinations
using RCall

import Base.show
import Gadfly.plot

global DEBUG = false #for debugging only
const DEBUGC = false #more detailed prints
global CHECKNET = false #for debugging only
global REDIRECT = false # changed for debugging to a file

export
HybridNetwork,
DataCF,
Quartet,
readTopology,
readTopologyLevel1,
tipLabels,
writeTopology,
deleteleaf!,
printEdges,
printNodes,
sorttaxa!,
readTrees2CF,
readTableCF,
readTableCF!,
readInputTrees,
summarizeDataCF,
snaq!,
readSnaqNetwork,
snaqDebug,
topologyMaxQPseudolik!,
topologyQPseudolik!,
rootatnode!,
rootonedge!,
directEdges!,
preorder!,
cladewiseorder!,
fittedQuartetCF,
plotNetGraphViz,
plot,
rotate!,
setLength!,
setGamma!,
mapAllelesCFtable,
deleteHybridThreshold!,
displayedTrees,
majorTree,
minorTreeAt,
displayedNetworkAt!,
hardwiredClusters,
hardwiredCluster,
hardwiredCluster!,
hardwiredClusterDistance,
treeEdgesBootstrap,
hybridDetection,
summarizeHFdf,
hybridBootstrapSupport,
bootsnaq,
readBootstrapTrees,
writeMultiTopology,
readMultiTopologyLevel1,
readMultiTopology,
hybridatnode!,
undirectedOtherNetworks,
phyloNetworklm,
PhyloNetworkLinearModel,
simulate,
TraitSimulation,
ParamsBM,
sharedPathMatrix,
ancestralStateReconstruction,
ReconstructedStates,
sigma2_estim,
mu_estim,
lambda_estim,
expectations,
expectationsPlot,
predint,
predintPlot,
parsimonyDiscrete

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
include("plotsRCall.jl")
include("compareNetworks.jl")
include("traits.jl")
include("parsimony.jl")

end #module
