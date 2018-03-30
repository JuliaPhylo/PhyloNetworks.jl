# Public Documentation

Documentation for `PhyloNetworks`'s public (exported) interface.

See [Internal Documentation](@ref) for documentation on internal functions.

```@meta
DocTestSetup = quote
    using PhyloNetworks
end
```

```@contents
Pages = ["public.md"]
```

## Index

```@index
Pages = ["public.md"]
```

## types

```@autodocs
Modules = [PhyloNetworks]
Private = false
Order   = [:type]
```

## utilities

```@docs
tipLabels
sorttaxa!
printEdges
printNodes
summarizeDataCF
directEdges!
preorder!
cladewiseorder!
rootatnode!
rootonedge!
hybridatnode!
setLength!
setGamma!
deleteleaf!
deleteHybridThreshold!
rotate!
getindex(::TraitSimulation, ::Symbol)
getNodeAges
pairwiseTaxonDistanceMatrix
```

## data and topology read/write

```@docs
readTopology
readTopologyLevel1
readInputTrees
readMultiTopology
readSnaqNetwork
readTrees2CF
readTableCF
readTableCF!
readBootstrapTrees
writeTopology
writeMultiTopology
mapAllelesCFtable
```

## network inference

```@docs
snaq!
snaqDebug
topologyMaxQPseudolik!
topologyQPseudolik!
fittedQuartetCF
bootsnaq
calibrateFromPairwiseDistances!
```
## network Comparisons

```@docs
majorTree
minorTreeAt
displayedTrees
displayedNetworkAt!
hardwiredClusters
hardwiredCluster
hardwiredClusterDistance
treeEdgesBootstrap
# hybridDetection
# summarizeHFdf
hybridBootstrapSupport
```

## continuous trait evolution

```@docs
simulate
shiftHybrid
getShiftEdgeNumber
getShiftValue
phyloNetworklm
sigma2_estim
mu_estim
lambda_estim
ancestralStateReconstruction
expectations
predint
expectationsPlot
predintPlot
descendenceMatrix
regressorShift
regressorHybrid
```

## discrete trait evolution

```@docs
parsimonyDiscrete
nStates
Q
P
randomTrait
```

```@meta
DocTestSetup = nothing
```
