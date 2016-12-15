# Public Documentation

Documentation for `PhyloNetworks`'s public (exported) interface.

See [Internal Documentation](@ref) for documentation on internal functions.

```@meta
DocTestSetup = quote
    using PhyloNetworks
end
```

## Contents

```@contents
Pages = ["public.md"]
```

## Index

```@index
Pages = ["public.md"]
```

## Types

```@docs
HybridNetwork
DataCF
Quartet
```

## Utilities

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
# plotNetGraphViz
plot
rotate!
```

## Data and Topology read/write

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

## Network inference

```@docs
snaq!
snaqDebug
topologyMaxQPseudolik!
topologyQPseudolik!
fittedQuartetCF
bootsnaq
```
## Network Comparisons

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

## Trait Evolution

```@docs
ParamsBM
simulate
phyloNetworklm
sigma2_estim
mu_estim
lambda_estim
ancestralStateReconstruction
expectations
predint
predintPlot
expectationsPlot
```

```@meta
DocTestSetup = nothing
```
