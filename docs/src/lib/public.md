# Public Documentation

Documentation for `PhyloNetworks`'s public (exported) interface.

See [Internal Documentation](@ref) for documentation on internal functions.

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
readBootstrapTrees
writeTopology
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