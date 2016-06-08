# `PhyloNetworks`

PhyloNetworks is a Julia package for the manipulation, visualization
and inference of phylogenetic networks.  SNaQ implements the
statistical inference method in
[Sol&iacute;s-Lemus and An&eacute; 2016](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896).
The procedure involves a numerical optimization of branch lengths and inheritance
probabilities and a heuristic search in the space of phylogenetic networks.

# PhyloNetworks.jl Documentation

```@contents
Pages = ["simpleJulia.md", "ticr_howtogetQuartetCFs.md"]
Depth = 2
```

## Functions

```@docs
HybridNetwork
DataCF
Quartet
readTopology
readTopologyLevel1
tipLabels
writeTopology
deleteleaf!
printEdges
printNodes
readTrees2CF
readTableCF
readInputTrees
summarizeDataCF
snaq!
readSnaqNetwork
snaqDebug
topologyMaxQPseudolik!
topologyQPseudolik!
rootatnode!
rootonedge!
directEdges!
preorder!
cladewiseorder!
fittedQuartetCF
plotNetGraphViz
plot
setLength!
setGamma!
mapAllelesCFtable
deleteHybridThreshold!
displayedTrees
majorTree
minorTreeAt
displayedNetworkAt!
hardwiredClusters
hardwiredCluster
hardwiredClusterDistance
treeEdgesBootstrap
hybridDetection
summarizeHFdf
hybridBootstrapSupport
bootsnaq
readBootstrapTrees
readMultiTopology
hybridatnode!
```

## Index

```@index
```