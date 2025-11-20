@deprecate getHeights getnodeheights false
@deprecate recursion_preorder  traversal_preorder  false
@deprecate recursion_postorder traversal_postorder false
@deprecate readTopology       readnewick
@deprecate readMultiTopology  readmultinewick
@deprecate readBootstrapTrees readmultinewick_files
@deprecate writeTopology      writenewick
@deprecate writeMultiTopology writemultinewick
Base.@deprecate_moved snaq! "SNaQ"
# minor name changes: camelCase to smallcase
@deprecate printEdges printedges
@deprecate printNodes printnodes
@deprecate tipLabels tiplabels
@deprecate setGamma! setgamma!
@deprecate directEdges! directedges!
@deprecate hardwiredClusters hardwiredclusters
@deprecate hardwiredCluster  hardwiredcluster
@deprecate hardwiredCluster! hardwiredcluster!
@deprecate hardwiredClusterDistance hardwiredclusterdistance
@deprecate biconnectedComponents biconnectedcomponents
@deprecate blobDecomposition! blobdecomposition!
@deprecate blobDecomposition  blobdecomposition
@deprecate deleteHybridThreshold! deletehybridthreshold!
@deprecate displayedNetworkAt! displayednetworkat!
@deprecate displayedTrees displayedtrees
@deprecate minorTreeAt    minortreeat
@deprecate majorTree      majortree
@deprecate sharedPathMatrix sharedpathmatrix
@deprecate descendenceMatrix descendencematrix
@deprecate getNodeAges getnodeages
@deprecate calibrateFromPairwiseDistances! calibratefrompairwisedistances!
@deprecate parsimonySoftwired parsimonysoftwired
