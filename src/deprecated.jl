@deprecate getHeights getnodeheights false
@deprecate recursion_preorder  traversal_preorder  false
@deprecate recursion_postorder traversal_postorder false
@deprecate readTopology       readnewick
@deprecate readMultiTopology  readmultinewick
@deprecate readBootstrapTrees readmultinewick_files
@deprecate writeTopology      writenewick
@deprecate writeMultiTopology writemultinewick
v16msg = """
    To roll back to PhyloNetworks v0.16, do this after typing ] to get into package mode:
    pkg> add PhyloNetworks#v0.16.4
    Then, also pin PhyloPlots to version v1.0.1:
    pkg> add PhyloPlots#v1.0.1
    See here for help about adding packages and checking versions:
    https://pkgdocs.julialang.org/v1/managing-packages/#Adding-packages
    """
@deprecate writeTableCF(obj...; kwargs...) writeTableCF()
writeTableCF() =
    error("""writeTableCF works with PhyloNetworks v0.16.
    Starting with PhyloNetworks v1 (and SNaQ v1), please use instead
    tablequartetCF(quartetlist::Vector{QuartetT} [, taxonnames];...)
    The resulting named tuple can be converted to a data frame like so for example:
    DataFrame(tablequartetCF(...), copycols=false)
    """ * "\n" * v16msg)
@deprecate snaq!(obj...; kwargs...) snaq!()
snaq!() =
    error("""snaq! was moved from PhyloNetworks.jl to SNaQ.jl. To use it, do this:
    pkg> add SNaQ
    Otherwise, snaq! is available in PhyloNetworks v0.16.
    """ * "\n" * v16msg)
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
@deprecate pairwiseTaxonDistanceMatrix pairwisetaxondistancematrix!
@deprecate sharedPathMatrix sharedpathmatrix
@deprecate descendenceMatrix descendencematrix
@deprecate getNodeAges getnodeages
@deprecate calibrateFromPairwiseDistances! calibratefrompairwisedistances!
@deprecate parsimonySoftwired parsimonysoftwired
