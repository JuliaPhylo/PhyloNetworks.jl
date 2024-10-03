# Network manipulation

The package contains a good number of utilities to manipulate phylogenetic
networks. Functions that are not exported are more likely to experience
breaking changes in future versions, but can be used by prefixing their
name with `PhyloNetworks.`.

Here is a list of the most useful functions that can handle networks of any level.
They typically assume a bicombining network, that is, a network in which
each hybrid node has exactly 2 parents (never more).

To traverse or learn something about a network, node or edge, see for example:
- [`tipLabels`](@ref),
  [`PhyloNetworks.descendants`](@ref) for the clade ("hardwired cluster") below an edge,
  [`PhyloNetworks.isdescendant`](@ref),
  [`PhyloNetworks.isconnected`](@ref),
  [`PhyloNetworks.getconnectingedge`](@ref)
- [`displayedTrees`](@ref), [`majorTree`](@ref),
  [`biconnectedComponents`](@ref), [`PhyloNetworks.blobInfo`](@ref)
- [`hardwiredCluster`](@ref), [`hardwiredClusters`](@ref)
- [`getroot`](@ref), [`isrootof`](@ref),
  [`isleaf`](@ref PhyloNetworks.isrootof), [`isexternal`](@ref PhyloNetworks.isrootof),
  [`isparentof`](@ref), [`ischildof`](@ref PhyloNetworks.isparentof),
  [`hassinglechild`](@ref),
- [`getchild`](@ref), [`getchildren`](@ref PhyloNetworks.getchild),
  [`getchildedge`](@ref PhyloNetworks.getchild)
  [`getparent`](@ref), [`getparents`](@ref PhyloNetworks.getparent),
  [`getparentminor`](@ref PhyloNetworks.getparent),
  [`getparentedge`](@ref PhyloNetworks.getparent),
  [`getparentedgeminor`](@ref PhyloNetworks.getparent),
  [`getpartneredge`](@ref)
- [`istimeconsistent`](@ref)

To modify a network, for example:
- [`rootonedge!`](@ref), [`rootatnode!`](@ref):
  very useful to root with an outgroup, and [`rotate!`](@ref) to improve plots
- [`deleteleaf!`](@ref),
  [`deleteaboveLSA!`](@ref) (the "least stable ancestor" may be different from the root)
- [`deleteHybridThreshold!`](@ref) to simplify a network by deleting edges with small γ's
- [`removedegree2nodes!`](@ref), [`shrink3cycles!`](@ref), [`shrink2cycles!`](@ref),
  [`PhyloNetworks.shrinkedge!`](@ref)
- [`PhyloNetworks.addleaf!`](@ref),
  [`PhyloNetworks.deletehybridedge!`](@ref),
  [`PhyloNetworks.addhybridedge!`](@ref)
- [`nni!`](@ref) (nearest neighbor interchange),
  [`PhyloNetworks.fliphybrid!`](@ref) to flip the direction of a hybrid edge
- [`PhyloNetworks.unzip_canonical!`](@ref) to "unzip" (or zip down) all
  reticulations, or [`PhyloNetworks.rezip_canonical!`](@ref) to undo.

To compare networks or compare nodes in a network, for example:
- [`hardwiredClusterDistance`](@ref): extends the Robinson-Foulds distance,
  it's a dissimilarity measure on networks
- [`pairwiseTaxonDistanceMatrix`](@ref) for *average* distances,
  [`getNodeAges`](@ref) if ultrametric network,
  [`getnodeheights`](@ref), [`getnodeheights_average`](@ref),
  [`vcv`](@ref)
