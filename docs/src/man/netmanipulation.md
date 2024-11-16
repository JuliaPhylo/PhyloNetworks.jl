```@setup network_getters
using PhyloNetworks
```

# Getting information about a network

The package contains many utilities to extract information about phylogenetic networks.
Functions that are not exported are more likely to experience
breaking changes in future versions, but can be used by prefixing their
name with `PhyloNetworks.` .

Here is a list of the most useful functions.
They typically assume a bicombining network, that is, a network in which
each hybrid node has exactly 2 parents (never more).

## information about the network overall

First, there are a variety of functions that we can use to learn about a network itself:
- [`tipLabels`](@ref) 
- [`getroot`](@ref) gives the index of the root in `net.node`,
- [`pairwiseTaxonDistanceMatrix`](@ref) for *average* distances
- [`vcv`](@ref) to compute the variance covariance matrix of a network under
  a Brownian Motion model
- [`istimeconsistent`](@ref) returns `true` or `false` if time consistent,
- [`getNodeAges`](@ref) if the network is time-consistent and ultrametric
- [`displayedTrees`](@ref) or [`majorTree`](@ref) to get the displayed trees
  or major tree, respectively
- [`hardwiredcluster`](@ref) or [`hardwiredclusters`](@ref) to learn about the
  clusters on a network for a given set of taxa
- [`biconnectedcomponents`](@ref) and [`PhyloNetworks.blobinfo`](@ref)
  give information about the blobs found within a network.
  Related utilities:
  [`PhyloNetworks.biconnectedcomponent_entrynodes`](@ref),
  [`PhyloNetworks.biconnectedcomponent_exitnodes`](@ref),
  [`blobdecomposition`](@ref).
- [`treeedgecomponents`](@ref) are the components after removing hybrid edges
- [`checkroot!`](@ref) to check that the graph is a valid semidirected
  network and its assigned root node is in a admissible rooting position

The following functions all compute distances from the root to each node.
Their differ in how they handle time inconsistency: when the distance from
the root to a node varies across multiple paths from the root to that node.
- [`getnodeheights`](@ref)
- [`getnodeheights_average`](@ref)
- [`getnodeheights_majortree`](@ref)

## information about nodes

To learn about the nodes related connected to a given node, one can use the
following functions:

- [`PhyloNetworks.isdescendant`](@ref) and [`PhyloNetworks.isconnected`](@ref) 
  can be used to learn about the relationship between two nodes
- [`hassinglechild`](@ref)
- [`getparent`](@ref) or [`getparents`](@ref PhyloNetworks.getparent) if the
  node is a hybrid
- [`getparentminor`](@ref PhyloNetworks.getparent)
- [`getchild`](@ref) or [`getchildren`](@ref PhyloNetworks.getchild) if the
  node has more than one child
- [`isrootof`](@ref)

To learn about the edges connected to a given node, one can use:
- [`getchildedge`](@ref PhyloNetworks.getchild)
- [`PhyloNetworks.getconnectingedge`](@ref)
- [`getparentedge`](@ref PhyloNetworks.getparent)
- [`getparentedgeminor`](@ref PhyloNetworks.getparent) for a hybrid node

## information about edges

To learn about the nodes related connected to a given edge, one can use the following functions: 

- [`PhyloNetworks.descendants`](@ref) for the clade ("hardwired cluster") below an edge
- [`getparent`](@ref)
- [`getchild`](@ref)
- [`getpartneredge`](@ref) gives the hybrid partner of an edge,
  if it is the parent edge of a hybrid node.
- [`isparentof`](@ref), [`ischildof`](@ref PhyloNetworks.isparentof) can
  inform whether an edge and node are connected

# Modifying a network

To modify some of the core components of a network:

- [`rootonedge!`](@ref), [`rootatnode!`](@ref):
  very useful to root with an outgroup, and [`rotate!`](@ref) to improve plots
- [`nni!`](@ref) to perform a semidirected nearest neighbor interchange
- [`PhyloNetworks.fliphybrid!`](@ref) to flip the direction of a hybrid edge
- [`PhyloNetworks.unzip_canonical!`](@ref) to "unzip" (or zip down) all
  reticulations, or [`PhyloNetworks.rezip_canonical!`](@ref) to undo.

To remove components from a network:

- [`deleteleaf!`](@ref)
- [`deleteaboveLSA!`](@ref) (the "least stable ancestor" may be different from the root)
- [`deleteHybridThreshold!`](@ref) to simplify a network by deleting edges with small γ's
- [`removedegree2nodes!`](@ref), [`shrink3cycles!`](@ref), [`shrink2cycles!`](@ref),
  [`PhyloNetworks.shrinkedge!`](@ref)
- [`PhyloNetworks.deletehybridedge!`](@ref)

To add components to a network:

- [`PhyloNetworks.addleaf!`](@ref)
- [`PhyloNetworks.addhybridedge!`](@ref)

# Comparing two networks

- [`hardwiredclusterdistance`](@ref): extends the Robinson-Foulds distance.
  It's a dissimilarity measure on networks: a dissimilarity of 0 does not
  guarantee that the 2 networks have the same topology in general.
  But it does if the networks are in some classes (e.g. trees, level-1,
  tree-child, and others).