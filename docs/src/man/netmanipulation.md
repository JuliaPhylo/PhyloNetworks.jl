```@setup network_getters
using PhyloNetworks
```

# Network Getter Functions 

The package contains many utilities to extract information about phylogenetic networks.
Functions that are not exported are more likely to experience
breaking changes in future versions, but can be used by prefixing their
name with `PhyloNetworks.` .

Here is a list of the most useful functions.
They typically assume a bicombining network, that is, a network in which
each hybrid node has exactly 2 parents (never more).

## Functions to get information about the network

First, there are a variety of functions that we can use to learn about a network itself:
- [`tipLabels`](@ref) 
- [`getroot`](@ref) gives the index of the root in `net.node`,
- [`pairwiseTaxonDistanceMatrix`](@ref) for *average* distances
- [`vcv`](@ref) to compute the variance covariance matrix of a network under Brownian Motion
- [`istimeconsistent`](@ref) returns `true` or `false` if time consistent,
- [`getNodeAges`](@ref) if the network is time-consistent and ultrametric
- [`displayedTrees`](@ref) or [`majorTree`](@ref) to get the displayed trees or major tree, respectively
- [`biconnectedComponents`](@ref) and [`PhyloNetworks.blobInfo`](@ref) give information about the blobs found within a network
- [`hardwiredCluster`](@ref) or [`hardwiredClusters`](@ref) to learn about the clusters on a network for a given set of taxa

The following functions all compute distances from the root to each node. 
- [`getnodeheights`](@ref)
- [`getnodeheights_average`](@ref)
- [`getnodeheights_majortree`](@ref)

When the network is time-consistent, they all give the same output:

``` @repl network_getters
#Read a time-consistent network
consistent_net = readTopology(
  "((A:2.5,#H1:1.5::0.4):0.25,(C:1.5,(B:1)#H1:0.5::0.6):1.25);" );

#These all return the same values because the network is time-consistent
heights = getnodeheights(consistent_net)
heights_average = getnodeheights_average(consistent_net);
heights_major = getnodeheights_majortree(consistent_net);
#All heights are the same
heights == heights_average == heights_major 

```
However, when the network is time-**in**consistent, each function handles the inconsistency differently:

``` @repl network_getters

#Read a time-inconsistent network
inconsistent_net = readTopology("((A:2.5,#H1:1.5::0.4):0.25,(C:1.5,(B:1)#H1:2.5::0.6):1.25);");

#getnodeheights(inconsistent_net) #throws an error
#The average heights across all paths to the time-inconsistent node
getnodeheights_average(inconsistent_net) 
getnodeheights_majortree(inconsistent_net) #Use the major tree heights at inconsistencies. 
```

## Functions to get information on nodes or edges

### Nodes

To learn about the nodes related connected to a given node, one can use the following functions: 
- [`PhyloNetworks.isdescendant`](@ref) and [`PhyloNetworks.isconnected`](@ref) can be used to learn about the relationship between two nodes
- [`hassinglechild`](@ref)
- [`getparent`](@ref) or [`getparents`](@ref PhyloNetworks.getparent) if the node is a hybrid
- [`getparentminor`](@ref PhyloNetworks.getparent)
- [`getchild`](@ref) or [`getchildren`](@ref PhyloNetworks.getchild) if the node has more than one child
- [`isrootof`](@ref)

To learn about the edges connected to a given node, one can use:
- [`getchildedge`](@ref PhyloNetworks.getchild)
- [`PhyloNetworks.getconnectingedge`](@ref)
- [`getparentedge`](@ref PhyloNetworks.getparent)
- [`getparentedgeminor`](@ref PhyloNetworks.getparent) for a hybrid node

### Edges

To learn about the nodes related connected to a given edge, one can use the following functions: 

- [`PhyloNetworks.descendants`](@ref) for the clade ("hardwired cluster") below an edge
- [`getparent`](@ref)
- [`getchild`](@ref)

- [`getpartneredge`](@ref) gives the hybrid partner of an edge, if it is the parent edge of a hybrid node.

- [`isparentof`](@ref), [`ischildof`](@ref PhyloNetworks.isparentof) can inform whether an edge and node are connected 

# Network Modification

To modify some of the core components of a network:
- [`rootonedge!`](@ref), [`rootatnode!`](@ref):
  very useful to root with an outgroup, and [`rotate!`](@ref) to improve plots
- [`nni!`](@ref) (nearest neighbor interchange)
- [`PhyloNetworks.fliphybrid!`](@ref) to flip the direction of a hybrid edge
- [`PhyloNetworks.unzip_canonical!`](@ref) to "unzip" (or zip down) all
  reticulations, or [`PhyloNetworks.rezip_canonical!`](@ref) to undo.

To remove components of the network:
- [`deleteleaf!`](@ref)
- [`deleteaboveLSA!`](@ref) (the "least stable ancestor" may be different from the root)
- [`deleteHybridThreshold!`](@ref) to simplify a network by deleting edges with small γ's
- [`removedegree2nodes!`](@ref), [`shrink3cycles!`](@ref), [`shrink2cycles!`](@ref),
  [`PhyloNetworks.shrinkedge!`](@ref)
- [`PhyloNetworks.deletehybridedge!`](@ref)

To add components to the network:

- [`PhyloNetworks.addleaf!`](@ref)
- [`PhyloNetworks.addhybridedge!`](@ref)



