```@setup network_getters
using PhyloNetworks
```

# Network manipulation

The package contains many utilities to extract information about phylogenetic networks, and to modify networks.
Functions that are not exported are more likely to experience
breaking changes in future versions, but can be used by prefixing their
name with `PhyloNetworks.` .

Below is a list of the most useful functions.
They typically assume a bicombining network, that is, a network in which
each hybrid node has exactly 2 parents (never more).

## Getting information on a network

### overall network information

- [`tiplabels`](@ref) for taxon labels
- [`getroot`](@ref) gives the root node
- [`leaststableancestor`](@ref) give the least stable ancestor (LSA) of
  all taxa (their MRCA on a tree), which may be below the root
- [`pairwisetaxondistancematrix`](@ref) for *average* distances
- [`vcv`](@ref) for the variance-covariance matrix between taxa under a
  Brownian Motion model along the network, and [`sharedpathmatrix`](@ref)
  for the variance-covariance between all nodes (not just leaves)
- [`istimeconsistent`](@ref): `true` or `false`, to know if for all nodes,
  the various paths from the root to that node have the same length
  (as expected if length was proportional to time)
- [`getnodeages`](@ref) assuming the network is time-consistent and ultrametric
- [`displayedtrees`](@ref) or [`majortree`](@ref) to get the displayed trees
  or major tree, respectively
- [`hardwiredclusters`](@ref) to get all clusters of taxa on a network
- [`checkroot!`](@ref) to check that the graph is a valid semidirected
  network with a root node in a admissible rooting position

on the network complexity:
- [`istreechild`](@ref)
- [`hashybridladder`](@ref) -- a hybrid ladder is also called a "stack":
  one hybrid parent of another
- [`isgalled`](@ref) in the sense of a galled *network*
- [`getlevel`](@ref)
- [`biconnectedcomponents`](@ref) and
  [`PhyloNetworks.process_biconnectedcomponents!`](@ref)
  to calculate (and store) the blobs in a network;
  also
  [`PhyloNetworks.biconnectedcomponent_entrynodes`](@ref),
  [`PhyloNetworks.biconnectedcomponent_exitnodes`](@ref),
  and [`blobdecomposition`](@ref)
- [`treeedgecomponents`](@ref) are the components after removing hybrid edges

The following functions all compute distances from the root to each node.
Their differ in how they handle time inconsistency: when the distance from
the root to a node varies across multiple paths from the root to that node.
- [`getnodeheights`](@ref)
- [`getnodeheights_average`](@ref)
- [`getnodeheights_majortree`](@ref)

### node information

To learn about nodes and how some might be related or connected, one can use:

- [`PhyloNetworks.isdescendant`](@ref) and [`PhyloNetworks.isconnected`](@ref) 
  can be used to learn about the relationship between two nodes
- [`hassinglechild`](@ref)
- [`getparent`](@ref) or [`getparents`](@ref PhyloNetworks.getparent) if the
  node is a hybrid
- [`getparentminor`](@ref PhyloNetworks.getparent)
- [`getchild`](@ref) or [`getchildren`](@ref PhyloNetworks.getchild) if the
  node has more than one child
- [`isrootof`](@ref)
- [`leaststableancestor_matrix`](@ref) for the LSA of all pairs of taxa
  (or MRCA on a tree: most recent common ancestor)

To learn about the edges connected to a given node, one can use:
- [`getchildedge`](@ref PhyloNetworks.getchild)
- [`PhyloNetworks.getconnectingedge`](@ref)
- [`getparentedge`](@ref PhyloNetworks.getparent)
- [`getparentedgeminor`](@ref PhyloNetworks.getparent) for a hybrid node

### edge information

To learn about the nodes related connected to a given edge, one can use the following functions: 

- [`PhyloNetworks.descendants`](@ref) for the clade ("hardwired cluster") below an edge
- [`getparent`](@ref)
- [`getchild`](@ref)
- [`getpartneredge`](@ref) gives the hybrid partner of an edge,
  if it is the parent edge of a hybrid node.
- [`isparentof`](@ref), [`ischildof`](@ref PhyloNetworks.isparentof) can
  inform whether an edge and node are connected
- [`hardwiredcluster`](@ref) to get the cluster of taxa below a particular edge

## Modifying a network

To modify some of the core components of a network:

- [`rootonedge!`](@ref), [`rootatnode!`](@ref):
  very useful to root with an outgroup, and [`rotate!`](@ref) to improve plots
- [`nni!`](@ref) to perform a semidirected nearest neighbor interchange
- [`PhyloNetworks.fliphybrid!`](@ref) to flip the direction of a hybrid edge
- [`PhyloNetworks.unzip_canonical!`](@ref) to "unzip" (or zip down) all
  reticulations, or [`PhyloNetworks.rezip_canonical!`](@ref) to undo.
- [`setlength!`](@ref) and [`setlengths!`](@ref) to change the length of one
  or more edges,
  and [`setgamma!`](@ref) to change a hybrid edge inheritance γ

To remove components from a network:

- [`deleteleaf!`](@ref)
- [`deleteaboveLSA!`](@ref): the "least stable ancestor" may be different
  from the root
- [`deletehybridthreshold!`](@ref) to simplify a network by deleting edges with small γ's
- [`PhyloNetworks.shrinkedge!`](@ref) to contract an edge
- [`removedegree2nodes!`](@ref) to suppress degree-2 nodes,
  [`suppressroot!`](@ref) to "unroot"
- [`shrink3cycles!`](@ref) and [`shrink2cycles!`](@ref) to contract "cycles"
  of 2 or 3 edges, which deletes 1 reticulation
- [`PhyloNetworks.deletehybridedge!`](@ref) to remove one hybrid edge
- [`treeofblobs`](@ref) to shrink blobs and get the strict tree-like part of
  the network

To add components to a network:

- [`PhyloNetworks.addleaf!`](@ref)
- [`PhyloNetworks.addhybridedge!`](@ref)

To modify some internal attributes, that don't affect the network topology
or edge parameters:

- [`nameinternalnodes!`](@ref)
- [`PhyloNetworks.resetnodenumbers!`](@ref) and
  [`PhyloNetworks.resetedgenumbers!`](@ref)

To calibrate a network (modify its edge lengths):

- [`calibratefrompairwisedistances!`](@ref). This documentation has little
  about calibration so far, but see this
  [tutorial](https://juliaphylo.github.io/networkPCM-tutorial/topic4-netcalibration.html)

## Comparing networks

- [`hardwiredclusterdistance`](@ref): extends the Robinson-Foulds distance.
  It's a dissimilarity measure on networks: a dissimilarity of 0 does not
  guarantee that the 2 networks have the same topology in general.
  But it does if the networks are in some classes (e.g. trees, level-1,
  tree-child, and others).
- μ-distances [`mudistance_rooted`](@ref) and [`mudistance_semidirected`](@ref):
  they also extend the Robinson-Foulds distance on trees.
  They are dissimilarities on general networks, but distances on some classes
  (including tree-child).

## Comparing taxa

- average pairwise distances on a network: [`pairwisetaxondistancematrix`](@ref),  
  and from data: [`PhyloNetworks.hammingdistancematrix`](@ref), that can be
  followed by [`PhyloNetworks.distancecorrection_JC!`](@ref)
- f2-distances expected from a network: [`expectedf2matrix`](@ref)
  which can be used to get expected f3 and expected f4 statistics:
  [`PhyloNetworks.expectedf3matrix`](@ref) and [`expectedf4table`](@ref).
