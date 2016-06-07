
<a id='PhyloNetworks-1'></a>

# `PhyloNetworks`


PhyloNetworks is a Julia package for the manipulation, visualization and inference of phylogenetic networks.  SNaQ implements the statistical inference method in [Sol&iacute;s-Lemus and An&eacute;](http://arxiv.org/pdf/1509.06075.pdf). The procedure involves a numerical optimization of branch lengths and inheritance probabilities and a heuristic search in the space of phylogenetic networks.


<a id='PhyloNetworks.jl-Documentation-1'></a>

# PhyloNetworks.jl Documentation

- [`PhyloNetworks`](index.md#PhyloNetworks-1)
- [PhyloNetworks.jl Documentation](index.md#PhyloNetworks.jl-Documentation-1)
    - [Functions](index.md#Functions-1)
    - [Index](index.md#Index-1)
- [Inserting data into the TICR pipeline at various stages](ticr_howtogetQuartetCFs.md#Inserting-data-into-the-TICR-pipeline-at-various-stages-1)
    - [To run MrBayes: You already have alignments](ticr_howtogetQuartetCFs.md#To-run-MrBayes:-You-already-have-alignments-1)
    - [To run BUCKy: You already have MrBayes output](ticr_howtogetQuartetCFs.md#To-run-BUCKy:-You-already-have-MrBayes-output-1)
- [Simple use of Julia types](simpleJulia.md#Simple-use-of-Julia-types-1)


<a id='Functions-1'></a>

## Functions

<a id='PhyloNetworks.HybridNetwork' href='#PhyloNetworks.HybridNetwork'>#</a>
**`PhyloNetworks.HybridNetwork`** &mdash; *Type*.



`HybridNetwork type` Explicit network or tree with the following attributes:

  * numTaxa
  * numNodes (total number of nodes)
  * numEdges
  * numHybrids (number of hybrid nodes)
  * edge (array of Edges)
  * node (array of Nodes)
  * root (index of root in vector 'node'. May be artificial, for printing and traversal purposes only.)
  * hybrid (array of Nodes: those are are hybrid nodes)
  * leaf (array of Nodes: those that are leaves)
  * loglik (negative log pseudolik after estimation)
  * isRooted (true or false)

<a id='PhyloNetworks.DataCF' href='#PhyloNetworks.DataCF'>#</a>
**`PhyloNetworks.DataCF`** &mdash; *Type*.



`DataCF type`

type that contains the following attributes:

  * quartet (vector of Quartets)
  * numQuartets
  * tree (vector of trees: empty if a table of CF was input instead of list of trees)
  * numTrees (-1 if a table CF was input instead of list of trees)
  * repSpecies (taxon names that were repeated in table of CF or input gene trees: used inside snaq for multiple alleles case)

The list of Quartet may be accessed with the attribute .quartet. If the input was a list of trees, the HybridNetwork's can be accessed with the attribute .tree. For example, if the DataCF object is named d, d.quartet[1] will show the first quartet and d.tree[1] will print the first input tree.

<a id='PhyloNetworks.Quartet' href='#PhyloNetworks.Quartet'>#</a>
**`PhyloNetworks.Quartet`** &mdash; *Type*.



`Quartet type`

type that saves the information on a given 4-taxon subset. It contains the following attributes:

  * number
  * taxon (vector of taxon names)
  * obsCF (vector of observed CF)
  * logPseudoLik
  * ngenes (number of gene trees used to compute the observed CF: -1 if unknown)
  * qnet (internal topological structure that saves the expCF after snaq estimation to emphasize that the expCF depend on a specific network, not the data)

<a id='PhyloNetworks.readTopology' href='#PhyloNetworks.readTopology'>#</a>
**`PhyloNetworks.readTopology`** &mdash; *Function*.



`readTopology(file name); readTopology(parenthetical description)`

function to read tree or network topology from parenthetical format. Input: text file or parenthetical format directly. The file name may not start with a left parenthesis, otherwise the file name itself would be interpreted as the parenthetical description.

<a id='PhyloNetworks.readTopologyLevel1' href='#PhyloNetworks.readTopologyLevel1'>#</a>
**`PhyloNetworks.readTopologyLevel1`** &mdash; *Function*.



`readTopologyLevel1(filename); readTopologyLevel1(parenthetical format)`

same as readTopology, reads a tree or network from parenthetical format, but this function enforces the necessary conditions for any starting topology in SNaQ: non-intersecting cycles, no polytomies, unrooted. It sets any missing branch length to 1.0.

<a id='PhyloNetworks.tipLabels' href='#PhyloNetworks.tipLabels'>#</a>
**`PhyloNetworks.tipLabels`** &mdash; *Function*.



`tipLabels(net::HybridNetwork)`

returns a vector of taxon names (at the leaves) from a HybridNetwork object

<a id='PhyloNetworks.writeTopology' href='#PhyloNetworks.writeTopology'>#</a>
**`PhyloNetworks.writeTopology`** &mdash; *Function*.



```
writeTopology(net)
writeTopology(net, filename)
```

write the parenthetical format of a HybridNetwork object, as a string or to a file. Optional arguments (default values):

  * di (false): write in format for Dendroscope
  * round (false): rounds branch lengths and heritabilities γ
  * digits (3): digits after the decimal place for rounding

If the current root placement is not admissible, other placements are tried. The network is updated with this new root placement, if successful.

<a id='PhyloNetworks.deleteleaf!' href='#PhyloNetworks.deleteleaf!'>#</a>
**`PhyloNetworks.deleteleaf!`** &mdash; *Function*.



`deleteleaf!(HybridNetwork,Int64; index=false, simplify=true)` `deleteleaf!(HybridNetwork,leafName::AbstractString; simplify=true)` `deleteleaf!(HybridNetwork,Node; simplify=true)`

Deletes a leaf node from the network, possibly from its name, number, or index in the network's array of nodes.

simplify: if true and if deleting the node results in 2 hybrid edges forming a cycle of k=2 nodes, then these hybrid edges are merged and simplified as a single tree edge.

The first version does **not** require that `node` is a leaf, so might be used to remove nodes of degree 2. The other versions do, and use the default simplify=true.

Warning: does **not** update attributes related to level-1 networks, such as inCycle, partition, gammaz, etc. Does not require branch lengths, and designed to work on networks of all levels.

<a id='PhyloNetworks.printEdges' href='#PhyloNetworks.printEdges'>#</a>
**`PhyloNetworks.printEdges`** &mdash; *Function*.



`printEdges(net::HybridNetwork)`

prints the information on the edges of net: edge number, node numbers of nodes attached to it, in which cycle it is contained (-1 if no cycle), can it contain root, is it an identifiable edge, length, is it hybrid, gamma value

<a id='PhyloNetworks.printNodes' href='#PhyloNetworks.printNodes'>#</a>
**`PhyloNetworks.printNodes`** &mdash; *Function*.



`printNodes(net::HybridNetwork)`

prints information on the nodes of net: node number, in which cycle it is contained (-1 if no cycle), is it hybrid, does it has hybrid edges, edges number attached to it

<a id='PhyloNetworks.readTrees2CF' href='#PhyloNetworks.readTrees2CF'>#</a>
**`PhyloNetworks.readTrees2CF`** &mdash; *Function*.



```
readTrees2CF(treefile)
readTrees2CF(vector of trees)
```

Read trees in parenthetical format from a file, or take a vector of trees already read, and calculate the proportion of these trees having a given quartet (concordance factor: CF), for all quartets or for a sample of quartets. Optional arguments include:

  * quartetfile: name of text file with list of 4-taxon subsets to be analyzed. If none is specified, the function will list all possible 4-taxon subsets.
  * whichQ="rand": to choose a random sample of 4-taxon subsets
  * numQ: size of random sample (ignored if whichQ is not set to "rand")
  * writeTab=false: does not write the observedCF to a table (default true)
  * CFfile: name of file to save the observedCF (default tableCF.txt)
  * writeQ=true: save intermediate files with the list of all 4-taxon subsets and chosen random sample (default false).

<a id='PhyloNetworks.readTableCF' href='#PhyloNetworks.readTableCF'>#</a>
**`PhyloNetworks.readTableCF`** &mdash; *Function*.



```
readTableCF(file)
readTableCF(data frame)
```

Read a file or DataFrame object containing a table of concordance factors (CF), with one row per 4-taxon set. The first 4 columns are assumed to give the labels of the 4 taxa in each set (tx1, tx2, tx3, tx4). Columns containing the CFs are assumed to be named 'CF12_34', 'CF13_24' and 'CF14_23', or else are assumed to be columns 5,6,7. If present, a column named 'ngenes' will be used to get the number of loci used to estimate the CFs for each 4-taxon set.

Optional arguments:

  * summaryfile: if specified, a summary file will be created with that name.
  * sep (for the second form only): to specify the type of separator in the file, with single quotes: sep=';'.

<a id='PhyloNetworks.readInputTrees' href='#PhyloNetworks.readInputTrees'>#</a>
**`PhyloNetworks.readInputTrees`** &mdash; *Function*.



`readInputTrees(file)`

function to read a text file with a list of trees in parenthetical format (one tree per line), it returns an array of HybridNetwork object.

<a id='PhyloNetworks.summarizeDataCF' href='#PhyloNetworks.summarizeDataCF'>#</a>
**`PhyloNetworks.summarizeDataCF`** &mdash; *Function*.



`summarizeDataCF(d::DataCF)`

function to summarize the information contained in a DataCF object. It has the following optional arguments: - filename: if provided, the summary will be saved in the filename, not to screen - pc (number between (0,1)): threshold of percentage of missing genes to identify 4-taxon subsets with fewer genes than the threshold

<a id='PhyloNetworks.snaq!' href='#PhyloNetworks.snaq!'>#</a>
**`PhyloNetworks.snaq!`** &mdash; *Function*.



`snaq!(T::HybridNetwork, d::DataCF)`

Estimate the network (or tree) to fit observed concordance factors (CFs) stored in a DataCF object, using maximum pseudo-likelihood. The search starts from topology `T`, which can be a tree or a network with no more than `hmax` hybrid nodes. The function name ends with ! because it modifies the CF data `d` by updating its attributes `expCF`: CFs expected under the network model. It does *not* modify `T`. The quartet pseudo-deviance is the negative log pseudo-likelihood, up to an additive constant, such that a perfect fit corresponds to a deviance of 0.0.

There are many optional arguments, including

  * hmax: maximum number of hybridizations allowed (default 1)
  * verbose: if true, it prints information about the numerical optimization
  * runs: number of independent starting points for the search (default 10)
  * outgroup: outgroup taxon to root the estimated topology at the very end
  * filename: root name for the output files. Default is "snaq". If empty (""),   files are *not* created, progress log goes to the screen only (standard out)   and the best network is returned.
  * seed: seed to replicate a given search

See also: `topologyMaxQPseudolik!` to optimize parameters on a fixed topology, and `topologyQPseudolik!` to get the deviance (pseudo log-likelihood up to a constant) of a fixed topology with fixed parameters.

Reference: Claudia Solís-Lemus and Cécile Ané (2016). Inferring phylogenetic networks with maximum pseudolikelihood under incomplete lineage sorting. [PLoS Genetics](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896) 12(3):e1005896

<a id='PhyloNetworks.readSnaqNetwork' href='#PhyloNetworks.readSnaqNetwork'>#</a>
**`PhyloNetworks.readSnaqNetwork`** &mdash; *Function*.



`readSnaqNetwork(output file)`

function to read the estimated network from an .out file generated by the snaq function

<a id='PhyloNetworks.snaqDebug' href='#PhyloNetworks.snaqDebug'>#</a>
**`PhyloNetworks.snaqDebug`** &mdash; *Function*.



`snaqDebug(currT::HybridNetwork,d::DataCF)`

function to replicate a given run that produces error according to the .err file generated by snaq. The same settings used in that run should be used in this function, specially the seed. See the readme file online for more details.

<a id='PhyloNetworks.topologyMaxQPseudolik!' href='#PhyloNetworks.topologyMaxQPseudolik!'>#</a>
**`PhyloNetworks.topologyMaxQPseudolik!`** &mdash; *Function*.



`topologyMaxQPseudolik!(net::HybridNetwork, d::DataCF)`

Estimate the branch lengths and inheritance probabilities (γ's) for a given network topology. The network is *not* modified, only the object `d` is, with updated expected concordance factors.

Ouput: new network, with optimized parameters (branch lengths and gammas). The maximized quartet pseudo-deviance is the negative log pseudo-likelihood, up to an additive constant, such that a perfect fit corresponds to a deviance of 0.0. This is also an attribute of the network, which can be accessed with `net.loglik`.

Optional arguments (default value):

  * verbose (false): if true, information on the numerical optimization is printed to screen
  * ftolRel (1e-5), ftolAbs (1e-6), xtolRel (1e-3), xtolAbs (1e-4):   absolute and relative tolerance values for the pseudo-deviance function   and the parameters

<a id='PhyloNetworks.topologyQPseudolik!' href='#PhyloNetworks.topologyQPseudolik!'>#</a>
**`PhyloNetworks.topologyQPseudolik!`** &mdash; *Function*.



`topologyQPseudolik!(net::HybridNetwork, d::DataCF)`

Calculate the quartet pseudo-deviance of a given network/tree for DataCF `d`. This is the negative log pseudo-likelihood, up to an additive constant, such that a perfect fit corresponds to a deviance of 0.0.

Be careful if the net object does not have all internal branch lengths specified because then the pseudolikelihood will be meaningless.

The loglik attribute of the network is undated, and `d` is updated with the expected concordance factors under the input network.

<a id='PhyloNetworks.rootatnode!' href='#PhyloNetworks.rootatnode!'>#</a>
**`PhyloNetworks.rootatnode!`** &mdash; *Function*.



`rootatnode!(HybridNetwork, nodeNumber::Int64; index=false::Bool)`

`rootatnode!(HybridNetwork, Node)`

`rootatnode!(HybridNetwork, nodeName::AbstractString)`

Roots the network/tree object at the node with name 'nodeName' or number 'nodeNumber' (by default) or with index 'nodeNumber' if index=true. Attributes isChild1 and containRoot are updated along the way. Use `plot(net, showNodeNumber=true, showEdgeLength=false)` to visualize and identify a node of interest.

Returns the network.

Warnings: - If the node is a leaf, the root will be placed along   the edge adjacent to the leaf, with a message. This might add a new node. - If the desired root placement is incompatible with one or more hybrids, then   * a RootMismatch error is thrown   * the input network will still have some attributes modified.

See also: `rootonedge!`.

<a id='PhyloNetworks.rootonedge!' href='#PhyloNetworks.rootonedge!'>#</a>
**`PhyloNetworks.rootonedge!`** &mdash; *Function*.



`rootonedge!(HybridNetwork, edgeNumber::Int64; index=false::Bool)`

`rootonedge!(HybridNetwork, Edge)`

Roots the network/tree object along an edge with number 'edgeNumber' (by default) or with index 'edgeNumber if index=true. Attributes isChild1 and containRoot are updated along the way.

This adds a new node and a new edge to the network. Use `plot(net, showEdgeNumber=true, showEdgeLength=false)` to visualize and identify an edge of interest.

See also: `rootatnode!`.

<a id='PhyloNetworks.directEdges!' href='#PhyloNetworks.directEdges!'>#</a>
**`PhyloNetworks.directEdges!`** &mdash; *Function*.



`directEdges!(net::HybridNetwork; checkMajor=true::Bool)`

Updates the edges' attribute `isChild1`, according to the root placement. Also updates edges' attribute `containRoot`, for other possible root placements compatible with the direction of existing hybrid edges. Relies on hybrid nodes having exactly 1 major hybrid parent edge, but checks for that if checkMajor=true.

Warning: Assumes that isChild1 is correct on hybrid edges (to avoid changing the identity of which nodes are hybrids and which are not).

Returns the network. Throws a 'RootMismatch' Exception if the root was found to conflict with the direction of any hybrid edge.

<a id='PhyloNetworks.preorder!' href='#PhyloNetworks.preorder!'>#</a>
**`PhyloNetworks.preorder!`** &mdash; *Function*.



`preorder!(net::HybridNetwork)`

Updates attribute net.nodes_changed in which the nodes are pre-ordered (also called topological sorting), such that each node is visited after its parent(s). The edges' direction needs to be correct before calling preorder!, using directEdges!

<a id='PhyloNetworks.cladewiseorder!' href='#PhyloNetworks.cladewiseorder!'>#</a>
**`PhyloNetworks.cladewiseorder!`** &mdash; *Function*.



`cladewiseorder!(net::HybridNetwork)`

Updates attribute net.cladewiseorder_nodeIndex. Used for plotting the network. In the major tree, all nodes in a given clade are consecutive. On a tree, this function also provides a pre-ordering of the nodes. The edges' direction needs to be correct before calling cladewiseorder!, using directEdges!

<a id='PhyloNetworks.fittedQuartetCF' href='#PhyloNetworks.fittedQuartetCF'>#</a>
**`PhyloNetworks.fittedQuartetCF`** &mdash; *Function*.



`fittedQuartetCF(d::DataCF)`

return a data frame with the observed and expected quartet concordance factors
after estimation of a network with snaq(T,d).
The format can be :wide (default) or :long.

* if wide, the output has one row per 4-taxon set, and each row has 10 columns: 4 columns
  for the taxon names, 3 columns for the observed CFs and 3 columns for the expected CF.
* if long, the output has one row per quartet, i.e. 3 rows per 4-taxon sets, and 7 columns:
  4 columns for the taxon names, one column to give the quartet resolution, one column for
  the observed CF and the last column for the expected CF.

<a id='PhyloNetworks.plotNetGraphViz' href='#PhyloNetworks.plotNetGraphViz'>#</a>
**`PhyloNetworks.plotNetGraphViz`** &mdash; *Function*.



`plotNetGraphViz(net::HybridNetwork)`

function to plot a HybridNetwork object. The plot will be saved in the working directory as a svg file. We are working on allowing other file formats, and to have the plot pop out in a window. This function has the following optional arguments: - imageName: name for plot file (default netImage) - mainTree: if true, only the underlying tree (with major hybrid edges) is plotted (default false) - width: width of image in inches (default 6) - height: height of image in inches (default 8) - vert: if true, plot displayed from top to bottom (default tre) - internalLabels: if true, prints number for internal nodes (default false) - fontSize: font size for taxon names (default 16.0) - hybridColor: color for hybrid edges (default green4) - unrooted: if true, prints the topology unrooted - nodeSeparation: minimum distance between nodes in inches (default 0.8) - labelAngle: angle for leaf label placement (default 180.0) - labelDistance: distance for leaf label placement (default 3.0) - includeGamma: if true, includes the gamma values in the plot (default false) - includeLength: if true, includes the branch lengths in the plot (default false

<a id='Gadfly.plot' href='#Gadfly.plot'>#</a>
**`Gadfly.plot`** &mdash; *Function*.



```
plot(net::HybridNetwork; useEdgeLength=false, mainTree=false, showTipLabel=true,
     showNodeNumber=false, showEdgeLength=false, showGamma=false, edgeColor=colorant"black",
     majorHybridEdgeColor=colorant"deepskyblue4", minorHybridEdgeColor=colorant"deepskyblue",
     showEdgeNumber=false, showIntNodeLabel=false, edgeLabel=[], nodeLabel=[])
```

Plots a network, from left to right.

  * useEdgeLength: if true, the tree edges and major hybrid edges are   drawn proportionally to their length. Minor hybrid edges are not, however.   Note that edge lengths in coalescent units may scale very poorly with time.
  * mainTree: if true, the minor hybrid edges are ommitted.
  * showTipLabel: if true, taxon labels are shown. You may need to zoom out to see them.
  * showNodeNumbers: if true, nodes are labelled with the number used internally.
  * showEdgeLength: if true, edges are labelled with their length (above)
  * showGamma: if true, hybrid edges are labelled with their heritability (below)
  * edgeColor: color for tree edges. black by default.
  * majorHybridEdgeColor: color for major hybrid edges
  * minorHybridEdgeColor: color for minor hybrid edges
  * showEdgeNumber: if true, edges are labelled with the number used internally.
  * showIntNodeLabel: if true, internal nodes are labelled with their names.   Useful for hybrid nodes, which do have tags like '#H1'.
  * edgeLabel: dataframe with two columns: the first with edge numbers, the second with labels   (like bootstrap values) to annotate edges. empty by default.
  * nodeLabel: dataframe with two columns: the first with node numbers, the second with labels   (like bootstrap values for hybrid relationships) to annotate nodes. empty by default.

Note that plot() actually modifies some (minor) attributes of the network, as it calls directEdges!, preorder! and cladewiseorder!.

If hybrid edges cross tree and major edges, you may choose to rotate some tree edges to eliminate crossing edges, using rotate!.

<a id='PhyloNetworks.setLength!' href='#PhyloNetworks.setLength!'>#</a>
**`PhyloNetworks.setLength!`** &mdash; *Function*.



`setLength!(Edge,new length)`

set a new length for an object Edge. The new length needs to be positive. For example, if you have a HybridNetwork object net, and do printEdges(net), you can see the list of Edges and their lengths. You can then change the length of the 3rd edge with setLength!(net.edge[3],1.2).

<a id='PhyloNetworks.setGamma!' href='#PhyloNetworks.setGamma!'>#</a>
**`PhyloNetworks.setGamma!`** &mdash; *Function*.



`setGamma!(Edge,new gamma)`

sets a gamma value for a hybrid edge (it has to be a hybrid edge) and new gamma needs to be (0,1). The function will automatically change the gamma value for the other hybrid edge to 1-gamma. For example, if you have a HybridNetwork object net, and do printEdges(net), you can see the list of Edges and their gammas. You can then change the length of the hybrid edge (assume it is in position 3) with setGamma!(net.edge[3],0.2). This will automatically set the gamma for the other hybrid edge to 0.8.

<a id='PhyloNetworks.mapAllelesCFtable' href='#PhyloNetworks.mapAllelesCFtable'>#</a>
**`PhyloNetworks.mapAllelesCFtable`** &mdash; *Function*.



`mapAllelesCFtable(mapping file, cf file)`

function that change the allele names in the CF table to species names. The new DataFrame object is returned. Optional argument: filename for the resulting CF table. If not specified, then no CF is saved as file.

<a id='PhyloNetworks.deleteHybridThreshold!' href='#PhyloNetworks.deleteHybridThreshold!'>#</a>
**`PhyloNetworks.deleteHybridThreshold!`** &mdash; *Function*.



`deleteHybridThreshold!(net::HybridNetwork,gamma::Float64)`

Deletes from a network all hybrid edges with heritability below a threshold gamma. Returns the network.

  * if gamma<0.5: deletes     minor hybrid edges with gamma value <  threshold
  * if gamma=0.5: deletes all minor hybrid edges (i.e gamma value <= threshold)

Warning: assumes correct isMajor attributes.

<a id='PhyloNetworks.displayedTrees' href='#PhyloNetworks.displayedTrees'>#</a>
**`PhyloNetworks.displayedTrees`** &mdash; *Function*.



`displayedTrees(net::HybridNetwork, gamma::Float64)`

Warning: assumes correct isMajor attributes.

Extracts all trees displayed in a network, following hybrid edges with heritability >= gamma threshold (or >0.5 if threshold=0.5) and ignoring any hybrid edge with heritability lower than gamma. Returns an array of trees, as HybridNetwork objects.

<a id='PhyloNetworks.majorTree' href='#PhyloNetworks.majorTree'>#</a>
**`PhyloNetworks.majorTree`** &mdash; *Function*.



`majorTree(net::HybridNetwork)`

Warning: assumes correct isMajor attributes.

Extracts the major tree displayed in a network, keeping the major edge and dropping the minor edge at each hybrid node. Returns a HybridNetwork object.

<a id='PhyloNetworks.minorTreeAt' href='#PhyloNetworks.minorTreeAt'>#</a>
**`PhyloNetworks.minorTreeAt`** &mdash; *Function*.



`minorTreeAt(net::HybridNetwork, hybindex::Int64)`

Warning: assumes correct isMajor attributes.

Extracts the tree displayed in the network, following the major hybrid edge at each hybrid node, except at the ith hybrid node (i=hybindex), where the minor hybrid edge is kept instead of the major hybrid edge.

<a id='PhyloNetworks.displayedNetworkAt!' href='#PhyloNetworks.displayedNetworkAt!'>#</a>
**`PhyloNetworks.displayedNetworkAt!`** &mdash; *Function*.



`displayedNetworkAt!(net::HybridNetwork, node::Node)`

Warning: assumes correct isMajor attributes.

Deletes all the minor hybrid edges, except at input node. The network is left with a single hybridization, and otherwise displays the same major tree as before.

<a id='PhyloNetworks.hardwiredClusters' href='#PhyloNetworks.hardwiredClusters'>#</a>
**`PhyloNetworks.hardwiredClusters`** &mdash; *Function*.



`hardwiredClusters(net::HybridNetwork, S::Union{Vector{ASCIIString},Vector{Int64}})`

Returns a matrix describing all the hardwired clusters in a network. Warnings: Clusters are rooted, so the root must be correct.           Allows for missing taxa, with entries all 0.

Each row corresponds to one internal edge, that is, external edges are excluded. If the root is a leaf node, the external edge to that leaf is included (first row). Both parent hybrid edges to a given hybrid node only contribute a single row (they share the same hardwired cluster).

  * first column: edge number
  * next columns: 0/1 values. 1=descendant of edge, 0=not a descendant, or missing taxon.
  * last column:  10/11 values. 10=tree edge, 11=hybrid edge

<a id='PhyloNetworks.hardwiredCluster' href='#PhyloNetworks.hardwiredCluster'>#</a>
**`PhyloNetworks.hardwiredCluster`** &mdash; *Function*.



```
hardwiredCluster(edge::Edge,taxa::Union{Vector{ASCIIString},Vector{Int64}})
hardwiredCluster!(v::Vector{Bool},edge::Edge,taxa::Union{Vector{ASCIIString},Vector{Int64}})
hardwiredCluster!(v::Vector{Bool},edge::Edge,taxa::Union{Vector{ASCIIString},Vector{Int64}},
                  visited::Vector{Int64})
```

Calculate the hardwired cluster of `node`, coded a vector of booleans: true for taxa that are descendent of nodes, false for other taxa (including missing taxa).

The node should belong in a rooted network for which isChild1 is up-to-date. Run directEdges! beforehand. This is very important, otherwise one might enter an infinite loop, and the function does not test for this.

visited: vector of node numbers, of all visited nodes.

**Examples: #"**

```
julia> net5 = "(A,((B,#H1),(((C,(E)#H2),(#H2,F)),(D)#H1)));" |> readTopology |> directEdges! ;

julia> taxa = net5 |> tipLabels # ABC EF D
6-element Array{ASCIIString,1}:
 "A"
 "B"
 "C"
 "E"
 "F"
 "D"

julia> hardwiredCluster(net5.edge[12], taxa) # descendants of 12th edge = CEF
6-element Array{Bool,1}:
 false
 false
  true
  true
  true
 false
```

<a id='PhyloNetworks.hardwiredClusterDistance' href='#PhyloNetworks.hardwiredClusterDistance'>#</a>
**`PhyloNetworks.hardwiredClusterDistance`** &mdash; *Function*.



`hardwiredClusterDistance(net1::HybridNetwork, net2::HybridNetwork, rooted::Bool)`

Takes 2 networks and returns their hardwired cluster distance, that is, the number of hardwired clusters found in one network and not in the other. Note that this is not a distance per se on the full space of hybrid networks: there are pairs of different networks for which this measure is 0. But it is a distance on some network subspaces.

If the 2 networks are trees, this is the Robinson-Foulds distance. If rooted=false, the trees are considered unrooted.

<a id='PhyloNetworks.treeEdgesBootstrap' href='#PhyloNetworks.treeEdgesBootstrap'>#</a>
**`PhyloNetworks.treeEdgesBootstrap`** &mdash; *Function*.



`treeEdgesBootstrap(net::Vector{HybridNetwork}, net0::HybridNetwork)`

read an array of bootstrap networks (net) and a reference network (net0), and calculates the bootstrap support of the tree edges in the reference network.

return a data frame with one row per tree edge and two columns: edge number, bootstrap support

<a id='PhyloNetworks.hybridDetection' href='#PhyloNetworks.hybridDetection'>#</a>
**`PhyloNetworks.hybridDetection`** &mdash; *Function*.



`hybridDetection(net::Vector{HybridNetwork}, net1::HybridNetwork, outgroup::AbstractString)`

function can only compare hybrid nodes in networks that have the same underlying major tree also, need to root all networks in the same place, and the root has to be compatible with the direction of the hybrid edges

it computes the rooted hardwired distance between networks, the root matters. input: vector of bootstrap networks (net), estimated network (net1), outgroup

returns

  * a matrix with one row per bootstrap network, and 2*number of hybrids in net1, column i corresponds to whether hybrid i (net1.hybrid[i]) is found in the bootstrap network, column 2i+1 corresponds to the estimated gamma on the bootstrap network (0.0 if hybrid not found)

  * list of discrepant trees (trees not matching the main tree in net1)

<a id='PhyloNetworks.summarizeHFdf' href='#PhyloNetworks.summarizeHFdf'>#</a>
**`PhyloNetworks.summarizeHFdf`** &mdash; *Function*.



`summarizeHFdf(HFmat::Matrix)`

function to summarize df output from hybridDetection input: HFdf (see hybridDetection) returns dataframe with one row per hybrid, and 5 columns:

  * hybrid index (order from estimated network, see hybridDetection),

  * number of bootstrap trees that match the underlying tree of estimated network, - number of bootstrap networks that have the hybrid

  * mean estimated gamma in the bootstrap networks that have the hybrid

  * sd estimated gamma in the bootstrap networks that have the hybrid also

last row has index -1, and the third column has the number of networks that have all hybrids (hybrid index, mean gamma, sd gamma are meaningless in this last row)

<a id='PhyloNetworks.hybridBootstrapSupport' href='#PhyloNetworks.hybridBootstrapSupport'>#</a>
**`PhyloNetworks.hybridBootstrapSupport`** &mdash; *Function*.



`hybridBootstrapSupport(boot_net::Vector{HybridNetwork}, ref_net::HybridNetwork; rooted=false)`

Match hybrid nodes in a reference network with those in an array of networks, like bootstrap networks. All networks must be fully resolved, and on the same taxon set. If `rooted=true`, all networks are assumed to have been properly rooted beforehand. Otherwise, the origin of each hybrid edge is considered as an unrooted bipartition (default).

Two hybrid edges in two networks are said to match if they share the same "hybrid" clade (or recipient) and the same "donor clade", which is a sister to the hybrid clade in the network. Since a hybrid clade has 2 parent edges, it is sister to two clades simultaneously: one is its major sister (following the major hybrid edge with γ>0.5) and one is its minor sister (following the major hybrid edge with γ<0.5).

To calculate these hybrid and sister clades at a given hybrid node, all other hybrid edges are first removed from the network. Then, the hybrid clade is the hardwired cluster (descendants) of either hybrid edge and major/minor clade is the hardwired cluster of the sibling edge of the major/minor hybrid parent. If `rooted=false`, sister clades are considered as bipartitions.

Output:

1. a "node" data frame with one row per clade and 9 columns giving:

  * **clade**: the clade's name, like the taxon name (if a hybrid is a single taxon) or      the hybrid tag (like 'H1') in the reference network
  * **node**: the node number in the reference network. NA if the clade is not in this network.
  * **hybridnode**: typically the same node number as above, except for hybrid clades in the      reference network. For those, the hybrid node number is listed here.
  * **edge**: number of the parent edge, parent to the node in column 2,      if found in the ref network. NA otherwise.
  * **BS_hybrid**: percentage of bootstrap networks in which the clade is found to be a hybrid clade.
  * **BS_sister**: percentage of bootstrap networks in which the clade is found to be sister to      some hybrid clade (sum of the next 2 columns)
  * **BS_major_sister**: percentage of bootstrap networks in which the clade is found to be the      major sister to some hybrid clade
  * **BS_minor_sister**: same as 7, but minor
  * **BS_hybrid_samesisters**: percentage of bootstrap networks in which the clade is found to be      a hybrid and with the same set of sister clades as in the reference network.      Applies to hybrid clades found in the reference network only, NA for all other clades.

1. an "edge" data frame with one row for each pair of clades, and 8 columns:

  * **edge**: hybrid edge number, if the edge appears in the reference network. NA otherwise.
  * **hybrid_clade**: name of the clade found to be a hybrid, descendent of 'edge'
  * **hybrid**: node number of that clade, if it appears in the reference network. NA otherwise.
  * **sister_clade**: name of the clade that is sister to 'edge', i.e. be sister to a hybrid
  * **sister**: node number of that clade, if in the ref network.
  * **BS_hybrid_edge**: percentage of bootstrap networks in which 'edge' is found to be a hybrid      edge, i.e. when the clade in the 'hybrid' column is found to be a hybrid and the clade in      the 'sister' column is one of its sisters.
  * **BS_major**: percentage of bootstrap networks in which 'edge' is found to be a major hybrid      edge, i.e. when 'hybrid' is found to be a hybrid clade and 'sister' is found to be its      major sister.
  * **BS_minor**: same as 7, but minor

1. a "clade" data frame to describe the make up of all clades found as hybrids or sisters,   starting with a column `taxa` that lists all taxa. All other columns correspond to a given   clade and contain true/false values. `true` means that a given taxon belongs in a given clade.   For a clade named `H1`, for instance, and if the data frame was named 'cla', the   list of taxa in this clade can be obtained with `cla[:taxa][cla[:H1]]`.

1. an array of gamma values, with one row for each bootstrap network and two columns (major/minor) for each hybrid   edge in the reference network. If this hybrid edge was found in the bootstrap network   (i.e. same hybrid and sister clades, after removal of all other hybrid nodes),   its bootstrap gamma value is recorded here. Otherwise, the gamma entry is 0.0.

1. an vector with the number of each hybrid edge in the reference network, in the same order   as for the columns in the array of gamma values above.

<a id='PhyloNetworks.bootsnaq' href='#PhyloNetworks.bootsnaq'>#</a>
**`PhyloNetworks.bootsnaq`** &mdash; *Function*.



```
bootsnaq(T::HybridNetwork, df::DataFrame)
bootsnaq(T::HybridNetwork, vector of tree lists)
```

Bootstrap analysis for SNaQ. Bootstrap data can be quartet concordance factors (CF), drawn from sampling uniformly in their credibility intervals, as given in the data frame `df`. Alternatively, bootstrap data can be gene trees sampled from a vector of tree lists: one list of bootstrap trees per locus (see `readBootstrapTrees` to generate this, from a file containing a list of bootstrap files: one per locus).

From each bootstrap replicate, a network is estimated with snaq!, with a search starting from topology `T`. Optional arguments include the following, with default values in parentheses:

  * hmax (1): max number of reticulations in the estimated networks
  * nrep (10): number of bootstrap replicates.
  * runs (10): number of independent optimization runs for each replicate
  * filename (bootsnaq): root name for output files
  * seed (0 to get a random seed from the clock): seed for random number generator

<a id='PhyloNetworks.readBootstrapTrees' href='#PhyloNetworks.readBootstrapTrees'>#</a>
**`PhyloNetworks.readBootstrapTrees`** &mdash; *Function*.



```
readBootstrapTrees(filename)
```

input: name of file containing the path/name to multiple bootstrap files, one per line. Each bootstrap file corresponds to bootstrap trees from a single gene.

output: vector of vectors of trees.

<a id='PhyloNetworks.readMultiTopology' href='#PhyloNetworks.readMultiTopology'>#</a>
**`PhyloNetworks.readMultiTopology`** &mdash; *Function*.



`readMultiTopology(file)`

Read a text file with a list of networks in parenthetical format (one per line). Crash if a network is broken over several lines. Return an array of HybridNetwork object.

<a id='PhyloNetworks.hybridatnode!' href='#PhyloNetworks.hybridatnode!'>#</a>
**`PhyloNetworks.hybridatnode!`** &mdash; *Function*.



`hybridatnode!(net::HybridNetwork, nodeNumber::Int64)`

Changes the hybrid in a cycle to the node with number `nodeNumber`. This node must be in one (and only one) cycle, otherwise an error will be thrown.

**Example #"**

```julia
julia> net = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
julia> plot(net, showNodeNumber=true)
julia> hybridatnode!(net, -4)
julia> plot(net)
```


<a id='Index-1'></a>

## Index

- [`PhyloNetworks.DataCF`](index.md#PhyloNetworks.DataCF)
- [`PhyloNetworks.HybridNetwork`](index.md#PhyloNetworks.HybridNetwork)
- [`PhyloNetworks.Quartet`](index.md#PhyloNetworks.Quartet)
- [`Gadfly.plot`](index.md#Gadfly.plot)
- [`PhyloNetworks.bootsnaq`](index.md#PhyloNetworks.bootsnaq)
- [`PhyloNetworks.cladewiseorder!`](index.md#PhyloNetworks.cladewiseorder!)
- [`PhyloNetworks.deleteHybridThreshold!`](index.md#PhyloNetworks.deleteHybridThreshold!)
- [`PhyloNetworks.deleteleaf!`](index.md#PhyloNetworks.deleteleaf!)
- [`PhyloNetworks.fittedQuartetCF`](index.md#PhyloNetworks.fittedQuartetCF)
- [`PhyloNetworks.directEdges!`](index.md#PhyloNetworks.directEdges!)
- [`PhyloNetworks.displayedNetworkAt!`](index.md#PhyloNetworks.displayedNetworkAt!)
- [`PhyloNetworks.displayedTrees`](index.md#PhyloNetworks.displayedTrees)
- [`PhyloNetworks.hardwiredCluster`](index.md#PhyloNetworks.hardwiredCluster)
- [`PhyloNetworks.hardwiredClusterDistance`](index.md#PhyloNetworks.hardwiredClusterDistance)
- [`PhyloNetworks.hardwiredClusters`](index.md#PhyloNetworks.hardwiredClusters)
- [`PhyloNetworks.hybridBootstrapSupport`](index.md#PhyloNetworks.hybridBootstrapSupport)
- [`PhyloNetworks.hybridDetection`](index.md#PhyloNetworks.hybridDetection)
- [`PhyloNetworks.hybridatnode!`](index.md#PhyloNetworks.hybridatnode!)
- [`PhyloNetworks.majorTree`](index.md#PhyloNetworks.majorTree)
- [`PhyloNetworks.mapAllelesCFtable`](index.md#PhyloNetworks.mapAllelesCFtable)
- [`PhyloNetworks.minorTreeAt`](index.md#PhyloNetworks.minorTreeAt)
- [`PhyloNetworks.plotNetGraphViz`](index.md#PhyloNetworks.plotNetGraphViz)
- [`PhyloNetworks.preorder!`](index.md#PhyloNetworks.preorder!)
- [`PhyloNetworks.printEdges`](index.md#PhyloNetworks.printEdges)
- [`PhyloNetworks.printNodes`](index.md#PhyloNetworks.printNodes)
- [`PhyloNetworks.readBootstrapTrees`](index.md#PhyloNetworks.readBootstrapTrees)
- [`PhyloNetworks.readInputTrees`](index.md#PhyloNetworks.readInputTrees)
- [`PhyloNetworks.readMultiTopology`](index.md#PhyloNetworks.readMultiTopology)
- [`PhyloNetworks.readSnaqNetwork`](index.md#PhyloNetworks.readSnaqNetwork)
- [`PhyloNetworks.readTableCF`](index.md#PhyloNetworks.readTableCF)
- [`PhyloNetworks.readTopology`](index.md#PhyloNetworks.readTopology)
- [`PhyloNetworks.readTopologyLevel1`](index.md#PhyloNetworks.readTopologyLevel1)
- [`PhyloNetworks.readTrees2CF`](index.md#PhyloNetworks.readTrees2CF)
- [`PhyloNetworks.rootatnode!`](index.md#PhyloNetworks.rootatnode!)
- [`PhyloNetworks.rootonedge!`](index.md#PhyloNetworks.rootonedge!)
- [`PhyloNetworks.setGamma!`](index.md#PhyloNetworks.setGamma!)
- [`PhyloNetworks.setLength!`](index.md#PhyloNetworks.setLength!)
- [`PhyloNetworks.snaq!`](index.md#PhyloNetworks.snaq!)
- [`PhyloNetworks.snaqDebug`](index.md#PhyloNetworks.snaqDebug)
- [`PhyloNetworks.summarizeDataCF`](index.md#PhyloNetworks.summarizeDataCF)
- [`PhyloNetworks.summarizeHFdf`](index.md#PhyloNetworks.summarizeHFdf)
- [`PhyloNetworks.tipLabels`](index.md#PhyloNetworks.tipLabels)
- [`PhyloNetworks.topologyMaxQPseudolik!`](index.md#PhyloNetworks.topologyMaxQPseudolik!)
- [`PhyloNetworks.topologyQPseudolik!`](index.md#PhyloNetworks.topologyQPseudolik!)
- [`PhyloNetworks.treeEdgesBootstrap`](index.md#PhyloNetworks.treeEdgesBootstrap)
- [`PhyloNetworks.writeTopology`](index.md#PhyloNetworks.writeTopology)

