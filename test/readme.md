### Tests functions
All in runtests.jl

old test files: checkout PhyloNetworks v0.9.1 or older to see those older files.

#### add hybridization

test_add2hyb.jl add a first hybrid, then a second hybrid that
makes a bad triangle, and the
functions should identify it

#### delete hybridization

test_deleteHybridizationUpdate.jl checks that all attributes are correctly updated after deleting a hybridization (1 and 2 hybrids)

#### read topology

test_5taxon_readTopology.jl runs all the tests for the eight 5taxon
networks by reading from parenthetical format and then updating

aux functions: test_functions_5taxon_read.jl

#### calculate exp CF

test_calculateExpCF.jl runs all the tests for Case G, bad diamond and
bad triangle for calculation of expCF

test_calculateExpCF2.jl computes the expCF for the n6 network

#### has Edge

test_hasEdge.jl tests if the attribute qnet.hasEdge is correctly
updated after extracting quartets for case G. It also checks if
net.ht, net.numht, qnet.indexht are correctly set for Case G.

#### parts of optBL

test_optBLparts.jl tests the parts of optBL separately to see if they
work.

#### parameters

test_parameters.jl get net.ht and net.numht for all the 5 taxon networks

#### Likelihood

test_correctLik.jl computes the pseudolik for a tree and a network with 1 hybrid and checks that it is correctly computed

#### partition

test_partition.jl (1 hybrid) and test_partition2.jl (2 hybrids) check if the attribute of partition is correctly set

#### network manipulations, comparisons, and plots

test_orderings_plot.jl tests functions in manipulateNet.jl and plotsGadfly.jl
- directEdges! rootatnode! rootonedge!
- preorder! cladewiseorder!
- plot rotate!

test_compareNetworks.jl tests:
 - deleteHybridEdge!
 - functions to extract displayed trees/subnetworks,
 - hardwiredClusterDistance (which uses functions above)

both run a few examplar tests by default (used by runtests.jl)
but can run more tests if one first defines doalltests = true
