# small tutorial of main functions in PhyloNetworks
Pkg.update() #check github website before
using PhyloNetworks #it can take a while to precompile everything

# read network in parenthetical format
# note: the first time a function is called, it takes longer because it compiles
net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
# net is a julia object of type "HybridNetwork", it can be a tree

# plot a network
plot(net)

# save a network plot as PDF
p=plot(net)
using Gadfly
draw(PDF("plot.pdf",4inch,4inch),p)

# write network in parenthetical format
writeTopology(net)

# root the network in a different edge
# first, identify the edge number where we want to root:
plot(net, showEdgeNumber=true, showEdgeLength=false)
# we want to root on the edge 1 (the one that leads to A):
rootonedge!(net,1)
# plot again to see root at A
plot(net)

# root the network in a different node
# first, identify the node number
plot(net, showNodeNumber=true, showEdgeLength=false)
# we want to root at node -3
rootatnode!(net,-3)
# plot again to see new root
plot(net)

# estimate network
# read in table of CF from TICR:
d = readTableCF("tableCFCI.csv")
# read in starting topology, use readTopologyLevel1
T = readTopologyLevel1("startTree.txt")
# estimate best network: hmax=1, runs=10
estNet1 = snaq!(T,d)
plot(estNet1)
# we can re-root at an outgroup
plot(estNet1, showEdgeNumber=true, showEdgeLength=false)
rootonedge!(estNet1,1)
plot(estNet1)
writeTopology(estNet1) #will only print branch lengths for internal edges

# read best network, and one per run, first is the best network overall
nets = readMultiTopology("snaq.out") #will change the name of the function later to readMultiTopologyLevel1
length(nets)
plot(nets[2]) #different to the best network, worse pseudolik


# estimate bigger network
estNet2 = snaq!(estNet1,d,hmax=2,runs=3)
plot(estNet2) #same network! recall the level-1 restriction

# bootstrap analysis
using DataFrames
df = readtable("tableCFCI.csv")
bootNet = bootsnaq(T,df,hmax=1,nrep=10,bestNet=estNet1, runs=1)
length(bootNet)
plot(bootNet[1])

# bootstrap support of tree edges:
df0,tree0 = treeEdgesBootstrap(bootNet,estNet1)
# df0 has the bootstrap support of tree edges by edge number, to see which edge is, plot:
# tree0 has the main underlying tree, plotted with bs support
plot(tree0, showEdgeLength=false, edgeLabel=df0)
plot(estNet1, showEdgeLength=false, edgeLabel=df0)
plot(estNet1, showEdgeLength=false, showGamma=false, edgeLabel=df0[df0[:proportion] .< 1.0, :])


f, fr, fd, fs, clade, gam, edgenum = hybridBootstrapSupport(bootNet, estNet1);
## - edgenum lists all the reticulation events, showing the number of the minor hybrid edge for each.
## - fs lists the bootstrap support (proportions) of recipient clades
## - clade lists recipient and donor clades
## - fd lists the bootstrap support (proportions) of donor clade
## - fs lists the bootstrap support (proportions) of sister clades
## - f lists the proportion of bootstrap networks in which a minor
##  hybrid matches that in the best network, in terms of both the
##  recipient and donor clades (for introgression), or in terms of
##  both the recipient clade and the set of parental clades {donor,
##  sibling} (for hybridization).
## - gam lists for each minor hybrid edge in a bootstrap network
##  matching a hybrid edge in the best network, its inheritance (γ)
##  value was extracted, γ=0 values are for bootstrap replicates that
##  did not have a match with the hybridization in the best network.
