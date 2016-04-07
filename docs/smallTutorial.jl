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

# read best network, and one per run
nets = readInputTrees("snaq.out") #will change the name of the function later to readMultiTopologyLevel1
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
# tree0 has the main underlying tree
plot(tree0, showEdgeNumber=true, showEdgeLength=false)

outgroup = "6" #to root all networks, should be compatible with hybrid edge directions
# count the bootstrap support for each hybrid in net0:
HFmat,discTrees = hybridDetection(bootNet,estNet1,outgroup)
# HFmat will have one row per bootstrap network, and number of columns depending on number of hybrids in estNet1
# if estNet1 had 3 hybrids for example, HFmat will have 6 columns:
# - first 3 columns indicate the presence (1.0) or absense (0.0) of each hybrid (column) for each bootstrap network (row)
# - the last 3 columns indicate the estimated gamma in the bootstrap network if the hybrid was found (0.0 if not found)
HFmat
# hybrid detection only makes sense if the underlying trees are the same, so discTrees has the list of trees that do not match the underlying tree in estNet1
length(discTrees)

# dfhyb has one row per hybrid, and 5 columns:
# - hybrid index
# - number of trees that match the underlying tree in estNet1 (same for all hybrids)
# - number of networks with that hybrid
# - mean estimated gamma among networks with the hybrid
# - sd estimated gamma among networks with the hybrid
# ** last row contains in 3er column the number of networks that have all
# same hybrids as estNet1 (hybrid index, mean gamma and sd gamma are meaningless)
dfhyb = summarizeHFdf(HFmat)
