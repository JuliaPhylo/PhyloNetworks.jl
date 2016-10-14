# problems found in readTopologyUpdate when trying to plot with John Spaw's functions
# will test here the problems
# Claudia July 2015

include("../src/types.jl")
include("../src/functions.jl")

# estNetworks form baseline and phylonet
n1 = "(((1,2))#H1,(6,(5,((4,3),#H1))));"
n2 = "((5)#H1,(6,((1,2),(3,(4,#H1)))));"
n3 = "((5,((3,4),((2,1))#H1)),(6,#H1));"
n4 = "(6,(5,((1,(2,(3)#H1)),(4,#H1))));"
n5 = "(6,((((3,4))#H1,(2,1)),(5,#H1)));"
n6 = "((((5,(3,4)))#H1,(1,2)),(6,#H1));"
n7 = "(((5)#H1,((2,1),(4,3))),(#H1,6));"
n8 = "(6,((5,((4,3))#H1),((1,2),#H1)));"
n9 = "(6,(((1,2),((3,4))#H1),(5,#H1)));"
n10 = "(6,((((4,3))#H1,(1,2)),(5,#H1)));"
n11 = "(6,(((2,1),((4,3))#H1),(5,#H1)));"
n12 = "(6,((5)#H1,((2,1),((4,3),#H1))));"
n13 = "(6,((((4,3))#H1,(1,2)),(#H1,5)));"
n14 = "(6,((5,((4,3))#H1),(#H1,(2,1))));"
n15 = "(6,(((2,1),((4,3))#H1),(5,#H1)));"
n16 = "(6,(5,(((4,3))#H1,(1,(#H1,2)))));"
n17 = "(6,(((1,2),((4)#H1,3)),(5,#H1)));"
n18 = "((5,((4,3),(1,(2)#H1))),(6,#H1));"
n19 = "((((5,(3,4)))#H1,6),((2,1),#H1));"
n20 = "(6,(((1,2),((4,3))#H1),(5,#H1)));"
n21 = "(6,((5,((3,4))#H1),((1,2),#H1)));"
n22 = "((5)#H1,(6,((2,1),(#H1,(4,3)))));"
n23 = "(6,(((1,2),((4)#H1,3)),(5,#H1)));"
n24 = "(6,((5,((3,4))#H1),((1,2),#H1)));"
n25 = "(6,((5,((4,3))#H1),((2,1),#H1)));"
n26 = "(((1,2))#H1,(6,(5,((3,4),#H1))));"
n27 = "(6,((5)#H1,((2,1),((4,3),#H1))));"
n28 = "(6,(((2,1),((3,4))#H1),(5,#H1)));"
n29 = "(6,((((4,3))#H1,(2,1)),(5,#H1)));"
n30 = "(6,(((1,2),(3,(4)#H1)),(#H1,5)));"

string = [n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30]

for s in string
    net = readTopologyUpdate(s);
    printEdges(net)
    printNodes(net)
    net.node[net.root].number
    net.node[net.root].leaf
    canBeRoot(net.node[net.root]) || error("root wrongly placed in $(s)")

    net2 = readTopologyUpdate(s,true);
    printEdges(net2)
    printNodes(net2)
    net2.node[net2.root].number
    net2.node[net2.root].leaf
    canBeRoot(net2.node[net2.root])  || error("root wrongly placed in $(s)")
end

println("NO ERRORS!")

## # --------------
## n31 = "(6,((5,(((1,(2,(3)#H7:9.380388137723473):1.215558354715711):1.533022718529031,#H7:0.8208307297991476):0.0,(4)#H9:-0.0):10.0):10.0,#H9:0.0):0.02572743545235116);"
## net = readTopologyUpdate(n31)
## printEdges(net)
## printNodes(net)
## net.node[net.root].number
## net.node[net.root].leaf
## canBeRoot(net.node[net.root]) || error("root wrongly placed in $(s)")
