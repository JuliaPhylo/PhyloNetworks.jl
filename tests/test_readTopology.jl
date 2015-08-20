# initial tests for readTopology
# reads only a tree and with taxon number, not taxon names
# Claudia October 2014
#########################################################

include("../src/types.jl")
include("../src/functions.jl")

using Base.Collections # for updateInCycle with priority queue

# good trees ---------------------------
f = open("prueba_tree.txt","w")
tree = "((1,2),(3,4));"
tree = "((11,22),(33,44));"
tree = "((Ant,Bear),(Cat,Dog));"
tree = "((Ant1,Bear2),(Cat3,Dog4));"
tree = "((1Ant,2Bear),(3Cat,4Dog));"
tree = "((1,2),3,4);"
tree = "(1,2,(3,4));"
tree = "(Ant,Bear,(Cat,Dog));"
tree = "(A,B,(C,D));"
tree = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
tree = "((((((1,2),3),4),5),6),7,8);" # yeast data tree0
write(f,tree)
close(f)

net = readTopology("prueba_tree.txt");
printEdges(net)
printNodes(net)
net.names

cleanAfterRead!(net)

# bad trees --------------------------
f = open("prueba_tree.txt","w")
tree = "(((1,2),(3,4)));" # extra parenthesis
tree = "((1,2,(3,4)));" # extra parenthesis
tree = "(((1,2),3,4));" # extra parenthesis
tree = "((1),2,3,4,5);" # not tree
tree = "((1,2),3,4,5);" # polytomy
tree = "((1,*),(3,4));" # not letter/number taxon name
tree = "(1,2::,(3,4));" #double :
tree = "(1,2:::,(3,4));" #triple :
tree = "(1,2:0.2:,(3,4));" #no 2nd :
tree = "(1,2:0.2:0.2:,(3,4));" #no 3rd :
tree = "(1,2:,(3,4));" #no 1st:
tree = "((1,2),(3,4))" # no ;
write(f,tree)
close(f)

net = readTopology("prueba_tree.txt")
printEdges(net)
printNodes(net)
net.names

cleanAfterRead!(net)


# trees with additional info -----------------
f = open("prueba_tree.txt","w")
tree = "((1:1.2,2:0.3),3:1.8,4);"
tree = "((1:1.2,2:0.3):1.2:0.6:0.3,3:1.8,4);"
tree = "((1:1.2,2:0.3):1.2::0.6,3:1.8,4);" # no bootstrap
tree = "((1:1.2,2:0.3):::0.6,3:1.8,4);" # no length nor bootstrap
tree = "((1:1.2,2:0.3):1.2:2.5:1.6,3:1.8,4);" # wrong value gamma
tree = "((1:1.2,2:0.3):1.2:2.5:,3:1.8,4);" # missing value gamma
write(f,tree)
close(f)

net = readTopology("prueba_tree.txt")
printEdges(net)
printNodes(net)
net.names

# beginning of networks ---------------------------
include("../src/types.jl")
include("../src/functions.jl")

using Base.Collections # for updateInCycle with priority queue


f = open("prueba_tree.txt","w")
tree = "(((3,4)Z#H1,1),(Z#H1,2));" # expand child, no gammas
tree = "((Z#H1,2),((3,4)Z#H1,1));" # expand child, no gammas, leaf read first
tree = "(((3,4)Z#H1:::0.9,1),(Z#H1:::0.2,2));" # gammas do not sum up to 1
tree = "(((3,4)Z#H1:5.0::0.9,1),(Z#H1:::0.2,2));"
tree = "(1,2,(3,4)A);"
tree = "(1,2,(3,4)A:0.8);"
tree = "((1,2),#H1,4);" #only one hybrid
tree = "((1,2),#H1,(3,#H2));" #two hybrids but different
tree = "((1,2)#H1,(3,4)#H1);" #error: both hybrid have childre
tree = "(((1,#H1),2),#H1,3);" # both hybrid leaves
tree = "(((1,4)#H1,2),#H1,(5,6)#H1);" # hybrid polytomy
write(f,tree)
close(f)

net = readTopology("prueba_tree.txt");
printEdges(net)
printNodes(net)
net.names
net.numHybrids

updateAllReadTopology!(net)

readTopologyUpdate("prueba_tree.txt")

#======================================================

include("../src/types.jl")
include("../src/functions.jl")

using Base.Collections # for updateInCycle with priority queue

# good trees ---------------------------
f = open("prueba_tree.txt","w")
tree = "((1,2),(3,4));"
tree = "((11,22),(33,44));"
tree = "((Ant,Bear),(Cat,Dog));"
tree = "((Ant1,Bear2),(Cat3,Dog4));"
tree = "((1Ant,2Bear),(3Cat,4Dog));"
tree = "((1,2),3,4);"
tree = "(1,2,(3,4));"
tree = "(Ant,Bear,(Cat,Dog));"
tree = "(A,B,(C,D));"
tree = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
tree = "((((((1,2),3),4),5),6),7,8);" # yeast data tree0
tree = "(Scer,((Smik,(Skud,Sbay)),Spar));" # yeast data astral output
tree = "(1,2,3,4,5);"
write(f,tree)
close(f)

net = readTopologyUpdate("prueba_tree.txt");
printEdges(net)
printNodes(net)
net.names

net2 = readTopologyUpdate(tree);
printEdges(net2)
printNodes(net2)
net2.names
