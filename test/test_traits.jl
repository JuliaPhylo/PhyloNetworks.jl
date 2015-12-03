 # test for trait evolution
 # Claudia November 2015

include("../src/types.jl")
include("../src/functions.jl")
include("../src/traits.jl")

tree = "(A,(B,(C,D)));"
tree= "(A,((B,#H1),(C,(D)#H1)));" #fails
net=readTopology(tree)
net=readTopologyLevel1(tree)
printEdges(net)
#[e.isChild1 for e in net.edge]
directEdges!(net)
#[e.isChild1 for e in net.edge]
postorder!(net)

sharedPathMatrix(net)
