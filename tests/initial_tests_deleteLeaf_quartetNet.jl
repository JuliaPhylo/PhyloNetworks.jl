# Initial tests to see if deleteLeaf works
# as expected with QuartetNetwork type
# Claudia November 2014

include("../bad_triangle_example.jl")
include("../case_f_example.jl")
qnet = QuartetNetwork(net);
printEdges(qnet)
printEdges(net)
qnet.hasEdge
updateHasEdge!(qnet,net);

# bad triangle
deleteLeaf!(qnet,qnet.node[8])
deleteLeaf!(qnet,qnet.node[7])
deleteLeaf!(qnet,qnet.node[5])
deleteLeaf!(qnet,qnet.node[6])
deleteLeaf!(qnet,qnet.node[1])

# bad diamond
deleteLeaf!(qnet,qnet.node[10])
deleteLeaf!(qnet,qnet.node[8])


# test extractQuartet

include("../bad_triangle_example.jl")
include("../case_f_example.jl")
qnet = extractQuartet(net,net.node[8],net.node[7],net.node[6],net.node[1]);
printEdges(qnet)
printEdges(net)
qnet.hasEdge
updateHasEdge!(qnet,net);

# bad triangle
deleteLeaf!(qnet,qnet.node[8])
deleteLeaf!(qnet,qnet.node[7])
deleteLeaf!(qnet,qnet.node[5])
deleteLeaf!(qnet,qnet.node[6])
deleteLeaf!(qnet,qnet.node[1])

# bad diamond
deleteLeaf!(qnet,qnet.node[10])
deleteLeaf!(qnet,qnet.node[8])
