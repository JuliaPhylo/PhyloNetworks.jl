# test functions for the 5taxon networks
# functions used in tests_5taxon_readTopology.jl
# fixit: it would be better that within each function, it would look
# for the needed edges/nodes starting from the hybrid node
# so that it would be a general function (not bounded to the specific indeces)
# Claudia November 2014

# Case C bad triangle II
function testCaseC(net::HybridNetwork)
    net.numHybrids != 1 ? error("network should have one hybrid node and it has $(net.numHybrids)") : nothing
    node=net.hybrid[1];
    node.k != 3 ? error("k diff than 3") : nothing
    (net.edge[5].inCycle != node.number || net.edge[6].inCycle != node.number || net.edge[7].inCycle != node.number) ? error("edges not incycle") : nothing
    (net.node[6].inCycle != node.number || net.node[7].inCycle != node.number || net.node[5].inCycle != node.number) ? error("nodes not incycle") : nothing
    net.edge[4].containRoot ? error("contain root wrong") : nothing
    (net.node[6].gammaz != exp(-1)*(1-0.1*(1-exp(-1))) || net.node[7].gammaz != exp(-1)*(1-0.9*(1-exp(-1))) || net.node[5].gammaz != exp(-2)*0.1*(1-exp(-1))*(1-exp(-1))*(-0.9)) ? error("gammaz updated wrong") : nothing
    all([!e.istIdentifiable for e in net.edge]) ? nothing : error("all edges should be not identifiable and some are not")
    node.isBadTriangleI ? error("thinks it is bad triangle I") : nothing
    !node.isBadTriangleII ? error("does not know it is bad triangle II") : nothing
    (net.edge[5].hybrid && net.edge[7].hybrid) ? nothing : error("hybrid edges not recognized as hybrids")
    [n.number for n in net.leaf] == [1,2,3,5,6] ? nothing : error("net.leaf is wrong")
end

# Case F bad diamond I
function testCaseF(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    net.numHybrids == 1 ? nothing : error("networks should have one hybrid node and it has $(net.numHybrids)")
    node=net.hybrid[1];
    node.k != 4 ? error("k diff than 4") : nothing
    (net.edge[3].inCycle != node.number || net.edge[4].inCycle != node.number || net.edge[5].inCycle != node.number || net.edge[7].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[3].inCycle  != node.number || net.node[4].inCycle  != node.number || net.node[6].inCycle  != node.number || net.node[7].inCycle  != node.number) ? error("nodes not correctly in cycle") : nothing
    !node.isBadDiamondI ? error("does not know it is bad diamond I") : nothing
    node.isBadDiamondII ? error("thinks it is bad diamond II") : nothing
    net.edge[2].containRoot ? error("edge can contain root") : nothing
    (!net.edge[3].hybrid || !net.edge[3].isMajor) ? error("edge 3 is not hybrid or major") : nothing
    (!net.edge[5].hybrid || net.edge[5].isMajor) ? error("edge 5 is not hybrid or is major") : nothing
    net.node[4].gammaz != net.edge[3].gamma*net.edge[4].z ? error("node 4 gammaz not correctly calculated") : nothing
    net.node[6].gammaz != net.edge[5].gamma*net.edge[7].z ? error("node 6 gammaz not correctly calculated") : nothing
    (net.edge[4].istIdentifiable || net.edge[7].istIdentifiable) ? error("edges identifiable that should not") : nothing
    !net.edge[8].istIdentifiable ? error("edge 8 not identifiable") : nothing
    net.visited[8] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    [n.number for n in net.leaf] == [1,2,5,7,8] ? nothing : error("net.leaf is wrong")
end

# Case G
function testCaseG(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    net.numHybrids == 1 ? nothing : error("networks should have one hybrid node and it has $(net.numHybrids)")
    node=net.hybrid[1];
    node.k != 4 ? error("k diff than 4") : nothing
    (net.edge[5].inCycle != node.number || net.edge[6].inCycle != node.number || net.edge[7].inCycle != node.number || net.edge[9].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[5].inCycle  != node.number || net.node[6].inCycle  != node.number || net.node[8].inCycle  != node.number || net.node[9].inCycle  != node.number) ? error("nodes not correctly in cycle") : nothing
    (node.isBadDiamondI || node.isBadDiamondII )? error("thinks it is bad diamond") : nothing
    net.edge[4].containRoot ? error("edge can contain root") : nothing
    (!net.edge[5].hybrid || !net.edge[5].isMajor) ? error("edge 5 is not hybrid or major") : nothing
    (!net.edge[7].hybrid || net.edge[7].isMajor) ? error("edge 7 is not hybrid or is major") : nothing
    net.node[5].gammaz != -1 ? error("node 5 gammaz should be -1") : nothing
    (!net.edge[6].istIdentifiable || !net.edge[3].istIdentifiable || !net.edge[9].istIdentifiable) ? error("edges identifiable as not identifiable") : nothing
    net.visited[6] = false;
    net.visited[3] = false;
    net.visited[9] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    [n.number for n in net.leaf] == [1,2,4,7,8] ? nothing : error("net.leaf is wrong")
end

# Case H
function testCaseH(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    net.numHybrids == 1 ? nothing : error("networks should have one hybrid node and it has $(net.numHybrids)")
    node=net.hybrid[1];
    node.k != 4 ? error("k diff than 4") : nothing
    (net.edge[4].inCycle != node.number || net.edge[5].inCycle != node.number || net.edge[7].inCycle != node.number || net.edge[9].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[4].inCycle  != node.number || net.node[6].inCycle  != node.number || net.node[8].inCycle  != node.number || net.node[10].inCycle  != node.number) ? error("nodes not correctly in cycle") : nothing
    (node.isBadDiamondI || node.isBadDiamondII) ? error("thinks it is bad diamond") : nothing
    net.edge[8].containRoot ? error("8 can contain root") : nothing
    (!net.edge[9].hybrid || !net.edge[9].isMajor) ? error("edge 9 is not hybrid or major") : nothing
    (!net.edge[4].hybrid || net.edge[4].isMajor) ? error("edge 4 is not hybrid or is major") : nothing
    node.gammaz != -1 ? error("hybrid node gammaz should be -1") : nothing
    (!net.edge[3].istIdentifiable || !net.edge[5].istIdentifiable || !net.edge[7].istIdentifiable) ? error("edge9,5,13not identifiable") : nothing
    net.visited[3] = false;
    net.visited[5] = false;
    net.visited[7] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    [n.number for n in net.leaf] == [1,2,4,5,6] ? nothing : error("net.leaf is wrong")
end


# Case J
function testCaseJ(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    net.numHybrids == 1 ? nothing : error("networks should have one hybrid node and it has $(net.numHybrids)")
    node=net.hybrid[1];
    node.k != 5 ? error("k diff than 5") : nothing
    (net.edge[2].inCycle != node.number || net.edge[4].inCycle != node.number || net.edge[6].inCycle != node.number || net.edge[10].inCycle != node.number || net.edge[8].inCycle != node.number) ? error("edges not correctly in cycle") : nothing
    (net.node[2].inCycle  != node.number || net.node[4].inCycle  != node.number || net.node[6].inCycle  != node.number || net.node[10].inCycle  != node.number || net.node[9].inCycle) != node.number ? error("nodes not correctly in cycle") : nothing
    net.edge[1].containRoot ? error("edge can contain root") : nothing
    (!net.edge[2].hybrid || !net.edge[2].isMajor) ? error("edge 2 is not hybrid or major") : nothing
    node.gammaz != -1 ? error("hybrid node gammaz should be -1") : nothing
    (!net.edge[4].istIdentifiable || !net.edge[6].istIdentifiable || !net.edge[10].istIdentifiable) ? error("edge9,5,10not identifiable") : nothing
    net.visited[4] = false;
    net.visited[6] = false;
    net.visited[10] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    [n.number for n in net.leaf] == [1,3,4,5,6] ? nothing : error("net.leaf is wrong")
end

# Case D bad triangle I
function testCaseD(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    net.numHybrids == 1 ? nothing : error("networks should have one hybrid node and it has $(net.numHybrids)")
    node = net.hybrid[1];
    node.k != 3 ? error("k diff than 3") : nothing
    (net.edge[4].inCycle != node.number || net.edge[5].inCycle != node.number || net.edge[7].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[4].inCycle  != node.number || net.node[7].inCycle  != node.number || net.node[6].inCycle  != node.number) ? error("nodes not correctly in cycle") : nothing
    !node.isBadTriangleI ? error("does not know it is bad triangle I") : nothing
    node.isBadTriangleII ? error("thinks it is bad triangle II") : nothing
    (net.edge[3].containRoot || net.edge[1].containRoot || net.edge[2].containRoot)  ? error("edges can contain root") : nothing
    (!net.edge[4].hybrid || !net.edge[4].isMajor) ? error("edge 4 is not hybrid or major") : nothing
    (!net.edge[5].hybrid || net.edge[5].isMajor) ? error("edge 5 is not hybrid or is major") : nothing
    node.gammaz != net.edge[5].gamma*net.edge[5].gamma*net.edge[7].z+(net.edge[4].gamma)*(net.edge[4].gamma)*net.edge[4].z ? error("hybrid node gammaz not correctly calculated") : nothing
    (net.edge[3].length != 0.0 || net.edge[3].istIdentifiable) ? error("edge 3 not correctly non identifiable length 0") : nothing
    (net.edge[5].length != 0.0 || net.edge[5].istIdentifiable) ? error("edge 5 not correctly non identifiable length 0") : nothing
    net.node[6].gammaz != net.edge[5].gamma*net.edge[7].z ? error("node 6 gammaz not correctly calculated") : nothing
    !net.edge[8].istIdentifiable ? error("edge 8 not identifiable") : nothing
    net.visited[8] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    [n.number for n in net.leaf] == [1,2,4,5,6] ? nothing : error("net.leaf is wrong")
end

# Case E bad triangle I
function testCaseE(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    net.numHybrids == 1 ? nothing : error("networks should have one hybrid node and it has $(net.numHybrids)")
    node = net.hybrid[1];
    node.k != 3 ? error("k diff than 3") : nothing
    (net.edge[4].inCycle != node.number || net.edge[6].inCycle != node.number || net.edge[7].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[4].inCycle  != node.number || net.node[7].inCycle  != node.number || net.node[6].inCycle  != node.number) ? error("nodes not correctly in cycle") : nothing
    !node.isBadTriangleI ? error("does not know it is bad triangle I") : nothing
    node.isBadTriangleII ? error("thinks it is bad triangle II") : nothing
    (net.edge[3].containRoot || net.edge[1].containRoot || net.edge[2].containRoot)  ? error("edges can contain root") : nothing
    (!net.edge[4].hybrid || !net.edge[4].isMajor) ? error("edge 4 is not hybrid or major") : nothing
    (!net.edge[7].hybrid || net.edge[7].isMajor) ? error("edge 7 is not hybrid or is major") : nothing
    node.gammaz != net.edge[4].gamma*net.edge[4].gamma*net.edge[6].z+(1-net.edge[4].gamma)*(1-net.edge[4].gamma)*net.edge[7].z ? error("hybrid node gammaz not correctly calculated") : nothing
    (net.edge[3].length != 0.0 || net.edge[3].istIdentifiable) ? error("edge 3 not correctly non identifiable length 0") : nothing
    (net.edge[4].length != 0.0 || net.edge[4].istIdentifiable) ? error("edge 7 not correctly non identifiable length 0") : nothing
    net.node[6].gammaz != net.edge[4].gamma*net.edge[6].z ? error("node 6 gammaz not correctly calculated") : nothing
    !net.edge[8].istIdentifiable ? error("edge 8 not identifiable") : nothing
    net.visited[8] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    [n.number for n in net.leaf] == [1,2,4,5,6] ? nothing : error("net.leaf is wrong")
end


# Case I bad diamond II
function testCaseI(net::HybridNetwork)
    net.numHybrids == 1 ? nothing : error("networks should have one hybrid node and it has $(net.numHybrids)")
    node = net.hybrid[1];
    node.isBadDiamondII ? nothing : error("does not know it is bad diamond II")
    node.isBadDiamondI ? error("thinks it is bad diamond I") : nothing
    node.k == 4 ? nothing : error("k should be 4")
    net.visited = [e.istIdentifiable for e in net.edge];
    edge4 = getIndexEdge(4,net);
    edge1 = getIndexEdge(1,net);
    edge2 = getIndexEdge(2,net);
    edge3 = getIndexEdge(3,net);
    edge9 = getIndexEdge(9,net);
    edge10 = getIndexEdge(10,net);
    edge6 = getIndexEdge(6,net);
    node1 = getIndexNode(-1,net);
    node2 = getIndexNode(-2,net);
    node5 = getIndexNode(-5,net);
    node3 = getIndexNode(3,net);
    (net.edge[edge4].inCycle != node.number || net.edge[edge9].inCycle != node.number || net.edge[edge6].inCycle != node.number || net.edge[edge9].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[node1].inCycle  != node.number || net.node[node2].inCycle  != node.number || net.node[node5].inCycle  != node.number || net.node[node3].inCycle  != node.number) ? error("nodes 1,5,11,12 not correctly in cycle") : nothing
    (net.edge[edge1].containRoot || net.edge[edge2].containRoot || net.edge[edge3].containRoot) ? error("edges can contain root and shouldn't") : nothing
    (!net.edge[edge4].hybrid || !net.edge[edge4].isMajor) ? error("edge 4 is not hybrid or major") : nothing
    (net.edge[edge4].length != 0 || net.edge[edge9].length != 0) ? error("edges should have length 0") : nothing
    (net.edge[edge9].istIdentifiable || net.edge[edge4].istIdentifiable) ? error("edge9,4 identifiable and should not") : nothing
    net.visited[edge3] = false;
    net.visited[edge6] = false;
    net.visited[edge10] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    [n.number for n in net.leaf] == [1,2,4,5,6] ? nothing : error("net.leaf is wrong")
end


# tree example
function testTree(net::HybridNetwork)
    !all([!e.hybrid for e in net.edge]) ? error("some edge is still hybrid") : nothing
    !all([!e.hybrid for e in net.node]) ? error("some node is still hybrid") : nothing
    !all([!e.hasHybEdge for e in net.node]) ? error("some node has hybrid edge") : nothing
    !all([e.isMajor for e in net.edge]) ? error("some edge is not major") : nothing
    !all([e.containRoot for e in net.edge]) ? error("some edge cannot contain root") : nothing
    edge9 = getIndexNumEdge(9,net);
    edge5 = getIndexNumEdge(5,net);
    (!net.edge[5].istIdentifiable || !net.edge[3].istIdentifiable) ? error("edge3,5 not identifiable") : nothing
    net.visited = [e.istIdentifiable for e in net.edge];
    net.visited[3] = false;
    net.visited[5] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    net.edge[2].length != 1.5 ? error("edge length for 2 is wrong") : nothing
    net.edge[4].length != 0.2 ? error("edge length for 4 is wrong") : nothing
    [n.number for n in net.leaf] == [1,2,4,6,7] ? nothing : error("net.leaf is wrong")
end
