# test functions for the 5taxon networks
# functions used in tests_5taxon.jl
# Claudia September 2014
# changed to new reparametrizations: bad diamondI,II and bad triangle I,II
# Claudia November 2014

# Case C bad triangle II
function testCaseC(net::HybridNetwork)
    n=searchHybridNode(net);
    node = n[1];
    node.k != 3 ? error("k diff than 3") : nothing
    node.isVeryBadTriangle ? nothing : error("does not know it is very bad triangle")
    node.isExtBadTriangle ? error("thinks it is extremely bad triangle") : nothing
    net.hasVeryBadTriangle ? nothing : error("net does not know it has very bad triangle")
    net.numBad == 0 ? nothing : error("net.numBad should be 0")
    net.numHybrids != 1 ? error("should have 1 hybrid, but net.numHybrids is $(net.numHybrids): $([n.number for n in net.hybrid])") : nothing
end

# Case F bad diamond
function testCaseF(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    n = searchHybridNode(net);
    node = n[1];
    node.k != 4 ? error("k diff than 4") : nothing
    edge5 = getIndexEdge(5,net);
    edge9 = getIndexEdge(9,net);
    edge11 = getIndexEdge(11,net);
    edge12 = getIndexEdge(12,net);
    edge15 = getIndexEdge(15,net);
    edge14 = getIndexEdge(14,net);
    node1 = getIndexNode(1,net);
    node5 = getIndexNode(5,net);
    node11 = getIndexNode(11,net);
    node12 = getIndexNode(12,net);
    (net.edge[edge5].inCycle != node.number || net.edge[edge11].inCycle != node.number || net.edge[edge12].inCycle != node.number || net.edge[edge15].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[node1].inCycle  != node.number || net.node[node5].inCycle  != node.number || net.node[node11].inCycle  != node.number || net.node[node12].inCycle  != node.number) ? error("nodes 1,5,11,12 not correctly in cycle") : nothing
    !node.isBadDiamondI ? error("does not know it is bad diamondI") : nothing
    node.isBadDiamondII ? error("thinks it is bad diamond II") : nothing
    net.edge[edge14].containRoot ? error("edge can contain root") : nothing
    (!net.edge[edge11].hybrid || !net.edge[edge11].isMajor) ? error("edge 11 is not hybrid or major") : nothing
    net.node[node12].gammaz != net.edge[edge15].gamma*net.edge[edge12].z ? error("node 12 gammaz not correctly calculated") : nothing
    net.node[1].gammaz != net.edge[edge11].gamma*net.edge[edge5].z ? error("node 11 gammaz not correctly calculated") : nothing
    !net.edge[edge9].istIdentifiable ? error("edge9 not identifiable") : nothing
    net.visited[edge9] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    net.numHybrids != 1 ? error("should have 1 hybrid, but net.numHybrids is $(net.numHybrids): $([n.number for n in net.hybrid])") : nothing
    net.numBad == 1 ? nothing : error("net.numBad should be 1")
end

# Case G
function testCaseG(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    n = searchHybridNode(net);
    node = n[1];
    node.k != 4 ? error("k diff than 4") : nothing
    edge9 = getIndexEdge(9,net);
    edge12 = getIndexEdge(12,net);
    edge15 = getIndexEdge(15,net);
    edge5 = getIndexEdge(5,net);
    edge13 = getIndexEdge(13,net);
    edge14 = getIndexEdge(14,net);
    node9 = getIndexNode(9,net);
    node5 = getIndexNode(5,net);
    node11 = getIndexNode(11,net);
    node12 = getIndexNode(12,net);
    (net.edge[edge9].inCycle != node.number || net.edge[edge12].inCycle != node.number || net.edge[edge13].inCycle != node.number || net.edge[edge15].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[node9].inCycle  != node.number || net.node[node5].inCycle  != node.number || net.node[node11].inCycle  != node.number || net.node[node12].inCycle  != node.number) ? error("nodes not correctly in cycle") : nothing
    (node.isBadDiamondI || node.isBadDiamondII) ? error("thinks it is bad diamond") : nothing
    net.edge[edge14].containRoot ? error("14 can contain root") : nothing
    (!net.edge[edge12].hybrid || !net.edge[edge12].isMajor) ? error("edge 12 is not hybrid or major") : nothing
    net.node[node12].gammaz != -1 ? error("node 12 gammaz should be -1") : nothing
    (!net.edge[edge9].istIdentifiable || !net.edge[edge5].istIdentifiable || !net.edge[edge13].istIdentifiable) ? error("edge9,5,13not identifiable") : nothing
    net.visited[edge9] = false;
    net.visited[edge5] = false;
    net.visited[edge13] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    net.numHybrids != 1 ? error("should have 1 hybrid, but net.numHybrids is $(net.numHybrids): $([n.number for n in net.hybrid])") : nothing
end

# Case H
function testCaseH(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    n = searchHybridNode(net);
    node = n[1];
    node.k != 4 ? error("k diff than 4") : nothing
    edge9 = getIndexEdge(9,net);
    edge14 = getIndexEdge(14,net);
    edge15 = getIndexEdge(15,net);
    edge13 = getIndexEdge(13,net);
    edge5 = getIndexEdge(5,net);
    edge8 = getIndexEdge(8,net);
    node9 = getIndexNode(9,net);
    node5 = getIndexNode(5,net);
    node11 = getIndexNode(11,net);
    node12 = getIndexNode(12,net);
    (net.edge[edge9].inCycle != node.number || net.edge[edge5].inCycle != node.number || net.edge[edge14].inCycle != node.number || net.edge[edge15].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[node9].inCycle  != node.number || net.node[node5].inCycle  != node.number || net.node[node11].inCycle  != node.number || net.node[node12].inCycle  != node.number) ? error("nodes 9,5,11,12 not correctly in cycle") : nothing
    (node.isBadDiamondI || node.isBadDiamondII )? error("thinks it is bad diamond") : nothing
    net.edge[edge8].containRoot ? error("8 can contain root") : nothing
    (!net.edge[edge14].hybrid || !net.edge[edge14].isMajor) ? error("edge 14 is not hybrid or major") : nothing
    net.node[node12].gammaz != -1 ? error("node 12 gammaz should be -1") : nothing
    (!net.edge[edge9].istIdentifiable || !net.edge[edge5].istIdentifiable || !net.edge[edge13].istIdentifiable) ? error("edge9,5,13not identifiable") : nothing
    net.visited[edge9] = false;
    net.visited[edge5] = false;
    net.visited[edge13] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    net.numHybrids != 1 ? error("should have 1 hybrid, but net.numHybrids is $(net.numHybrids): $([n.number for n in net.hybrid])") : nothing
end


# Case J
function testCaseJ(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    n = searchHybridNode(net);
    node = n[1];
    node.k != 5 ? error("k diff than 5") : nothing
    edge9 = getIndexEdge(9,net);
    edge10 = getIndexEdge(10,net);
    edge5 = getIndexEdge(5,net);
    edge15 = getIndexEdge(15,net);
    edge6 = getIndexEdge(6,net);
    edge14 = getIndexEdge(14,net);
    node9 = getIndexNode(9,net);
    node1 = getIndexNode(1,net);
    node5 = getIndexNode(5,net);
    node11 = getIndexNode(11,net);
    node12 = getIndexNode(12,net);
    (net.edge[edge9].inCycle != node.number || net.edge[edge10].inCycle != node.number || net.edge[edge6].inCycle != node.number || net.edge[edge15].inCycle != node.number || net.edge[edge5].inCycle != node.number) ? error("edges not correctly in cycle") : nothing
    (net.node[node9].inCycle  != node.number || net.node[node5].inCycle  != node.number || net.node[node11].inCycle  != node.number || net.node[node12].inCycle  != node.number || net.node[node1].inCycle) != node.number ? error("nodes 9,5,11,12,1 not correctly in cycle") : nothing
    net.edge[edge14].containRoot ? error("14 can contain root") : nothing
    (!net.edge[edge6].hybrid || !net.edge[edge6].isMajor) ? error("edge 6 is not hybrid or major") : nothing
    net.node[node12].gammaz != -1 ? error("node 12 gammaz should be -1") : nothing
    (!net.edge[edge9].istIdentifiable || !net.edge[edge5].istIdentifiable || !net.edge[edge10].istIdentifiable) ? error("edge9,5,10not identifiable") : nothing
    net.visited[edge9] = false;
    net.visited[edge5] = false;
    net.visited[edge10] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    net.numHybrids != 1 ? error("should have 1 hybrid, but net.numHybrids is $(net.numHybrids): $([n.number for n in net.hybrid])") : nothing
end

# Case D bad triangle I
function testCaseD(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    n = searchHybridNode(net);
    node = n[1];
    node.k != 3 ? error("k diff than 3") : nothing
    node.isVeryBadTriangle ? nothing : error("does not know it is very bad triangle")
    node.isExtBadTriangle ? error("thinks it is extremely bad triangle") : nothing
    net.hasVeryBadTriangle ? nothing : error("net does not know it has very bad triangle")
    net.numBad == 0 ? nothing : error("net.numBad should be 0")
    net.numHybrids != 1 ? error("should have 1 hybrid, but net.numHybrids is $(net.numHybrids): $([n.number for n in net.hybrid])") : nothing
end

# Case E bad triangle I
function testCaseE(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    n = searchHybridNode(net);
    node = n[1];
    node.k != 3 ? error("k diff than 3") : nothing
    node.isVeryBadTriangle ? nothing : error("does not know it is very bad triangle")
    node.isExtBadTriangle ? error("thinks it is extremely bad triangle") : nothing
    net.hasVeryBadTriangle ? nothing : error("net does not know it has very bad triangle")
    net.numBad == 0 ? nothing : error("net.numBad should be 0")
    net.numHybrids != 1 ? error("should have 1 hybrid, but net.numHybrids is $(net.numHybrids): $([n.number for n in net.hybrid])") : nothing
end


# Case I bad diamond II
function testCaseI(net::HybridNetwork)
    n = searchHybridNode(net);
    node = n[1];
    node.isBadDiamondII ? nothing : error("does not know it is bad diamond II")
    !node.isBadDiamondI ? nothing : error("thinks it is bad diamond I")
    net.numHybrids != 1 ? error("should have 1 hybrid, but net.numHybrids is $(net.numHybrids): $([n.number for n in net.hybrid])") : nothing
    net.visited = [e.istIdentifiable for e in net.edge];
    edge14 = getIndexEdge(14,net);
    edge10 = getIndexEdge(10,net);
    edge5 = getIndexEdge(5,net);
    edge15 = getIndexEdge(15,net);
    edge9 = getIndexEdge(9,net);
    edge11 = getIndexEdge(11,net);
    edge8 = getIndexEdge(8,net);
    node1 = getIndexNode(1,net);
    node12 = getIndexNode(12,net);
    node5 = getIndexNode(5,net);
    node11 = getIndexNode(11,net);
    (net.edge[edge5].inCycle != node.number || net.edge[edge9].inCycle != node.number || net.edge[edge11].inCycle != node.number || net.edge[edge15].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[node1].inCycle  != node.number || net.node[node12].inCycle  != node.number || net.node[node5].inCycle  != node.number || net.node[node11].inCycle  != node.number) ? error("nodes 1,5,11,12 not correctly in cycle") : nothing
    (net.edge[edge14].containRoot || net.edge[edge10].containRoot || net.edge[edge8].containRoot) ? error("edges can contain root and shouldn't") : nothing
    (!net.edge[edge9].hybrid || !net.edge[edge9].isMajor) ? error("edge 4 is not hybrid or major") : nothing
    net.edge[edge14].length != 0 ? error("edges should have length 0") : nothing
    net.edge[edge14].istIdentifiable ? error("edge14 identifiable and should not") : nothing
    net.visited[edge9] = false;
    net.visited[edge5] = false;
    net.visited[edge11] = false;
    net.visited[edge15] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    net.numBad == 0 ? nothing : error("net.numBad should be 0")
end


# tree example
# warning: if added hybridization (bad diamond/triangle) and then delete,
#          original edge lengths cannot be recovered
function testTree(net::HybridNetwork)
    !all([!e.hybrid for e in net.edge]) ? error("some edge is still hybrid") : nothing
    !all([!e.hybrid for e in net.node]) ? error("some node is still hybrid") : nothing
    !all([!e.hasHybEdge for e in net.node]) ? error("some node has hybrid edge") : nothing
    !all([e.isMajor for e in net.edge]) ? error("some edge is not major") : nothing
    !all([e.containRoot for e in net.edge]) ? error("some edge cannot contain root") : nothing
    edge9 = getIndexEdge(9,net);
    edge5 = getIndexEdge(5,net);
    (!net.edge[edge9].istIdentifiable || !net.edge[edge5].istIdentifiable) ? error("edge9,5 not identifiable") : nothing
    net.visited = [e.istIdentifiable for e in net.edge];
    net.visited[edge9] = false;
    net.visited[edge5] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
    isempty(net.hybrid) ? nothing : error("something in net.hybrid")
    net.numHybrids != 0 ? error("should have 0 hybrid, but net.numHybrids is $(net.numHybrids)") : nothing
    #edge11 = getIndexEdge(11,net);
    #edge12 = getIndexEdge(12,net);
    #net.edge[edge11].length != 1.5 ? error("edge length for 11 is wrong") : nothing
    #net.edge[edge12].length != 0.2 ? error("edge length for 12 is wrong") : nothing
end
