# test functions for the 5taxon networks
# functions used in tests_5taxon.jl
# Claudia September 2014

# Case C
function testCaseC(net::HybridNetwork)
    node=searchHybridNode(net);
    node.k != 3 ? error("k diff than 3") : nothing
    edge9 = getIndexNumEdge(9,net);
    edge5 = getIndexNumEdge(5,net);
    edge13 = getIndexNumEdge(13,net);
    edge15 = getIndexNumEdge(15,net);
    edge12 = getIndexNumEdge(12,net);
    edge14 = getIndexNumEdge(14,net);
    node5 = getIndexNumNode(5,net);
    node11 = getIndexNumNode(11,net);
    node12 = getIndexNumNode(12,net);
    net.visited = [e.istIdentifiable for e in net.edge];
    (net.edge[edge9].inCycle != 11 || net.edge[edge12].inCycle != 11 || net.edge[edge15].inCycle != 11) ? error("edges not incycle") : nothing
    (net.node[node5].inCycle != 11 || net.node[node11].inCycle != 11 || net.node[node12].inCycle != 11) ? error("nodes not incycle") : nothing
    net.edge[edge14].containRoot ? error("contain root wrong") : nothing
    (net.edge[edge13].istIdentifiable || net.edge[edge13].length != 0.0) ? error("gammaz updated wrong") : nothing
    node.isBadTriangle ? error("thinks it is bad triangle") : nothing
    (!net.edge[edge9].istIdentifiable || !net.edge[edge5].istIdentifiable) ? error("edge9 or 15 not identifiable") : nothing
    net.visited[edge9] = false;
    net.visited[edge5] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
end

# Case F bad diamond
function testCaseF(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    node=searchHybridNode(net);
    node.k != 4 ? error("k diff than 4") : nothing
    edge5 = getIndexNumEdge(5,net);
    edge9 = getIndexNumEdge(9,net);
    edge11 = getIndexNumEdge(11,net);
    edge12 = getIndexNumEdge(12,net);
    edge15 = getIndexNumEdge(15,net);
    edge14 = getIndexNumEdge(14,net);
    node1 = getIndexNumNode(1,net);
    node5 = getIndexNumNode(5,net);
    node11 = getIndexNumNode(11,net);
    node12 = getIndexNumNode(12,net);
    (net.edge[edge5].inCycle != node.number || net.edge[edge11].inCycle != node.number || net.edge[edge12].inCycle != node.number || net.edge[edge15].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[node1].inCycle  != node.number || net.node[node5].inCycle  != node.number || net.node[node11].inCycle  != node.number || net.node[node12].inCycle  != node.number) ? error("nodes 1,5,11,12 not correctly in cycle") : nothing
    !node.isBadDiamond ? error("does not know it is bad diamond") : nothing
    net.edge[edge14].containRoot ? error("edge can contain root") : nothing
    (!net.edge[edge11].hybrid || !net.edge[edge11].isMajor) ? error("edge 11 is not hybrid or major") : nothing
    net.node[node12].gammaz != net.edge[edge15].gamma*net.edge[edge12].z ? error("node 12 gammaz not correctly calculated") : nothing
    !net.edge[edge9].istIdentifiable ? error("edge9 not identifiable") : nothing
    net.visited[edge9] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
end

# Case G
function testCaseG(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    node=searchHybridNode(net);
    node.k != 4 ? error("k diff than 4") : nothing
    edge9 = getIndexNumEdge(9,net);
    edge12 = getIndexNumEdge(12,net);
    edge15 = getIndexNumEdge(15,net);
    edge5 = getIndexNumEdge(5,net);
    edge13 = getIndexNumEdge(13,net);
    edge14 = getIndexNumEdge(14,net);
    node9 = getIndexNumNode(9,net);
    node5 = getIndexNumNode(5,net);
    node11 = getIndexNumNode(11,net);
    node12 = getIndexNumNode(12,net);
    (net.edge[edge9].inCycle != node.number || net.edge[edge12].inCycle != node.number || net.edge[edge13].inCycle != node.number || net.edge[edge15].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[node9].inCycle  != node.number || net.node[node5].inCycle  != node.number || net.node[node11].inCycle  != node.number || net.node[node12].inCycle  != node.number) ? error("nodes not correctly in cycle") : nothing
    node.isBadDiamond ? error("thinks it is bad diamond") : nothing
    net.edge[edge14].containRoot ? error("14 can contain root") : nothing
    (!net.edge[edge12].hybrid || !net.edge[edge12].isMajor) ? error("edge 12 is not hybrid or major") : nothing
    net.node[node12].gammaz != -1 ? error("node 12 gammaz should be -1") : nothing
    (!net.edge[edge9].istIdentifiable || !net.edge[edge5].istIdentifiable || !net.edge[edge13].istIdentifiable) ? error("edge9,5,13not identifiable") : nothing
    net.visited[edge9] = false;
    net.visited[edge5] = false;
    net.visited[edge13] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
end

# Case H
function testCaseH(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    node=searchHybridNode(net);
    node.k != 4 ? error("k diff than 4") : nothing
    edge9 = getIndexNumEdge(9,net);
    edge14 = getIndexNumEdge(14,net);
    edge15 = getIndexNumEdge(15,net);
    edge13 = getIndexNumEdge(13,net);
    edge5 = getIndexNumEdge(5,net);
    edge8 = getIndexNumEdge(8,net);
    node9 = getIndexNumNode(9,net);
    node5 = getIndexNumNode(5,net);
    node11 = getIndexNumNode(11,net);
    node12 = getIndexNumNode(12,net);
    (net.edge[edge9].inCycle != node.number || net.edge[edge5].inCycle != node.number || net.edge[edge14].inCycle != node.number || net.edge[edge15].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[node9].inCycle  != node.number || net.node[node5].inCycle  != node.number || net.node[node11].inCycle  != node.number || net.node[node12].inCycle  != node.number) ? error("nodes 9,5,11,12 not correctly in cycle") : nothing
    node.isBadDiamond ? error("thinks it is bad diamond") : nothing
    net.edge[edge8].containRoot ? error("8 can contain root") : nothing
    (!net.edge[edge14].hybrid || !net.edge[edge14].isMajor) ? error("edge 14 is not hybrid or major") : nothing
    net.node[node12].gammaz != -1 ? error("node 12 gammaz should be -1") : nothing
    (!net.edge[edge9].istIdentifiable || !net.edge[edge5].istIdentifiable || !net.edge[edge13].istIdentifiable) ? error("edge9,5,13not identifiable") : nothing
    net.visited[edge9] = false;
    net.visited[edge5] = false;
    net.visited[edge13] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
end


# Case J
function testCaseJ(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    node=searchHybridNode(net);
    node.k != 5 ? error("k diff than 5") : nothing
    edge9 = getIndexNumEdge(9,net);
    edge10 = getIndexNumEdge(10,net);
    edge5 = getIndexNumEdge(5,net);
    edge15 = getIndexNumEdge(15,net);
    edge6 = getIndexNumEdge(6,net);
    edge14 = getIndexNumEdge(14,net);
    node9 = getIndexNumNode(9,net);
    node1 = getIndexNumNode(1,net);
    node5 = getIndexNumNode(5,net);
    node11 = getIndexNumNode(11,net);
    node12 = getIndexNumNode(12,net);
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
end

# Case D bad triangle
function testCaseD(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    node = searchHybridNode(net);
    edge = getHybridEdge(node);
    node.k != 3 ? error("k diff than 3") : nothing
    edge6 = getIndexNumEdge(6,net);
    edge9 = getIndexNumEdge(9,net);
    edge11 = getIndexNumEdge(11,net);
    edge14 = getIndexNumEdge(14,net);
    edge15 = getIndexNumEdge(15,net);
    edge12 = getIndexNumEdge(12,net);
    edge5 = getIndexNumEdge(5,net);
    node5 = getIndexNumNode(5,net);
    node11 = getIndexNumNode(11,net);
    node12 = getIndexNumNode(12,net);
    (net.edge[edge12].inCycle != node.number || net.edge[edge15].inCycle != node.number || net.edge[edge5].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[node5].inCycle  != node.number || net.node[node11].inCycle  != node.number || net.node[node12].inCycle  != node.number) ? error("nodes 5,11,12 not correctly in cycle") : nothing
    !node.isBadTriangle ? error("does not know it is bad triangle") : nothing
     (net.edge[edge14].containRoot || net.edge[edge6].containRoot || net.edge[edge11].containRoot)  ? error("14,6,11 can contain root") : nothing
    (!net.edge[edge5].hybrid || !net.edge[edge5].isMajor) ? error("edge 5 is not hybrid or major") : nothing
    net.node[node11].gammaz != edge.gamma*edge.gamma*net.edge[edge12].z+(1-edge.gamma)*(1-edge.gamma)*net.edge[edge5].z ? error("node 11 gammaz not correctly calculated") : nothing
    (net.edge[edge14].length != 0.0 || net.edge[edge14].istIdentifiable) ? error("edge 14 not correctly non identifiable length 0") : nothing
    (net.edge[edge15].length != 0.0 || net.edge[edge15].istIdentifiable) ? error("edge 15 not correctly non identifiable length 0") : nothing
    net.node[node12].gammaz != edge.gamma*net.edge[edge12].z ? error("node 12 gammaz not correctly calculated") : nothing
    !net.edge[edge9].istIdentifiable ? error("edge9 not identifiable") : nothing
    net.visited[edge9] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
end

# Case E
function testCaseE(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    node = searchHybridNode(net);
    edge = getHybridEdge(node);
    node.k != 3 ? error("k diff than 3") : nothing
    edge8 = getIndexNumEdge(8,net);
    edge10 = getIndexNumEdge(10,net);
    edge14 = getIndexNumEdge(14,net);
    edge13 = getIndexNumEdge(13,net);
    edge15 = getIndexNumEdge(15,net);
    edge9 = getIndexNumEdge(9,net);
    edge5 = getIndexNumEdge(5,net);
    node5 = getIndexNumNode(5,net);
    node11 = getIndexNumNode(11,net);
    node12 = getIndexNumNode(12,net);
    (net.edge[edge9].inCycle != node.number || net.edge[edge5].inCycle != node.number || net.edge[edge15].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[node5].inCycle  != node.number || net.node[node11].inCycle  != node.number || net.node[node12].inCycle  != node.number) ? error("nodes 5,11,12 not correctly in cycle") : nothing
    node.isBadTriangle ? error("thinks it is bad triangle") : nothing
     (net.edge[edge14].containRoot || net.edge[edge8].containRoot || net.edge[edge10].containRoot)  ? error("14,8,10 can contain root") : nothing
    (!net.edge[edge9].hybrid || !net.edge[edge9].isMajor) ? error("edge 9 is not hybrid or major") : nothing
    net.node[node11].gammaz != -1 ? error("node 11 gammaz not correctly calculated") : nothing
    (net.edge[edge13].length != 0.0 || net.edge[edge13].istIdentifiable) ? error("edge 13 not correctly non identifiable length 0") : nothing
    (net.edge[edge14].length != 0.0 || net.edge[edge14].istIdentifiable) ? error("edge 14 not correctly non identifiable length 0") : nothing
    (net.edge[edge15].length != 0.0 || net.edge[edge15].istIdentifiable) ? error("edge 15 not correctly non identifiable length 0") : nothing
    net.node[node12].gammaz != -1 ? error("node 12 gammaz not correctly calculated") : nothing
    (!net.edge[edge9].istIdentifiable || !net.edge[edge5].istIdentifiable) ? error("edge9 or 5 not identifiable") : nothing
    net.visited[edge9] = false;
    net.visited[edge5] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
end


# Case I
function testCaseI(net::HybridNetwork)
    net.visited = [e.istIdentifiable for e in net.edge];
    node=searchHybridNode(net);
    node.k != 4 ? error("k diff than 4") : nothing
    edge5 = getIndexNumEdge(5,net);
    edge9 = getIndexNumEdge(9,net);
    edge10 = getIndexNumEdge(10,net);
    edge14 = getIndexNumEdge(14,net);
    edge15 = getIndexNumEdge(15,net);
    edge11 = getIndexNumEdge(11,net);
    edge8 = getIndexNumEdge(8,net);
    node1 = getIndexNumNode(1,net);
    node5 = getIndexNumNode(5,net);
    node11 = getIndexNumNode(11,net);
    node12 = getIndexNumNode(12,net);
    (net.edge[edge5].inCycle != node.number || net.edge[edge9].inCycle != node.number || net.edge[edge15].inCycle != node.number || net.edge[edge11].inCycle != node.number ) ? error("edges not correctly in cycle") : nothing
    (net.node[node1].inCycle  != node.number || net.node[node5].inCycle  != node.number || net.node[node11].inCycle  != node.number || net.node[node12].inCycle  != node.number) ? error("nodes 1,5,11,12 not correctly in cycle") : nothing
    node.isBadDiamond ? error("thinks it is bad diamond") : nothing
    (net.edge[edge14].containRoot || net.edge[edge8].containRoot || net.edge[edge10].containRoot) ? error("14,8,10 can contain root") : nothing
    (!net.edge[edge9].hybrid || !net.edge[edge9].isMajor) ? error("edge 9 is not hybrid or major") : nothing
    net.node[node12].gammaz != -1 ? error("node 12 gammaz not correctly calculated") : nothing
    (net.edge[edge15].length != 0.0 || net.edge[edge15].istIdentifiable) ? error("edge 15 not non identifiable or length not 0") : nothing
     (!net.edge[edge9].istIdentifiable || !net.edge[edge5].istIdentifiable || !net.edge[edge11].istIdentifiable || !net.edge[edge14].istIdentifiable) ? error("edge9,5,11,14 not identifiable") : nothing
    net.visited[edge9] = false;
    net.visited[edge5] = false;
    net.visited[edge11] = false;
    net.visited[edge14] = false;
    !all([!id for id in net.visited]) ? error("edges not identifiable as identifiable") : nothing
end
