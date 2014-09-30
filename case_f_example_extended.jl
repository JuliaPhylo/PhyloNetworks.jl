# example "Case f" (plot in ipad)
# this is the bad diamond case
# Claudia August 2014
#
# extended with a second hybridization below
#
# in julia: include("case_f_example_extended.jl")

include("types.jl")
include("functions.jl")
using Base.Collections # for updateInCycle with priority queue

ed1=Edge(1,0.6,true,0.7);
ed2=Edge(2,0.7,true,0.3);
ed3=Edge(3,0.9);
ed4=Edge(4,0.1);
ed5=Edge(5,0.2);
ed6=Edge(6,0.1);
ed7=Edge(7,0.1);
ed8=Edge(8,0.1);
ed9=Edge(9,0.1);
ed10=Edge(10,0.1);
ed11=Edge(11,0.1);
ed12=Edge(12,0.1);
ed13=Edge(13,0.1);
ed14=Edge(14,0.1);
ed15=Edge(15,0.1);
ed16=Edge(16,0.1);

n1=Node(1,false,false,[ed1,ed5,ed6]);
n2=Node(2,false,true,[ed1,ed2,ed3]);
n3=Node(3,false,false,[ed2,ed4,ed7]);
n4=Node(4,false,false,[ed3,ed11,ed12]);
n5=Node(5,false,false,[ed4,ed5,ed9]);
n6=Node(6,true, false,[ed6]);
n7=Node(7,true, false,[ed7]);
n8=Node(8,true, false,[ed8]);
n9=Node(9,false,false,[ed8,ed9,ed10]);
n10=Node(10,true, false,[ed10]);
n11=Node(11,true, false,[ed11]);
n12=Node(12,false, false,[ed12,ed13,ed14]);
n13=Node(13,false, false,[ed13,ed15,ed16]);
n14=Node(14,true, false,[ed14]);
n15=Node(15,true, false,[ed15]);
n16=Node(16,true, false,[ed16]);


setNode!(ed1,n1);
setNode!(ed1,n2);
setNode!(ed2,n2);
setNode!(ed2,n3);
setNode!(ed5,[n5,n1]);
setNode!(ed4,[n3,n5]);
setNode!(ed3,n2);
setNode!(ed3,n4);
setNode!(ed6,[n1,n6]);
setNode!(ed7,[n3,n7]);
setNode!(ed8,[n8,n9]);
setNode!(ed9,[n5,n9]);
setNode!(ed10,[n9,n10]);
setNode!(ed11,[n4,n11]);
setNode!(ed12,[n4,n12]);
setNode!(ed13,[n12,n13]);
setNode!(ed14,[n12,n14]);
setNode!(ed15,[n13,n15]);
setNode!(ed16,[n13,n16]);

#ed1.inCycle=2;
#ed2.inCycle=2;
#ed4.inCycle=2;
#ed5.inCycle=2;

net=HybridNetwork([n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16],[ed1,ed2,ed3,ed4,ed5,ed6,ed7,ed8,ed9,ed10,ed11,ed12,ed13,ed14,ed15,ed16]);

# we add the new hybridization from edge 14 to edge 13

# we need a different chooseEdgesGamma to control the edges1,2 chosen
# change index1, index2 according to the specific case
# warning: not real chooseEdgesGamma, the real one is in functions.jl
function chooseEdgesGamma(net::HybridNetwork)
warn("function chooseEdgesGamma is deterministic")
    index1 = 14;
    index2 = 13;
    edge1 = net.edge[index1];
    edge2 = net.edge[index2];
    gamma = rand()*0.5;
    return edge1, edge2, gamma
end


printEdges(net)
printNodes(net)
success,hybrid,flag,nocycle,flag2,flag3 = addHybridizationUpdate!(net);
printEdges(net)
printNodes(net)


hybnodes=searchHybridNode(net);

# now we want to update the first hybridization
success,hybrid,flag,nocycle,flag2,flag3 = updateAllNewHybrid!(hybnodes[1],net)

nocycle1, edgesIncycle1, nodesIncycle1 = identifyInCycle(net,hybnodes[1]);
nocycle2, edgesIncycle2, nodesIncycle2 = identifyInCycle(net,hybnodes[2]);

[e.number for e in nodesIncycle1]
[e.number for e in nodesIncycle2]
[e.number for e in edgesIncycle1]
[e.number for e in edgesIncycle2]

# only done to check, it changes istIdentifiable of an edge in the bad triangle
#flag1 = identifyGammaz(net,hybnodes[1]);
#flag2 = identifyGammaz(net,hybnodes[2]);

edgesRoot1 = identifyContainRoot(net,hybnodes[1]);
edgesRoot2 = identifyContainRoot(net,hybnodes[2]);

[e.number for e in edgesRoot1]
[e.number for e in edgesRoot2]
