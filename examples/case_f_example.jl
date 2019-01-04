# example "Case f" (plot in ipad)
# this is the bad diamond case
# Claudia August 2014
#
# in julia: include("case_f_example.jl")

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

n1=Node(1,false,false,[ed1,ed5,ed6]);
n2=Node(2,false,true,[ed1,ed2,ed3]);
n3=Node(3,false,false,[ed2,ed4,ed7]);
n4=Node(4,true,false,[ed3]);
n5=Node(5,false,false,[ed4,ed5,ed9]);
n6=Node(6,true, false,[ed6]);
n7=Node(7,true, false,[ed7]);
n8=Node(8,true, false,[ed8]);
n9=Node(9,false,false,[ed8,ed9,ed10]);
n10=Node(10,true, false,[ed10]);

setNode!(ed1,n1)
setNode!(ed1,n2)
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

#ed1.inCycle=2;
#ed2.inCycle=2;
#ed4.inCycle=2;
#ed5.inCycle=2;

net=HybridNetwork([n1,n2,n3,n4,n5,n6,n7,n8,n9,n10],[ed1,ed2,ed3,ed4,ed5,ed6,ed7,ed8,ed9,ed10]);
node=searchHybridNode(net);
net.names=["1","2","3","4","5","6","7","8","9","10"]

flag, nocycle,edges, nodes = updateInCycle!(net,node[1]);
flag2, edges2 = updateContainRoot!(net,node[1]);
flag3, edges3 = updateGammaz!(net,node[1]);

#printEdges(net)

#deleteHybrid!(node[1],net,false)

#changeDirectionUpdate!(node[1],net);

