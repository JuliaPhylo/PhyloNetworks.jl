#Graph sandbox
#John Spaw 6-8-15


#Creating edges
e1 = Edge(1,5.0)
e2 = Edge(2,1.0)
e3 = Edge(3,7.0)

#Creating nodes
n1 = Node(1, false, false, [e1, e2]);
n2 = Node(2, true, false, [e1]);
n3 = Node(3, false, false, [e2, e3]);
n4 = Node(4, true, false, [e3]);

#Setting nodes with edges
setNode!(e1,[n1,n2])
setNode!(e2,[n1,n3])
setNode!(e3,[n3,n4])

net=HybridNetwork([n1,n2,n3,n4],[e1,e2,e3])

printNodes(net)
printEdges(net)