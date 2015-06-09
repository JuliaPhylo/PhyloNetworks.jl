#Graphs for testing algorithms
#John Spaw 6-9-15
#Contains numerous examples of CF trees and networks used for testing algorithms
################################################################################


#####################################################
#Simple tree with 4 nodes, 3 edges, 0 hybrid nodes
#####################################################
function create_g1()
	#Edges
	e1 = Edge(1,5.0)
	e2 = Edge(2,1.0)
	e3 = Edge(3,7.0)
	println("Edges created")
	
	#Nodes
	n1 = Node(1, false, false, [e1, e2]);
	n2 = Node(2, true, false, [e1]);
	n3 = Node(3, false, false, [e2, e3]);
	n4 = Node(4, true, false, [e3]);
	println("Nodes created")

	#Nodes --> edges
	setNode!(e1,[n1,n2])
	setNode!(e2,[n1,n3])
	setNode!(e3,[n3,n4])
	
	g1=HybridNetwork([n1,n2,n3,n4],[e1,e2,e3])
	println("CF network created")
	
	return g1;
end
		#This graph has been tested and works correctly

#####################################################
