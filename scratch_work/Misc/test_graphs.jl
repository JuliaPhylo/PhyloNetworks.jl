#Graphs for testing algorithms
#John Spaw 6-9-15
#Contains numerous examples of CF trees and networks used for testing algorithms
################################################################################

#Steps for creating CF network

#Creating edges: example_edge = Edge(int label, weight/length, hybridBool, .. , .. , gamma)
#Creating nodes: example_node = Node(int label, leafBool, hybridBool, gamma, edge array)
#Linking nodes with edges: setNode!(example_edge, array of nodes)
#Create CF network: example_network = HybridNetwork(array of nodes, array of edges)

################################################################################

#NOTE: Drawings of all graphs are included on physical notesheet titled "test_graphs.jl"

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

			#Notes:
			#This graph has been successfully converted into LightGraphs type

#####################################################
#Simple Network with 4 nodes, 3 edges, 0 hybrid nodes
#####################################################
function create_g2()
	#Edges

	e1 = Edge(1,1.0)
	e2 = Edge(2,1.0)
	h3 = Edge(3,1.0)
	h4 = Edge(4,1.0)
	println("Edges created")

	#e1 = Edge(1,1.0,false)
	#e2 = Edge(2,1.0,false)
	#h3 = Edge(3,1.0,true)
	#h4 = Edge(4,1.0,true)
	#println("Edges created")

	#Nodes
	n1 = Node(1,false,false,[e1,e2])
	n2 = Node(2,false,true,[e1,h3])
	n3 = Node(3,false,true,[e2,h4])
	n4 = Node(4,false,false,[h3,h4])
	println("Nodes created")

	#Nodes --> edges
	setNode!(e1,[n1,n2])
	setNode!(e2,[n1,n3])
	setNode!(h3,[n2,n4])
	setNode!(h4,[n3,n4])
	println("Nodes and edges linked")

	g2=HybridNetwork([n1,n2,n3,n4],[e1,e2,h3,h4])
	println("CF network with hybrid edges/nodes created")

	return g2
end

#######################################################################
#Medium Network with 8 nodes, 7 edges, 2 hybrid edges, 1 hybrid node
#######################################################################

function create_g3()

  #Edge creation
  e1 = Edge(1,1.0)
  e2 = Edge(2,1.0)
  e3 = Edge(3,1.0)
  h4 = Edge(4,1.0)
  h4.gamma = 0.9
  h5 = Edge(5,1.0)
  h5.isMajor = false
  h5.gamma = 0.1
  e6 = Edge(6,1.0)
  e7 = Edge(7,1.0)
  e8 = Edge(8,1.0)
  e9 = Edge(9,1.0)
  println("Edges created")

  #Node creation
  n1 = Node(1,false,false,[e1,e2])
  n2 = Node(2,false,false,[e1,e3,h4])
  n3 = Node(3,false,false,[e2,h5,e6])
  n4 = Node(4,true,false,[e3])
  n5 = Node(5,false,true,[h4,h5,e7])
  n6 = Node(6,true,false,[e6])
  n7 = Node(7,false,false,[e7,e8,e9])
  n8 = Node(8,true,false,[e8])
  n9 = Node(9,true,false,[e9])
  println("Nodes created")

  #Nodes --> Edges
  setNode!(e1,[n1,n2])
  setNode!(e2,[n1,n3])
  setNode!(e3,[n2,n4])
  setNode!(h4,[n2,n5])
  setNode!(h5,[n3,n5])
  setNode!(e6,[n3,n6])
  setNode!(e7,[n5,n7])
  setNode!(e8,[n7,n8])
  setNode!(e9,[n7,n9])
  println("Nodes and edges attached")

  g3 = HybridNetwork([n1,n2,n3,n4,n5,n6,n7,n8,n9],[e1,e2,e3,h4,h5,e6,e7,e8,e9])
  println("Network successfully created!")

  return g3
end







