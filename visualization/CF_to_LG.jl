#using  LightGraphs

#####################################################
#TEST VALUES
#Creating a graph to test the function with

#Creating edges
e1 = Edge(1,5.0)
e2 = Edge(2,1.0)
e3 = Edge(3,7.0)
println("Edges created")

#Creating nodes
n1 = Node(1, false, false, [e1, e2]);
n2 = Node(2, true, false, [e1]);
n3 = Node(3, false, false, [e2, e3]);
n4 = Node(4, true, false, [e3]);
println("Nodes created")

#Setting nodes with edges
setNode!(e1,[n1,n2])
setNode!(e2,[n1,n3])
setNode!(e3,[n3,n4])

net=HybridNetwork([n1,n2,n3,n4],[e1,e2,e3])
println("CF network created")
#####################################################

function CF_to_Light_Graph(net)
  n = net.numNodes
  println("The number of nodes in the network is $(n)")

  lgNet = LightGraphs.Graph(n)
  println(lgNet)

  eCont = net.edge

  println(" ")

  for i in eCont
    nodesAttached = i.node
    n1=nodesAttached[1].number
    n2=nodesAttached[2].number
    println("There is an edge between $(n1) and $(n2)")
    add_edge!(lgNet,n1,n2)
    println("Added an edge between $(n1) and $(n2) in the lgNet")
    println(" ")
  end

  return lgNet
end


CF_to_Light_Graph(net)




