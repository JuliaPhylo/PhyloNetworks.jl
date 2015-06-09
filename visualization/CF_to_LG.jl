include("../types.jl")
include("../functions.jl")
using  LightGraphs

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






