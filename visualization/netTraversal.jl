#John Spaw
#Network Edge Traversal Function
#6-25-15

#This function will start at the root of a network and traverse each edge a single time
#Will eventually be used to append edges to a .dot file while being able to keep track of which node is a parent of a given edge

function traverseEdges(net, node, parentEdge=false)
  nNum = node.number
  println("On node $nNum")
  if node.number == net.root
    println("Node is the root")
    for edge in node.edge
      eNum = edge.number
      println("Checking edge $eNum")
      println("Append edge to dot file")
      if edge.node[1] == node
        node = edge.node[2]
      else
        node = edge.node[1]
      end
      nnN = node.number
      println("New node is $nnN")
      traverseEdges(net,node,edge)
      print()
      end #if else
    end #for
  end #if
  if node.leaf
    println("Node is a leaf")
  end
  if (node.number != net.root) && ~(node.leaf)
    for edge in node.edge
      eNum = edge.number
      println("Checking edge $eNum")
      println("Append edge to dot file")
      println("Node is internal")
      if edge != parentEdge
        println("Append edge to dot file")
      end
      if (~(edge.hybrid)) || ((edge.hybrid) && edge.isMajor)
          if edge.node[1] == node
            node = edge.node[2]
          else
            node = edge.node[1]
          end #if else
          traverseEdges(node,edge)
      end #if
    end #for
  end #if
end #function
