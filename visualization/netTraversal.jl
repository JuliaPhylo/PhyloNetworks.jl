#John Spaw
#Network Edge Traversal Function
#6-25-15

#This function will start at the root of a network and traverse each edge a single time
#Will eventually be used to append edges to a .dot file while being able to keep track of which node is a parent of a given edge

function traverseEdges(net, node; parentEdge=null::any)
  nNum = node.number
  println("On node $nNum")

  #Case 1: Node is the root
  if node.number == net.root
    println("Node is the root")
    for edge in node.edge                      #Iterates through all edges associated with the node
      eNum = edge.number                       #Stores the edge number for use with $ syntax
      println("Checking edge $eNum")
      println("Append edge to dot file")       #This will later be replaced with a different action (writing the edge to the dot file)
      if edge.node[1] == node                  #Decides which node is the parent and which is the child
        node = edge.node[2]
      else
        node = edge.node[1]
      end #if else
      nnN = node.number                        #Simply used for progress/debugging statements... will be removed in the end
      println("New node is $nnN")
      traverseEdges(net,node,edge)             #Recursively call the function on the child node (as new parent node) using the edge as the new parent edge
    end #for
  end #if

  #Case 2: Node is a leaf
  if node.leaf                                 #Don't need to do anything if the node is a leaf
    println("Node is a leaf")                  #There is probably a better statement to use here... pass/break?
  end

  #Case 3: Node is internal
  if (node.number != net.root) && ~(node.leaf)
    for edge in node.edge                      #Iterate through edges associated with the node
      eNum = edge.number                       #Needed for $ notation in progress statements... will remove later
      println("Checking edge $eNum")
      println("Node is internal")
      if edge != parentEdge
        println("Append edge to dot file")     #Action statement here... this is where we will append the edge to the .dot file
      end
      if (~(edge.hybrid)) || ((edge.hybrid) && edge.isMajor)   #Continue with the function if the edge is either a tree edge or a major hybrid
          if edge.node[1] == node                              #Determining which node is the parent/child
            node = edge.node[2]
          else
            node = edge.node[1]
          end #if else
          traverseEdges(net,node,edge)                             #Recursively call function with new child node (as parent) and previous edge as parent edge
      end #if
    end #for
  end #if
end #function
