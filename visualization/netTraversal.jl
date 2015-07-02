#John Spaw
#Network Edge Traversal Function
#6-25-15

#This function will start at the root of a network and traverse each edge a single time
#Will eventually be used to append edges to a .dot file while being able to keep track of which node is a parent of a given edge

dummy = Edge(-1,1.0)

function traverseEdges(net::HybridNetwork, node::Node, parentEdge=dummy::Edge)
  #nNum = node.number
  #println("On node $nNum")

  #Case 1: Node is the root
  if node.number == net.root
    println("Current node is the root")
    for edge in node.edge                      #Iterates through all edges associated with the node
      eNum = edge.number                       #Stores the edge number for use with $ syntax
      #println("Checking edge $eNum")
      println("Appending edge $eNum")       #This will later be replaced with a different action (writing the edge to the dot file)
      println("Found an acceptable edge to traverse... edge number $eNum")
      println("Setting new node...")
      println("Current node is $(node.number)")
      println("Node[1] = $((edge.node[1]).number)")
      println("Node[2] = $((edge.node[2]).number)")
      if edge.node[1] == node                  #Decides which node is the parent and which is the child
        newnode = edge.node[2]
      else
        newnode = edge.node[1]
      end #if else
      nnN = newnode.number                        #Simply used for progress/debugging statements... will be removed in the end
      println("About to call function with new node as $nnN")
      traverseEdges(net,newnode,edge)             #Recursively call the function on the child node (as new parent node) using the edge as the new parent edge
    end #for
  end #if

  #Case 2: Node is a leaf
  if node.leaf                                 #Don't need to do anything if the node is a leaf
    println("Node is a leaf")                  #There is probably a better statement to use here... pass/break?
  end

  #Case 3: Node is internal
  if (node.number != net.root) && ~(node.leaf)
    println("Function called with new node $(node.number)")
    println("Node is internal")
    println("Checking edges attached to node $(node.number)")
    for edge in node.edge                      #Iterate through edges associated with the node
      eNum = edge.number                       #Needed for $ notation in progress statements... will remove later
      println("Current edge is $eNum")
      println("Parent edge is $(parentEdge.number)")
      if edge != parentEdge
        println("$(edge.number) != $(parentEdge.number)")
        println("Appending edge $eNum")        #Action statement here... this is where we will append the edge to the .dot file
      end
      if (((~(edge.hybrid)) || ((edge.hybrid) && edge.isMajor)) && edge != parentEdge)  #Continue with the function if the edge is either a tree edge or a major hybrid
          if edge.gamma > 0.5


          println("Found an acceptable edge to traverse... edge number $eNum")
          println("Setting new node...")
          println("Current node is $(node.number)")
          println("Node[1] = $((edge.node[1]).number)")
          println("Node[2] = $((edge.node[2]).number)")
          if edge.node[1] == node                                                       #Determining which node is the parent/child
            println("Case A Succeeded: Node[1] is current node")
            newnode = edge.node[2]
            println("Set new node = Node[2] = $((edge.node[2]).number) ")
          elseif edge.node[2] == node
            println("Case B succeeded: Node[2] is current node")
            newnode = edge.node[1]
            println("Set new node = Node[1] = $((edge.node[1]).number)")
          else
            println("something is wrong")
          end #if else
          println("About to call function with new node $(newnode.number)")
          traverseEdges(net,newnode,edge)                             #Recursively call function with new child node (as parent) and previous edge as parent edge

          end #if gamma
      end #if
    end #for
  end #if
end #function
