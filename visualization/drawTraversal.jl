function traverseEdges(net::HybridNetwork,
                       node::Node,
                       mainTree::Bool,
                       dotIo,
                       gammaThreshold,
                       parentEdge=dummy::Edge,
                       hybridColor="green4"::String,
                       layoutStyle="dot"::String,
                       labelAngle= 180.0::FloatingPoint,
                       labelDistance= 3.0::FloatingPoint
                       )
  #*************************************************************************************************************************
  #Case 1: Node is the root
  if node.number == net.root
    #println("Node number $(node.number) is the root")
    for edge in node.edge                      #Iterates through all edges associated with the node
      eNum = edge.number                       #Stores the edge number for use with $ syntax
      #println("Appending edge $eNum")
      parent = node                            #Sets current node as parent node
      if (edge.node)[1] == node
        child = (edge.node)[2]                 #Finds other nodes attached to the edge and makes it the child
      else
        child = (edge.node)[1]
      end

      i = edge
      gma = round(i.gamma,3)        #Gamma value for an edge i (default to 1.0 for tree edges)
      hThickness = gma*4
      edgeNum = i.number
      node1 = parent;               #Parent node
      node1Num = node1.number       #Necessary for $ notation
      node2 = child;                #Child node
      node2Num = node2.number       #Necessary for $ notation


      #Creates image for the entire network including all hybridization events
          if node2.hybrid
              #Applies a thicker edge AND a gamma label on the edge to the dominant edge
              if gma < gammaThreshold
                write(dotIo,"     $node1Num -- $node2Num
                                  [color=$(hybridColor)]
                                  [penwidth=$hThickness]")
                if layoutStyle == "neato"
                  write(dotIo,"   [weight=$(i.length)]")
                end
                write(dotIo,"
                                  [headlabel=\" &gamma; = $gma\"]
                                  [labeldistance = 3.5]
                                  [labelangle=45.0]; \n")
              elseif mainTree == false                                 #Skip this step if you only want mainTree
                write(dotIo,"     $node1Num -- $node2Num
                                  [color=$(hybridColor)]")
                if layoutStyle == "neato"
                  write(dotIo,"   [weight=$(i.length)]")
                end
                write(dotIo,"
                                  [penwidth=$hThickness]; \n")
        end #if gma > 0.5 elseif mainTree == false
          else
              write(dotIo,"     $node1Num -- $node2Num [penwidth=4]; \n")
        end #if node2.hybrid else (net)

      if edge.node[1] == node                  #Decides which node is the parent and which is the child
        newnode = edge.node[2]
      else
        newnode = edge.node[1]
      end #if edge.node[1] == node else
      nnN = newnode.number                        #Simply used for progress/debugging statements... will be removed in the end
      traverseEdges(net,newnode,mainTree,dotIo,gammaThreshold,edge,hybridColor,layoutStyle,labelAngle,labelDistance)
    end #for edge in node.edge
  end #if node.number == net.root

  #*************************************************************************************************************************

  #Case 2: Node is a leaf
  if node.leaf                                                #Don't need to do anything if the node is a leaf
    #println("Node $(node.number) is a leaf")                  #There is probably a better statement to use here... pass/break?
  end #if node.leaf

  #*************************************************************************************************************************

  #Case 3: Node is internal
  if (node.number != net.root) && ~(node.leaf)
    #println("Node $(node.number) is internal")
    for edge in node.edge                      #Iterate through edges associated with the node
      eNum = edge.number                       #Needed for $ notation in progress statements... will remove later
      #println("Checking edge number $eNum")
      if edge != parentEdge

        parent = node
        if (edge.node)[1] == node
          child = (edge.node)[2]
        else
          child = (edge.node)[1]
        end

        if (edge.isMajor) || node.hybrid == false
        #println("Appending edge $eNum")        #Action statement here... this is where we will append the edge to the .dot file

        i = edge
        gma = round(i.gamma,3)                 #Gamma value for an edge i (default to 1.0 for tree edges)
        hThickness = gma*4
        edgeNum = i.number
        node1 = parent;            #Parent node
        node1Num = node1.number       #Necessary for $ notation
        node2 = child;            #Child node
        node2Num = node2.number       #Necessary for $ notation



        #Creates image for the entire network including all hybridization events
            if node2.hybrid
                #Applies a thicker edge AND a gamma label on the edge to the dominant edge
                if gma < gammaThreshold
                  write(dotIo,"   $node1Num -- $node2Num
                                  [color=$(hybridColor)]
                                  [penwidth=$hThickness]")
                   if layoutStyle == "neato"
                      write(dotIo,"   [len=$(i.length)]")
                   end
                  write(dotIo,"
                                  [label=\" &gamma; = $gma\"]
                                  [labeldistance = 3.5]
                                  [labelangle=45.0]; \n")
                elseif mainTree == false
                  write(dotIo,"   $node1Num -- $node2Num
                                  [color=$(hybridColor)]")
                    if layoutStyle == "neato"
                      write(dotIo,"   [weight=$(i.length)]")
                    end
                  write(dotIo,"
                                  [penwidth=$hThickness]; \n")
                end #if elseif mainTree == false
            elseif child.leaf
                write(dotIo,"     $node1Num -- $node2Num
                                  [headlabel=$(child.name)]")
                if layoutStyle == "neato"
                  write(dotIo,"   [weight=$(i.length)]")
                end
                write(dotIo,"
                                  [labeldistance=$(labelDistance)]
                                  [labelangle=$(labelAngle)]
                                  [penwidth=4]; \n")
            else
                write(dotIo,"     $node1Num -- $node2Num")
                if layoutStyle == "neato"
                  write(dotIo,"   [weight=$(i.length)]")
                end
                write(dotIo,"
                                  [penwidth=4]; \n")
              end #if else

      end
      if (((~(edge.hybrid)) || ((edge.hybrid) && edge.isMajor)) && edge != parentEdge)  #Continue with the function if the edge is either a tree edge or a major hybrid
          if edge.gamma > 0.5
            if edge.node[1] == node                                                       #Determining which node is the parent/child
              newnode = edge.node[2]
            elseif edge.node[2] == node
              newnode = edge.node[1]
            else
              println("something is wrong")
            end #if else
            traverseEdges(net,newnode,mainTree,dotIo,gammaThreshold,edge,hybridColor,layoutStyle, labelAngle, labelDistance)                             #Recursively call function with new child node (as parent) and previous edge as parent edge
          end #if gamma
        end #if
      end #if
    end #for
  end #if
end #function
