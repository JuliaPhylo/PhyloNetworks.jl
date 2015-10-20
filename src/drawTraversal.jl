function traverseEdges(net::HybridNetwork,
                       node::Node,
                       mainTree::Bool,
                       dotIo,
                       gammaThreshold,
                       parentEdge=dummy::Edge,
                       hybridColor="green4"::String,
                       layoutStyle="dot"::String,
                       labelAngle= 180.0::FloatingPoint,
                       labelDistance= 3.0::FloatingPoint,
                       includeGamma=false::Bool,
                       includeLength=false::Bool
                       )
  #*************************************************************************************************************************
  #Case 1: Node is the root
  if node.number == net.root
    #println("Node number $(node.number) is the root")
    for edge in node.edge                      #Iterates through all edges associated with the node
      eNum = edge.number                       #Stores the edge number for use with $ syntax
      #println("Appending edge $eNum")
      if (edge.node)[1] == node
        child = (edge.node)[2]                 #Finds other node attached to the edge and makes it the child
      else
        child = (edge.node)[1]
      end

      i = edge
      gma = round(i.gamma,3)        #Gamma value for an edge i (default to 1.0 for tree edges)
      hThickness = gma*4
      edgeNum = i.number
      node1Num = node.number        # parent node number. Necessary for $ notation
      node2Num = child.number       # child  node number. Necessary for $ notation


      #Creates image for the entire network including all hybridization events
          if child.hybrid
              #Applies a thicker edge AND a gamma label on the edge to the dominant edge
              if gma > gammaThreshold
                write(dotIo,"     $node1Num -- $node2Num
                                  [color=$(hybridColor)]
                                  [penwidth=$hThickness]")
                if mainTree == true && includeGamma == true
                  if includeLength == true
                    write(dotIo,"
                                  [label=\" &gamma; = $gma \n l = $(i.length)\"]; \n")
                  else
                    write(dotIo,"
                                  [label=\" &gamma; = $gma\"]; \n")
                  end #if includeLength == false else
                end  #mainTree == true && includeGamma == true
              elseif mainTree == false                                 #Skip this step if you only want mainTree
                write(dotIo,"     $node1Num -- $node2Num
                                  [color=$(hybridColor)]
                                  [penwidth=$hThickness]")
                if includeGamma == true
                  write(dotIo,"
                                  [label=\" &gamma; = $gma\"];\n")
                end  #if includeGamma == true
              end #if gma > 0.5 elseif mainTree == false
          else
              write(dotIo,"     $node1Num -- $node2Num [penwidth=4]")
              if includeLength == true
                write(dotIo," [label=\" l = $(i.length)\"]; \n")
              else
                write(dotIo,"; \n")
              end
        end #if child.hybrid else (net)

      nnN = child.number                        #Simply used for progress/debugging statements... will be removed in the end
      traverseEdges(net,child,mainTree,dotIo,gammaThreshold,edge,hybridColor,layoutStyle,labelAngle,labelDistance,includeGamma,includeLength)
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

        if (edge.node)[1] == node
          child = (edge.node)[2]
        elseif (edge.node)[2] == node
          child = (edge.node)[1]
        else
          warn("something is wrong")
        end

        if (edge.isMajor) || node.hybrid == false
        #println("Appending edge $eNum")        #Action statement here... this is where we will append the edge to the .dot file

        i = edge
        gma = round(i.gamma,3)                 #Gamma value for an edge i (default to 1.0 for tree edges)
        hThickness = gma*4
        edgeNum = i.number
        node1Num = node.number       # parent node number. Necessary for $ notation
        node2Num = child.number      # child  node number. Necessary for $ notation
        DEBUG && println("$(i.number)     $(i.length)")



        #Creates image for the entire network including all hybridization events
            if child.hybrid
                #Applies a thicker edge AND a gamma label on the edge to the dominant edge
                if gma > gammaThreshold
                  write(dotIo,"   $node1Num -- $node2Num
                                  [color=$(hybridColor)]
                                  [penwidth=$hThickness]")
                  if mainTree == true && includeGamma == true
                      if includeLength == true
                        write(dotIo,"
                                  [label=\" &gamma; = $gma \n l = $(i.length)\"]; \n")
                      else
                        write(dotIo,"
                                  [label=\" &gamma; = $gma\"]; \n")
                      end
                  else #if mainTree == true && include Gamma == true else
                    if includeLength == true
                      write(dotIo,"
                                  [label=\" l = $(i.length)\"]; \n")
                    end
                  end #if mainTree == true
                elseif mainTree == false
                  write(dotIo,"   $node1Num -- $node2Num
                                  [color=$(hybridColor)]")
                  if includeGamma == true
                    if includeLength == true
                      write(dotIo,"[label=\" &gamma; = $gma \n l = $(i.length)\"] \n")
                    else
                      write(dotIo,"
                                    [label=\" &gamma; = $gma\"]")
                    end
                end #if includeGamma == true
                  write(dotIo,"
                                  [penwidth=$hThickness]; \n")
                end #if child.hybrid elseif mainTree == false
            elseif child.leaf
                write(dotIo,"     $node1Num -- $node2Num
                                  [headlabel=$(child.name)]
                                  [labeldistance=$(labelDistance)]
                                  [labelangle=$(labelAngle)]
                                  [penwidth=4]")
                if includeLength == true
                  write(dotIo," [label=\" l = $(i.length)\"]; \n")
                else
                  write(dotIo,"; \n")
                end
            else
                write(dotIo,"     $node1Num -- $node2Num
                                  [penwidth=4]")
                if includeLength == true
                  write(dotIo," [label=\" l = $(i.length)\"]; \n")
                else
                  write(dotIo,"; \n")
                end
            end #if else

      end
      if (((~(edge.hybrid)) || ((edge.hybrid) && edge.isMajor)) && edge != parentEdge)  # recursive call if the edge is either a tree edge or a major hybrid
          if edge.gamma > 0.5 # fixit: would be best to remove this, instead check that 
            traverseEdges(net,child,mainTree,dotIo,gammaThreshold,edge,hybridColor,layoutStyle, labelAngle, labelDistance,includeGamma,includeLength)
          end #if gamma
        end #if
      end #if
    end #for
  end #if
end #function
