function traverseEdges(net::HybridNetwork,
                       node::Node,
                       mainTree::Bool,
                       dotIo,
                       gammaThreshold,
                       parentEdge=dummy::Edge,
                       hybridColor="green4"::AbstractString,
                       layoutStyle="dot"::AbstractString,
                       labelAngle= 180.0::AbstractFloat,
                       labelDistance= 3.0::AbstractFloat,
                       includeGamma=false::Bool,
                       includeLength=false::Bool
                       )
  # *********************************************
  # Case 1: Node is the root
  if node.number == net.node[net.root].number
    #println("Node number $(node.number) is the root")
    for edge in node.edge
      if (edge.node)[1] == node
        child = (edge.node)[2]      # Finds other node attached to the edge and makes it the child
      else
        child = (edge.node)[1]
      end

      gma = round(edge.gamma,3)        # Gamma for edge i (default to 1.0 for tree edges)
      hThickness = gma*4
      edgeNum = edge.number         # Stores edge number for use with $ syntax
      node1Num = node.number        # parent node number. Necessary for $ notation
      node2Num = child.number       # child  node number. Necessary for $ notation


      if child.hybrid
          #Applies a thicker edge AND a gamma label on the edge to the dominant edge
          if gma > gammaThreshold
            write(dotIo,"     $node1Num -- $node2Num\n\t[color=$(hybridColor)]\n")
            write(dotIo,"\t[penwidth=$hThickness]")
            if mainTree == true && includeGamma == true
              if includeLength == true
                write(dotIo,"\n\t[label=\" &gamma; = $gma \n l = $(edge.length)\"]; \n")
              else
                write(dotIo,"\n\t[label=\" &gamma; = $gma\"]; \n")
              end
            end
          elseif mainTree == false                                 #Skip this step if you only want mainTree
            write(dotIo,"     $node1Num -- $node2Num\n\t[color=$(hybridColor)]\n")
            write(dotIo,"\t[penwidth=$hThickness]")
            if includeGamma == true
              write(dotIo,"\n\t[label=\" &gamma; = $gma\"];\n")
            end
          end #if gma > 0.5 elseif mainTree == false

      elseif child.leaf
          write(dotIo,"     $node1Num -- $node2Num\n\t[headlabel=$(child.name)]\n")
          write(dotIo,"\t[labeldistance=$(labelDistance)]\n\t[labelangle=$(labelAngle)]\n\t[penwidth=4]")
          if includeLength == true
            write(dotIo,"\n\t[label=\" l = $(edge.length)\"]; \n")
          else
            write(dotIo,"; \n")
          end

      else # child is an internal tree node
          write(dotIo,"     $node1Num -- $node2Num\n\t[penwidth=4]")
          if includeLength == true
            write(dotIo,"\n\t[label=\" l = $(edge.length)\"]; \n")
          else
            write(dotIo,"; \n")
          end
      end

      traverseEdges(net,child,mainTree,dotIo,gammaThreshold,edge,hybridColor,layoutStyle,labelAngle,labelDistance,includeGamma,includeLength)
    end # loop over edges
  end # if node is root

  # *********************************************

  # Case 2: Node is a leaf
  if node.leaf                                  # do nothing
    #println("Node $(node.number) is a leaf")   # better statement to use here... pass/break?
  end #if node.leaf

  # *********************************************

  # Case 3: Node is internal
  if (node.number != net.node[net.root].number) && ~(node.leaf)
    #println("Node $(node.number) is internal")
    for edge in node.edge                      #Iterate through edges associated with the node
      if edge != parentEdge

        if (edge.node)[1] == node
          child = (edge.node)[2]
        elseif (edge.node)[2] == node
          child = (edge.node)[1]
        else
          warn("something is wrong")
        end

        if (edge.isMajor) || node.hybrid == false

        gma = round(edge.gamma,3)    # Gamma for edge (default to 1.0 for tree edges)
        hThickness = gma*4
        edgeNum = edge.number        # edge number. Needed for $ notation
        node1Num = node.number       # parent node number. Necessary for $ notation
        node2Num = child.number      # child  node number. Necessary for $ notation
        DEBUG && println("$(edge.number)     $(edge.length)")

        # Creates image for the entire network including all hybridization events
        if child.hybrid
            # thicker edge AND a gamma label on major hybrid edge
            if gma > gammaThreshold
              write(dotIo,"   $node1Num -- $node2Num\n\t[color=$(hybridColor)]\n\t[penwidth=$hThickness]")
              if mainTree == true && includeGamma == true
                  if includeLength == true
                    write(dotIo,"\n\t[label=\" &gamma; = $gma \n l = $(edge.length)\"]; \n")
                  else
                    write(dotIo,"\n\t[label=\" &gamma; = $gma\"]; \n")
                  end
              else 
                if includeLength == true
                  write(dotIo,"\n\t[label=\" l = $(edge.length)\"]; \n")
                end
              end
            elseif mainTree == false
              write(dotIo,"   $node1Num -- $node2Num\n\t[color=$(hybridColor)]")
              if includeGamma == true
                if includeLength == true
                  write(dotIo,"\n\t[label=\" &gamma; = $gma \n l = $(edge.length)\"] \n")
                else
                  write(dotIo,"\n\t[label=\" &gamma; = $gma\"]")
                end
              end #if includeGamma == true
              write(dotIo,"\n\t[penwidth=$hThickness]; \n")
            end #if child.hybrid elseif mainTree == false

        elseif child.leaf
            write(dotIo,"     $node1Num -- $node2Num\n\t[headlabel=$(child.name)]\n")
            write(dotIo,"\t[labeldistance=$(labelDistance)]\n\t[labelangle=$(labelAngle)]\n\t[penwidth=4]")
            if includeLength == true
              write(dotIo,"\n\t[label=\" l = $(edge.length)\"]; \n")
            else
              write(dotIo,"; \n")
            end

        else  # child is an internal tree node
            write(dotIo,"     $node1Num -- $node2Num\n\t[penwidth=4]")
            if includeLength == true
              write(dotIo,"\n\t[label=\" l = $(edge.length)\"]; \n")
            else
              write(dotIo,"; \n")
            end
        end
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
