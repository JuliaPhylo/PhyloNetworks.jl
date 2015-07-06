function drawEdge(edge::Edge, parent::Node, child::Node, dotIo, gammaThreshold=0.5::FloatingPoint, mainTree=false::Bool)
  i = edge
  gma = i.gamma                 #Gamma value for an edge i (default to 1.0 for tree edges)
    hThickness = gma*4
    edgeNum = i.number
    node1 = parent;            #Parent node
    node1Num = node1.number       #Necessary for $ notation
    node2 = child;            #Child node
    node2Num = node2.number       #Necessary for $ notation

  #Creates image for underlying tree structure according the gamma threshold
  #mainTree is a bool type (optional) function parameter that decides whether the image will be the network or underlying tree structure (defaults as net)
  if mainTree
    if node2.hybrid
      if i.gamma > gammaThreshold
        write(dotIo,"     $node1Num -- $node2Num
                          [color=blue]
                          [penwidth=4]
                          [taillabel=\" &gamma; = $gma\"]
                          [labeldistance = 5.0]
                          [labelangle=45.0]; \n")
      end #if
    else
      write(dotIo,"     $node1Num -- $node2Num
                        [penwidth=4]; \n")
    end #if else
  #Creates image for the entire network including all hybridization events
  else
      if node2.hybrid
          #Applies a thicker edge AND a gamma label on the edge to the dominant edge
          if gma > 0.5
            write(dotIo,"     $node1Num -- $node2Num
                              [color=blue]
                              [penwidth=$hThickness]
                              [taillabel=\" &gamma; = $gma\"]
                              [labeldistance = 3.5]
                              [labelangle=60.0]; \n")
          else
            write(dotIo,"     $node1Num -- $node2Num
                              [color=red]
                              [penwidth=$hThickness]; \n")
          end #if else
      else
          write(dotIo,"     $node1Num -- $node2Num [penwidth=4]; \n")
        end #if else
      end #if else
end #function

