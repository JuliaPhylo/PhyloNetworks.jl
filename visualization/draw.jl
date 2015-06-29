#John Spaw
#drawCF takes a CF network object and converts it into .dot file, which is then converted into a .svg image file

function drawCF(graph::Network; gammaThreshold=0.5::FloatingPoint, mainTree=false::Bool,imageName="netImage")

  dotIo = open("visualization/drawCF.dot","w+")

  #Creates a list of all node numbers that are leaves... will be used for ranking of nodes later on
  leafNodes = Int64[]
  for i in graph.node;
    num = i.number
    if i.leaf
      push!(leafNodes,i.number)
    end
  end

  netRoot = graph.root;

  #********************************************************************************************************************

  #Writes initial preample lines for .dot file

  println("Creating preamble statement")
  write(dotIo,"Graph { \n")
  write(dotIo,"labelloc=b \n")                  #Ensures that labels do not overlap each other
  write(dotIo,"    ratio=\"fill\"; \n")         #Fits graph to the full image size
  write(dotIo,"    size=\"8,5\"; \n")           #Changes the size of the entire graph
  write(dotIo,"    node [shape = point] \n")    #Sets the shape of the nodes

  write(dotIo, "    rank=max $netRoot \n     subgraph    { ")   #Places root node at top of tree
  leafArraySize = countnz(leafNodes)
  for i in leafNodes                            #Groups leaf nodes so they are all placed at bottom of tree
    if i != leafNodes[leafArraySize]
      write(dotIo,"$i , ")
    else
      write(dotIo,"$i")
    end
  end
  write(dotIo," } \n")

  for i in leafNodes
    write(dotIo,"    $i [label=$i] [shape = circle] [height = 0.1] \n")        #Applies labels to only leaf nodes
  end
  write(dotIo," \n")

  println("Preamble written!")
  #********************************************************************************************************************

  #Draws all edges and nodes

  println("Drawing nodes and edges")
  for i in graph.edge
    gma = i.gamma
    edgeNum = i.number
    node1 = i.node[1];
    node1Num = node1.number
    node2 = i.node[2];
    node2Num = node2.number

  #Creates image for underlying tree structure according the gamma threshold
  if mainTree
    if node2.hybrid
      if i.gamma > gammaThreshold
        write(dotIo,"     $node1Num -- $node2Num [color=blue] [penwidth=4] [taillabel=\" &gamma; = $gma\"] [labeldistance = 5.0] [labelangle=45.0]; \n")
      end #if
    else
      write(dotIo,"     $node1Num -- $node2Num [penwidth=4]; \n")
      end #if else
  #Creates image for the entire network including all hybridization events
  else
      if node2.hybrid
          if gma > gammaThreshold
            write(dotIo,"     $node1Num -- $node2Num [color=blue] [penwidth=4] [taillabel=\" &gamma; = $gma\"] [labeldistance = 5.0]; \n")
          else
            write(dotIo,"     $node1Num -- $node2Num [color=blue] [penwidth=2]; \n")
          end #if else
        else
          write(dotIo,"     $node1Num -- $node2Num [penwidth=4]; \n")
        end #if else
      end #if else
    end #if else

  println("All nodes and edges drawn")

  #********************************************************************************************************************

  #Writes closing lines to .dot file and closes the IO stream
  write(dotIo,"}")
  close(dotIo)
  println("Final lines written and IO stream has been closed")

  #********************************************************************************************************************

  #Converts .dot file into .svg image
  print("Exporting .dot file as .svg")
  dotExport("visualization/drawCF.dot",filename="$imageName")

end
