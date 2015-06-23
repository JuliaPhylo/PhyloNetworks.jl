#drawCF takes a CF network object and converts it into .dot file, which is then converted into a .svg image file

function drawCF(graph)
  dotIo = open("visualization/drawCF.dot","w+")

  #Creates a list of all node numbers that are leaves... will be used for ranking of nodes
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
  write(dotIo,"    ratio=\"fill\"; \n")
  write(dotIo,"    size=\"8,5\"; \n")
  write(dotIo,"    node [shape = point] \n")

  #Places root node at top of tree
  write(dotIo, "    rank=max $netRoot \n     subgraph    {
                                        ")
  #Groups leaf nodes so they are all placed at bottom of tree
  leafArraySize = countnz(leafNodes)
  for i in leafNodes
    if i != leafNodes[leafArraySize]
      write(dotIo,"$i , ")
    else
      write(dotIo,"$i")
    end
  end

  write(dotIo," } \n")

  for i in leafNodes
    write(dotIo,"    $i [xlabel=$i] \n")
  end
  write(dotIo," \n")

  #********************************************************************************************************************

  println("Preamble written successfully")

  #Draws all edges and nodes
  println("Drawing nodes and edges")
  for i in graph.edge
    edgeNum = i.number
    node1 = i.node[1];
    node1Num = node1.number
    node2 = i.node[2];
    node2Num = node2.number

    if node2.hybrid
      write(dotIo,"     $node1Num -- $node2Num [color=blue] [penwidth=4]; \n")
    else
      write(dotIo,"     $node1Num -- $node2Num [penwidth=4]; \n")
    end
  end

  println("All nodes and edges drawn")

  #********************************************************************************************************************

  #Writes closing lines to .dot file and closes the IO stream
  write(dotIo,"}")
  close(dotIo)
  println("Final lines written and IO stream has been closed")

  #********************************************************************************************************************

  #Converts .dot file into .svg image
  print("Exporting .dot file as .svg")
  dotExport("visualization/drawCF.dot")

end
