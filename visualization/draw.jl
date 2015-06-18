#Function that will write to .dot file

#Function that will create .dot file

function drawCF(graph)
  dotIo = open("testdraw.dot","w+")
  nodes = Int64[]
  edges = Int64[]
  leafNodes = Int64[]
  hybridNodes = Int64[]
  hybrideEdges = Int64[]

  #Creates a list of all node numbers that are leaves... will be used for ranking of nodes
  for i in graph.node;
    num = i.number
    if i.leaf
      push!(leafNodes,i.number)
    end
  end

  netRoot = graph.root;

  #Writes initial preample lines for .dot file
  println("Creating preamble statement")
  write(dotIo,"Graph { \n")
  write(dotIo,"    node [shape = point] \n")

  #Places all leaf nodes in a 'subgraph' which ensures they are plotted on the same level
  write(dotIo,"    rank=max $netRoot")    #Places the tree root at the top of the plot
  write(dotIo,"\n")
  write(dotIo,"    subgraph { \n")
  write(dotIo,"        ")
  leafArraySize = countnz(leafNodes)
  for i in leafNodes
    if i != leafNodes[leafArraySize]
      write(dotIo,"$i")
      write(dotIo,", ")
    else
      write(dotIo,"$i")
    end
  end
  write(dotIo,"; \n")
  write(dotIo,"} \n")


  println("Preamble written successfully")

  #Draws all edges and nodes
  println("Drawing nodes and edges")
  for i in graph.edge
    edgeNum = i.number
    node1 = i.node[1];
    node1Num = node1.number
    node2 = i.node[2];
    node2Num = node2.number
    write(dotIo,"    $node1Num")
    write(dotIo," -- ")
    write(dotIo,"$node2Num; \n")
  end

  println("All nodes and edges drawn")

  #Writes closing lines to .dot file and closes the IO stream
  write(dotIo,"}")
  close(dotIo)
  println("Final lines written and IO stream has been closed")

  #Converts .dot file into .svg image
  print("Exporting .dot file as .svg")
  dotExport("testdraw.dot")
end
