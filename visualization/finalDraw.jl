#John Spaw
#drawCF takes a CF network object and converts it into .dot file, which is then converted into a .svg image file

include("visualization/drawTraversal.jl")

function plotNet(graph::Network; gammaThreshold=0.5::FloatingPoint, mainTree=false::Bool,imageName="netImage",
                 width=8::Number,height=5::Number)

  #Argument Breakdown
      #graph: Network object you are trying to visualize
      #gammaThreshold: Lower bound for gamma value.
                      #When ploting networks, gamma values above threshold will be bolded.
                      #When plotting trees, edges below threshold are ignored.
      #mainTree: If true, function will plot underlying tree structure. If false, function will plot entire network.
      #imageName: Names the image being output as imageName.svg
      #width: Maximum width of image (in inches)
      #height: Maximum height of image (in inches)


  #IO stream for writing to .dot file
  dotIo = open("visualization/drawCF.dot","w+")

  #Creates a list of all node numbers that are leaves... will be used for ranking of nodes later on
  leafNodes = Int64[]
  for i in graph.node;
    num = i.number
    if i.leaf
      push!(leafNodes,i.number)
    end
  end

  #Assigns the integer for the network root
  netRoot = graph.root;

  #********************************************************************************************************************

  #Writes initial preample lines for .dot file
  println("Creating preamble statement")
  write(dotIo,"Graph { \n")
  write(dotIo,"labelloc=b \n")                                  #Ensures that labels do not overlap each other (DOUBLE CHECK THIS)
  write(dotIo,"    ratio=\"expand\"; \n")                         #Fits graph to the full image size             (TEST OTHER RATIO OPTIONS)
  write(dotIo,"    size=\"$width ,$height\"; \n")                 #Changes the size of the entire graph
  write(dotIo,"    node [shape = point] \n")                    #Sets the shape of the nodes
  write(dotIo, "    rank=max $netRoot \n     subgraph    { ")   #Places root node at top of tree

  leafArraySize = countnz(leafNodes)                 #Probably redundant... could delete later on
  for i in leafNodes                                 #Groups leaf nodes so they are all placed at bottom of tree
    if i != leafNodes[leafArraySize]                 #First appends each leaf node (except the last) to .dot file followed by a comma
      write(dotIo,"$i , ")
    else                                             #Appends final leaf node to .dot file WITHOUT a comma
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
  #Drawing edges and nodes

  #Establishing edge/node variables... necessary for calling values using $ notation when appending strings to a file
  println("Drawing nodes and edges")

  dummy = Edge(-1,1.0)

  traverseEdges(graph, (graph.node)[1], mainTree, dotIo, dummy)

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
