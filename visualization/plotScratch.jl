#John Spaw
#drawCF takes a CF network object and converts it into .dot file, which is then converted into a .svg image file

function plotPhylonet(graph::Network;                           #Network object you are trying to plot
                 gammaThreshold=0.5::FloatingPoint,        #Gamma Threshold for extracting underlying tree structure
                 mainTree=false::Bool,                     #When true, the function will plot only the underlying tree structure
                 imageName="netImage",                     #Name for the file to be output
                 width=6::Number,                          #Maximum width of image in inches
                 height=8::Number,                         #Maximum height of image in inches
                 vert=true::Bool,                          #When true, function will display heirarchy from top to bottom. Otherwise, it will be from left to right
                 internalLabels=false::Bool,               #When true, all nodes will have labels (including internal nodes)
                 fontSize=16.0::FloatingPoint,             #Font size for labels in points
                 layoutStyle="dot"::String,
                 hybridColor="green4"::String
                 )

  #IO stream for writing to .dot file
  dotIo = open("visualization/drawCF.dot","w+")

  #Creates a list of all node numbers and all leaf numbers
  #Leaf numbers will be used for ranking process later on
  allNodes = Int64[]
  leafNodes = Int64[]
  for i in graph.node;
    num = i.number
    push!(allNodes,i.number)
    if i.leaf
      push!(leafNodes,i.number)
    end
  end

  netRoot = graph.root; #Assigns the integer for the network root... necessary for $ notation when appending strings

  global rootNode = graph.node[net.root]
  #for i in graph.node
  #  println("checking to see if $(i.number) is the root")
  #  if i.number == netRoot
  #    rootNode = i;
  #    println("Root found!")
  #  end
  #end

  println("Root node found... it is $(rootNode.number)")

  #********************************************************************************************************************

  #Writes initial preample lines for .dot file
  println("Creating preamble statement")
  write(dotIo,"Graph { \n")
  if vert == false
    write(dotIo,"     rankdir=LR; \n")
  end
  write(dotIo,"labelloc=b \n")                                    #Ensures that labels do not overlap each other (DOUBLE CHECK THIS)
  write(dotIo,"    ratio=\"fill\"; \n")                           #Fits graph to the full image size             (TEST OTHER RATIO OPTIONS)
  write(dotIo,"    size=\"$width ,$height\"; \n")                 #Changes the size of the entire graph
  write(dotIo,"    nodesep=0.5; \n")
  write(dotIo,"    edge [fontsize=$fontSize]; \n")
  write(dotIo,"    node [shape = point] \n")                      #Sets the shape of the nodes
  write(dotIo,"    rank=max $(rootNode.number) \n     subgraph    { ")      #Places root node at top of tree

  leafArraySize = countnz(leafNodes)                              #Probably redundant... could delete later on
  for i in leafNodes                                              #Groups leaf nodes so they are all placed at bottom of tree
    if i != leafNodes[leafArraySize]                              #First appends each leaf node (except the last) to .dot file followed by a comma
      write(dotIo,"$i , ")
    else                                                          #Appends final leaf node to .dot file WITHOUT a comma
      write(dotIo,"$i")
    end
  end
  write(dotIo," } \n")

  #Choosing which nodes have labels
  if internalLabels
    for i in allNodes
      write(dotIo,"    $i [xlabel=$i] [shape = point] \n")
    end #for i in allNodes
  else
    for i in leafNodes
      write(dotIo,"    $i [label=$i] [shape = point] \n")      #Applies labels to only leaf nodes... this will be removed once xlabels have been fixed
    end #for i in leafNodes
    write(dotIo," \n")
  end #if internalLabels else

  println("Preamble written!")
  #********************************************************************************************************************
  #Drawing edges and nodes

  println("Drawing nodes and edges")

  dummy = Edge(-1,1.0)      #dummy edge used to begin recursion because net.root does not have a parent edge

  #Traverse the network using a pseudo-depth-first method
  #traverseEdges is a recursive function that traverses each edge only once and appends it to the dot file
  traverseEdges(graph, rootNode, mainTree, dotIo, gammaThreshold, dummy,hybridColor)

  #********************************************************************************************************************

  #Writes closing lines to .dot file and closes the IO stream
  write(dotIo,"}")
  close(dotIo)
  println("Final lines written and IO stream has been closed")

  #********************************************************************************************************************

  #Converts .dot file into .svg image
  print("Exporting .dot file as .svg")
  generalExport("visualization/drawCF.dot",filename="$imageName",layoutEngine=layoutStyle)

end

function plotPhylonet(netString::String;
                 gammaThreshold=0.5::FloatingPoint,        #Gamma Threshold for extracting underlying tree structure
                 mainTree=false::Bool,                     #When true, the function will plot only the underlying tree structure
                 imageName="netImage",                     #Name for the file to be output
                 width=6::Number,                          #Maximum width of image in inches
                 height=8::Number,                         #Maximum height of image in inches
                 vert=true::Bool,                          #When true, function will display heirarchy from top to bottom. Otherwise, it will be from left to right
                 internalLabels=false::Bool,               #When true, all nodes will have labels (including internal nodes)
                 fontSize=16.0::FloatingPoint              #Font size for labels in points
                 )

  f = open("string_Net.txt","w")
  write(f,netString)
  close(f)
  net = readTopologyUpdate("string_Net.txt");
  plotNet(net,
                 gammaThreshold=0.5::FloatingPoint,        #Gamma Threshold for extracting underlying tree structure
                 mainTree=false::Bool,                     #When true, the function will plot only the underlying tree structure
                 imageName="netImage",                     #Name for the file to be output
                 width=6::Number,                          #Maximum width of image in inches
                 height=8::Number,                         #Maximum height of image in inches
                 vert=true::Bool,                          #When true, function will display heirarchy from top to bottom. Otherwise, it will be from left to right
                 internalLabels=false::Bool,               #When true, all nodes will have labels (including internal nodes)
                 fontSize=16.0::FloatingPoint,              #Font size for labels in points
                 layoutStyle="dot"::String
                 )
end