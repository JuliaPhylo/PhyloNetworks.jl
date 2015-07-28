#John Spaw
#drawCF takes a CF network object and converts it into .dot file, which is then converted into a .svg image file

function plotPhylonet(graph::Network;                      #Network object you are trying to plot
                 gammaThreshold=0.5::FloatingPoint,        #Gamma Threshold for extracting underlying tree structure
                 mainTree=false::Bool,                     #When true, the function will plot only the underlying tree structure
                 imageName="netImage",                     #Name for the file to be output
                 width=6::Number,                          #Width of image in inches
                 height=8::Number,                         #Height of image in inches
                 vert=true::Bool,                          #When true, function will display heirarchy from top to bottom. Otherwise, it will be from left to right
                 internalLabels=false::Bool,               #When true, all nodes will have labels (including internal nodes)
                 fontSize=16.0::FloatingPoint,             #Font size for labels in points
                 layoutStyle="dot"::String,                #Chooses the layout engine for placing nodes and edges (dot,neato,circo,twopi)
                 hybridColor="green4"::String,             #Sets color for hybrid edges
                 forcedLeaf=true::Bool,                    #When true, places all leaf nodes on the same line
                 unrooted=false::Bool,                     #Defaults to neato engine
                 nodeSeparation=0.8::FloatingPoint,        #Sets the minimum distance between nodes in inches
                 edgeStyle="false"::String,                 #Sets the style of edges used. Options include "line", "ortho", "composite" (which uses both lines and curved splines), "curved"
                 labelAngle= 180.0::FloatingPoint,         #Sets the angle for leaf label placement
                 labelDistance= 3.0::FloatingPoint,        #Sets the distance for leaf label placement
                 includeGamma=false::Bool,                   #When true, gamma labels are included on hybrid edges
                 includeLength=false::Bool

                 )

  #IO stream for writing to .dot file
  dotIo = open("drawCF.dot","w+")

  #Creates a list of all node numbers and all leaf numbers
  #Leaf numbers will be used for ranking process later on
  allNodes = Int64[]
  leafNodes = Int64[]

  #Change the default engine if unrooted is selected
  if unrooted
    layoutStyle = "neato"
  end

  #Create list of nodes in order to append labels later on
  if internalLabels
    for i in graph.node;
      num = i.number
      push!(allNodes,i.number)
    end #for i in graph.node
  end #if internalLabels

  #Create array of leaf nodes for use in ranking and label creation
  for i in graph.leaf
    #if i != graph.node[graph.root]
     push!(leafNodes,i.number)
    #end
  end


  netRoot = graph.root;
  global rootNode = graph.node[graph.root]

  #********************************************************************************************************************

  #Writes initial preample lines for .dot file
  write(dotIo,"Graph Phylogeny { \n")
  if vert == false
    write(dotIo,"     rankdir=LR; \n")
  end
  write(dotIo,"    labelloc=b \n")                                #Ensures that labels do not overlap each other (DOUBLE CHECK THIS)
  write(dotIo,"    ratio=\"fill\"; \n")                           #Fits graph to the full image size             (TEST OTHER RATIO OPTIONS)
  write(dotIo,"    size=\"$width ,$height !\"; \n")               #Changes the size of the entire graph
  write(dotIo,"    nodesep=$(nodeSeparation); \n")                #Minimum distance between any two nodes
  write(dotIo,"    splines=$(edgeStyle); \n")                     #Edge style argument (see function argument declaration)
  write(dotIo,"    edge [fontsize=$fontSize]; \n")                #Fontsize argument
  write(dotIo,"    node [shape = point] \n")                      #Sets the shape of the nodes
  write(dotIo,"    rank=max $(rootNode.number) \n")               #Guarantees placement of root node at top of hierarchy (dot engine only)

  if unrooted
    write(dotIo,"  mode = KK; \n")
    write(dotIo,"  epsilon = 0.01; \n")
  end

  #Ranking leaves together at the bottom of the tree
  if forcedLeaf
    write(dotIo,"  subgraph    { ")                                 #Places root node at top of tree
    leafArraySize = countnz(leafNodes)                              #Probably redundant... could delete later on
    for i in leafNodes                                              #Groups leaf nodes so they are all placed at bottom of tree
      if i != leafNodes[leafArraySize]                              #First appends each leaf node (except the last) to .dot file followed by a comma
        write(dotIo,"$i , ")
      else                                                          #Appends final leaf node to .dot file WITHOUT a comma
        write(dotIo,"$i")
      end
    end
    write(dotIo," } \n")
  end

  #Choosing which nodes have labels
  if internalLabels
    for i in allNodes
      write(dotIo,"    $i [xlabel=$i] [shape = point] \n")           #Places node label on internal nodes and leaf nodes
    end #for i in allNodes
  else
    for i in leafNodes
      write(dotIo,"    $i [label=$i] [shape = point] \n")            #Applies labels to only leaf nodes... this will be removed once xlabels have been fixed
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
  traverseEdges(graph, rootNode, mainTree, dotIo, gammaThreshold, dummy,hybridColor,layoutStyle,labelAngle,labelDistance,includeGamma,includeLength)

  #********************************************************************************************************************

  #Writes closing lines to .dot file and closes the IO stream
  write(dotIo,"}")
  close(dotIo)
  println("Final lines written and IO stream has been closed")

  #********************************************************************************************************************

  #Converts .dot file into .svg image
  print("Exporting .dot file as .svg")
  generalExport("drawCF.dot",filename="$imageName",layoutEngine=layoutStyle)

end

function plotPhylonet(netString::String;
                 gammaThreshold=0.5::FloatingPoint,        #Gamma Threshold for extracting underlying tree structure
                 mainTree=false::Bool,                     #When true, the function will plot only the underlying tree structure
                 imageName="netImage",                     #Name for the file to be output
                 width=6::Number,                          #Maximum width of image in inches
                 height=8::Number,                         #Maximum height of image in inches
                 vert=true::Bool,                          #When true, function will display heirarchy from top to bottom. Otherwise, it will be from left to right
                 internalLabels=false::Bool,               #When true, all nodes will have labels (including internal nodes)
                 fontSize=16.0::FloatingPoint,             #Font size for labels in points
                 layoutStyle="dot"::String,                #Chooses the layout engine for placing nodes and edges (dot,neato,circo,twopi)
                 hybridColor="green4"::String,             #Sets color for hybrid edges
                 forcedLeaf=true::Bool,                    #When true, places all leaf nodes on the same line
                 unrooted=false::Bool,                     #Defaults to neato engine
                 nodeSeparation=0.5::FloatingPoint,        #Sets the minimum distance between nodes in inches
                 edgeStyle="false"::String,                #Sets the style of edges used. Options include "line", "ortho", "composite" (which uses both lines and curved splines), "curved"
                 labelAngle= 180.0::FloatingPoint,         #Sets the angle for leaf label placement
                 labelDistance= 3.0::FloatingPoint,        #Sets the distance for leaf label placement
                 includeGamma=false::Bool,                 #When true, gamma labels are included on hybrid edges
                 includeLength=false::Bool                 #When true, edge length labels are included
                 )


  net = readTopologyUpdate(netString,true);
  plotPhylonet(net,
                 gammaThreshold=gammaThreshold,
                 mainTree=mainTree,
                 imageName=imageName,
                 width=width,
                 height=height,
                 vert=vert,
                 internalLabels=internalLabels,
                 fontSize=fontSize,
                 layoutStyle=layoutStyle,
                 hybridColor=hybridColor,
                 forcedLeaf=forcedLeaf,
                 unrooted=unrooted,
                 nodeSeparation=nodeSeparation,
                 edgeStyle=edgeStyle,
                 labelAngle=labelAngle,
                 labelDistance=labelDistance,
                 includeGamma=includeGamma,
                 includeLength=includeLength
                 )
end
