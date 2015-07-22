#John Spaw
#Contains function that will convert a .dot file into a .svg image

using  GraphViz
#Converts a .dot file
function generalExport(file;filename="genImage"::String,layoutEngine="dot")
  dot = open(file,"r") do io Graph(io) end
  GraphViz.layout!(dot,engine=layoutEngine)
  open("visualization/$filename.svg","w") do f
    GraphViz.writemime(f, MIME"image/svg+xml"(),dot)
  end #do
  print("File saved")
end


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
                       includeGamma=true::Bool
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
              if gma > gammaThreshold
                write(dotIo,"     $node1Num -- $node2Num
                                  [color=$(hybridColor)]
                                  [penwidth=$hThickness]")
                if mainTree == true && includeGamma == true
                  write(dotIo,"
                                  [headlabel=\" &gamma; = $gma\"]
                                  [labeldistance = 3.5]
                                  [labelangle=45.0]; \n")
                end  #mainTree == true && includeGamma == true
              elseif mainTree == false                                 #Skip this step if you only want mainTree
                write(dotIo,"     $node1Num -- $node2Num
                                  [color=$(hybridColor)]
                                  [penwidth=$hThickness]")
                if includeGamma == true
                  write(dotIo,"
                                  [headlabel=\" &gamma; = $gma\"]
                                  [labeldistance = 3.5]
                                  [labelangle=45.0];\n")
                end  #if includeGamma == true
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
      traverseEdges(net,newnode,mainTree,dotIo,gammaThreshold,edge,hybridColor,layoutStyle,labelAngle,labelDistance,includeGamma)
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
                if gma > gammaThreshold
                  write(dotIo,"   $node1Num -- $node2Num
                                  [color=$(hybridColor)]
                                  [penwidth=$hThickness]")
                   if mainTree == true && includeGamma == true
                      write(dotIo,"
                                  [label=\" &gamma; = $gma\"]
                                  [labeldistance = 3.5]
                                  [labelangle=45.0]; \n")
                   end #if mainTree == true
                elseif mainTree == false
                  write(dotIo,"   $node1Num -- $node2Num
                                  [color=$(hybridColor)]")
                  if includeGamma == true
                    write(dotIo,"
                                  [label=\" &gamma; = $gma\"]
                                  [labeldistance = 3.5]
                                  [labelangle=45.0]")
                  end #if includeGamma == true
                  write(dotIo,"
                                  [penwidth=$hThickness]; \n")
                end #if node2.hybrid elseif mainTree == false
            elseif child.leaf
                write(dotIo,"     $node1Num -- $node2Num
                                  [headlabel=$(child.name)]
                                  [labeldistance=$(labelDistance)]
                                  [labelangle=$(labelAngle)]
                                  [penwidth=4]; \n")
            else
                write(dotIo,"     $node1Num -- $node2Num
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
            traverseEdges(net,newnode,mainTree,dotIo,gammaThreshold,edge,hybridColor,layoutStyle, labelAngle, labelDistance,includeGamma)                             #Recursively call function with new child node (as parent) and previous edge as parent edge
          end #if gamma
        end #if
      end #if
    end #for
  end #if
end #function

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
                 edgeStyle="truep"::String,                 #Sets the style of edges used. Options include "line", "ortho", "composite" (which uses both lines and curved splines), "curved"
                 labelAngle= 180.0::FloatingPoint,         #Sets the angle for leaf label placement
                 labelDistance= 3.0::FloatingPoint,        #Sets the distance for leaf label placement
                 includeGamma=true::Bool                   #When true, gamma labels are included on hybrid edges
                 )

  #IO stream for writing to .dot file
  dotIo = open("visualization/drawCF.dot","w+")

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
  traverseEdges(graph, rootNode, mainTree, dotIo, gammaThreshold, dummy,hybridColor,layoutStyle,labelAngle,labelDistance,includeGamma)

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
                 fontSize=16.0::FloatingPoint,             #Font size for labels in points
                 layoutStyle="dot"::String,                #Chooses the layout engine for placing nodes and edges (dot,neato,circo,twopi)
                 hybridColor="green4"::String,             #Sets color for hybrid edges
                 forcedLeaf=true::Bool,                    #When true, places all leaf nodes on the same line
                 unrooted=false::Bool,                     #Defaults to neato engine
                 nodeSeparation=0.5::FloatingPoint,        #Sets the minimum distance between nodes in inches
                 edgeStyle="true"::String,                 #Sets the style of edges used. Options include "line", "ortho", "composite" (which uses both lines and curved splines), "curved"
                 labelAngle= 180.0::FloatingPoint,         #Sets the angle for leaf label placement
                 labelDistance= 3.0::FloatingPoint,         #Sets the distance for leaf label placement
                 includeGamma=true::Bool                   #When true, gamma labels are included on hybrid edges
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
                 includeGamma=includeGamma
                 )
end
