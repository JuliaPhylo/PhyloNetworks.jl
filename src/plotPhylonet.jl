# initial function: John Spaw, 7/2015
# takes a network object and converts it into .dot file, which is then converted into a .svg image file
"""
`plotPhylonet(net::HybridNetwork)`

function to plot a HybridNetwork object. The plot will be saved in the working directory as a svg file. We are working on allowing other file formats, and to have the plot pop out in a window.
This function has the following optional arguments:
- imageName: name for plot file (default netImage)
- mainTree: if true, only the underlying tree (with major hybrid edges) is plotted (default false)
- width: width of image in inches (default 6)
- height: height of image in inches (default 8)
- vert: if true, plot displayed from top to bottom (default tre)
- internalLabels: if true, prints number for internal nodes (default false)
- fontSize: font size for taxon names (default 16.0)
- hybridColor: color for hybrid edges (default green4)
- unrooted: if true, prints the topology unrooted
- nodeSeparation: minimum distance between nodes in inches (default 0.8)
- labelAngle: angle for leaf label placement (default 180.0)
- labelDistance: distance for leaf label placement (default 3.0)
- includeGamma: if true, includes the gamma values in the plot (default false)
- includeLength: if true, includes the branch lengths in the plot (default false
"""
function plotPhylonet(graph::Network;    # Network object to plot
  imageName="netImage",                  # Name for output files
  gammaThreshold=0.5::FloatingPoint,     # threshold for ...
  mainTree=false::Bool,                  # When true, only the underlying tree will be plotted
  width=6::Number,                       # Width of image in inches
  height=8::Number,                      # Height of image in inches
  vert=true::Bool,                       # When true, hierarchy displayed from top to bottom. Otherwise, left to right
  internalLabels=false::Bool,            # When true, all nodes will have labels (including internal nodes)
  fontSize=16.0::FloatingPoint,          # Font size for labels in points
  layoutStyle="dot"::AbstractString,             # layout engine for placing nodes and edges (dot,neato,circo,twopi)
  hybridColor="green4"::AbstractString,          # color for hybrid edges
  forcedLeaf=true::Bool,                 # When true, places all leaf nodes on the same line
  unrooted=false::Bool,                  # if true, enforces layoutStyle to 'neato'
  nodeSeparation=0.8::FloatingPoint,     # minimum distance between nodes in inches
  edgeStyle="true"::AbstractString,              # style of edges. Options: "line", "ortho", "composite" (which uses both lines and curved splines), "curved"
  labelAngle= 180.0::FloatingPoint,      # angle for leaf label placement
  labelDistance= 3.0::FloatingPoint,     # distance for leaf label placement
  includeGamma=false::Bool,              # When true, gamma labels are displayed on hybrid edges
  includeLength=false::Bool
  )

  # check that the root is correcly placed (Claudia, Aug2015)
  if(!isTree(graph))
      if(!graph.cleaned)
          DEBUG && println("net not cleaned inside plotPhylonet, need to run updateCR")
          for(n in graph.hybrid)
              flag,edges = updateContainRoot!(graph,n)
              flag || error("hybrid node $(n.hybrid) has conflicting containRoot")
          end
      end
      checkRootPlace!(graph)
  end

  #IO stream for writing to .dot file
  dotIo = open("$imageName.dot","w+")

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

  rootNode = graph.node[graph.root]

  #********************************************************************************************************************

  #Writes initial preample lines for .dot file
  write(dotIo,"Graph Phylogeny { \n")
  if vert == false
    write(dotIo,"     rankdir=LR; \n")
  end
  write(dotIo,"    labelloc=b \n")                     # Ensures that labels do not overlap each other (DOUBLE CHECK THIS)
  write(dotIo,"    ratio=\"fill\"; \n")                # Fits graph to the full image size             (TEST OTHER RATIO OPTIONS)
  write(dotIo,"    size=\"$width ,$height !\"; \n")    # size of the entire graph
  write(dotIo,"    nodesep=$(nodeSeparation); \n")     # Minimum distance between any two nodes
  write(dotIo,"    splines=$(edgeStyle); \n")          # Edge style argument
  write(dotIo,"    edge [fontsize=$fontSize]; \n")
  write(dotIo,"    node [shape = point] \n")
  write(dotIo,"    rank=max $(rootNode.number) \n")    # Guarantees placement of root node at top of hierarchy (dot engine only)

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

  DEBUG && println("Preamble written!")
  #********************************************************************************************************************
  #Drawing edges and nodes

  DEBUG && println("Drawing nodes and edges")

  dummy = Edge(-1,1.0)      #dummy edge used to begin recursion because net.root does not have a parent edge

  #Traverse the network using a pseudo-depth-first method
  #traverseEdges is a recursive function that traverses each edge only once and appends it to the dot file
  traverseEdges(graph, rootNode, mainTree, dotIo, gammaThreshold, dummy,hybridColor,layoutStyle,labelAngle,labelDistance,includeGamma,includeLength)

  #********************************************************************************************************************

  #Writes closing lines to .dot file and closes the IO stream
  write(dotIo,"}")
  close(dotIo)
  DEBUG && println("Final lines written and IO stream has been closed")

  #********************************************************************************************************************

  #Converts .dot file into .svg image
  print("Exporting $(imageName).dot to image in $(imageName).svg\n")
  generalExport("$imageName.dot",filename="$imageName",layoutEngine=layoutStyle)

end

function plotPhylonet(netString::AbstractString;
  imageName="netImage",                  # Name for output files
  gammaThreshold=0.5::FloatingPoint,     # threshold for ...
  mainTree=false::Bool,                  # When true, only the underlying tree will be plotted
  width=6::Number,                       # Width of image in inches
  height=8::Number,                      # Height of image in inches
  vert=true::Bool,                       # When true, hierarchy displayed from top to bottom. Otherwise, left to right
  internalLabels=false::Bool,            # When true, all nodes will have labels (including internal nodes)
  fontSize=16.0::FloatingPoint,          # Font size for labels in points
  layoutStyle="dot"::AbstractString,             # layout engine for placing nodes and edges (dot,neato,circo,twopi)
  hybridColor="green4"::AbstractString,          # color for hybrid edges
  forcedLeaf=true::Bool,                 # When true, places all leaf nodes on the same line
  unrooted=false::Bool,                  # if true, enforces layoutStyle to 'neato'
  nodeSeparation=0.8::FloatingPoint,     # minimum distance between nodes in inches. Default was 0.5 before, here only.
  edgeStyle="true"::AbstractString,              # style of edges. Options: "line", "ortho", "composite" (which uses both lines and curved splines), "curved". Default was "false", here only.
  labelAngle= 180.0::FloatingPoint,      # angle for leaf label placement
  labelDistance= 3.0::FloatingPoint,     # distance for leaf label placement
  includeGamma=false::Bool,              # When true, gamma labels are displayed on hybrid edges
  includeLength=false::Bool
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
