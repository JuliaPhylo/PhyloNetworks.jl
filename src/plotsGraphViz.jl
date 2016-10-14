# initial function: John Spaw, 7/2015
# takes a network object and converts it into .dot file, which is then converted into a .svg image file
"""
`plotNetGraphViz(net::HybridNetwork)`

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
function plotNetGraphViz(graph::Network; # Network object to plot
  imageName="netImage",                  # Name for output files
  gammaThreshold=0.5::AbstractFloat,     # threshold for ...
  mainTree=false::Bool,                  # When true, plots underlying tree only
  width=6::Number,                       # Width of image in inches
  height=8::Number,                      # Height of image in inches
  vert=true::Bool,                       # top to bottom if true, left to right if false
  internalLabels=false::Bool,            # When true, all nodes labeled (internal nodes too)
  fontSize=16.0::AbstractFloat,          # Font size for labels in points
  layoutStyle="dot"::AbstractString,     # layout engine (dot,neato,circo,twopi)
  hybridColor="green4"::AbstractString,  # color for hybrid edges
  forcedLeaf=true::Bool,                 # When true, places all leaf nodes on the same line
  unrooted=false::Bool,                  # if true, enforces layoutStyle to 'neato'
  nodeSeparation=0.8::AbstractFloat,     # min distance between nodes in inches
  edgeStyle="true"::AbstractString,      # options: "line", "ortho", "composite" (uses both lines and curved splines), "curved"
  labelAngle= 180.0::AbstractFloat,      # angle for leaf label placement
  labelDistance= 3.0::AbstractFloat,     # distance for leaf label placement
  includeGamma=false::Bool,              # When true, gamma labels are displayed on hybrid edges
  includeLength=false::Bool
  )

  # check that the root is correcly placed (Claudia, Aug2015)
  if(!isTree(graph))
      if(!graph.cleaned)
          DEBUG && println("net not cleaned inside plotNetGraphViz, need to run updateCR")
          for n in graph.hybrid
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
    end
  end

  #Create array of leaf nodes for use in ranking and label creation
  for i in graph.leaf
    #if i != graph.node[graph.root]
     push!(leafNodes,i.number)
    #end
  end

  rootNode = graph.node[graph.root]

  #*******************************************
  #Writes initial preample lines for .dot file
  write(dotIo,"Graph Phylogeny { \n")
  if vert == false
    write(dotIo,"     rankdir=LR; \n")
  end
  write(dotIo,"    labelloc=b \n")         # avoid overlapped labels (DOUBLE CHECK THIS)
  write(dotIo,"    ratio=\"fill\"; \n")    # Fits graph to full image size (TEST OTHER RATIO OPTIONS)
  write(dotIo,"    size=\"$width ,$height !\"; \n") # size of entire graph
  write(dotIo,"    nodesep=$(nodeSeparation); \n")  # Min distance between any two nodes
  write(dotIo,"    splines=$(edgeStyle); \n")
  write(dotIo,"    edge [fontsize=$fontSize]; \n")
  write(dotIo,"    node [shape = point] \n")
  write(dotIo,"    rank=max $(rootNode.number) \n") # to place root at top of hierarchy (dot engine only)

  if unrooted
    write(dotIo,"  mode = KK; \n")
    write(dotIo,"  epsilon = 0.01; \n")
  end

  #Ranking leaves together at the bottom of the tree
  if forcedLeaf
    write(dotIo,"  subgraph    { ")    #Places root node at top of tree
    leafArraySize = countnz(leafNodes) #Probably redundant... could delete later on
    for i in leafNodes                 #Groups leaves so they are all placed at bottom of tree
      if i != leafNodes[leafArraySize] # appends each leaf (except last) to .dot file followed by a comma
        write(dotIo,"$i , ")
      else                             # appends final leaf WITHOUT a comma
        write(dotIo,"$i")
      end
    end
    write(dotIo," } \n")
  end

  #Choosing which nodes have labels
  if internalLabels
    for i in allNodes
      write(dotIo,"    $i [xlabel=$i] [shape = point] \n") # node labels, internal and leaves
    end
  else
    for i in leafNodes
      write(dotIo,"    $i [label=$i] [shape = point] \n")  # labels only leaf nodes... this will be removed once xlabels have been fixed
    end
    write(dotIo," \n")
  end

  DEBUG && println("Preamble written!")

  #************************************
  #Drawing edges and nodes

  DEBUG && println("Drawing nodes and edges")

  dummy = Edge(-1,1.0) # to begin recursion because root has no parent edge

  #Traverse the network using a pseudo-depth-first method
  writeEdgesToDotFile(graph, rootNode, mainTree, dotIo, gammaThreshold, dummy,hybridColor,layoutStyle,labelAngle,labelDistance,includeGamma,includeLength)

  #************************************
  #Writes closing lines to .dot file and closes the IO stream
  write(dotIo,"}")
  close(dotIo)
  DEBUG && println("Final lines written and IO stream has been closed")

  #************************************
  #Converts .dot file into .svg image
  print("Exporting $(imageName).dot to image in $(imageName).svg\n")
  generalExport("$imageName.dot",filename="$imageName",layoutEngine=layoutStyle)

end

function plotNetGraphViz(netString::AbstractString;
  imageName="netImage",                 # Name for output files
  gammaThreshold=0.5::AbstractFloat,    # threshold for ...
  mainTree=false::Bool,                 # When true, plots underlying tree only
  width=6::Number,                      # Width of image in inches
  height=8::Number,                     # Height of image in inches
  vert=true::Bool,                      # top to bottom if true, left to right if false
  internalLabels=false::Bool,           # When true, all nodes labeled (internal nodes too)
  fontSize=16.0::AbstractFloat,         # Font size for labels in points
  layoutStyle="dot"::AbstractString,    # layout engine (dot,neato,circo,twopi)
  hybridColor="green4"::AbstractString, # color for hybrid edges
  forcedLeaf=true::Bool,                # When true, places all leaf nodes on the same line
  unrooted=false::Bool,                 # if true, enforces layoutStyle to 'neato'
  nodeSeparation=0.8::AbstractFloat,    # min distance between nodes, inches. Default was 0.5 before.
  edgeStyle="true"::AbstractString,     # options: "line", "ortho", "composite" (uses both lines and curved splines), "curved". Default was "false", here only.
  labelAngle= 180.0::AbstractFloat,     # angle for leaf label placement
  labelDistance= 3.0::AbstractFloat,    # distance for leaf label placement
  includeGamma=false::Bool,             # When true, gamma labels are displayed on hybrid edges
  includeLength=false::Bool
  )

  net = readTopologyUpdate(netString,true);
  plotNetGraphViz(net,
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



# recursive function that traverses each edge only once
# and appends each edge's info to the dot file
function writeEdgesToDotFile(net::HybridNetwork,
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
    # println("root node number:"); @show node.number
    for edge in node.edge
      child = (edge.node[1]==node ? edge.node[2]  : edge.node[1])
      gma = round(edge.gamma,3)  # Gamma for edge i (default 1.0 for tree edges)
      hThickness = gma*4
      edgeNum = edge.number      # for easier use with $ syntax
      node1Num = node.number     # parent node number. for $ notation
      node2Num = child.number    # child  node number. for $ notation

      if child.hybrid
          # thicker edge AND gamma label on edge to dominant edge
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
          elseif mainTree == false  # Skip this step if you only want mainTree
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

      writeEdgesToDotFile(net,child,mainTree,dotIo,gammaThreshold,edge,hybridColor,layoutStyle,labelAngle,labelDistance,includeGamma,includeLength)
    end # loop over edges
  end # if node is root

  # *********************************************

  # Case 2: Node is a leaf
  if node.leaf                                  # do nothing
    #println("Node $(node.number) is a leaf")   # better statement to use here... pass/break?
  end

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
      if (((~(edge.hybrid)) || ((edge.hybrid) && edge.isMajor)) && edge != parentEdge)  # recursive call if edge is a tree edge or a major hybrid
          if edge.gamma > 0.5 # fixit: would be best to remove this, instead check that
            writeEdgesToDotFile(net,child,mainTree,dotIo,gammaThreshold,edge,hybridColor,layoutStyle, labelAngle, labelDistance,includeGamma,includeLength)
          end #if gamma
        end #if
      end #if
    end #for
  end #if
end #function

# Converts a .dot file into a .svg image
function generalExport(file;filename="genImage"::AbstractString,layoutEngine="dot")
  dot = open(file,"r") do io Graph(io) end
  GraphViz.layout!(dot,engine=layoutEngine)
  open("$filename.svg","w") do f
    GraphViz.writemime(f, MIME"image/svg+xml"(),dot)
  end #do
  DEBUG && print("File saved")
end
