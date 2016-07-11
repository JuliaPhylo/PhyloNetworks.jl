"""
    plot(net::HybridNetwork; useEdgeLength=false, mainTree=false, showTipLabel=true,
         showNodeNumber=false, showEdgeLength=false, showGamma=false, edgeColor=colorant"black",
         majorHybridEdgeColor=colorant"deepskyblue4", minorHybridEdgeColor=colorant"deepskyblue",
         showEdgeNumber=false, showIntNodeLabel=false, edgeLabel=[], nodeLabel=[])

Plots a network, from left to right.

- useEdgeLength: if true, the tree edges and major hybrid edges are
  drawn proportionally to their length. Minor hybrid edges are not, however.
  Note that edge lengths in coalescent units may scale very poorly with time.
- mainTree: if true, the minor hybrid edges are ommitted.
- showTipLabel: if true, taxon labels are shown. You may need to zoom out to see them.
- showNodeNumber: if true, nodes are labelled with the number used internally.
- showEdgeLength: if true, edges are labelled with their length (above)
- showGamma: if true, hybrid edges are labelled with their heritability (below)
- edgeColor: color for tree edges. black by default.
- majorHybridEdgeColor: color for major hybrid edges
- minorHybridEdgeColor: color for minor hybrid edges
- showEdgeNumber: if true, edges are labelled with the number used internally.
- showIntNodeLabel: if true, internal nodes are labelled with their names.
  Useful for hybrid nodes, which do have tags like '#H1'.
- edgeLabel: dataframe with two columns: the first with edge numbers, the second with labels
  (like bootstrap values) to annotate edges. empty by default.
- nodeLabel: dataframe with two columns: the first with node numbers, the second with labels
  (like bootstrap values for hybrid relationships) to annotate nodes. empty by default.

Note that `plot` actually modifies some (minor) attributes of the network,
as it calls `directEdges!`, `preorder!` and `cladewiseorder!`.

If hybrid edges cross tree and major edges, you may choose to rotate some tree
edges to eliminate crossing edges, using `rotate!`.
"""
function Gadfly.plot(net::HybridNetwork; useEdgeLength=false::Bool,
        mainTree=false::Bool, showTipLabel=true::Bool, showNodeNumber=false::Bool,
        showEdgeLength=false::Bool, showGamma=false::Bool,
        edgeColor=colorant"black"::ColorTypes.Colorant,
        majorHybridEdgeColor=colorant"deepskyblue4"::ColorTypes.Colorant,
        minorHybridEdgeColor=colorant"deepskyblue"::ColorTypes.Colorant,
        showEdgeNumber=false::Bool, showIntNodeLabel=false::Bool,
        edgeLabel=DataFrame()::DataFrame, nodeLabel=DataFrame()::DataFrame)

    try
        directEdges!(net)   # to update isChild1
    catch e
        if isa(e, RootMismatch)
            e.msg *= "\nPlease change the root, perhaps using rootatnode! or rootatedge!"
        end
        rethrow(e)
    end
    preorder!(net)       # to update net.nodes_changed: true pre-ordering
    cladewiseorder!(net) # to update cladewiseorder_nodeIndex: cladewise on major tree

    !net.node[net.root].leaf ||
        warn("the root is leaf $(net.node[net.root].name): the plot will look weird...")

    # determine y for each node = y of its parent edge: post-order traversal
    # also [yB,yE] for each internal node: range of y's of all children nodes
    ymin = 1.0; ymax = Float64(net.numTaxa);
    node_y  = zeros(Float64, net.numNodes) # order: in net.nodes, *!not in nodes_changed!*
    node_yB = zeros(Float64,net.numNodes) # min (B=begin) and max (E=end)
    node_yE = zeros(Float64,net.numNodes) #   of at children's nodes
    nexty = ymax # first tips at the top, last at bottom
    for i=length(net.node):-1:1 # post-order traversal in major tree
        ni = net.cladewiseorder_nodeIndex[i]
        if net.node[ni].leaf
            node_y[ni] = nexty
            nexty -= 1.0
        else
            node_yB[ni]=ymax; node_yE[ni]=ymin;
            for (e in net.node[ni].edge)
                if net.node[ni] == (e.isChild1 ? e.node[2] : e.node[1]) # if e = child of node
                    if (!e.isMajor) continue; end
                    yy = node_y[getIndex(getOtherNode(e, net.node[ni]), net)]
                    yy!=0 || error("oops, child has not been visited and its y value is 0.")
                    node_yB[ni] = min(node_yB[ni], yy)
                    node_yE[ni] = max(node_yE[ni], yy)
                end
                node_y[ni] = (node_yB[ni]+node_yE[ni])/2
            end
        end
    end

    # setting branch lengths for plotting
    elenCalculate = !useEdgeLength
    if (useEdgeLength)
        allBLmissing = true; nonBLmissing = true;
        for (e in net.edge)
            if (nonBLmissing && e.length==-1.0) nonBLmissing=false; end
            if (allBLmissing && e.length!=-1.0) allBLmissing=false; end
        end
        if (allBLmissing)
            println("All edge lengths are missing, won't be used for plotting.")
            elenCalculate = true
        end
        if (!nonBLmissing && !allBLmissing) # not all, but some are missing
            println("At least one non-missing edge length: plotting any NA length as 1.0")
        end
    end
    elen = Float64[] # edge lengths to be used for plotting. same order as net.edge.
    if (elenCalculate)
        # setting elen such that the age of each node = 1 + age of oldest child
        # (including minor hybrid edges): need true post-ordering.
        # calculating node ages first, elen will be calculated later.
        elen     = zeros(Float64,net.numEdges)
        node_age = zeros(Float64,net.numNodes)
        for i=length(net.node):-1:1 # post-order traversal
            if (net.nodes_changed[i].leaf) continue; end
            ni = getIndex(net.nodes_changed[i], net)
            for (e in net.nodes_changed[i].edge) # loop over children only
                if net.nodes_changed[i] == (e.isChild1 ? e.node[2] : e.node[1])
                    node_age[ni] = max(node_age[ni], 1 +
                     node_age[getIndex(getOtherNode(e, net.nodes_changed[i]), net)])
                end
            end
        end
    else
        for (e in net.edge)
            push!(elen, (e.length==-1.0 ? 1.0 : e.length))
        end
    end

    # determine xB,xE for each edge: pre-order traversal, uses branch lengths
    # then x and yB,yE for each node: x=xE of parent edge
    xmin = 1.0; xmax=xmin
    node_x  = zeros(Float64,net.numNodes) # order: in net.nodes, *!not in nodes_changed!*
    edge_xB = zeros(Float64,net.numEdges) # min (B=begin) and max (E=end)
    edge_xE = zeros(Float64,net.numEdges) # xE-xB = edge length
    edge_yE = zeros(Float64,net.numEdges) # yE of edge = y of child node
    node_x[net.root] = xmin # root node: x=xmin=0
    for i=2:length(net.node)              # true pre-order, skipping the root (i=1)
        ni = getIndex(net.nodes_changed[i], net)
        ei = 0 # index of major parent edge of current node
        for (e in net.nodes_changed[i].edge)
            if (e.isMajor && net.nodes_changed[i] == e.node[e.isChild1 ? 1 : 2]) # major parent edge
                ei = getIndex(e,net)
                break
            end
        end
        ei>0 || error("oops, could not find major parent edge of node number $ni.")
        edge_yE[ei] = node_y[ni]
        pni = getIndex(getOtherNode(net.edge[ei], net.nodes_changed[i]), net) # parent node index
        edge_xB[ei] = node_x[pni]
        if (elenCalculate)
            elen[ei] = node_age[pni] - node_age[ni]
        end
        edge_xE[ei] = edge_xB[ei] + elen[ei]
        node_x[ni] = edge_xE[ei]
        xmax = max(xmax, edge_xE[ei])
    end
    edge_yB = copy(edge_yE) # true for tree and major edges
    for i=1:net.numEdges
        if (!net.edge[i].isMajor) # minor hybrid edges
            cni = getIndex((net.edge[i].isChild1 ? net.edge[i].node[1]:  net.edge[i].node[2]), net)
            pni = getIndex((net.edge[i].isChild1 ? net.edge[i].node[2]:  net.edge[i].node[1]), net)
            # indices of child and parent nodes
            edge_xB[i] = node_x[pni]
            edge_xE[i] = node_x[cni]
            edge_yB[i] = node_y[pni]
            edge_yE[i] = node_y[cni]
            #@show i; @show net.edge[i]; @show pni; @show net.node[pni]; @show cni; @show net.node[cni]
        end
    end

    #@show node_x;  @show node_yB; @show node_y;  @show node_yE
    #@show edge_xB; @show edge_xE; @show edge_yB; @show edge_yE

    mylayers = Layer[] # gadfly layers
    # one layers for each edge
    for i=1:net.numEdges
        if (!mainTree || net.edge[i].isMajor)
            col = edgeColor
            if net.edge[i].hybrid
              if (net.edge[i].isMajor) col = majorHybridEdgeColor;
              else col = minorHybridEdgeColor; end
            end
            push!(mylayers,
              layer(x = [edge_xB[i],edge_xE[i]],
                    y = [edge_yB[i],edge_yE[i]], Geom.line,
                    Theme(default_color=col))[1])
        end
    end
    # one layer for each (vertical) clade
    for i=1:net.numNodes
        if (net.node[i].leaf) continue; end
        push!(mylayers,
              layer(y = [node_yB[i],node_yE[i]],
                    x = [node_x[i], node_x[i]], Geom.line)[1])
    end
    # check data frame for node annotations
    labelnodes = size(nodeLabel,1)>0
    if (labelnodes && (size(nodeLabel,2)<2 || !(eltype(nodeLabel[:,1]) <: Integer)))
        warn("nodeLabel should have 2+ columns, the first one giving the node numbers (Integer)")
        labelnodes = false
    end
    if labelnodes # remove rows with no node number, check if at least one row remains
        nodeLabel = nodeLabel[~isna(nodeLabel[1]),:]
        labelnodes = size(nodeLabel,1)>0
    end
    if labelnodes
      tmp = setdiff(nodeLabel[1], [e.number for e in net.node])
      if length(tmp)>0
        msg = "Some node numbers in the nodeLabel data frame are not found in the network:\n"
        for (a in tmp) msg *= string(" ",a); end
        warn(msg)
      end
    end
    # data frame to place tip names and node annotations (labels)
    if (showTipLabel || showNodeNumber || showIntNodeLabel || labelnodes)
      # white dot beyond tip labels/name and root node label to force enough zoom out
      expfac = 0.1
      push!(mylayers, layer(x=[xmin-(xmax-xmin)*expfac,xmax+(xmax-xmin)*expfac],
                            y=[ymin,ymax+(ymax-ymin)*expfac],
               Geom.point, Theme(default_color=colorant"white"))[1])
      nrows = (showNodeNumber || showIntNodeLabel || labelnodes ? net.numNodes : net.numTaxa)
      ndf = DataFrame([ASCIIString,ASCIIString,ASCIIString,Bool,Float64,Float64], # column types, column names, nrows
               [symbol("name"),symbol("num"),symbol("lab"),symbol("lea"),symbol("x"),symbol("y")], nrows)
      j=1
      for i=1:net.numNodes
        if (net.node[i].leaf  || showNodeNumber || showIntNodeLabel || labelnodes)
            ndf[j,:name] = net.node[i].name
            ndf[j,:num] = string(net.node[i].number)
            if (labelnodes)
              jn = findfirst(nodeLabel[:,1],net.node[i].number)
              ndf[j,:lab] = (jn==0 || isna(nodeLabel[jn,2]) ? "" :  # node label not in table or NA
                (eltype(nodeLabel[:,2])<:AbstractFloat ?
                  @sprintf("%0.3g",nodeLabel[jn,2]) : string(nodeLabel[jn,2])))
            end
            ndf[j,:lea] = net.node[i].leaf # use this later to remove #H? labels
            ndf[j,:y] = node_y[i]
            ndf[j,:x] = node_x[i]
            j += 1
        end
      end
      # @show ndf
      if (showTipLabel)
        push!(mylayers, layer(ndf[ndf[:lea], [:x,:y,:name]], y="y", x="x", label="name",
            Geom.label(position=:right ;hide_overlaps=true))[1])
      end
      if (showIntNodeLabel)
        push!(mylayers, layer(ndf[(!ndf[:lea]), [:x,:y,:name]], y="y", x="x", label="name",
            Geom.label(position=:above ;hide_overlaps=true))[1])
      end
      if (showNodeNumber)
        push!(mylayers, layer(ndf, y="y", x="x", label="num",
            Geom.label(position=:dynamic ;hide_overlaps=true))[1])
      end
      if labelnodes
        push!(mylayers, layer(ndf[:,[:x,:y,:lab]], y="y", x="x", label="lab",
            Geom.label(position=:left ;hide_overlaps=false))[1])
      end
    end
    # data frame for edge annotations.
    nrows = net.numEdges - (mainTree ? net.numHybrids : 0)
    edf = DataFrame([ASCIIString,ASCIIString,ASCIIString,ASCIIString,Bool,Bool,Float64,Float64],
                  [symbol("len"),symbol("gam"),symbol("num"),symbol("lab"),
                   symbol("hyb"),symbol("min"),symbol("x"),symbol("y")], nrows)
    labeledges = size(edgeLabel,1)>0
    if (labeledges && (size(edgeLabel,2)<2 || !(eltype(edgeLabel[:,1]) <: Integer)))
        warn("edgeLabel should have 2+ columns, the first one giving the edge numbers (Integer)")
        labeledges = false
    end
    if labeledges # remove rows with no edge number and check if at least one remains
        edgeLabel = edgeLabel[~isna(edgeLabel[1]),:]
        labeledges = size(edgeLabel,1)>0
    end
    if labeledges
      tmp = setdiff(edgeLabel[1], [e.number for e in net.edge])
      if length(tmp)>0
        msg = "Some edge numbers in the edgeLabel data frame are not found in the network:\n"
        for (a in tmp) msg *= string(" ",a); end
        warn(msg)
      end
    end
    j=1
    for (i = 1:length(net.edge))
        if (!mainTree || !net.edge[i].hybrid || net.edge[i].isMajor)
            edf[j,:len] = (net.edge[i].length==-1.0 ? "" : @sprintf("%0.3g",net.edge[i].length))
            # @sprintf("%c=%0.3g",'γ',net.edge[i].length)
            edf[j,:gam] = (net.edge[i].gamma==-1.0  ? "" : @sprintf("%0.3g",net.edge[i].gamma))
            edf[j,:num] = string(net.edge[i].number)
            if (labeledges)
              je = findfirst(edgeLabel[:,1],net.edge[i].number)
              edf[j,:lab] = (je==0 || isna(edgeLabel[je,2]) ? "" :  # edge label not found in table
                (eltype(edgeLabel[:,2])<:AbstractFloat ?
                  @sprintf("%0.3g",edgeLabel[je,2]) : string(edgeLabel[je,2])))
            end
            edf[j,:hyb] = net.edge[i].hybrid
            edf[j,:min] = !net.edge[i].isMajor
            edf[j,:y] = (edge_yB[i] + edge_yE[i])/2
            edf[j,:x] = (edge_xB[i] + edge_xE[i])/2
            j += 1
        end
    end
    # @show edf
        if labeledges
            push!(mylayers, layer(edf[:,[:x,:y,:lab]], y="y", x="x", label="lab",
                  Geom.label(position=:above ;hide_overlaps=false))[1])
        end
        if (showEdgeLength)
            push!(mylayers, layer(edf[:,[:x,:y,:len]], y="y", x="x", label="len",
                  Geom.label(position=:below ;hide_overlaps=false))[1])
        end
        if (showGamma && net.numHybrids>0)
            if !mainTree
            push!(mylayers, layer(edf[(edf[:hyb]) & (edf[:min]), [:x,:y,:gam]], y="y", x="x",label="gam",
                  Geom.label(position=:below ;hide_overlaps=true),
                  Theme(point_label_color=minorHybridEdgeColor))[1])
            end
            push!(mylayers, layer(edf[(edf[:hyb]) & (!edf[:min]),[:x,:y,:gam]], y="y", x="x",label="gam",
                  Geom.label(position=:below ;hide_overlaps=true),
                  Theme(point_label_color=majorHybridEdgeColor))[1])
        end
        if (showEdgeNumber)
            push!(mylayers, layer(edf[:,[:x,:y,:num]], y="y", x="x", label="num",
                  Geom.label(position=:dynamic ;hide_overlaps=false))[1])
        end

    plot(mylayers, Guide.xlabel("time"), Guide.ylabel(nothing),
         Guide.xticks(ticks=:auto, label=false), # ticks=[xmin,xmin,xmax,xmax*1.1],
         Guide.yticks(ticks=:auto, label=false), # ticks=[ymin,ymax],
         Theme(default_color=edgeColor,grid_color=colorant"white"))
end
