function Gadfly.plot(net::HybridNetwork; useEdgeLength=false::Bool,
        mainTree=false::Bool, showTipLabel=true::Bool, showNodeNumber=false::Bool,
        showEdgeLength=true::Bool, showGamma=true::Bool,
        edgeColor=colorant"black"::ColorTypes.Colorant,
        majorHybridEdgeColor=colorant"blue1"::ColorTypes.Colorant,
        minorHybridEdgeColor=colorant"blue1"::ColorTypes.Colorant,
        preorder=true::Bool)

    if (preorder)
        directEdges!(net)
        preorder!(net)
    else
       (length(net.node)==length(net.nodes_changed)      &&
        length(net.node)==length(net.preorder_nodeIndex) &&
        length(net.node)==length(net.preorder_edgeIndex)    ) ||
       error("please run preorder! first.")
    end
    # determine y for each node = y of its parent edge: post-order traversal
    # also [yB,yE] for each internal node: range of y's of all children nodes
    ymin = 1.0; ymax = Float64(net.numTaxa);
    node_y  = zeros(Float64, net.numNodes) # order: in net.nodes, *!not in nodes_changed!*
    node_yB = zeros(Float64,net.numNodes) # min (B=begin) and max (E=end)
    node_yE = zeros(Float64,net.numNodes) #   of at children's nodes
    nexty = ymin
    for i=length(net.nodes_changed):-1:1 # post-order traversal
        ni = net.preorder_nodeIndex[i]
        if net.nodes_changed[i].leaf
            node_y[ni] = nexty
            nexty += 1.0
        else
            node_yB[ni]=ymax; node_yE[ni]=ymin;
            for (e in net.nodes_changed[i].edge)
                if net.nodes_changed[i] == (e.isChild1 ? e.node[2] : e.node[1]) # if e = child of node
                    if (!e.isMajor) continue; end
                    yy = node_y[getIndex(getOtherNode(e, net.nodes_changed[i]), net)]
                    if (yy==0) println("oops, error with yy=0"); end
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
        if (!nonBLmissing)
            println("At least one non-missing edge length: plotting any NA length as 1.0")
        end
    end
    elen = Float64[] # edge lengths to be used for plotting. same order as net.edge.
    if (elenCalculate)
        # setting elen such that the age of each node = 1 + age of oldest child
        # calculating node ages first, elen will be calculated later.
        elen     = zeros(Float64,net.numEdges)
        node_age = zeros(Float64,net.numNodes)
        for i=length(net.nodes_changed):-1:1 # post-order traversal
            if (net.nodes_changed[i].leaf) continue; end
            ni = net.preorder_nodeIndex[i]
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
    xmin = 0.0; xmax=xmin
    node_x  = zeros(Float64,net.numNodes) # order: in net.nodes, *!not in nodes_changed!*
    edge_xB = zeros(Float64,net.numEdges) # min (B=begin) and max (E=end)
    edge_xE = zeros(Float64,net.numEdges) # xE-xB = edge length
    edge_yE = zeros(Float64,net.numEdges) # yE of edge = y of child node
    node_x[net.root] = xmin # root node: x=xmin=0
    for i=2:length(net.nodes_changed)  # pre-order traversal, skipping the root (i=1)
        ei = net.preorder_edgeIndex[i] # major parent edge of current node
        edge_yE[ei] = node_y[net.preorder_nodeIndex[i]]
        pni = getIndex(getOtherNode(net.edge[ei], net.nodes_changed[i]), net) # parent node index
        edge_xB[ei] = node_x[pni]
        if (elenCalculate)
            elen[ei] = node_age[pni] - node_age[net.preorder_nodeIndex[i]]
        end
        edge_xE[ei] = edge_xB[ei] + elen[ei]
        node_x[net.preorder_nodeIndex[i]] = edge_xE[ei]
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
    # data frame to place tip labels
    if (showTipLabel || showNodeNumber)
      # white dot beyond tip labels to force enough zoom out
      push!(mylayers, layer(x=[xmax*1.1], y=[ymax*1.1],
               Geom.point, Theme(default_color=colorant"white"))[1])
      nrows = (showNodeNumber ? net.numNodes : net.numTaxa)
      ndf = DataFrame([ASCIIString,ASCIIString,Bool,Float64,Float64], # column types, column names, nrows
               [symbol("lab"),symbol("num"),symbol("lea"),symbol("x"),symbol("y")], nrows)
      j=1
      for i=1:net.numNodes
        if (net.node[i].leaf  || showNodeNumber)
            ndf[j,:lab] = net.node[i].name
            ndf[j,:num] = string(net.node[i].number)
            ndf[j,:lea] = net.node[i].leaf # use this later to remove #H? labels
            ndf[j,:y] = node_y[i]
            ndf[j,:x] = node_x[i]
            j += 1
        end
      end
      if (showTipLabel)
        push!(mylayers, layer(ndf, y="y", x="x", label="lab",
            Geom.label(position=:right ;hide_overlaps=true))[1])
      end
      if (showNodeNumber)
        push!(mylayers, layer(ndf, y="y", x="x", label="num",
            Geom.label(position=:dynamic ;hide_overlaps=true))[1])
      end
    end
    # data frame for edge annotations
    if (showEdgeLength || showGamma)
        edf = DataFrame([ASCIIString,ASCIIString,ASCIIString,Bool,Bool,Float64,Float64],
                  [symbol("len"),symbol("gam"),symbol("num"),
                   symbol("hyb"),symbol("min"),symbol("x"),symbol("y")], net.numEdges)

        for (i = 1:length(net.edge))
            edf[i,:len] = string(net.edge[i].length)
            edf[i,:gam] = string(net.edge[i].gamma)
            edf[i,:num] = string(net.edge[i].number)
            edf[i,:hyb] = net.edge[i].hybrid
            edf[i,:min] = !net.edge[i].isMajor
            edf[i,:y] = (edge_yB[i] + edge_yE[i])/2
            edf[i,:x] = (edge_xB[i] + edge_xE[i])/2
        end
        if (showEdgeLength)
            push!(mylayers, layer(edf, y="y", x="x", label="len",
                  Geom.label(position=:above ;hide_overlaps=true))[1])
        end
        if (showGamma)
            push!(mylayers, layer(edf[edf[:hyb],:], y="y", x="x", label="gam",
                  Geom.label(position=:below ;hide_overlaps=true),
                  Theme(point_label_color=minorHybridEdgeColor))[1])
        end
        #@show edf
    end

    plot(mylayers, Guide.xlabel("time"), Guide.ylabel(nothing),
         Guide.xticks(ticks=:auto, label=false), # ticks=[xmin,xmin,xmax,xmax*1.1],
         Guide.yticks(ticks=:auto, label=false), # ticks=[ymin,ymax],
         Theme(default_color=edgeColor,grid_color=colorant"white"))
end
