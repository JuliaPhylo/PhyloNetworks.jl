import Gadfly.plot

function Gadfly.plot(net::HybridNetwork)
    preorder!(net)
    # fixit: add option checkPreorder to bypass the preordering.
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
                    # fixit: what if e = minor hybrid edge: how to draw minor edges?
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
    allBLmissing = true
    for (e in net.edge)
        if (e.length != -1.0)
            println("At least one non-missing edge length: plotting any NA length as 1.0")
            allBLmissing = false
            break
        end
    end
    elen = Float64[] # edge lengths to be used for plotting. same order as net.edge.
    if (allBLmissing)
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
    edge_y  = zeros(Float64,net.numEdges) # y of edge = y of child node
    node_x[net.root] = xmin # root node: x=xmin=0
    for i=2:length(net.nodes_changed)  # pre-order traversal, skipping the root (i=1)
        ei = net.preorder_edgeIndex[i] # major parent edge of current node
        edge_y[ei] = node_y[net.preorder_nodeIndex[i]]
        pni = getIndex(getOtherNode(net.edge[ei], net.nodes_changed[i]), net) # parent node index
        edge_xB[ei] = node_x[pni]
        if (allBLmissing)
            elen[ei] = node_age[pni] - node_age[net.preorder_nodeIndex[i]]
        end
        edge_xE[ei] = edge_xB[ei] + elen[ei]
        node_x[net.preorder_nodeIndex[i]] = edge_xE[ei]
        xmax = max(xmax, edge_xE[ei])
    end
    # @show node_yE;  @show node_y; @show node_yB

    edgecol = colorant"black"
    mylayers = Layer[] # gadfly layers
    # one layers for each edge
    for i=1:net.numEdges
        col = edgecol
        if (!net.edge[i].isMajor) # minor hybrid edge
            col = colorant"green"
            cni = getIndex((net.edge[i].isChild1 ? net.edge[i].node[1]:  net.edge[i].node[2]), net)
            pni = getIndex((net.edge[i].isChild1 ? net.edge[i].node[2]:  net.edge[i].node[1]), net)
            # indices of child and parent nodes
            xE
            push!(mylayers,
                  layer(x = [node_x[pni],node_x[cni]],
                        y = [node_y[pni],node_y[cni]], Geom.line,
                        Theme(default_color=col))[1])
            continue;
        end
        if net.edge[i].hybrid col = colorant"blue"; end
        push!(mylayers,
              layer(x = [edge_xB[i],edge_xE[i]],
                    y = [edge_y[i], edge_y[i]], Geom.line,
                    Theme(default_color=col))[1])
    end
    # one layer for each (vertical) clade
    for i=1:net.numNodes
        if (net.node[i].leaf) continue; end
        push!(mylayers,
              layer(y = [node_yB[i],node_yE[i]],
                    x = [node_x[i], node_x[i]], Geom.line)[1])
    end
    # data frame to place tip labels
    tipDF = DataFrame([ASCIIString,Float64,Float64], # column types, column names, nrows
                      [symbol("lab"),symbol("x"),symbol("y")],net.numTaxa)
    j=1
    for i=1:net.numNodes
        if (net.node[i].leaf)
            tipDF[j,:lab] = net.node[i].name
            tipDF[j,:y] = node_y[i]
            tipDF[j,:x] = node_x[i]
            j += 1
        end
    end
    push!(mylayers,
         layer(tipDF, y="y", x="x", label="lab",
            Geom.label(position=:right ;hide_overlaps=false))[1])

    plot(mylayers, Guide.xlabel("time"), Guide.ylabel(nothing),
         Guide.xticks(ticks=[xmin,xmax,xmax*1.1], label=false),
         Guide.yticks(ticks=[ymin,ymax], label=false),
         Theme(default_color=edgecol))
end
