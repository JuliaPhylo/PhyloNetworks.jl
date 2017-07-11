"""
`plot(net::HybridNetwork, method::Symbol)`

Plot a network using R graphics.
`method` should be `:RCall` or `:ape`.

optional arguments, shared with the Gadfly-based plot function:
- useEdgeLength: if true, the tree edges and major hybrid edges are
  drawn proportionally to their length. Minor hybrid edges are not, however.
  Note that edge lengths in coalescent units may scale very poorly with time.
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

optional arguments specific to this function:
- xlim, ylim: array of 2 values
- tipOffset: to offset tip labels

Note that `plot` actually modifies some (minor) attributes of the network,
as it calls `directEdges!`, `preorder!` and `cladewiseorder!`.

If hybrid edges cross tree and major edges, you may choose to rotate some tree
edges to eliminate crossing edges, using `rotate!`.
"""

function plot(net::HybridNetwork, method::Symbol; useEdgeLength=false::Bool,
    mainTree=false::Bool, showTipLabel=true::Bool, showNodeNumber=false::Bool,
    showEdgeLength=false::Bool, showGamma=false::Bool,
    edgeColor="black"::String,
    majorHybridEdgeColor="deepskyblue4"::String,
    minorHybridEdgeColor="deepskyblue"::String,
    showEdgeNumber=false::Bool, showIntNodeLabel=false::Bool,
    edgeLabel=DataFrame()::DataFrame, nodeLabel=DataFrame()::DataFrame,
    xlim=Float64[]::Array{Float64,1}, ylim=Float64[]::Array{Float64,1},
    tipOffset=0.0::Float64)

    (edge_xB, edge_xE, edge_yB, edge_yE, node_x, node_y, node_yB, node_yE,
     xmin, xmax, ymin, ymax) = getEdgeNodeCoordinates(net, useEdgeLength)
    labelnodes, nodeLabel = checkNodeDataFrame(net, nodeLabel)
    ndf = prepareNodeDataFrame(net, nodeLabel, showNodeNumber,
            showIntNodeLabel, labelnodes, node_x, node_y)
    if (showTipLabel || showNodeNumber || showIntNodeLabel || labelnodes)
        expfac = 0.1  # force 10% more space to show tip/node/root name
        expfacy = 0.5 # additive expansion for y axis
        xmin -= (xmax-xmin)*expfac
        xmax += (xmax-xmin)*expfac
        ymin -= expfacy
        ymax += expfacy
    end
    xmax += tipOffset
    if length(xlim)==2
        xmin=xlim[1]; xmax=xlim[2]
    end
    if length(ylim)==2
        ymin=ylim[1]; xmax=ylim[2]
    end
    leaves = [n.leaf for n in net.node]
    eCol = fill(edgeColor, length(net.edge))
    eCol[ [ e.hybrid  for e in net.edge] ] = majorHybridEdgeColor
    eCol[ [!e.isMajor for e in net.edge] ] = minorHybridEdgeColor
    
    R"""
    plot($(node_x[leaves]), $(node_y[leaves]), type='n',
         xlim=c($xmin,$xmax), ylim=c($ymin,$ymax),
         axes=FALSE, xlab='', ylab='')
    segments($edge_xB, $edge_yB, $edge_xE, $edge_yE, col=$eCol)
    segments($node_x, $node_yB, $node_x, $node_yE, col=$edgeColor)
    """
    if showTipLabel
    R"""
    text($(node_x[leaves])+$tipOffset, $(node_y[leaves]),
         $(tipLabels(net)), adj=0, font=3)
    """
    end
    if showIntNodeLabel
    R"""
    text($(ndf[!ndf[:lea],:x]), $(ndf[!ndf[:lea],:y]),
         $(ndf[!ndf[:lea],:name]), adj=c(.5,0))
    """
    end
    if showNodeNumber
    R"""
    text($(ndf[:x]), $(ndf[:y]), $(ndf[:num]), adj=1)
    """
    end
    if labelnodes
    R"""
    text($(ndf[:x]), $(ndf[:y]), $(ndf[:lab]), adj=1)
    """
    end
    labeledges, edf = prepareEdgeDataFrame(net, edgeLabel, mainTree,
                        edge_xB, edge_xE, edge_yB, edge_yE)
    if labeledges
    R"""
    text($(edf[:x]), $(edf[:y]), $(edf[:lab]), adj=c(.5,0))
    """
    end
    if showEdgeLength
    R"""
    text($(edf[:x]), $(edf[:y]), $(edf[:len]), adj=c(.5,1))
    """
    end
    if (showGamma && net.numHybrids>0)
    im = edf[:hyb] .& edf[:min]
    iM = edf[:hyb] .& .!edf[:min]
    R"""
    text($(edf[im,:x]), $(edf[im,:y]),$(edf[im,:gam]),
         adj=c(.5,1), col=$minorHybridEdgeColor)
    text($(edf[iM,:x]), $(edf[iM,:y]),$(edf[iM,:gam]),
         adj=c(.5,1), col=$majorHybridEdgeColor)
    """
    end
    if showEdgeNumber
    R"""
    text($(edf[:x]), $(edf[:y]), $(edf[:num]), adj=c(.5,0))
    """
    end
    return (xmin, xmax, ymin, ymax, node_x, node_y, node_yB, node_yE,
      edge_xB, edge_xE, edge_yB, edge_yE, ndf, edf)
end
