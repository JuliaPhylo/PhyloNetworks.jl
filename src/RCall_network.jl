using PhyloNetworks, RCall
import PhyloNetworks.getOtherNode
import PhyloNetworks.getIndex

net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")


"""
generate taxa labels for network
"""

function generateLeaves(net::HybridNetwork, tipLabels::Array{String}, leaves::Array{Bool})
    for n in net.node
        if n.leaf == true
            push!(tipLabels, n.name)
            push!(leaves, true)
        else
            push!(tipLabels, "")
            push!(leaves, false)
        end
    end
end

"""
color hybrid edges based on inheritance values
"""

function netColors(net::HybridNetwork, edgeColors::Array{String})    
    for e in net.edge
        if e.hybrid == true && e.gamma >= 0.5
            push!(edgeColors, "blue")
        elseif e.hybrid == true && e.gamma < 0.5
            push!(edgeColors, "lightblue")
        else
            push!(edgeColors, "black")
        end
    end
end

"""
generate gamma inheritance value labels for hybrid edges
"""

function generateGamma(net::HybridNetwork, gammaLabel::Array{Any})
    for e in net.edge
        if e.hybrid == true
            push!(gammaLabel, e.gamma)
        else
            push!(gammaLabel, "")
        end
    end
end

"""
plot hybrid network using r graphics
"""

function RPlotNetworks(net::HybridNetwork)
    
    directEdges!(net)
    preorder!(net)
    cladewiseorder!(net)
    
    tipLabels = Array{String}(0)
    leaves = Array{Bool}(0)
    generateLeaves(net, tipLabels, leaves)
    
    edgeColors = Array{String}(0)
    netColors(net, edgeColors)
    
    gammaLabel = Array{Any}(0)
    generateGamma(net, gammaLabel)
    
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
            for e in net.node[ni].edge
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
    
    useEdgeLength=false
    # setting branch lengths for plotting
    elenCalculate = !useEdgeLength
    if (useEdgeLength)
        allBLmissing = true; nonBLmissing = true;
        for e in net.edge
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
        elen = zeros(Float64,net.numEdges)
        node_age = zeros(Float64,net.numNodes)
        for i=length(net.node):-1:1 # post-order traversal
            if (net.nodes_changed[i].leaf) continue; end
            ni = getIndex(net.nodes_changed[i], net)
            for e in net.nodes_changed[i].edge # loop over children only
                if net.nodes_changed[i] == (e.isChild1 ? e.node[2] : e.node[1])
                    node_age[ni] = max(node_age[ni], 1 +
                    node_age[getIndex(getOtherNode(e, net.nodes_changed[i]), net)])
                end
            end
        end
    else
        for e in net.edge
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
        for e in net.nodes_changed[i].edge
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
    
    R"""
    plot($node_x[$leaves], $node_y[$leaves], pch=46, xlim=c(1,5), ylim=c(1,4), axes=FALSE, xlab='', ylab='', text($node_x, $node_y, $tipLabels, pos=4, offset=0.0))
    segments($edge_xB, $edge_yB, $edge_xE, $edge_yE, lwd=2, col=$edgeColors, text(($edge_xB+$edge_xE)/2, ($edge_yB+$edge_yE)/2, $gammaLabel, pos=3, offset=0.5))
    segments($node_x, $node_yB, $node_x, $node_yE)
    """
end

RPlotNetworks(net)
