net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")

"""
change numbers of internal nodes to predictable, non-negative integers
"""

function phyloNumbers(tree::HybridNetwork)  
    lnum = 1
    inum = length(tree.leaf) + 1 #ensures internal node numbers do not overlap with leaf numbers
    #root will be labeled ntips +1, consistent with ape "phylo" class
    for n in tree.nodes_changed #traverse tree in topological order excluding leaves
        if n.leaf
            n.number = lnum
            lnum += 1
        else
            n.number = inum #change node number
            inum += 1
        end
    end
end

"""
colors hybrid edges based on inheritance values
"""

function treeColors(net::HybridNetwork, tree::HybridNetwork, edgeColors::Array{String,1})
    for n in net.nodes_changed #traverse tree in topological order excluding leaves
        if !n.leaf
            for e in n.edge #iterate over each edge attatched to a node
                if e.node[e.isChild1 ? 1 : 2] != n #exclude parent edges
                    if e.hybrid == true && e.gamma >= 0.5 && issubset(e.node[e.isChild1 ? 2 : 1].number, [n.number for n in tree.node]) == true
                        push!(edgeColors, "blue")
                    elseif e.hybrid == true && e.gamma < 0.5 && issubset(e.node[e.isChild1 ? 2 : 1].number, [n.number for n in tree.node]) == true
                        push!(edgeColors, "lightblue")
                    elseif issubset(e.node[e.isChild1 ? 2 : 1].number, [n.number for n in tree.node]) == true
                        push!(edgeColors, "black")
                    else 
                        continue
                    end
                end
            end
        end
    end
    return edgeColors
end

"""
generates edge matrix from parent and child nodes
"""

function generateEdge(tree::HybridNetwork) 
    
    edge = Matrix{Int}(length(tree.edge), 2) #initialize edge matrix
    
    #fill edge matrix
    i = 1 #row index for edge matrix
    for n in tree.nodes_changed #traverse tree in topological order excluding leaves
        if !n.leaf
            for e in n.edge #iterate over each edge attatched to a node
                if e.node[e.isChild1 ? 1 : 2] != n #exclude node that is parent of current edge
                    edge[i,1] = e.node[e.isChild1 ? 2 : 1].number #push parent node to first column
                    edge[i,2] = e.node[e.isChild1 ? 1 : 2].number #push child node to second column
                    i += 1 #increase index value
                end
            end
        end
    end
    return edge
end

"""
plots Julia PhyloNetworks::HybridNetwork object as R S3 phylo object
""" 

function RPlot(net::HybridNetwork)
    
    directEdges!(net)
    preorder!(net) #organize nodes for preorder traversal
    
    #tree = majorTree(net)
    trees = displayedTrees(net, 0.0)
    tree = trees[2]
    
    preorder!(tree) #organize nodes for preorder traversal
    
    edgeColors = Array{String}(0)
    
    treeColors(net, tree, edgeColors) #fill edgeColors vector
    
    phyloNumbers(tree)
    ntips = length(tree.leaf) #number of leaves for Nnode
    totalnodes = length(tree.node) #total number of nodes for edge matrix and Nnode
    Nnode = totalnodes - ntips #integer for number of internal nodes

    tipLabel = [node.name for node in tree.leaf] #tiplabel vector

    edge = generateEdge(tree)

    #push Nnode integer, edge matrix, tiplabel vector to RObject $phy and assign "phylo" class attribute and plot
    R"""
    library(ape)
    phy = list(Nnode = $Nnode, edge = $edge, tip.label = $tipLabel)
    class(phy) = "phylo"
    plot(phy, edge.color = $edgeColors)
    """
end

RPlot(net)

