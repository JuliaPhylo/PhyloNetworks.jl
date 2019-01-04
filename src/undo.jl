# functions to undo incycle, containRoot, gammaz
# originally in functions.jl
# Claudia March 2015


# ---------------------------------------- undo update of new hybridization --------------------------------


# function to undo updateInCycle which returns an array
# of edges/nodes changed
function undoInCycle!(edges::Array{Edge,1},nodes::Array{Node,1})
    for e in edges
        e.inCycle = -1;
    end
    for n in nodes
        n.inCycle = -1;
    end
end

# function to undo updateContainRoot (which returns an array
# of edges changed) and gives value bool, by default true
function undoContainRoot!(edges::Array{Edge,1}, bool::Bool)
    for e in edges
        e.containRoot = bool
    end
end

undoContainRoot!(edges::Vector{Edge}) = undoContainRoot!(edges,true)

# function to undo updateGammaz which returns an array
# of edges changed
# it only changes the status of istIdentifiable to true
function undoistIdentifiable!(edges::Array{Edge,1})
    for e in edges
        !e.istIdentifiable ? e.istIdentifiable = true : e.istIdentifiable = false;
    end
end


"""
    undoGammaz!(node, network)

Undo `updateGammaz!` for the 2 cases: bad diamond I,II.
`node` should be a hybrid node.
Set length to edges that were not identifiable and
change edges' `gammaz` attribute to -1.0.
Recalculate branch lengths in terms of `gammaz`.  
*warning*: needs to know `incycle` attributes
"""
function undoGammaz!(node::Node, net::HybridNetwork)
    node.hybrid || error("cannot undo gammaz if starting node is not hybrid")
    if(node.isBadDiamondI)
        edge_maj, edge_min, tree_edge2 = hybridEdges(node);
        other_maj = getOtherNode(edge_maj,node);
        other_min = getOtherNode(edge_min,node);
        edgebla,tree_edge_incycle1,tree_edge = hybridEdges(other_min);
        edgebla,tree_edge_incycle2,tree_edge = hybridEdges(other_maj);
        other_min.gammaz != -1 || error("bad diamond I in node $(node.number) but no gammaz updated correctly")
        setLength!(tree_edge_incycle1,-log(1-other_min.gammaz))
        other_maj.gammaz != -1 || error("bad diamond I in node $(node.number) but no gammaz updated correctly")
        setLength!(tree_edge_incycle2,-log(1-other_maj.gammaz))
        if approxEq(other_maj.gammaz,0.0) && approxEq(other_min.gammaz,0.0)
            setGamma!(edge_maj,0.0, true) # gamma could be anything if both gammaz are 0.0, but will set to 0.0
            setLength!(edge_maj,0.0)
            setLength!(edge_min,0.0)
        else
            setGamma!(edge_maj,other_maj.gammaz / (other_maj.gammaz+other_min.gammaz), true)
        end
        other_min.gammaz = -1.0
        other_maj.gammaz = -1.0
        tree_edge_incycle1.istIdentifiable = true;
        tree_edge_incycle2.istIdentifiable = true;
        edge_maj.istIdentifiable = true;
        edge_min.istIdentifiable = true;
        node.isBadDiamondI = false
        net.numBad -= 1
    elseif(node.isBadDiamondII)
        edge_maj, edge_min, tree_edge2 = hybridEdges(node);
        tree_edge2.istIdentifiable = true
        node.isBadDiamondII = false
    elseif(node.isBadTriangle)
        edge_maj, edge_min, tree_edge2 = hybridEdges(node);
        tree_edge2.istIdentifiable = true
        node.isBadTriangle = false
    elseif(node.isVeryBadTriangle || node.isExtBadTriangle)
        node.isVeryBadTriangle = false
        node.isExtBadTriangle = false
        net.hasVeryBadTriangle = false
    else
        edge_maj, edge_min, tree_edge2 = hybridEdges(node);
        edge_maj.istIdentifiable = isEdgeIdentifiable(edge_maj)
        edge_min.istIdentifiable = isEdgeIdentifiable(edge_min)
    end
end
