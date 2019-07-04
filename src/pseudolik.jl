# functions to calculate the pseudolik function
# originally in functions.jl
# Claudia March 2015


# -------------------- delete Leaf-------------------------------------


# aux function to make a hybrid node a tree node
# used in deleteLeaf
# input: hybrid node
function makeNodeTree!(net::Network, hybrid::Node)
    hybrid.hybrid || error("cannot make node $(hybrid.number) tree node because it already is")
    @debug "we make node $(hybrid.number) a tree node, but it can still have hybrid edges pointing at it, hasHybEdge set to false"
    hybrid.gammaz = -1
    hybrid.isBadDiamondI = false
    hybrid.isBadDiamondII = false
    hybrid.isBadTriangle = false
    hybrid.k = -1
    removeHybrid!(net,hybrid)
    hybrid.hybrid = false
    hybrid.hasHybEdge = false
    hybrid.name = ""
end

# aux function to make a hybrid edge tree edge
# used in deleteLeaf
# input: edge and hybrid node it pointed to
function makeEdgeTree!(edge::Edge, node::Node)
    edge.hybrid || error("cannot make edge $(edge.number) tree because it is tree already")
    node.hybrid || error("need the hybrid node at which edge $(edge.number) is pointing to, node $(node.number) is tree node")
    @debug "we make edge $(edge.number) a tree edge, but it will still point to hybrid node $(node.number)"
    edge.hybrid = false
    edge.isMajor = true
    edge.gamma = 1.0
    edge.inCycle = -1 #warn:changed recently because I believe it does not affect the parameters function for bad diamond I
    other = getOtherNode(edge,node)
    hyb = false
    for e in other.edge #is it attached to other hybrid?
        if(e.hybrid)
            hyb = true
        end
    end
    other.hasHybEdge = hyb
    edge.istIdentifiable = isEdgeIdentifiable(edge)
end

# function to delete internal nodes until there is no
# more internal node with only two edges
# calls deleteIntLeaf! inside a while
# returns middle node (see ipad drawing)
# negative=true allows negative branch lengths
function deleteIntLeafWhile!(net::Network, middle::Node, leaf::Node, negative::Bool)
    while(size(middle.edge,1) == 2)
        middle = deleteIntLeaf!(net,middle,leaf, negative)
    end
    return middle
end

deleteIntLeafWhile!(net::Network, middle::Node, leaf::Node) = deleteIntLeafWhile!(net, middle, leaf, false)

# function to delete internal nodes until there is no
# more internal node with only two edges
# calls deleteIntLeaf! inside a while
# negative=true allows negative branch lengths
function deleteIntLeafWhile!(net::Network, leafedge::Edge, leaf::Node, negative::Bool)
    middle = getOtherNode(leafedge, leaf)
    while(size(middle.edge,1) == 2)
        middle = deleteIntLeaf!(net,middle,leaf, negative)
    end
end

deleteIntLeafWhile!(net::Network, leafedge::Edge, leaf::Node) = deleteIntLeafWhile!(net, leafedge, leaf, false)

# function to delete an internal node in an external edge
# input: network, internal node and leaf
# returns the new middle
# note: leafedge is the edge that survives
# negative=true allows negative branch lengths
function deleteIntLeaf!(net::Network, middle::Node, leaf::Node, negative::Bool)
    #println("calling deleteIntLeaf for middle $(middle.number) and leaf $(leaf.number)")
    if(size(middle.edge,1) == 2)
        if(isEqual(getOtherNode(middle.edge[1],middle),leaf))
            leafedge = middle.edge[1]
            otheredge = middle.edge[2]
        elseif(isEqual(getOtherNode(middle.edge[2],middle),leaf))
            leafedge = middle.edge[2]
            otheredge = middle.edge[1]
        else
            error("leaf $(leaf.number) is not attached to internal node $(middle.number) by any edge")
        end
        othernode = getOtherNode(otheredge,middle)
        setLength!(leafedge,leafedge.length + otheredge.length, negative)
        removeNode!(middle,leafedge)
        removeEdge!(othernode,otheredge)
        setNode!(leafedge,othernode)
        setEdge!(othernode,leafedge)
        deleteNode!(net,middle)
        deleteEdge!(net,otheredge)
        return othernode
    else
        @error "internal node $(middle.number) does not have two edges only, it has $(size(middle.edge,1))"
    end
end

deleteIntLeaf!(net::Network, middle::Node, leaf::Node) = deleteIntLeaf!(net, middle, leaf, false)

# function to delete an internal node in an external edge
# input: network, edge (from leaf.edge) to move in that direction and leaf
# returns the new middle
# note: leafedge is the edge that survives
# negative=true allows negative branch lengths
function deleteIntLeaf!(net::Network, leafedge::Edge, leaf::Node, negative::Bool)
    middle = getOtherNode(leafedge,leaf)
    if(size(middle.edge,1) == 2)
        if(isEqual(getOtherNode(middle.edge[1],middle),leaf))
            otheredge = middle.edge[2]
        elseif(isEqual(getOtherNode(middle.edge[2],middle),leaf))
            otheredge = middle.edge[1]
        else
            error("leaf $(leaf.number) is not attached to internal node $(middle.number) by any edge")
        end
        othernode = getOtherNode(otheredge,middle)
        setLength!(leafedge,leafedge.length + otheredge.length, negative)
        removeNode!(middle,leafedge)
        removeEdge!(othernode,otheredge)
        setNode!(leafedge,othernode)
        setEdge!(othernode,leafedge)
        deleteNode!(net,middle)
        deleteEdge!(net,otheredge)
        return othernode
    else
        @error "internal node $(middle.number) does not have two edges only, it has $(size(middle.edge,1))"
    end
end

deleteIntLeaf!(net::Network, leafedge::Edge, leaf::Node) = deleteIntLeaf!(net, leafedge, leaf, false)

# fixit: still missing to test hybrid-> bad triangle II and
#        normal case (not bad triangle/diamond) of hasHybEdge
#        and normal case, delete not hybrid leaf bad triangle II
function deleteLeaf!(net::Network, leaf::Node)
    leaf.leaf || error("node $(leaf.number) is not a leaf, cannot delete it")
    isNodeNumIn(leaf,net.leaf) || error("node $(leaf.number) is not in net.leaf, cannot delete it")
    size(leaf.edge,1) == 1 || error("strange leaf $(leaf.number) with $(size(leaf.edge,1)) edges instead of 1")
    other = getOtherNode(leaf.edge[1],leaf);
    DEBUGC && @debug "leaf is $(leaf.number) and other is $(other.number)"
    if(other.hybrid)
        DEBUGC && @debug "entra al caso other is hybrid node, other is $(other.number)"
        edge1,edge2 = hybridEdges(other,leaf.edge[1]);
        @debug "edge1 $(edge1.number), edge2 $(edge2.number)"
        (edge1.hybrid && edge2.hybrid) || error("hybrid node $(other.node) does not have two hybrid edges, they are tree edges: $(edge1.number), $(edge2.number)")
        other1 = getOtherNode(edge1,other);
        other2 = getOtherNode(edge2,other);
        @debug "in deleteLeaf, other hybrid node, other1 $(other1.number), other2 $(other2.number)"
        removeEdge!(other1,edge1)
        removeEdge!(other2,edge2)
        deleteEdge!(net,edge1)
        deleteEdge!(net,edge2)
        deleteEdge!(net,leaf.edge[1])
        deleteNode!(net,other)
        deleteNode!(net,leaf)
        if(size(other1.edge,1) == 2  && isNodeNumIn(other1,net.node)) # need to delete internal nodes with only 2 edges (one of them external edge)
            node = getOtherNode(other1.edge[1],other1)
            leaf1 = node.leaf ? node : getOtherNode(other1.edge[2],other1)
            DEBUGC && @debug "other1 si tiene solo dos edges, y leaf1 es $(leaf1.number)"
            if(leaf1.leaf)
                deleteIntLeafWhile!(net,other1,leaf1);
            end
        end
        if(size(other2.edge,1) == 2 && isNodeNumIn(other2,net.node))
            node = getOtherNode(other2.edge[1],other2)
            leaf1 = node.leaf ? node : getOtherNode(other2.edge[2],other2)
            DEBUGC && @debug "other2 si tiene solo dos edges, y leaf1 es $(leaf1.number)"
            if(leaf1.leaf)
                deleteIntLeafWhile!(net,other2,leaf1);
            end
        end
        if(length(other1.edge) == 1 && isNodeNumIn(other1,net.node)) #internal node with only one edge
            !other1.leaf || error("node $(other1.number) was attached no hybrid edge, cannot be a leaf")
            removeWeirdNodes!(net,other1)
        end
        if(length(other2.edge) == 1 && isNodeNumIn(other2,net.node))
            !other2.leaf || error("node $(other2.number) was attached no hybrid edge, cannot be a leaf")
            removeWeirdNodes!(net,other2)
        end
    else
        if(other.hasHybEdge)
            @debug "entro a caso tree node has hyb edge"
            edge1,edge2 = hybridEdges(other,leaf.edge[1]);
            (edge1.hybrid || edge2.hybrid) || error("node $(other.number) has hybrid edge attribute true, but the edges $(edge1.number), $(edge2.number) are not hybrid (and the third edge has a leaf $(leaf.number)")
            other1 = getOtherNode(edge1,other);
            other2 = getOtherNode(edge2,other);
            # warning: if something changed for other1, also need to do for other2
            if(other1.hybrid)
                if(other1.isBadDiamondI)
                    edgemaj,edgemin,edgebla = hybridEdges(other1)
                    edge4 = isEqual(edge1,edgemaj) ? edgemin : edgemaj
                    other3 = getOtherNode(edge4,other1)
                    ind = isEqual(getOtherNode(edgemaj,other1),other3) ? 1 : 2
                    edgebla,edge3,edge5 = hybridEdges(other3)
                    leaf5 = getOtherNode(edge5,other3)
                    DEBUGC && @debug "edge2 is $(edge2.number) and is identifiable $(edge2.istIdentifiable)"
                    edge2.fromBadDiamondI = true # to keep track of edges from bad diamondI
                    removeNode!(other,edge2)
                    removeEdge!(other1,edge1)
                    setNode!(edge2,other1)
                    setEdge!(other1,edge2)
                    other3.gammaz != -1 || error("hybrid node $(other1.number) is bad diamond, but for node $(other3.number), gammaz is not well updated, it is $(other3.gammaz)")
                    DEBUGC && @debug "entro a cambiar length en edge $(edge2.number) con gammaz $(other3.gammaz)"
                    setLength!(edge2,-log(1-other3.gammaz))
                    edge2.number = parse(Int,string(string(other1.number),string(ind)))
                    DEBUGC && @debug "edge2 is $(edge2.number) should be not identifiable now $(edge2.istIdentifiable)"
                    makeEdgeTree!(edge4,other1)
                    makeNodeTree!(net,other1)
                    removeEdge!(other2,edge3)
                    removeEdge!(other3,edge3)
                    deleteNode!(net,other)
                    deleteEdge!(net,edge1)
                    deleteEdge!(net,edge3)
                    removeNode!(other3,edge4)
                    removeEdge!(leaf5,edge5)
                    setNode!(edge4,leaf5)
                    setEdge!(leaf5,edge4)
                    setLength!(edge4,edge5.length)
                    deleteNode!(net,other3)
                    deleteEdge!(net,edge5)
                    other1.hasHybEdge = false
                    @debug begin printEdges(net); "printed edged" end
                else
                    removeEdge!(other,leaf.edge[1])
                    edgebla,edgebla,treeedge = hybridEdges(other1)
                    if(getOtherNode(treeedge,other1).leaf)
                        setLength!(edge1,edge2.length+edge1.length)
                        removeEdge!(other2,edge2)
                        removeNode!(other,edge1)
                        setNode!(edge1,other2)
                        setEdge!(other2,edge1)
                        deleteEdge!(net,edge2)
                        deleteNode!(net,other)
                    end
                end
                deleteNode!(net,leaf)
                deleteEdge!(net,leaf.edge[1])
            elseif(other2.hybrid)
                if(other2.isBadDiamondI)
                    edgemaj,edgemin,edgebla = hybridEdges(other2)
                    edge4 = isEqual(edge2,edgemaj) ? edgemin : edgemaj
                    other3 = getOtherNode(edge4,other2)
                    ind = isEqual(getOtherNode(edgemaj,other2),other3) ? 1 : 2
                    edgebla,edge3,edge5 = hybridEdges(other3)
                    leaf5 = getOtherNode(edge5,other3)
                    DEBUGC && @debug "edge1 is $(edge1.number) and is identifiable $(edge1.istIdentifiable)"
                    edge1.fromBadDiamondI = true # to keep track of edges from bad diamondI
                    removeNode!(other,edge1)
                    removeEdge!(other2,edge2)
                    setNode!(edge1,other2)
                    setEdge!(other2,edge1)
                    other3.gammaz != -1 || error("hybrid node $(other2.number) is bad diamond, but for node $(other3.number), gammaz is not well updated, it is $(other3.gammaz)")
                    setLength!(edge1,-log(1-other3.gammaz))
                    edge1.number = parse(Int,string(string(other2.number),string(ind)))
                    DEBUGC && @debug "edge1 is $(edge1.number) should be not identifiable now $(edge1.istIdentifiable)"
                    makeEdgeTree!(edge4,other2)
                    makeNodeTree!(net,other2)
                    removeEdge!(other1,edge3)
                    removeEdge!(other3,edge3)
                    deleteNode!(net,other)
                    deleteEdge!(net,edge2)
                    deleteEdge!(net,edge3)
                    removeNode!(other3,edge4)
                    removeEdge!(leaf5,edge5)
                    setNode!(edge4,leaf5)
                    setEdge!(leaf5,edge4)
                    setLength!(edge4,edge5.length)
                    deleteNode!(net,other3)
                    deleteEdge!(net,edge5)
                else
                    removeEdge!(other,leaf.edge[1])
                    edgebla,edgebla,treeedge = hybridEdges(other2)
                    if(getOtherNode(treeedge,other2).leaf)
                        setLength!(edge2,edge2.length+edge1.length)
                        removeEdge!(other1,edge1)
                        removeNode!(other,edge2)
                        setNode!(edge2,other1)
                        setEdge!(other1,edge2)
                        deleteEdge!(net,edge1)
                        deleteNode!(net,other)
                    end
                end
                deleteNode!(net,leaf)
                deleteEdge!(net,leaf.edge[1])
            else
                error("node $(other.number) has hybrid edge, but neither of the other nodes $(other1.number), $(other2.number) are hybrid")
            end
        else # other is tree node without hybrid edges
            DEBUGC && @debug "entra al caso other is tree node no hyb edges"
            if(length(other.edge) == 2)
                edge1 = isEqual(other.edge[1],leaf.edge[1]) ? other.edge[2] : other.edge[1]
                other1 = getOtherNode(edge1,other)
                removeEdge!(other1,edge1)
                removeNode!(other1,edge1)
                deleteNode!(net,other)
                deleteNode!(net,leaf)
                deleteEdge!(net,leaf.edge[1])
                deleteEdge!(net,edge1)
                DEBUGC && @debug "other1 is $(other1.number)"
                if(size(other1.edge,1) == 2  && isNodeNumIn(other1,net.node)) # need to delete internal nodes with only 2 edges (one of them external edge)
                    node = getOtherNode(other1.edge[1],other1)
                    leaf1 = node.leaf ? node : getOtherNode(other1.edge[2],other1)
                    DEBUGC && @debug "other1 si tiene solo dos edges, y leaf1 es $(leaf1.number)"
                    if(leaf1.leaf)
                        other1 = deleteIntLeafWhile!(net,other1,leaf1);
                    end
                end
                if(other1.hybrid && length(other1.edge) == 2)
                    DEBUGC && @debug "entra a tratar de arreglar a other1 $(other1.number)"
                    removeWeirdNodes!(net,other1)
                end
                if(length(other1.edge) == 1 && isNodeNumIn(other1,net.node)) #internal node with only one edge
                    DEBUGC && @debug "entra a tratar de arreglar a other1 $(other1.number)"
                    !other1.leaf || error("node $(other1.number) was attached no leaf edge, cannot be a leaf")
                    removeWeirdNodes!(net,other1)
                end
            else
                edge1,edge2 = hybridEdges(other,leaf.edge[1]);
                other1 = getOtherNode(edge1,other);
                other2 = getOtherNode(edge2,other);
                removeEdge!(other,leaf.edge[1])
                deleteNode!(net,leaf)
                deleteEdge!(net,leaf.edge[1])
                if(other1.leaf || other2.leaf)
                    (!other1.leaf || !other2.leaf) || error("just deleted a leaf $(leaf.number) and its two attached nodes are leaves also $(other1.number), $(other2.number)")
                    newleaf = other1.leaf ? other1 : other2
                    middle = other
                    DEBUGC && @debug "middle is $(middle.number), middle.hybrid $(middle.hybrid), middle.hasHybEdge $(middle.hasHybEdge)"
                    middle = deleteIntLeafWhile!(net,middle,newleaf)
                    DEBUGC && @debug "middle is $(middle.number), middle.hybrid $(middle.hybrid), middle.hasHybEdge $(middle.hasHybEdge)"
                    if(middle.hybrid)
                        edges = hybridEdges(middle)
                        edges[1].istIdentifiable = false
                        edges[2].istIdentifiable = false
                    end
                end
            end
        end
    end
end

"""
`deleteLeaf!(net::HybridNetwork, leaf::AbstractString)`
`deleteLeaf!(net::Network, leaf::Node)`

Deletes the leaf taxon from the network. The leaf argument is the name of the taxon to delete.

Warnings:

- requires a level-1 network with up-to-date attributes for snaq! (e.g. non-missing branch lengths, gammaz, etc.)
- does not care where the root is and does not update it to a sensible location if the root
  is affected by the leaf removal.
- does not merge edges, i.e. does not remove all nodes of degree 2. Within snaq!, this
  is used to extract quartets and to keep track of which
  edge lengths in the original network map to the quartet network.
"""
function deleteLeaf!(net::Network, leaf::AbstractString)
    found = false
    for l in net.leaf
        if(l.name == leaf)
            found = true
            deleteLeaf!(net,l)
            break
        end
    end
    found || error("cannot delete leaf $(leaf) because it was not found in the network")
end


# -------------------- extract quartet ---------------------------------------

# function to update hasEdge attribute in a
# QuartetNetwork after leaves deleted
# with deleteLeaf!
function updateHasEdge!(qnet::QuartetNetwork, net::HybridNetwork)
    #@warn "function to compare edges depends on edges number being unique"
    edges = Bool[]
    h = Bool[]
    hz = Bool[]
    for e in net.edge
        if(e.istIdentifiable)
            if(isEdgeNumIn(e,qnet.edge) && qnet.edge[getIndexEdge(e.number,qnet)].istIdentifiable)
                #println("found identifiable edge $(e.number) in net and qnet")
                push!(edges,true)
            else
                #println("not found identifiable edge $(e.number) in qnet but identifiable in net")
                push!(edges,false)
            end
        end
        if e.hybrid && !e.isMajor
            node = e.node[e.isChild1 ? 1 : 2]
            node.hybrid || error("strange thing, hybrid edge $(e.number) pointing at tree node $(node.number)")
            if(!node.isBadDiamondI)
                #println("found hybrid edge $(e.number) in net and qnet")
                push!(h,isNodeNumIn(node,qnet.hybrid))
            else
                if(isNodeNumIn(node,qnet.hybrid))
                    push!(hz,true)
                    push!(hz,true)
                else
                    ind1 = parse(Int,string(string(node.number),"1"))
                    ind2 = parse(Int,string(string(node.number),"2"))
                    found1 = true
                    found2 = true
                    try
                        getIndexEdge(ind1,qnet)
                    catch
                        found1 = false
                    end
                    try
                        getIndexEdge(ind2,qnet)
                    catch
                        found2 = false
                    end
                    if(!found1 && !found2)
                        push!(hz,false)
                        push!(hz,false)
                    elseif(found1 && !found2)
                        index = getIndexEdge(ind1,qnet)
                        if(!qnet.edge[index].istIdentifiable && all((n->!n.leaf), qnet.edge[index].node))
                            push!(hz,true)
                        else
                            push!(hz,false)
                        end
                        push!(hz,false)
                    elseif(found2 && !found1)
                        push!(hz,false)
                        index = getIndexEdge(ind2,qnet)
                        if(!qnet.edge[index].istIdentifiable && all((n->!n.leaf), qnet.edge[index].node))
                            push!(hz,true)
                        else
                            push!(hz,false)
                        end
                    else
                        index1 = getIndexEdge(ind1,qnet)
                        index2 = getIndexEdge(ind2,qnet)
                        if(!qnet.edge[index1].istIdentifiable && all((n->!n.leaf), qnet.edge[index1].node) && !qnet.edge[index2].istIdentifiable && all((n->!n.leaf), qnet.edge[index2].node))
                            error("strange qnet when net has node $(node.number) Bad Diamond I: qnet should have only one of the gammaz if it does not have node, but it has two")
                        end
                    end
                end
            end
        end
    end # for in net.edge
    #println("edges $(edges), h $(h), hz $(hz)")
    qnet.hasEdge = vcat(h,edges,hz)
end


# function to extract a quartet from a network
# input: QuartetNetwork (already created from HybridNetwork)
#        quartet: array with the 4 leaf nodes to keep
# return: QuartetNetwork with only 4 tips
# it updates qnet.hasEdge and qnet.indexht
function extractQuartet(net::HybridNetwork,quartet::Array{Node,1})
    size(quartet,1) == 4 || error("quartet array should have 4 nodes, it has $(size(quartet,1))")
    (quartet[1].leaf && quartet[2].leaf && quartet[3].leaf && quartet[4].leaf) || error("all four nodes to keep when extracting the quartet should be leaves: $([q.number for q in quartet])")
    qnet = QuartetNetwork(net) # fixit: try to re-use memory? quartetTaxon has not changed for instance.
    leaves = copy(qnet.leaf)
    for n in leaves
        if(!isNodeNumIn(n,quartet))
            DEBUGC && @debug "delete leaf $(n.number)"
            deleteLeaf!(qnet,n)
            DEBUGC && printEdges(qnet)
        end
    end
    DEBUGC && @debug "deletion of leaves successful"
    return qnet
end

# function to extract a quartet from a Quartet object
# it calls the previous extractQuartet
# returns qnet (check: maybe not needed later) and assigns
# quartet.qnet = qnet
function extractQuartet!(net::HybridNetwork, quartet::Quartet)
    list = Node[]
    for q in quartet.taxon
        try
            getIndexNode(getIndex(q,net.names),net)
        catch
            error("taxon $(q) not in network")
        end
        push!(list, net.node[getIndexNode(getIndex(q,net.names),net)])
    end
    qnet = extractQuartet(net,list)
    @debug "EXTRACT: extracted quartet $(quartet.taxon)"
    redundantCycle!(qnet) #removes no leaves, cleans external edges
    updateHasEdge!(qnet,net)
    parameters!(qnet,net)
    qnet.quartetTaxon = quartet.taxon
    quartet.qnet = qnet
    #return qnet
end


# function to extract all quartets from net according
# to the array of quartets of a Data object
# it updates expCF, hasEdge, indexht
function extractQuartet!(net::HybridNetwork, quartet::Vector{Quartet})
    @debug "EXTRACT: begins extract quartets for network"
    for q in quartet
        extractQuartet!(net,q)
        qnet = deepcopy(q.qnet); #there is a reason not to mess up with the original q.qnet, i believe to keep ht consistent
        calculateExpCFAll!(qnet);
        q.qnet.expCF = qnet.expCF
    end
end

extractQuartet!(net::HybridNetwork, d::DataCF) = extractQuartet!(net, d.quartet)

# function to check if there are potential redundant cycles in net
# return the flag (true=redundant cycle found) and the hybrid node for the redundant cycle
function hasRedundantCycle(net::Network)
    length(net.hybrid) == net.numHybrids || error("found net with length net.hybrid of $(length(net.hybrid)) and net.numHybrids of $(net.numHybrids)")
    if net.numHybrids > 0
        for h in net.hybrid
            k = sum([(n.inCycle == h.number && length(n.edge) == 3) ? 1 : 0 for n in net.node])
            if k < 2
                return true,h
            end
        end
    end
    return false,nothing
end


# function to delete redundante cycles on all hybrid nodes in net
function redundantCycle!(net::Network)
    length(net.hybrid) == net.numHybrids || error("found net with length net.hybrid of $(length(net.hybrid)) and net.numHybrids of $(net.numHybrids)")
    if(length(net.hybrid) > 0)
        redCycle, node = hasRedundantCycle(net)
        while(redCycle)
            !isa(node,Nothing) || error("redundant cycle found, but the hybrid node is set to nothing")
            redundantCycle!(net,node)
            DEBUGC && @debug "after redundante cycle for hybrid node $(n.number)"
            DEBUGC && printEdges(net)
            DEBUGC && printNodes(net)
            redCycle, node = hasRedundantCycle(net)
        end
    end
    cleanExtEdges!(net)
end


# function to delete redundant cycles (k=1,k=2), see ipad notes
# warning: should not modify the quartet too much because we need to extract ht afterwards
function redundantCycle!(net::Network,n::Node)
    n.hybrid || error("cannot clean a cycle on a tree node $(n.number)")
    edges = hybridEdges(n)
    edges[1].hybrid && edges[2].hybrid || error("hybrid node $(n.number) does not have two hybrid edges $(edges[1].number), $(edges[2].number)")
    @debug "search for redundantCycle for node $(n.number) with edges are $([e.number for e in edges])"
    other = getOtherNode(edges[1],n)
    if(length(other.edge) == 2)
        e = edges[1]
        while(length(other.edge) == 2)
            ind = isEqual(e,other.edge[1]) ? 2 : 1
            e = other.edge[ind]
            other = getOtherNode(e,other)
        end
        if(isEqual(n,other))
            n1 = getOtherNode(edges[1],n)
            n2 = getOtherNode(edges[2],n)
            @debug "redundant cycle found! n1 $(n1.number), n2 $(n2.number)"
            deleteIntLeafWhile!(net,n1,n)
            edge = n.edge[1].hybrid ? n.edge[1] : n.edge[2]
            @debug "edge is $(edge.number), should be the first (or only) edge in hybrid node $(n.number)"
            if(isEqual(edge.node[1],edge.node[2]))
                @debug "entra a q son iguales los nodes de edge"
                n3 = getOtherNode(edges[3],n)
                @debug "edges[3] is $(edges[3].number), n3 is $(n3.number)"
                removeEdge!(n3,edges[3])
                deleteNode!(net,n)
                deleteEdge!(net,edge)
                deleteEdge!(net,edges[3])
                if(!n3.leaf && length(n3.edge)==1)
                    removeNoLeafWhile!(net,n3);
                end
            end
        end
    end
end

# function to delete an internal node with only one edge
function removeNoLeaf!(net::Network,n::Node)
    !n.leaf || error("node $(n.number) is a leaf, so we cannot remove it")
    length(n.edge) == 1 || error("node $(n.number) has $(length(n.edge)) edges (not 1), so we do not have to remove it")
    node = getOtherNode(n.edge[1],n)
    removeEdge!(node,n.edge[1])
    deleteNode!(net,n)
    deleteEdge!(net,n.edge[1])
    return node
end

# function to do a while for removeNoLeaf
# warning: we do not want to remove edges because we need to know
# exactly which edges in qnet are affected by changes in net.ht
function removeNoLeafWhile!(net::Network,n::Node)
    while(!n.leaf && length(n.edge)==1)
        n = removeNoLeaf!(net,n)
    end
    return n
end

# Function to check that external edges do not have nodes with only two edges
function cleanExtEdges!(net::Network)
#    if(net.numHybrids > 0)
        for l in net.leaf
            length(l.edge) == 1 || error("leaf $(l.number) with $(length(l.edge)) edges, not 1")
            deleteIntLeafWhile!(net,getOtherNode(l.edge[1],l),l)
        end
#    end
end

# function to delete a hybrid node that has only two edges: the hybrid ones
# it fixes other1,other2 for internal nodes with only 2 edges (one external),
# returns other1, other2 if they are internal nodes with only one edge
# returns array of nodes (to avoid treating o1,o2 as a tuple later)
function removeLoneHybrid!(net::Network, n::Node)
    n.hybrid || error("cannot remove a lone hybrid if it is not hybrid: $(n.number)")
    length(n.edge) == 2 || error("hybrid node $(n.number) should have only 2 edges to be deleted, it has $(length(n.edge))")
    (n.edge[1].hybrid && n.edge[2].hybrid) || error("both edges have to be hybrid and they are not both: $(n.edge[1].number), $(n.edge[2].number)")
    other1 = getOtherNode(n.edge[1],n)
    other2 = getOtherNode(n.edge[2],n)
    removeEdge!(other1,n.edge[1])
    removeEdge!(other2,n.edge[2])
    deleteEdge!(net,n.edge[1])
    deleteEdge!(net,n.edge[2])
    deleteNode!(net,n)
    if(length(other1.edge) == 2 && isNodeNumIn(other1,net.node))
        leaf1 = getOtherNode(other1.edge[1],other1)
        leaf = leaf1.leaf ? leaf1 : getOtherNode(other1.edge[2],other1)
        if(leaf.leaf)
            deleteIntLeafWhile!(net,other1,leaf)
        end
    end
    if(length(other2.edge) == 2 && isNodeNumIn(other2,net.node))
        leaf1 = getOtherNode(other2.edge[1],other2)
        leaf = leaf1.leaf ? leaf1 : getOtherNode(other2.edge[2],other2)
        if(leaf.leaf)
            deleteIntLeafWhile!(net,other2,leaf)
        end
    end
    o1 = false
    o2 = false
    if(length(other1.edge) == 1)
        !other1.leaf || error("node $(other1.number) should not be a leaf because it was attached to hybrid edge")
        o1 = true
    end
    if(length(other2.edge) == 1)
        !other2.leaf || error("node $(other2.number) should not be a leaf because it was attached to hybrid edge")
        o2 = true
    end
    if(o1 && !o2)
        return [other1]
    elseif(!o1 && o2)
        return [other2]
    elseif(o1 && o2)
        return [other1,other2]
    else
        return nothing
    end
end

function removeWeirdNode!(net::Network,n::Node)
    other = nothing
    @debug "calling removeWeirdNode in node $(n.number)"
    if(length(n.edge) == 1 && !n.leaf)
        @debug "node $(n.number) is not a leaf and has one edge"
        n = removeNoLeafWhile!(net,n)
    end
    if(n.hybrid && length(n.edge) == 2)
        @debug "node $(n.number) is hybrid with only two edges"
        (n.edge[1].hybrid && n.edge[2].hybrid) || error("two edges for this hybrid node $(n.number) must be hybrid, and they are not")
        other = removeLoneHybrid!(net,n) #returns array of nodes (length 1 or 2) or nothing
    end
    return other
end

# function to remove internal nodes with only one edge
# keep deleting those nodes, and hybrid nodes without descendants
function removeWeirdNodes!(net::Network, n::Node)
    @debug "calling removeWeirdNodes in node $(n.number)"
    list = Node[]
    push!(list,n)
    while !isempty(list)
        n = pop!(list)
        if length(n.edge) == 1 && !n.leaf
            n = removeWeirdNode!(net,n)
        elseif n.hybrid && length(n.edge) == 2
            (n.edge[1].hybrid && n.edge[2].hybrid) || error("hybrid node $(n.number) only has two edges and they must be hybrids")
            n = removeWeirdNode!(net,n)
        else
            @error "calling removeWeirdNode on normal node $(n.number)"
            n = nothing
        end
        if n !== nothing && !isa(n,Vector{Nothing})
            @debug "typeof n $(typeof(n))"
            for l in n
                @debug "typeof l $(typeof(l))"
                push!(list,l)
            end
        end
    end
end
# ------------------------------- calculate expCF -------------------------------------

# ------- Identify Quartet

# function to identify the type of hybridization of a given hybrid node
# in a quartet network
# sets node.k to the updated count
# sets the node.typeHyb (1,2,3,4,5 see ipad notes) and pushes into qnet.typeHyb
# sets node.prev = other to keep track of the "other" node
# fixit: I think qnet.typeHyb is never used
function identifyQuartet!(qnet::QuartetNetwork, node::Node)
    node.hybrid || error("cannot identify the hybridization around node $(node.number) because it is not hybrid node.")
    k = sum([(n.inCycle == node.number && size(n.edge,1) == 3) ? 1 : 0 for n in qnet.node])
    node.k = k
    if k < 2
        @debug begin
            printEdges(qnet)
            printNodes(qnet)
            "printed edges and nodes"
        end
        error("strange quartet network with a hybrid node $(node.number) but no cycle")
    elseif(k == 2)
        other = qnet.node[getIndex(true, [(n.inCycle == node.number && size(n.edge,1) == 3 && !isEqual(n,node)) for n in qnet.node])]
        edgebla,edgebla,edge1 = hybridEdges(node)
        edgebla,edgebla,edge2 = hybridEdges(other)
        if(getOtherNode(edge1,node).leaf || getOtherNode(edge2,other).leaf)
            node.typeHyb = 1
            node.prev = other
            push!(qnet.typeHyb,1)
        else
            node.typeHyb = 3
            node.prev = other
            push!(qnet.typeHyb,3)
        end
    elseif(k == 3)
        #println("identifyQuartet: entra a k=3")
        edge1,edge2,edge3 = hybridEdges(node)
        #println("node is $(node.number), edge3 is $(edge3.number)")
        if(getOtherNode(edge3,node).leaf)
            node.typeHyb = 2
            push!(qnet.typeHyb,2)
            for n in qnet.node
                if(n.inCycle == node.number && size(n.edge,1) == 3 && !isEqual(n,node))
                    #println("found a node $(n.number)")
                    edge1,edge2,edge3 = hybridEdges(n)
                    #println("edge1,egde2,edge3 are $(edge1.number), $(edge2.number), $(edge3.number)")
                    if(getOtherNode(edge3,n).leaf)
                        node.prev = n
                        break
                    end
                end
            end
        else
            node.typeHyb = 4
            push!(qnet.typeHyb,4)
            for n in qnet.node
                if(n.inCycle == node.number && size(n.edge,1) == 3 && !isEqual(n,node))
                    node.prev = n
                    break
                end
            end
        end
    elseif(k == 4)
        node.typeHyb = 5
        push!(qnet.typeHyb,5)
        node.prev = nothing
    else
        @debug begin
            printEdges(qnet)
            printNodes(qnet)
            "printed edges and nodes"
        end
        error("strange quartet network with $(k) nodes in cycle, maximum should be 4")
    end
    DEBUGC && @debug "qnet identified as type $(node.typeHyb)"
end

# function to identify the Quartet network as
# 1 (equivalent to tree), 2 (minor CF different)
function identifyQuartet!(qnet::QuartetNetwork)
    #if(qnet.which == -1)
        if(qnet.numHybrids == 0)
            qnet.which = 1
        elseif(qnet.numHybrids == 1)
            identifyQuartet!(qnet, qnet.hybrid[1])
            if(qnet.typeHyb[1] == 5)
                qnet.which = 2
            else
                qnet.which = 1
            end
        else
            for n in qnet.hybrid
                cleanUpNode!(qnet,n)
                identifyQuartet!(qnet, n)
            end
            if(all((n->(n.typeHyb != 5)), qnet.hybrid))
                qnet.which = 1
            else
                if(all((n->(n.typeHyb == 5 || n.typeHyb == 1)), qnet.hybrid))
                    #@warn "hybridization of type 5 found in quartet network along with other hybridizations of type 1. there is the possibility of overlapping cycles."
                    qnet.which = 2
                else
                    @debug begin
                        printEdges(qnet)
                        printNodes(qnet)
                        "warning: found in the same quartet, two hybridizations with overlapping cycles: types of hybridizations are $([n.typeHyb for n in qnet.hybrid]), maybe this will cause problems if the hyb do not become all but one type 1"
                    end
                    qnet.which = 2
                end
            end
        end
    #end
end

# function that will get rid of internal nodes with only
# two edges for all three directions of a node
function cleanUpNode!(net::Network,node::Node)
    edge1,edge2,edge3 = hybridEdges(node)
    deleteIntLeafWhile!(net,getOtherNode(edge1,node),node)
    deleteIntLeafWhile!(net,getOtherNode(edge2,node),node)
    deleteIntLeafWhile!(net,getOtherNode(edge3,node),node)
end

# ----------------------- Eliminate Hybridizations

# aux function to eliminate hybridizations
# in quartet
# input: quartet network,
#        node correspond to the hybrid node
#        internal: true if loop is in internal edge
function eliminateLoop!(qnet::QuartetNetwork, node::Node, internal::Bool)
    node.hybrid || error("cannot eliminate loop around node $(node.number) since it is not hybrid")
    edge1,edge2,edge3 = hybridEdges(node)
    deleteIntLeafWhile!(qnet, edge1, node)
    deleteIntLeafWhile!(qnet, edge2, node)
    other = getOtherNode(edge1,node)
    isEqual(getOtherNode(edge2,node),other) || error("node $(node.number) and other $(other.number) are not the two nodes in a cycle with k=2")
    other.number == node.prev.number || error("something strange, other node $(other.number) should be the same as the stored in node.prev $(node.prev.number)")
    removeEdge!(node,edge2)
    removeEdge!(other,edge2)
    deleteEdge!(qnet,edge2)
    if(internal)
        setLength!(edge1, -log(1-edge1.gamma*edge1.gamma*edge1.z-edge2.gamma*edge2.gamma*edge2.z))
    else
        leaf = getOtherNode(edge3,node)
        #if(leaf.leaf)
            deleteIntLeafWhile!(qnet,node,leaf)
        #else
        #    edge1,edge2,edge3 = hybridEdges(other)
        #    deleteIntLeafWhile!(qnet,other,getOtherNode(edge3,other))
        #end
    end
end

# aux function to identify intermediate edge between two nodes
function intermediateEdge(node1::Node,node2::Node)
    edge = nothing
    for e in node1.edge
        if(isEqual(getOtherNode(e,node1),node2))
            edge = e
        end
    end
    if(isa(edge, Nothing))
        error("nodes $(node1.number), $(node2.number) are not connected by an edge")
    end
    return edge
end

# function to eliminate a triangle hybridization
# input: quartet network,
#        node, other nodes in the hybridization
#        case: 1 (global case 2),2 (global case 4), 1 (global case 5)
# warning: special treatment for bad diamond II
function eliminateTriangle!(qnet::QuartetNetwork, node::Node, other::Node, case::Integer)
    #println("start eliminateTriangle----")
    node.hybrid || error("cannot eliminate triangle around node $(node.number) since it is not hybrid")
    #println("hybrid node is $(node.number), with edges $([e.number for e in node.edge]), with gammas $([e.gamma for e in node.edge])")
    edgemaj, edgemin, treeedge = hybridEdges(node)
    isa(edgemaj,Nothing) ? error("edge maj is nothing for node $(node.number), other $(other.number) and taxon $(qnet.quartetTaxon), $(printEdges(qnet))") : nothing
    isa(edgemin,Nothing) ? error("edge min is nothing for node $(node.number), other $(other.number) and taxon $(qnet.quartetTaxon), $(printEdges(qnet))") : nothing
    deleteIntLeafWhile!(qnet, edgemaj, node)
    deleteIntLeafWhile!(qnet, edgemin, node)
    if(isEqual(getOtherNode(edgemaj,node),other))
        hybedge = edgemaj
        otheredge = edgemin
    elseif(isEqual(getOtherNode(edgemin,node),other))
        hybedge = edgemin
        otheredge = edgemaj
    else
        error("node $(node.number) and other node $(other.number) are not connected by an edge")
    end
    #println("hybedge is $(hybedge.number), otheredge is $(otheredge.number)")
    middle = qnet.node[getIndex(true, [(n.inCycle == node.number && size(n.edge,1) == 3 && !isEqual(n,other) && !isEqual(n,node)) for n in qnet.node])]
    #println("middle node is $(middle.number) in eliminateTriangle")
    ind = getIndex(true,[(e.inCycle == node.number && !isEqual(getOtherNode(e,middle),node)) for e in middle.edge])
    edge = middle.edge[ind]
    #println("edge is $(edge.number) with length $(edge.length) in eliminateTriangle, will do deleteIntLeaf from middle through edge")
    deleteIntLeafWhile!(qnet,edge,middle)
    #println("after deleteIntLeaf, edge $(edge.number) has length $(edge.length)")
    isEqual(getOtherNode(edge,middle),other) || error("middle node $(middle.number) and other node $(other.number) are not connected by an edge")
    if(case == 1)
        setLength!(edge,-log(1 - hybedge.gamma*edge.z))
        #println("Case 1: edge $(edge.number) length is $(edge.length) after updating")
        removeEdge!(middle,otheredge)
        removeEdge!(other,hybedge)
        removeEdge!(node,treeedge)
        removeNode!(node,treeedge)
        setEdge!(other, treeedge)
        setNode!(treeedge, other)
        deleteEdge!(qnet,otheredge)
        deleteEdge!(qnet,hybedge)
        deleteNode!(qnet,node)
    elseif(case == 2)
        (hybedge.hybrid && otheredge.hybrid) || error("hybedge $(hybedge.number) and otheredge $(otheredge.number) should by hybrid edges in eliminateTriangle Case 2")
        setLength!(hybedge, -log(otheredge.gamma*otheredge.gamma*otheredge.y + hybedge.gamma*otheredge.gamma*(3-edge.y) + hybedge.gamma*hybedge.gamma*hybedge.y), true)
        #println("Case 2: edge $(edge.number) length is $(edge.length) after updating")
        removeEdge!(middle,otheredge)
        removeEdge!(node,otheredge)
        deleteEdge!(qnet,otheredge)
        deleteIntLeafWhile!(qnet, node, getOtherNode(treeedge,node),true)
        deleteIntLeafWhile!(qnet, edge, other)
    else
        error("unknown case $(case), should be 1 or 2")
    end
    #println("end eliminateTriangle ---")
end

# function to polish quartet with hybridization type 5
# (2 minor CF different) and calculate the expCF
# CF calculated in the order 12|34, 13|24, 14|23 of the qnet.quartet.taxon
function quartetType5!(qnet::QuartetNetwork, node::Node)
    (node.hybrid && node.typeHyb == 5) || error("cannot polish the quartet type 5 hybridization since either the node is not hybrid: $(!node.hybrid) or it has type $(node.typeHyb), different than 5")
    edge1,edge2,edge3 = hybridEdges(node);
    if(!node.isBadDiamondI)
        deleteIntLeafWhile!(qnet,edge1,node)
        deleteIntLeafWhile!(qnet,edge2,node)
    end
    other1 = getOtherNode(edge1,node);
    other2 = getOtherNode(edge2,node);
    edgebla,edge5, edgetree1 = hybridEdges(other1);
    edgebla,edge6, edgetree2 = hybridEdges(other2);
    if(!node.isBadDiamondI)
        deleteIntLeafWhile!(qnet,edge5,other1)
        deleteIntLeafWhile!(qnet,edge6,other2)
    end
    if(node.isBadDiamondI)
        (other1.gammaz != -1 && other2.gammaz != -1) || error("node $(node.number) is bad diamond I but gammaz are -1")
        @debug "it will calculate the expCF in a bad diamond I case with gammaz: $(other1.gammaz) and $(other2.gammaz)"
        cf1 = (1 + 2*other1.gammaz - other2.gammaz)/3
        cf2 = (1 + 2*other2.gammaz - other1.gammaz)/3
        cf3 = (1 - other1.gammaz - other2.gammaz)/3
    else
        cf1 = edge1.gamma*(1-2/3*edge5.y) + edge2.gamma*1/3*edge6.y
        cf2 = edge1.gamma*1/3*edge5.y + edge2.gamma*(1-2/3*edge6.y)
        cf3 = edge1.gamma*1/3*edge5.y + edge2.gamma*1/3*edge6.y
    end
    #println("cf1,cf2,cf3: $(cf1),$(cf2),$(cf3)")
    leaf1 = getOtherNode(edge3,node)
    if(isa(edgetree1,Nothing))
        println("node $(node.number), edge3 $(edge3.number), other1 $(other1.number), leaf1 $(leaf1.number), other2 $(other2.number)")
        println("edge1 $(edge1.number), edge2 $(edge2.number), edge5 $(edge5.number), edge6 $(edge6.number)")
        printEdges(qnet)
        printNodes(qnet)
    end
    leaf2 = getOtherNode(edgetree1,other1)
    leaf3 = getOtherNode(edgetree2, other2)
    leaf4 = qnet.leaf[getIndex(true,[(!isEqual(n,leaf1) && !isEqual(n,leaf2) && !isEqual(n,leaf3)) for n in qnet.leaf])]
    #println("leaf1 is $(leaf1.number)")
    #println("leaf2 is $(leaf2.number)")
    #println("leaf3 is $(leaf3.number)")
    #println("leaf4 is $(leaf4.number)")
    tx = whichLeaves(qnet,qnet.quartetTaxon[1],qnet.quartetTaxon[2], leaf1,leaf2,leaf3,leaf4)
    #println("tx is $(tx)")
    if(tx == (1,2) || tx == (2,1) || tx == (3,4) || tx == (4,3))
        qnet.expCF[1] = cf1
        tx = whichLeaves(qnet,qnet.quartetTaxon[1],qnet.quartetTaxon[3], leaf1,leaf2,leaf3,leaf4)
        if(tx == (1,3) || tx == (3,1) || tx == (2,4) || tx == (4,2))
            qnet.expCF[2] = cf2
            qnet.expCF[3] = cf3
        elseif(tx == (1,4) || tx == (4,1) || tx == (3,2) || tx == (2,3))
            qnet.expCF[2] = cf3
            qnet.expCF[3] = cf2
        else
            error("strange quartet network, could not find which leaves correspond to taxon1, taxon3")
        end
    elseif(tx == (1,3) || tx == (3,1) || tx == (2,4) || tx == (4,2))
        qnet.expCF[1] = cf2
        tx = whichLeaves(qnet,qnet.quartetTaxon[1],qnet.quartetTaxon[3], leaf1,leaf2,leaf3,leaf4)
        if(tx == (1,2) || tx == (2,1) || tx == (3,4) || tx == (4,3))
            qnet.expCF[2] = cf1
            qnet.expCF[3] = cf3
        elseif(tx == (1,4) || tx == (4,1) || tx == (3,2) || tx == (2,3))
            qnet.expCF[2] = cf3
            qnet.expCF[3] = cf1
        else
            error("strange quartet network, could not find which leaves correspond to taxon1, taxon3")
        end
    elseif(tx == (1,4) || tx == (4,1) || tx == (3,2) || tx == (2,3))
        qnet.expCF[1] = cf3
        tx = whichLeaves(qnet,qnet.quartetTaxon[1],qnet.quartetTaxon[3], leaf1,leaf2,leaf3,leaf4)
        if(tx == (1,3) || tx == (3,1) || tx == (2,4) || tx == (4,2))
            qnet.expCF[2] = cf2
            qnet.expCF[3] = cf1
        elseif(tx == (1,2) || tx == (2,1) || tx == (3,4) || tx == (4,3))
            qnet.expCF[2] = cf1
            qnet.expCF[3] = cf2
        else
            error("strange quartet network, could not find which leaves correspond to taxon1, taxon3")
        end
    else
        error("strange quartet network, could not find which leaves correspond to taxon1, taxon2")
    end
    if(!approxEq(sum(qnet.expCF),1.))
        error("strange quartet network with hybridization in node $(node.number) of type 5: expCF do not add up to 1")
    end
end

# function to eliminate a hybridization around a given
# hybrid node
function eliminateHybridization!(qnet::QuartetNetwork, node::Node)
    node.hybrid || error("cannot eliminate hybridization around node $(node.number) since it is not hybrid node")
    if(node.typeHyb == 1)
        eliminateLoop!(qnet,node,false)
    elseif(node.typeHyb == 3)
        eliminateLoop!(qnet,node,true)
    elseif(node.typeHyb == 4)
        #println("node is $(node.number), other node is $(node.prev.number)")
        eliminateTriangle!(qnet,node,node.prev,2)
    elseif(node.typeHyb == 2)
        #println("node is $(node.number), other node is $(node.prev.number)")
        eliminateTriangle!(qnet,node,node.prev,1)
    elseif(node.typeHyb != 5)
        error("node type of hybridization should be 1,2,3,4 or 5, but for node $(node.number), it is $(node.typeHyb)")
    end
end

# aux function to eliminate all internal nodes with only
# two edges in a quartet network with qnet.which=1
# eliminate internal nodes in every direction
function internalLength!(qnet::QuartetNetwork)
    if(qnet.which == 1)
        try
            getIndex(true,[size(n.edge,1) == 3 for n in qnet.node])
        catch
            printEdges(qnet)
            printNodes(qnet)
            error("not found internal node in qnet with 3 edges")
        end
        node = qnet.node[getIndex(true,[size(n.edge,1) == 3 for n in qnet.node])]
        try
            getIndex(true,[size(n.edge,1) == 3 && !isEqual(n,node) for n in qnet.node])
        catch
            println("first node found with 3 edges $(node.number)")
            printEdges(qnet)
            printNodes(qnet)
            error("not found another internal node in qnet with 3 edges")
        end
        node2 = qnet.node[getIndex(true,[size(n.edge,1) == 3 && !isEqual(n,node) for n in qnet.node])]
        for e in node.edge
            deleteIntLeafWhile!(qnet,e,node,true)
        end
        for e in node2.edge
            deleteIntLeafWhile!(qnet,e,node2,true)
        end
        edge = nothing
        for e in node.edge
            if(!getOtherNode(e,node).leaf)
                edge = e
            end
        end
        !isa(edge,Nothing) || error("cannot find internal edge attached to node $(node.number) in qnet")
        isEqual(getOtherNode(edge,node),node2) || error("strange internal edge $(edge.number) found in qnet, should have as nodes $(node.number), $(node2.number)")
        qnet.t1 = edge.length
    end
end

# function to eliminate hybridizations in a quartet network
# first step to later calculate expCF
# input: quartet network
# fixit: need to add a loop, eliminate type1 until there are none, and then identify again
function eliminateHybridization!(qnet::QuartetNetwork)
    qnet.which != -1 || error("qnet which has to be updated by now to 1 or 2, and it is $(qnet.which)")
    if(qnet.numHybrids == 1)
        eliminateHybridization!(qnet,qnet.hybrid[1])
    elseif(qnet.numHybrids > 1)
        #eliminate in order: first type1 only
        DEBUGC && @debug "starting eliminateHyb for more than one hybrid with types $([n.typeHyb for n in qnet.hybrid])"
        while(qnet.numHybrids > 0 && any([n.typeHyb == 1 for n in qnet.hybrid]))
            hybrids = copy(qnet.hybrid)
            for n in hybrids
                if(n.typeHyb == 1) #only delete type 1 hybridizations (non identifiable ones)
                    !isa(n.prev,Nothing) || error("hybrid node $(n.number) is type 1 hybridization, prev should be automatically set")
                    eliminateHybridization!(qnet,n)
                end
            end
            qnet.typeHyb = Int[]
            if(qnet.numHybrids > 0)
                DEBUGC && @debug "need to identify hybridizations again after deleting type 1 hybridizations"
                identifyQuartet!(qnet)
            end
        end
        DEBUGC && @debug "now types are $([n.typeHyb for n in qnet.hybrid])"
        hybrids = copy(qnet.hybrid)
        for n in hybrids
            eliminateHybridization!(qnet,n)
        end
    end
    if(qnet.which == 1)
        internalLength!(qnet)
    end
end

# ------------------------- update qnet.formula

# function to identify to which of the 4 leaves
# two taxa correspond. this is to identify later
# if the two taxa correspond to major/minor CF
# input: qnet, taxon1,taxon2
# returns leaf for taxon1, leaf for taxon2 (i.e. 12)
# warning: assumes that the numbers for the taxon in the output.csv table are the names
function whichLeaves(qnet::QuartetNetwork, taxon1::String, taxon2::String, leaf1::Node, leaf2::Node, leaf3::Node, leaf4::Node)
    if(taxon1 == qnet.names[leaf1.number])
        if(taxon2 == qnet.names[leaf2.number])
            return 1,2
        elseif(taxon2 == qnet.names[leaf3.number])
            return 1,3
        elseif(taxon2 == qnet.names[leaf4.number])
            return 1,4
        else
            error("taxon2 $(taxon2) is not one of the three remaining leaves: $(leaf2.number), $(leaf3.number), $(leaf4.number)")
        end
    elseif(taxon1 == qnet.names[leaf2.number])
        if(taxon2 == qnet.names[leaf1.number])
            return 2,1
        elseif(taxon2 == qnet.names[leaf3.number])
            return 2,3
        elseif(taxon2 == qnet.names[leaf4.number])
            return 2,4
        else
            error("taxon2 $(taxon2) is not one of the three remaining leaves: $(leaf1.number), $(leaf3.number), $(leaf4.number)")
        end
    elseif(taxon1 == qnet.names[leaf3.number])
        if(taxon2 == qnet.names[leaf1.number])
            return 3,1
        elseif(taxon2 == qnet.names[leaf2.number])
            return 3,2
        elseif(taxon2 == qnet.names[leaf4.number])
            return 3,4
        else
            error("taxon2 $(taxon2) is not one of the three remaining leaves: $(leaf2.number), $(leaf1.number), $(leaf4.number)")
        end
    elseif(taxon1 == qnet.names[leaf4.number])
        if(taxon2 == qnet.names[leaf2.number])
            return 4,2
        elseif(taxon2 == qnet.names[leaf3.number])
            return 4,3
        elseif(taxon2 == qnet.names[leaf1.number])
            return 4,1
        else
            error("taxon2 $(taxon2) is not one of the three remaining leaves: $(leaf2.number), $(leaf3.number), $(leaf4.number)")
        end
    else
        error("taxon1: $(taxon1) is not one of the 4 leaves in quartet: $(leaf1.number), $(leaf2.number), $(leaf3.number), $(leaf4.number)")
    end
end

# function to update the attribute split in qnet
# qnet.leaf is a vector [x,y,z,w] and qnet.split is a vector
# [1,1,2,2] that says in which side of the split is each leaf
# warning: it needs to be run after eliminating hybridization and uniting
# internal edge
function updateSplit!(qnet::QuartetNetwork)
    if(qnet.which == 1)
        if(qnet.split == [-1,-1,-1,-1])
            qnet.split = [2,2,2,2]
            middle = qnet.node[getIndex(true,[size(n.edge,1) == 3 for n in qnet.node])]
            leaf1 = middle.edge[getIndex(true,[getOtherNode(e,middle).leaf for e in middle.edge])]
            leaf2 = middle.edge[getIndex(true,[(getOtherNode(e,middle).leaf && !isEqual(leaf1,e)) for e in middle.edge])]
            leaf1 = getOtherNode(leaf1,middle) #leaf1 was edge, now it is node
            leaf2 = getOtherNode(leaf2,middle)
            ind1 = getIndex(leaf1,qnet.leaf)
            ind2 = getIndex(leaf2,qnet.leaf)
            qnet.split[ind1] = 1
            qnet.split[ind2] = 1
        end
    elseif(qnet.which == -1)
        error("cannot update split in quartet network if it has not been identified, and eliminated hybridizations")
    end
end

# function to know which formula (minor/major) to compute
# for qnet.expCF[1,2,3] depending on the order of taxa in
# qnet.quartet
# warning: needs qnet.split updated already
function updateFormula!(qnet::QuartetNetwork)
    if(qnet.which == 1)
        if(qnet.formula == [-1,-1,-1])
            qnet.split != [-1,-1,-1,-1] || error("cannot update qnet.formula if qnet.split is not updated: $(qnet.split)")
            qnet.formula = [2,2,2]
            for i in 2:4
                size(qnet.leaf,1) == 4 || error("strange quartet with $(size(qnet.leaf,1)) leaves instead of 4")
                tx1,tx2 = whichLeaves(qnet,qnet.quartetTaxon[1],qnet.quartetTaxon[i], qnet.leaf[1], qnet.leaf[2], qnet.leaf[3], qnet.leaf[4]) # index of leaf in qnet.leaf
                if(qnet.split[tx1] == qnet.split[tx2])
                    qnet.formula[i-1] = 1
                    break
                end
            end
        end
    elseif(qnet.which != 2)
        error("qnet.which should be updated already to 1 or 2, and it is $(qnet.which)")
    end
end


# --------------- calculate exp CF ----------------------


# function to calculate expCF for a quartet network
# warning: needs qnet.formula and qnet.t1 already updated
function calculateExpCF!(qnet::QuartetNetwork)
    if(qnet.which == 1)
        if(qnet.formula != [-1,-1,-1] && qnet.t1 != -1)
            for i in 1:3
                qnet.expCF[i] = qnet.formula[i] == 1 ? 1-2/3*exp(-qnet.t1) : 1/3*exp(-qnet.t1)
            end
        else
            error("quartet network needs to have updated formula and t1 before computing the expCF")
        end
    elseif(qnet.which == 2)
        if(qnet.numHybrids == 1)
            if(qnet.hybrid[1].typeHyb == 5)
                quartetType5!(qnet,qnet.hybrid[1])
            else
                error("strange quartet network type $(qnet.which) with one hybrid node $(qnet.hybrid[1].number) but it is not type 5, it is type $(qnet.hybrid[1].typeHyb)")
            end
        else
            error("quartet network with type $(qnet.which) but with $(qnet.numHybrids) hybrid nodes. all hybridizations type 1 (not identifiable) have been eliminated already, so there should only be one hybrid.")
        end
    end
end

# function to compute all the process of calculating the expCF
# for a given qnet
function calculateExpCFAll!(qnet::QuartetNetwork)
    identifyQuartet!(qnet)
    eliminateHybridization!(qnet)
    updateSplit!(qnet)
    updateFormula!(qnet)
    calculateExpCF!(qnet)
end

# function to calculate expCF for all the quartets in data
# after extractQuartet(net,data) that updates quartet.qnet
# warning: only updates expCF for quartet.changed=true
function calculateExpCFAll!(data::DataCF)
    !all((q->(q.qnet.numTaxa != 0)), data.quartet) ? error("qnet in quartets on data are not correctly updated with extractQuartet") : nothing
    #@warn "assume the numbers for the taxon read from the observed CF table match the numbers given to the taxon when creating the object network"
    for q in data.quartet
        if(q.qnet.changed)
            qnet = deepcopy(q.qnet);
            calculateExpCFAll!(qnet);
            q.qnet.expCF = qnet.expCF
        end
    end
end


# function to calculate expCF for all the quartets in data
# after extractQuartet(net,data) that updates quartet.qnet
# first updates the edge lengths according to x
# warning: assumes qnet.indexht is updated already
# warning: only updates expCF for quartet.qnet.changed=true
function calculateExpCFAll!(data::DataCF, x::Vector{Float64},net::HybridNetwork)
    !all((q->(q.qnet.numTaxa != 0)), data.quartet) ? error("qnet in quartets on data are not correctly updated with extractQuartet") : nothing
    #println("calculateExpCFAll in x: $(x) with net.ht $(net.ht)")
    for q in data.quartet
        update!(q.qnet,x,net)
        if(q.qnet.changed)
            #println("enters to recalculate expCF for some quartet")
            qnet = deepcopy(q.qnet);
            calculateExpCFAll!(qnet);
            q.qnet.expCF = qnet.expCF
        end
    end
end

# function to simply calculate the pseudolik of a given network
"""
`topologyQPseudolik!(net::HybridNetwork, d::DataCF)`

Calculate the quartet pseudo-deviance of a given network/tree for
DataCF `d`. This is the negative log pseudo-likelihood,
up to an additive constant, such that a perfect fit corresponds to a deviance of 0.0.

Be careful if the net object does
not have all internal branch lengths specified because then the
pseudolikelihood will be meaningless.

The loglik attribute of the network is undated, and `d` is updated with the expected
concordance factors under the input network.
"""
function topologyQPseudolik!(net0::HybridNetwork,d::DataCF; verbose=false::Bool)
    for ed in net0.edge
      !ed.hybrid || (ed.gamma >= 0.0) ||
        error("hybrid edge has missing γ value. Cannot compute quartet pseudo-likelihood.\nTry `topologyMaxQPseudolik!` instead, to estimate these γ's.")
    end
    missingBL = any([e.length < 0.0 for e in net0.edge]) # at least one BL was missing
    net = readTopologyUpdate(writeTopologyLevel1(net0))  # update level-1 attributes. Changes <0 BL into 1.0
    if(!isempty(d.repSpecies))
      expandLeaves!(d.repSpecies, net)
      net = readTopologyLevel1(writeTopologyLevel1(net)) # dirty fix to multiple alleles problem with expandLeaves
    end
    missingBL && any([(e.length == 1.0 && e.istIdentifiable) for e in net.edge]) &&
      @warn "identifiable edges lengths were originally missing, so assigned default value of 1.0"
    try
        checkNet(net)
    catch
        error("starting topology not a level 1 network")
    end
    extractQuartet!(net,d) # quartets are all updated: hasEdge, expCF, indexht
    all((q->(q.qnet.numTaxa != 0)), d.quartet) || error("qnet in quartets on data are not correctly updated with extractQuartet")
    for q in d.quartet
        if verbose println("computing expCF for quartet $(q.taxon)") # to stdout
        else @debug        "computing expCF for quartet $(q.taxon)"; end # to logger if debug turned on by user 
        qnet = deepcopy(q.qnet);
        calculateExpCFAll!(qnet);
        q.qnet.expCF = qnet.expCF
        if verbose println("$(qnet.expCF)") # to stdout
        else @debug        "$(qnet.expCF)"; end # to logger
    end
    val = logPseudoLik(d)
    if verbose println("$the value of pseudolikelihood is $(val)") # to stdout
    else @debug        "$the value of pseudolikelihood is $(val)"; end # to logger
    net0.loglik = val
    return val
end


# ---------------------------- Pseudolik for a quartet -------------------------

# function to calculate the log pseudolikelihood function for a single
# quartet
# sum_i=1,2,3 (obsCF_i)*log(expCF_i/obsCF_i)
# warning: assumes that quartet.qnet is already updated with extractQuartet and
#          calculateExpCF
function logPseudoLik(quartet::Quartet)
    sum(quartet.qnet.expCF) != 0.0 || error("expCF not updated for quartet $(quartet.number)")
    #@debug "quartet= $(quartet.taxon), obsCF = $(quartet.obsCF), expCF = $(quartet.qnet.expCF)"
    suma = 0
    for i in 1:3
        if(quartet.qnet.expCF[i] < 0)
            @debug "found expCF negative $(quartet.qnet.expCF[i]), will set loglik=-1.e15"
            suma += -1.e15
        else
            suma += quartet.obsCF[i] == 0 ? 0.0 : 100*quartet.obsCF[i]*log(quartet.qnet.expCF[i]/quartet.obsCF[i])
            # WARNING: 100 should be replaced by -2*ngenes to get the deviance.
            # below: negative sign used below in logPseudoLik() when summing up across 4-taxon sets
        end
    end
    ## to account for missing data:
    ## if(quartet.ngenes > 0)
    ##     suma = quartet.ngenes*suma
    ## end
    quartet.logPseudoLik = suma
    return suma
end

# function to calculate the -log pseudolikelihood function for array of
# quartets
# sum_q [sum_i=1,2,3 (obsCF_i)*log(expCF_i)]
# warning: assumes that quartet.qnet is already updated with extractQuartet and
#          calculateExpCF for all quartets
function logPseudoLik(quartet::Array{Quartet,1})
    suma = 0
    for q in quartet
        suma += logPseudoLik(q)
    end
    return -suma
end

logPseudoLik(d::DataCF) = logPseudoLik(d.quartet)

