# functions to delete hybridization
# originally in functions.jl
# Claudia MArch 2015

# --------------------------------- delete hybridization -------------------------------

# function to identify inCycle (with priority queue)
# based on updateInCycle as it is the exact same code
# only without changing the incycle attribute
# only returning array of edges/nodes affected by the hybrid
# used when attempting to delete
# input: hybrid node around which we want to identify inCycle
# returns tuple: nocycle, array of edges changed, array of nodes changed
# check: is this traversal much faster than a simple loop over
#        all edges/nodes and check if incycle==hybrid.number?
function identifyInCycle(net::Network,node::Node)
    node.hybrid || error("node $(node.number) is not hybrid, cannot identifyInCycle")
    start = node;
    hybedge = getHybridEdge(node);
    last = getOtherNode(hybedge,node);
    dist = 0;
    queue = PriorityQueue();
    path = Node[];
    net.edges_changed = Edge[];
    net.nodes_changed = Node[];
    push!(net.edges_changed,hybedge);
    push!(net.nodes_changed,node);
    found = false;
    net.visited = [false for i = 1:size(net.node,1)];
    enqueue!(queue,node,dist);
    while(!found)
        if(isempty(queue))
            return true, net.edges_changed, net.nodes_changed
        else
            curr = dequeue!(queue);
            if(isEqual(curr,last))
                found = true;
                push!(path,curr);
            else
                if(!net.visited[getIndex(curr,net)])
                    net.visited[getIndex(curr,net)] = true;
                    if(isEqual(curr,start))
                        for e in curr.edge
                            if(!e.hybrid || e.isMajor)
                                other = getOtherNode(e,curr);
                                other.prev = curr;
                                dist = dist+1;
                                enqueue!(queue,other,dist);
                            end
                        end
                    else
                        for e in curr.edge
                            if(!e.hybrid || e.isMajor)
                                other = getOtherNode(e,curr);
                                if(!other.leaf && !net.visited[getIndex(other,net)])
                                    other.prev = curr;
                                    dist = dist+1;
                                    enqueue!(queue,other,dist);
                                end
                            end
                        end
                    end
                end
            end
        end
    end # end while
    curr = pop!(path);
    while(!isEqual(curr, start))
        if(curr.inCycle == start.number)
            push!(net.nodes_changed, curr);
            edge = getConnectingEdge(curr,curr.prev);
            if(edge.inCycle == start.number)
                push!(net.edges_changed, edge);
            end
            curr = curr.prev;
        end
    end
    return false, net.edges_changed, net.nodes_changed
end



# aux function to traverse the network
# similar to traverseContainRoot but only
# identifying the edges that would be changed by a given
# hybridization
# warning: it does not go accross hybrid node/edge,
#          nor tree node with minor hybrid edge
function traverseIdentifyRoot(node::Node, edge::Edge, edges_changed::Array{Edge,1})
    if(!node.leaf && !node.hybrid)
        for e in node.edge
            if(!isEqual(edge,e) && e.isMajor && !e.hybrid)
                other = getOtherNode(e,node);
                push!(edges_changed, e);
                if(!other.hybrid)
                    #if(!other.hasHybEdge)
                        traverseIdentifyRoot(other,e, edges_changed);
                    #else
                    #    if(hybridEdges(other)[1].isMajor)
                    #        traverseIdentifyRoot(other,e, edges_changed);
                    #    end
                    #end
                end
            end
        end
    end
end


# function to identify containRoot
# depending on a hybrid node on the network
# input: hybrid node (can come from searchHybridNode)
# return array of edges affected by the hybrid node
function identifyContainRoot(net::HybridNetwork, node::Node)
    node.hybrid || error("node $(node.number) is not hybrid, cannot identify containRoot")
    net.edges_changed = Edge[];
    for e in node.edge
        if(!e.hybrid)
            other = getOtherNode(e,node);
            push!(net.edges_changed,e);
            traverseIdentifyRoot(other,e, net.edges_changed);
        end
    end
    return net.edges_changed
end

# function to undo the effect of a hybridization
# and then delete it
# input: network, hybrid node, random flag
#        random = true, deletes one hybrid egde at random
#        (minor with prob 1-gamma, major with prob gamma)
#        random = false, deletes the minor edge always
# warning: it uses the gamma of the hybrid edges even if
#          it is not identifiable like in bad diamond I (assumes undone by now)
# blacklist = true: add the edge as a bad choice to put a hybridization (not fully tested)
function deleteHybridizationUpdate!(net::HybridNetwork, hybrid::Node, random::Bool, blacklist::Bool)
    hybrid.hybrid || error("node $(hybrid.number) is not hybrid, so we cannot delete hybridization event around it")
    @debug "MOVE: delete hybridization on hybrid node $(hybrid.number)"
    nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,hybrid);
    !nocycle || error("the hybrid node $(hybrid.number) does not create a cycle")
    edgesRoot = identifyContainRoot(net,hybrid);
    edges = hybridEdges(hybrid);
    undoGammaz!(hybrid,net);
    othermaj = getOtherNode(edges[1],hybrid)
    edgesmaj = hybridEdges(othermaj)
    @debug "edgesmaj[3] $(edgesmaj[3].number) is the one to check if containRoot=false already: $(edgesmaj[3].containRoot)"
    if(edgesmaj[3].containRoot) #if containRoot=true, then we need to undo
        push!(edgesRoot, edges[1]) ## add hybrid edges to edgesRoot to undo containRoot
        push!(edgesRoot, edges[2])
        undoContainRoot!(edgesRoot);
    end
    @debug begin
        msg = ""
        if edges[1].gamma < 0.5
            msg = "strange major hybrid edge $(edges[1].number) with gamma $(edges[1].gamma) less than 0.5"
        end
        if edges[1].gamma == 1.0
            msg = "strange major hybrid edge $(edges[1].number) with gamma $(edges[1].gamma) equal to 1.0"
        end
        msg
    end
    limit = edges[1].gamma
    if(random)
        minor = rand() < limit ? false : true
    else
        minor = true;
    end
    deleteHybrid!(hybrid,net,minor, blacklist)
    undoInCycle!(edgesInCycle, nodesInCycle); #moved after deleteHybrid to mantain who is incycle when deleteEdge and look for partition
    undoPartition!(net,hybrid, edgesInCycle)
end

deleteHybridizationUpdate!(net::HybridNetwork, hybrid::Node) = deleteHybridizationUpdate!(net, hybrid, true, false)

# function to delete a hybridization event
# input: hybrid node and network
#        minor: true (deletes minor edge), false (deletes major)
# warning: it is meant after undoing the effect of the
#          hybridization in deleteHybridizationUpdate!
#          by itself, it leaves things as is
# branch lengths of -1.0 are interpreted as missing.
function deleteHybrid!(node::Node,net::HybridNetwork,minor::Bool, blacklist::Bool)
    node.hybrid || error("node $(node.number) has to be hybrid for deleteHybrid")
    if(minor)
        hybedge1,hybedge2,treeedge1 = hybridEdges(node);
        other1 = getOtherNode(hybedge1,node);
        other2 = getOtherNode(hybedge2,node);
        other3 =  getOtherNode(treeedge1,node);
        if(hybedge1.number > treeedge1.number)
            setLength!(treeedge1, addBL(treeedge1.length, hybedge1.length));
            removeNode!(node,treeedge1);
            setNode!(treeedge1,other1);
            setEdge!(other1,treeedge1);
            removeEdge!(other1, hybedge1);
            deleteEdge!(net,hybedge1);
            #treeedge1.containRoot = (!treeedge1.containRoot || !hybedge1.containRoot) ? false : true #causes problems if hybrid.CR=false
            if(blacklist)
                println("put in blacklist edge $(treeedge1.number)")
                push!(net.blacklist, treeedge1.number)
            end
        else
            makeEdgeTree!(hybedge1,node)
            other1.hasHybEdge = false;
            setLength!(hybedge1, addBL(hybedge1.length, treeedge1.length));
            removeNode!(node,hybedge1);
            setNode!(hybedge1,other3);
            setEdge!(other3,hybedge1);
            removeEdge!(other3,treeedge1);
            deleteEdge!(net,treeedge1);
            hybedge1.containRoot = (!treeedge1.containRoot || !hybedge1.containRoot) ? false : true
            if(blacklist)
                println("put in blacklist edge $(hybedge1.number)")
                push!(net.blacklist, hybedge1.number)
            end
        end
        hybindex = getIndex(true,[e.hybrid for e in other2.edge]);
        if(hybindex == 1)
            treeedge1 = other2.edge[2];
            treeedge2 = other2.edge[3];
        elseif(hybindex == 2)
            treeedge1 = other2.edge[1];
            treeedge2 = other2.edge[3];
        elseif(hybindex == 3)
            treeedge1 = other2.edge[1];
            treeedge2 = other2.edge[2];
        else
            error("strange node has more than three edges")
        end
        treenode1 = getOtherNode(treeedge1,other2);
        treenode2 = getOtherNode(treeedge2,other2);
        if(abs(treeedge1.number) > abs(treeedge2.number))
            setLength!(treeedge2, addBL(treeedge2.length, treeedge1.length));
            removeNode!(other2,treeedge2);
            setNode!(treeedge2,treenode1);
            setEdge!(treenode1,treeedge2);
            removeEdge!(treenode1,treeedge1);
            deleteEdge!(net,treeedge1);
            treeedge2.containRoot = (!treeedge1.containRoot || !treeedge2.containRoot) ? false : true
            if(blacklist)
                println("put in blacklist edge $(treeedge2.number)")
                push!(net.blacklist, treeedge2.number)
            end
        else
            setLength!(treeedge1, addBL(treeedge2.length, treeedge1.length));
            removeNode!(other2,treeedge1);
            setNode!(treeedge1,treenode2);
            setEdge!(treenode2,treeedge1);
            removeEdge!(treenode2,treeedge2);
            deleteEdge!(net,treeedge2);
            treeedge1.containRoot = (!treeedge1.containRoot || !treeedge2.containRoot) ? false : true
            if(blacklist)
                println("put in blacklist edge $(treeedge1.number)")
                push!(net.blacklist, treeedge1.number)
            end
        end
        #removeHybrid!(net,node);
        deleteNode!(net,node);
        deleteNode!(net,other2);
        deleteEdge!(net,hybedge2);
    else
        hybedge1,hybedge2,treeedge1 = hybridEdges(node);
        other1 = getOtherNode(hybedge1,node);
        other2 = getOtherNode(hybedge2,node);
        setLength!(treeedge1, addBL(treeedge1.length, hybedge2.length))
        removeEdge!(other2,hybedge2)
        removeNode!(node,treeedge1)
        setEdge!(other2,treeedge1)
        setNode!(treeedge1,other2)
        #removeHybrid!(net,node)
        deleteNode!(net,node)
        deleteEdge!(net,hybedge1)
        deleteEdge!(net,hybedge2)
        removeEdge!(other1,hybedge1)
        size(other1.edge,1) == 2 || error("strange node $(other1.number) had 4 edges")
        if(abs(other1.edge[1].number) < abs(other1.edge[2].number))
            edge = other1.edge[1]
            otheredge = other1.edge[2]
        else
            edge = other1.edge[2]
            otheredge = other1.edge[1]
        end
        setLength!(other1.edge[1], addBL(other1.edge[1].length, other1.edge[2].length))
        other3 =  getOtherNode(otheredge,other1);
        removeNode!(other1,edge)
        removeEdge!(other3,otheredge)
        setEdge!(other3,edge)
        setNode!(edge,other3)
        deleteNode!(net,other1)
        deleteEdge!(net,otheredge)
    end
end

deleteHybrid!(node::Node,net::HybridNetwork,minor::Bool) = deleteHybrid!(node,net,minor, false)

"""
    deleteHybridEdge!(net::HybridNetwork, edge::Edge, keepNodes=false)

Deletes a hybrid edge from a network. The network does not have to be of level 1,
and may contain some polytomies. Updates branch lengths, allowing for missing values.
Returns the network.

At each of the 2 junctions, the child edge is retained (below the hybrid node).
If `keepNodes` is true, all nodes are retained during edge removal.

Warnings:

- if `keepNodes` is true: partner hybrid parent edge has its γ value unchanged
- if the parent of `edge` is the root and if `keepNodes` is false, the root
  is moved to keep the network unrooted with a root of degree two.
- does **not** update containRoot (could be implemented later)
- does **not** update attributes needed for snaq! (like containRoot, inCycle, edge.z, edge.y etc.)
"""
function deleteHybridEdge!(net::HybridNetwork, edge::Edge, keepNodes=false::Bool)
    edge.hybrid || error("edge $(edge.number) has to be hybrid for deleteHybridEdge!")
    n1 = getChild(edge)  # child  of edge, to be deleted
    n1.hybrid || error("child node $(n1.number) of hybrid edge $(edge.number) should be a hybrid.")
    n2 = getParent(edge)  # parent of edge, to be deleted too.
    # next: keep hybrid node n1 if it has 4+ edges (2 parents and 2+ children).
    #       2 or 1 edges should never occur.
    if length(n1.edge) < 3
        error("node $(n1.number) has $length(n1.edge) edges instead of 3+");
    elseif length(n1.edge) == 3 && !keepNodes
        pe = nothing # will be other parent (hybrid) edge of n1
        ce = nothing # will be child edge of n1, to be merged with pe
        for e in n1.edge
            if e.hybrid && e!=edge && n1==getChild(e) pe = e; end
            if !e.hybrid || n1==getParent(e)  ce = e; end
            # ce could be a hybrid edge (if so, not level-1 network). If tree edge, its isChild1 may be outdated
        end
        pn = getParent(pe); # parent node of n1, other than n2
        pn ≢ n2 || error("k=2 cycle: 2 hybrid edges from node $(n2.number) to node $(n1.number)")
        atRoot = (net.node[net.root] ≡ n1) # n1 should not be root, but if so, pn will be new root
        # next: replace ce by pe+ce, remove n1 and pe from network.
        ce.length = addBL(ce.length, pe.length)
        removeNode!(n1,ce) # ce now has 1 single node cn
        setNode!(ce,pn)    # ce now has 2 nodes in this order: cn, pn
        ce.isChild1 = true
        setEdge!(pn,ce)
        removeEdge!(pn,pe)
        # if (pe.number<ce.number) ce.number = pe.number; end # bad to match edges between networks
        deleteEdge!(net,pe,part=false) # decreases net.numEdges   by 1
        deleteNode!(net,n1) # decreases net.numHybrids by 1, numNodes too.
        # warning: containRoot could be updated in ce and down the tree.
        if atRoot
            i = findfirst(x -> x===pn, net.node)
            i !== nothing || error("node $(pn.number) not in net!")
            net.root = i
        end
    else # n1 has 4+ edges (polytomy) or 3 edges but we want to keep it anyway:
        # keep n1 but detach it from 'edge', set its remaining parent to major tree edge
        pe = getPartner(edge, n1) # partner edge: keep it this time
        if !pe.isMajor pe.isMajor=true; end
        pe.hybrid = false
        # note: pe.gamma *not* set to 1.0 here
        removeEdge!(n1,edge) # does not update n1.hybrid at this time
        removeHybrid!(net,n1) # removes n1 from net.hybrid, updates net.numHybrids
        n1.hybrid = false
    end

    # next: keep n2 if it has 4+ edges (or if keepNodes). 2 or 1 edges should never occur.
    #       If root, would have no parent: treat network as unrooted and change the root.
    if length(n2.edge) < 2
        error("node $(n2.number) (parent of hybrid edge $(edge.number) to be deleted) has 1 edge only!")
    elseif length(n2.edge) == 2 && !keepNodes
        # n2=root or degree-2 node: remove n2 and both of its edges (including 'edge')
        pe = (edge ≡ n2.edge[1] ? n2.edge[2] : n2.edge[1]) #  <--edge-- n2 ---pe--- pn
        pn = (  n2 ≡ pe.node[1] ? pe.node[2] : pe.node[1])
        if net.node[net.root] ≡ n2 # if n2 was root, new root = pn
            net.root = findfirst(x -> x===pn, net.node)
        end
        # remove n2 and pe
        removeEdge!(pn,pe)
        deleteEdge!(net,pe,part=false)
        deleteNode!(net,n2)
    elseif length(n2.edge) == 3 && !keepNodes && !n2.hybrid # delete n2, merge its 2 parent edges
        # if n2 is hybrid: its 2 parents cannot be merged: direction conflict
        ce, pe = [e for e in n2.edge if e != edge]
        # ce will be kept. pe will be folded into new ce = pe+ce
        switch = false
        if getOtherNode(pe, n2).leaf # difficulty: isChild1 may be outdated for tree edges
            !getOtherNode(ce, n2).leaf ||
                error("root $(n2.number) connected to 1 hybrid and 2 leaves: edges $(pe.number) and $(ce.number).")
            switch = true
        elseif n2 ≡ getChild(ce) && n2 ≡ getParent(pe)
            switch = true
        end
        if switch # to ensure correct direction isChild1 for new edges if the original was up-to-date
            # but no error if original isChild1 was outdated
            # and to ensure that pn can be the new root if needed.
            ce,pe = pe,ce
        end
        atRoot = (net.node[net.root] ≡ n2) # if n2=root, new root will be 'pn' = other node of pe
        # next: replace ce by pe+ce, remove n2 and pe from network.
        pn = getOtherNode(pe,n2) # parent node of n2 if n2 not root. Otherwise, pn will be new root.
        ce.length = addBL(ce.length, pe.length)
        removeNode!(n2,ce) # ce now has 1 single node cn
        setNode!(ce,pn)    # ce now has 2 nodes in this order: cn, pn
        ce.isChild1 = true
        setEdge!(pn,ce)
        removeEdge!(pn,pe)
        # if (pe.number<ce.number) ce.number = pe.number; end # bad to match edges between networks
        deleteEdge!(net,pe,part=false)
        deleteNode!(net,n2)
        if (atRoot)
            try
                net.root = getIndex(pn,net)
            catch e
                if isa(e, ErrorException) error("node $(pn.number) not in net!"); end
            end
        end
    elseif length(n2.edge) == 3 && !keepNodes # n2 is a hybrid node then. "edge" is its single child. Not level-1
        @warn "case not fully implemented yet: hybrid node $(n2.number) will be kept"
        removeEdge!(n2,edge)
        # fixit: mark n2 as leaf, but do not add to net.leaf?
        # fixit: deleteleaf!(net, n2; simplify) to delete the hybrid leaf & its parents edges??
        #        but delete "edge" before (not below), because need to simplify again recursively
    else # n2 has 4+ edges (polytomy) or keepNodes: keep n2 but detach it from 'edge'
        removeEdge!(n2,edge)
    end
    # finally: remove hybrid 'edge' from network
    deleteEdge!(net,edge,part=false)
    return net
end

# function to update net.partition after deleting a hybrid node
# needs a list of the edges in cycle
function undoPartition!(net::HybridNetwork, hybrid::Node, edgesInCycle::Vector{Edge})
    hybrid.hybrid || error("node $(hybrid.number) is not hybrid, and we need hybrid node inside deleteHybUpdate for undoPartition")
    if(net.numHybrids == 0)
        net.partition = Partition[]
    else
        cycles = Int[]
        edges = Edge[]
        N = length(net.partition)
        i = 1
        while(i <= N)
            @debug "hybrid number is $(hybrid.number) and partition is $([e.number for e in net.partition[i].edges]), with cycle $(net.partition[i].cycle)"
            if(in(hybrid.number,net.partition[i].cycle))
                @debug "hybrid number matches with partition.cycle"
                p = splice!(net.partition,i)
                @debug "after splice, p partition has edges $([e.number for e in p.edges]) and cycle $(p.cycle)"
                ind = getIndex(hybrid.number,p.cycle)
                deleteat!(p.cycle,ind) #get rid of that hybrid number
                cycles = vcat(cycles,p.cycle)
                edges = vcat(edges,p.edges)
                @debug "edges is $([e.number for e in edges]) and cycles is $(cycles)"
                N = length(net.partition)
            else
                i += 1
            end
        end
        for e in edgesInCycle
            @debug "edge in cycle is $(e.number)"
            if(isEdgeNumIn(e,net.edge)) #only include edge if still in net
                @debug "edge is in net still"
                push!(edges,e)
            end
        end
        newPartition = Partition(unique(cycles),edges)
        @debug "new partition with cycle $(newPartition.cycle), edges $([e.number for e in newPartition.edges])"
        push!(net.partition,newPartition)
    end
end

