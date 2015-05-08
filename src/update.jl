# functions to update incycle, containRoot, gammaz
# originally in functions.jl
# Claudia March 2015
#####################

# --------------------------------------- update incycle, root, gammaz -------------------------------------------

# function to update inCycle (with priority queue) after becoming part of a network
# based on program 3 CS367 with priority queue
# expected to be much faster than the other two udpateInCycle (queue and recursive)
# input: hybrid node around which we want to update inCycle
# needs module "Base.Collections"
# returns tuple: flag, nocycle, array of edges changed, array of nodes changed
#   flag: false if cycle intersects existing cycle or number of nodes in cycle < 3
#         (there is the possibility of returning edges in intersection: path)
#         true if cycle does not intersect existing cycle
#   nocycle: true if there is no cycle (instead of error). it is used in addHybridization
# calculates also the number of nodes in the cycle and put as hybrid node attribute "k"
# warning: it is not checking if hybrid node or minor hybrid edge
#          were already part of a cycle (inCycle!= -1)
#          But it is checking so for the other edges in cycle
# warning: it needs extra things: visited attribute, prev attribute
#          unlike updateInCycle recursive, but it is expected
#          to be much faster
function updateInCycle!(net::HybridNetwork,node::Node)
    if(node.hybrid)
        start = node;
        node.inCycle = node.number;
        node.k = 1;
        hybedge = getHybridEdge(node);
        hybedge.inCycle = node.number;
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
        println("start is $(start.number), last is $(last.number)")
        while(!found)
            if(isempty(queue))
                return false, true, net.edges_changed, net.nodes_changed
            else
                #println("at this moment, queue is $([n.number for n in queue])")
                curr = dequeue!(queue);
                println("curr $(curr.number)")
                if(isEqual(curr,last))
                    println("encuentra a la last $(curr.number)")
                    found = true;
                    push!(path,curr);
                else
                    if(!net.visited[getIndex(curr,net)])
                        println("curr not visited $(curr.number)")
                        net.visited[getIndex(curr,net)] = true;
                        if(isEqual(curr,start))
                            println("curr is start")
                            for(e in curr.edge)
                                if(!e.hybrid || e.isMajor)
                                    other = getOtherNode(e,curr);
                                    println("other is $(other.number), pushed to queue")
                                    other.prev = curr;
                                    dist = dist+1;
                                    enqueue!(queue,other,dist);
                                end
                            end
                        else
                            for(e in curr.edge)
                                if(!e.hybrid || e.isMajor)
                                    other = getOtherNode(e,curr);
                                    println("other is $(other.number)")
                                    if(!other.leaf && !net.visited[getIndex(other,net)])
                                        println("dice que other $(other.number) no es leaf ni visited, lo mete a queue")
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
        println("after while, path has $([n.number for n in path])")
        curr = pop!(path);
        println("path should be empty here: $(isempty(path))")
        while(!isEqual(curr, start))
            println("curr es $(curr.number)")
            if(curr.inCycle!= -1)
                push!(path,curr);
                curr = curr.prev;
            else
                curr.inCycle = start.number;
                push!(net.nodes_changed, curr);
                node.k  =  node.k + 1;
                edge = getConnectingEdge(curr,curr.prev);
                edge.inCycle = start.number;
                push!(net.edges_changed, edge);
                curr = curr.prev;
            end
        end
        if(!isempty(path) || node.k<3)
            warn("new cycle intersects existing cycle")
            return false, false, net.edges_changed, net.nodes_changed
        else
            return true, false, net.edges_changed, net.nodes_changed
        end
    else
        error("node is not hybrid")
    end
end


# aux function to traverse the network for updateContainRoot
# it changes the containRoot argument to false
# of all the edges visited
# changed to recursive after Cecile's idea
# warning: it does not go accross hybrid node, minor hybrid edge
# fixit: should have exclamation sign
function traverseContainRoot(node::Node, edge::Edge, edges_changed::Array{Edge,1})
    if(!node.leaf && !node.hybrid)
        for(e in node.edge)
            if(!isEqual(edge,e) && e.isMajor)
                other = getOtherNode(e,node);
                if(e.containRoot) # only considered changed those that were true and not hybrid
                    e.containRoot = false;
                    push!(edges_changed, e);
                end
                traverseContainRoot(other,e, edges_changed);
            end
        end
    end
end


# function to update containRoot (by default true)
# depending on the network
# input: hybrid node (can come from searchHybridNode)
# return flag, array of edges changed
#        flag: false if the set of edges to place the root is empty
function updateContainRoot!(net::HybridNetwork, node::Node)
    node.hybrid || error("node $(node.number )is not hybrid, cannot update containRoot")
    net.edges_changed = Edge[];
    for (e in node.edge)
        if(!e.hybrid)
            other = getOtherNode(e,node);
            e.containRoot = false;
            push!(net.edges_changed,e);
            traverseContainRoot(other,e, net.edges_changed);
        end
    end
    if(all([!e.containRoot for e in net.edge]))
        return false,net.edges_changed
    else
        return true,net.edges_changed
    end
end

# function to identify if the network is one of the pathological cases
# see ipad notes: k = 0 (nonidentifiable), k = 1 (nonidentifiable
# "extreme bad triangle", "very" bad triangle I,II) k = 2 (bad diamond
# I,II) also checks if hybrid node has leaf child, in which case,
# major edge is non identifiable
# input: hybrid node around which to check (can come from searchHybridNode)
#        updates gammaz with whatever
# edge lengths are originally in the network
#        allow = true, returns true always, used when reading topology
# returns net.hasVeryBadTriangle, array of edges changed (istIdentifiable, except hybrid edges)
#         false if the network has extremely/very bad triangles
# warning: needs to have updateInCycle already done as it needs inCycle, and k
# check: assume any tree node that has hybrid Edge has only
# one tree edge in cycle (true?)
# updates net.numBad attribute when found a bad diamond I
function updateGammaz!(net::HybridNetwork, node::Node, allow::Bool)
    node.hybrid || error("node $(node.number) is not hybrid, cannot updategammaz")
    node.k != -1 || error("udpate in cycle should have been run before: node.k not -1")
    node.isExtBadTriangle = false
    node.isVeryBadTriangle = false
    node.isBadTriangle = false
    node.isBadDiamondI = false
    node.isBadDiamondII = false
    net.edges_changed = Edge[];
    edge_maj, edge_min, tree_edge2 = hybridEdges(node);
    other_maj = getOtherNode(edge_maj,node);
    other_min = getOtherNode(edge_min,node);
    node.k > 2 || return false, []
    if(node.k == 4) # could be bad diamond I,II
        net.numTaxa >= 5 || return false, []
        edgebla,edge_min2,tree_edge3 = hybridEdges(other_min);
        edgebla,edge_maj2,tree_edge1 = hybridEdges(other_maj);
        other_min2 = getOtherNode(edge_min2,other_min);
        isLeaf1 = getOtherNode(tree_edge1,other_maj);
        isLeaf2 = getOtherNode(tree_edge2,node);
        isLeaf3 = getOtherNode(tree_edge3,other_min);
        tree_edge4 = nothing;
        for(e in other_min2.edge)
            if(isa(tree_edge4,Nothing) && e.inCycle == -1 && !e.hybrid)
                tree_edge4 = e;
            end
        end
        if(isEqual(other_min2,getOtherNode(edge_maj2,other_maj)) && isLeaf1.leaf && isLeaf2.leaf && isLeaf3.leaf) # bad diamond I
            warn("bad diamond I found")
            net.numBad += 1
            node.isBadDiamondI = true;
            other_min.gammaz = edge_min.gamma*edge_min2.z;
            other_maj.gammaz = edge_maj.gamma*edge_maj2.z;
            edge_min2.istIdentifiable = false;
            edge_maj2.istIdentifiable = false;
            edge_maj.istIdentifiable = false;
            edge_min.istIdentifiable = false;
            push!(net.edges_changed,edge_min2);
            push!(net.edges_changed,edge_min);
            push!(net.edges_changed,edge_maj2);
            push!(net.edges_changed,edge_maj);
        elseif(isEqual(other_min2,getOtherNode(edge_maj2,other_maj)) && isLeaf1.leaf && !isLeaf2.leaf && isLeaf3.leaf && getOtherNode(tree_edge4,other_min2).leaf) # bad diamond II
            warn("bad diamond II found")
            node.isBadDiamondII = true;
            setLength!(edge_maj,edge_maj.length+tree_edge2.length)
            setLength!(tree_edge2,0.0)
            push!(net.edges_changed,tree_edge2)
            tree_edge2.istIdentifiable = false
        end
    elseif(node.k == 3) # could be extreme/very bad triangle or just bad triangle
        if(net.numTaxa <= 5)
            warn("extremely or very bad triangle found")
            node.isVeryBadTriangle = true
            net.hasVeryBadTriangle = true
        elseif(net.numTaxa >= 6)
            edgebla,tree_edge_incycle,tree_edge1 = hybridEdges(other_min);
            edgebla,edgebla,tree_edge3 = hybridEdges(other_maj);
            isLeaf1 = getOtherNode(tree_edge1,other_min);
            isLeaf2 = getOtherNode(tree_edge2,node);
            isLeaf3 = getOtherNode(tree_edge3,other_maj);
            if(isLeaf1.leaf || isLeaf2.leaf || isLeaf3.leaf)
                if(sum([l.leaf?1:0 for l in [isLeaf1,isLeaf2,isLeaf3]]) >= 2)
                    warn("extremely bad triangle found")
                    node.isExtBadTriangle = true;
                    net.hasVeryBadTriangle = true
                elseif(sum([l.leaf?1:0 for l in [isLeaf1,isLeaf2,isLeaf3]]) == 1)
                    warn("bad triangle I or II found")
                    node.isVeryBadTriangle = true;
                    net.hasVeryBadTriangle = true
                end
            else
                node.isBadTriangle = true
                setLength!(edge_maj,edge_maj.length+tree_edge2.length)
                setLength!(tree_edge2,0.0)
                tree_edge2.istIdentifiable = false
                push!(net.edges_changed, tree_edge2);
            end
        end
    end #ends the search for bad things
    if(node.k > 3 && !node.isBadDiamondI && !node.isBadDiamondII)
        #println("si entra el ultimo if de k>3 y no bad diamondI,II")
        edgebla,tree_edge_incycle,tree_edge1 = hybridEdges(other_min);
        if(!tree_edge_incycle.istIdentifiable)
            tree_edge_incycle.istIdentifiable = true;
            push!(net.edges_changed,tree_edge_incycle);
        end
        edge_maj.istIdentifiable = isEdgeIdentifiable(edge_maj)
        edge_min.istIdentifiable = isEdgeIdentifiable(edge_min)
    end
    if(allow)
        return true, net.edges_changed
    else
        return !net.hasVeryBadTriangle, net.edges_changed
    end
end

updateGammaz!(net::HybridNetwork, node::Node) = updateGammaz!(net, node, false)

#function to check if edge should be identifiable
#it is not only if followed by leaf, or if a newly converted tree edge
function isEdgeIdentifiable(edge::Edge)
    if(edge.hybrid)
        node = edge.node[edge.isChild1 ? 1 : 2]
        #println("is edge $(edge.number) identifiable, node $(node.number)")
        node.hybrid || error("hybrid edge $(edge.number) pointing at tree node $(node.number)")
        major,minor,tree = hybridEdges(node)
        #println("major $(major.number), minor $(minor.number), tree $(tree.number)")
        if(getOtherNode(tree,node).leaf)
            return false
        else
            return true
        end
    else
        if(all([!edge.node[1].leaf,!edge.node[2].leaf]))
            if(!edge.node[1].hybrid && !edge.node[2].hybrid && !edge.fromBadDiamondI)
                return true
            elseif(edge.node[1].hybrid || edge.node[2].hybrid)
                ind = edge.node[1].hybrid ? 1 : 2
                if(!edge.node[ind].isBadDiamondII && !edge.node[ind].isBadTriangle)
                    return true
                else
                    return false
                end
            end
        else
            return false
        end
    end
end
