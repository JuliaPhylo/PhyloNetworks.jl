# functions for moves: change direction, move origin, move target, nni
# originally functions.jl
# Claudia March 2015

# fixit: in moveOrigin/moveTarget when updating partitions with
# splice, looking for same partition twice (before and after splice),
# maybe there is a more efficient way to handle this

# -------------------------- change direction of minor hybrid edge ---------------------------------

# aux function to transform tree edge into hybrid edge
# input: new edge, hybrid node (needs to be attached to new edge)
#        new gamma
# swtichHyb=true if called in hybridatnode
function makeEdgeHybrid!(edge::Edge,node::Node,gamma::Float64; switchHyb=false::Bool)
    !edge.hybrid || error("edge $(edge.number) already hybrid, cannot make it hybrid")
    node.hybrid || error("to make edge $(edge.number) hybrid, you need to give the hybrid node it is going to point to and node $(node.number) is not hybrid")
    #println("estamos en make edge hybrid en edge $(edge.number) y node $(node.number)")
    #println("vamos a hacer hashybedge true para $(getOtherNode(edge,node).number)")
    #println("$(getOtherNode(edge,node).hasHybEdge) debe ser true")
    size(edge.node,1) == 2 || error("strange edge $(edge.number) has $(size(edge.node,1)) nodes instead of 2")
    if(isEqual(edge.node[1],node))
        edge.isChild1 = true
    elseif(isEqual(edge.node[2],node))
        edge.isChild1 = false
    else
        error("node $(node.number) is not attached to edge $(edge.number)")
    end
    edge.hybrid = true
    getOtherNode(edge,node).hasHybEdge = true
    if switchHyb
        edge.gamma = gamma;
        edge.isMajor = (gamma>=0.5) ? true : false
    else
        setGamma!(edge,gamma,false)
    end
    edge.istIdentifiable = isEdgeIdentifiable(edge)
    edge.containRoot = false
end

# aux function to exchange who is the hybrid node
# input: current hybrid, new hybrid
# returns false if there is no need to updategammaz after
#         true if there is need to updategammaz after
function exchangeHybridNode!(net::HybridNetwork, current::Node,new::Node)
    (current.hybrid && !new.hybrid) || error("either current node $(current.number) is not hybrid: current.hybrid $(current.hybrid) or new node $(new.number) is already hybrid: new.hybrid $(new.hybrid)")
    #println("find cycle for current node $(current.number)")
    nocycle,edgesInCycle,nodesInCycle = identifyInCycle(net,current)
    !nocycle || error("strange here: change direction on hybrid node $(current.number) that does not have a cycle to begin with")
    #println("edges in cycle for node $(current.number) are $([e.number for e in edgesInCycle]) and nodes in cycle $([n.number for n in nodesInCycle])")
    for e in edgesInCycle
        e.inCycle = new.number
    end
    for n in nodesInCycle
        n.inCycle = new.number
    end
    update = true
    new.hybrid = true
    new.k = current.k
    removeHybrid!(net,current)
    pushHybrid!(net,new)
    current.hybrid = false
    current.k = -1
    if(current.isBadDiamondI || current.isBadDiamondII)
        current.isBadDiamondI = false
        current.isBadDiamondII = false
        update = false
    elseif(current.isBadTriangle)
        current.isBadTriangle = false
        update = true
    end
    return update
end


# function to change the direction of a hybrid edge
# input: hybrid node, network, isminor=true: changes dir of minor edge
# warning: it assumes that undoGammaz has been run before
#          and that the new hybridization has been
#          approved by the containRoot criterion
# warning: it needs incycle attributes
# returns flag of whether it needs update gammaz, and alreadyNoRoot = !edgemin2.containRoot
function changeDirection!(node::Node, net::HybridNetwork, isminor::Bool)
    node.hybrid || error("node $(node.number) is not hybrid, so we cannot change the direction of the hybrid edge")
    major,minor,tree = hybridEdges(node);
    othermin = getOtherNode(minor,node);
    othermaj = getOtherNode(major,node);
    othertree = getOtherNode(tree,node);
    if(isminor)
        edgebla,edgemin1,edgemin2 = hybridEdges(othermin);
        nodemin1 = getOtherNode(edgemin1,othermin)
        nodemin2 = getOtherNode(edgemin2,othermin)
        removeEdge!(nodemin1,edgemin1)
        removeEdge!(nodemin2,edgemin2)
        removeEdge!(othermaj,major)
        removeEdge!(othertree,tree)
        removeNode!(nodemin1,edgemin1)
        removeNode!(nodemin2,edgemin2)
        removeNode!(othermaj,major)
        removeNode!(othertree,tree)
        setNode!(edgemin1,othermaj)
        setNode!(major,nodemin1)
        setNode!(edgemin2,othertree)
        setNode!(tree,nodemin2)
        setEdge!(othermaj,edgemin1)
        setEdge!(nodemin1,major)
        setEdge!(othertree,edgemin2)
        setEdge!(nodemin2,tree)
        tmajor = major.length
        ttree = tree.length
        tedmin1 = edgemin1.length
        tedmin2 = edgemin2.length
        setLength!(major,tedmin1)
        setLength!(tree,tedmin2)
        setLength!(edgemin1,tmajor)
        setLength!(edgemin2,ttree)
        # -- update partition
        indexPedgemin = whichPartition(net,edgemin2,node.number)
        indexPtree = whichPartition(net,tree,node.number)
        ind = getIndex(edgemin2,net.partition[indexPedgemin].edges)
        deleteat!(net.partition[indexPedgemin].edges,ind)
        ind = getIndex(tree,net.partition[indexPtree].edges)
        deleteat!(net.partition[indexPtree].edges,ind)
        push!(net.partition[indexPtree].edges,edgemin2)
        push!(net.partition[indexPedgemin].edges,tree)
        # --
        @debug "edgemin2 is $(edgemin2.number) and its containRoot is $(edgemin2.containRoot)"
        alreadyNoRoot = !edgemin2.containRoot
        @debug "edgemin2 is $(edgemin2.number) and its containRoot is $(edgemin2.containRoot), alreadyNoRoot $(alreadyNoRoot)"
    else
        edgebla,edgemaj1,edgemaj2 = hybridEdges(othermaj);
        nodemaj1 = getOtherNode(edgemaj1,othermaj)
        nodemaj2 = getOtherNode(edgemaj2,othermaj)
        removeEdge!(nodemaj1,edgemaj1)
        removeEdge!(nodemaj2,edgemaj2)
        removeEdge!(othermin,minor)
        removeEdge!(othertree,tree)
        removeNode!(nodemaj1,edgemaj1)
        removeNode!(nodemaj2,edgemaj2)
        removeNode!(othermin,minor)
        removeNode!(othertree,tree)
        setNode!(edgemaj1,othermin)
        setNode!(minor,nodemaj1)
        setNode!(edgemaj2,othertree)
        setNode!(tree,nodemaj2)
        setEdge!(othermin,edgemaj1)
        setEdge!(nodemaj1,minor)
        setEdge!(othertree,edgemaj2)
        setEdge!(nodemaj2,tree)
        tminor = minor.length
        ttree = tree.length
        tedmaj1 = edgemaj1.length
        tedmaj2 = edgemaj2.length
        setLength!(minor,tedmaj1)
        setLength!(tree,tedmaj2)
        setLength!(edgemaj1,tminor)
        setLength!(edgemaj2,ttree)
        # -- update partition
        indexPedgemaj = whichPartition(net,edgemaj2,node.number)
        indexPtree = whichPartition(net,tree,node.number)
        ind = getIndex(edgemaj2,net.partition[indexPedgemaj].edges)
        deleteat!(net.partition[indexPedgemaj].edges,ind)
        ind = getIndex(tree,net.partition[indexPtree].edges)
        deleteat!(net.partition[indexPtree].edges,ind)
        push!(net.partition[indexPtree].edges,edgemaj2)
        push!(net.partition[indexPedgemaj].edges,tree)
        # --
        @debug "edgemaj2 is $(edgemaj2.number) and its containRoot is $(edgemaj2.containRoot)"
        alreadyNoRoot = !edgemaj2.containRoot
        @debug "edgemaj2 is $(edgemaj2.number) and its containRoot is $(edgemaj2.containRoot), alreadyNoRoot $(alreadyNoRoot)"
    end
    if(node.isBadDiamondI || node.isBadDiamondII)
        node.isBadDiamondI = false
        node.isBadDiamondII = false
        return false, alreadyNoRoot
    elseif(node.isBadTriangle)
        return false, alreadyNoRoot
    end
    return true, alreadyNoRoot
end

changeDirection!(node::Node, net::HybridNetwork) = changeDirection!(node, net, true)

# function to change the direction of the minor hybrid edge
# and update necessary stepts before and after
# input: hybrid node and network, random=false, changes dir of minor edge always
# returns success
# if success=false, it undoes the change direction
function changeDirectionUpdate!(net::HybridNetwork,node::Node, random::Bool)
    global CHECKNET
    node.hybrid || error("cannot change the direction of minor hybrid edge since node $(node.number) is not hybrid")
    undoGammaz!(node,net)
    edgesRoot = identifyContainRoot(net,node);
    #undoContainRoot!(edgesRoot); now inside updateContainRootChangeDir
    if(random)
        edges = hybridEdges(node)
        edges[1].hybrid || error("hybrid node $(node.number) has as major edge a tree edge $(edges[1].number)")
        edges[1].gamma >= 0.5 || error("major hybrid edge $(edges[1].number) has gamma less than 0.5: $(edges[1].gamma)")
        minor = rand() < edges[1].gamma
    else
        minor = true
    end
    @debug "MOVE: change direction around hybrid node $(node.number) for minor $(minor)------- "
    update, alreadyNoRoot = changeDirection!(node,net,minor)
    if(node.k == 4 && update)
        flag2, edgesgammaz = updateGammaz!(net,node)
    else
        flag2 = true
    end
    if(flag2)
#        flag3,edgesroot = updateContainRoot!(net,node); now in updateContainRootChangeDir
        flag3,edgesroot = updateContainRootChangeDir!(net,node,edgesRoot, alreadyNoRoot)
        if(flag3)
            parameters!(net);
            @debug "MOVE: change direction around hybrid node $(node.number) SUCCESSFUL"
            return true
        else
            @debug "MOVE: change direction around hybrid node $(node.number) FAILED because of containRoot"
            CHECKNET && checkNet(net)
            node.isBadDiamondI || undoGammaz!(node,net)
            alreadyNoRoot || undoContainRoot!(edgesroot) #only undo edgesroot if there were changed: alreadyNoRoot=true => no change
            update,a = changeDirection!(node,net,minor);
            alreadyNoRoot == a || error("change direction alreadyNoRoot should match when done the first time and when undoing")
            if(node.k == 4 && update)
                flag2, edgesgammaz = updateGammaz!(net,node)
            else
                flag2 = true
            end
            flag2 || error("when undoing change direction, we should be able to update gammaz again")
            alreadyNoRoot || undoContainRoot!(edgesRoot,false);
            CHECKNET && checkNet(net)
            return false
        end
    else
        @debug "MOVE: change direction around hybrid node $(node.number) FAILED because of gammaz"
        CHECKNET && checkNet(net)
        node.isBadDiamondI || undoGammaz!(node,net)
        update,a = changeDirection!(node,net,minor);
        if(node.k == 4 && update)
            flag2, edgesgammaz = updateGammaz!(net,node)
        else
            flag2 = true
        end
        flag2 || error("when undoing change direction, we should be able to update gammaz again")
        #undoContainRoot!(edgesRoot,false); #redo containRoot as before NO NEED ANYMORE
        CHECKNET && checkNet(net)
        return false
    end
end

changeDirectionUpdate!(net::HybridNetwork,node::Node) = changeDirectionUpdate!(net,node, false)


# function that will take care of the case when one cycle is above another cycle
# edgesRoot come from identifyContainRoot before changeDirection
# this function should be run after changeDirection which gives alreadyNoRoot
# returns flag2, edgesroot
function updateContainRootChangeDir!(net::HybridNetwork,node::Node,edgesRoot::Vector{Edge}, alreadyNoRoot::Bool)
    node.hybrid || error("cannot update contain root on node $(node.number) because it is not hybrid")
    @debug "updating contain root for hybrid node $(node.number), with alreadyNoRoot $(alreadyNoRoot)"
    if(!alreadyNoRoot) #only update root if new descendats were not forbidden to carry root already
        undoContainRoot!(edgesRoot);
        flag3,edgesroot = updateContainRoot!(net,node);
        @debug "alreadyNoRoot is $(alreadyNoRoot), undoing containRoot for $([e.number for e in edgesRoot]), updating containRoot for $([e.number for e in edgesroot]), with flag3 $(flag3)"
        return flag3,edgesroot
    else
        flag3,edgesroot = updateContainRoot!(net,node);
        return flag3,edgesroot
    end
end


# ------------------------- move origin of hybrid edge ---------------------------------


# function to choose an edge to move origin of hybrid
# but for a nni-like move: only to neighbors
# input: hybridnode for the move
# will choose one of 4 neighbor edges at random and see if
# it is suitable. if none is suitable, stops
# vector a is the vector of possible edges: comes from getNeighborsOrigin/Target
# returns: success (bool), edge, ind
function chooseEdgeOriginTarget!(net::HybridNetwork, neighbor::Vector{Edge}, node::Node)
    length(neighbor) < 5 || error("aux vector a should have only 4 entries: $([n.number for n in neighbor])")
    while(!isempty(neighbor))
        ind = 0
        while(ind == 0 || ind > length(neighbor))
            ind = round(Integer,rand()*length(neighbor));
        end
        #println("ind es $(ind), neighbor edge $(neighbor[ind].number)")
        if(!neighbor[ind].hybrid && (neighbor[ind].inCycle == -1 || neighbor[ind].inCycle == node.number))
            return true, neighbor[ind], ind
        else
            deleteat!(neighbor,ind)
        end
    end
    return false, nothing, nothing
end



# function to move the origin of a hybrid edge
# warning: it does not do any update
# input: necessary nodes/edges, see moveOriginUpdate and
#        ipad notes
# warning: it changes the branch lengths of newedge, tree1, tree2 to match the
#          optimum branch lengths in the corresponding other edge (see ipad notes)
# needed tree1 and tree2 as parameters because we need to do undogammaz before moveOrigin
# returns flag=true if newedge was incycle before, to be able to undo if needed (newedgeincycle)
function moveOrigin!(net::HybridNetwork,node::Node,othermin::Node,tree1::Edge, tree2::Edge,newedge::Edge, undo::Bool, newedgeincycle::Bool)
    node.hybrid || error("cannot move origin of hybridization because node $(node.number) is not hybrid")
    size(newedge.node,1) == 2 || error("strange edge $(newedge.number) that has $(size(newedge.node,1)) nodes instead of 2")
    in(newedge,net.edge) || error("newedge $(newedge.number) not in net.edge")
    @debug "othermin $(othermin.number) with edges $([e.number for e in othermin.edge])"
    @debug "tree1 is $(tree1.number), tree2 is $(tree2.number)"
    if(tree1.inCycle == node.number && tree2.inCycle == -1)
        treej = tree1
        treei = tree2
    elseif(tree2.inCycle == node.number && tree1.inCycle == -1)
        treej = tree2
        treei = tree1
    else
        error("tree1 edge $(tree1.number) and tree2 edge $(tree2.number) should be one in cycle and one not incycle")
    end
    @debug "MOVE: newedge is $(newedge.number)"
    @debug "treei is $(treei.number), treej is $(treej.number) (the one incycle)"
    otheri = getOtherNode(treei,othermin);
    otherj = getOtherNode(treej,othermin);
    @debug "otheri is $(otheri.number), otherj is $(otherj.number)"
    node1 = newedge.node[1]; # not waste of memory, needed step
    node2 = newedge.node[2];
    @debug "node1 is $(node1.number), node2 is $(node2.number)"
    neighbor = false
    from_otheri = false
    from_otherj = false
    if(isEqual(otheri,node1))
        n1 = node1
        n2 = node2
        neighbor = true
        from_otheri = true
    elseif(isEqual(otheri,node2))
        n1 = node2
        n2 = node1
        neighbor = true
        from_otheri = true
    elseif(isEqual(otherj,node1))
        n1 = node1
        n2 = node2
        neighbor = true
        from_otherj = true
    elseif(isEqual(otherj,node2))
        n1 = node2
        n2 = node1
        neighbor = true
        from_otherj = true
    end
    @debug "neighbor $(neighbor), from otheri $(otheri.number) $(from_otheri), from otherj $(otherj.number) $(from_otherj), n1 $(n1.number), n2 $(n2.number)"
    if(neighbor && from_otheri)
        ## println("leaving n1 $(n1.number) as it is")
        ##println("removing n2 $(n2.number) from newedge $(newedge.number) and viceversa")
        removeEdge!(n2,newedge)
        removeNode!(n2,newedge)
        #println("removing otherj $(otherj.number) from treej $(treej.number) and viceversa")
        removeEdge!(otherj,treej)
        removeNode!(otherj,treej)
        #println("setting newedge $(newedge.number) to otherj $(otherj.number) and viceversa")
        setNode!(newedge,otherj)
        setEdge!(otherj,newedge)
        #println("setting treej edge $(treej.number) to n2 $(n2.number) and viceversa")
        setNode!(treej,n2)
        setEdge!(n2,treej)
    elseif(neighbor && from_otherj)
        #println("leaving n1 $(n1.number) as it is")
        #println("removing n2 $(n2.number) from newedge $(newedge.number) and viceversa")
        removeEdge!(n2,newedge)
        removeNode!(n2,newedge)
        #println("removing otheri $(otheri.number) from treei edge $(treei.number) and viceversa")
        removeEdge!(otheri,treei)
        removeNode!(otheri,treei)
        #println("setting newedge $(newedge.number) to otheri $(otheri.number) and viceversa")
        setNode!(newedge,otheri)
        setEdge!(otheri,newedge)
        #println("setting treei edge $(treei.number) to n2 $(n2.number) and viceversa")
        setNode!(treei,n2)
        setEdge!(n2,treei)
    else
        error("move origin to edge $(newedge.number) not neighbor to tree1 edge $(tree1.number) nor tree2 edge $(tree2.number), function not debugged!")
        other1 = getOtherNode(tree1,othermin);
        other2 = getOtherNode(tree2,othermin);
        removeEdge!(other1,tree1)
        removeNode!(other1,tree1)
        removeEdge!(node1,newedge)
        removeNode!(node1,newedge)
        setEdge!(other1,newedge)
        setNode!(newedge,other1)
        setEdge!(node1,tree1)
        setNode!(tree1,node1)
        removeEdge!(other2,tree2)
        removeNode!(other2,tree2)
        removeEdge!(node2,newedge)
        removeNode!(node2,newedge)
        setEdge!(other2,newedge)
        setNode!(newedge,other2)
        setEdge!(node2,tree2)
        setNode!(tree2,node2)
    end
    ti = treei.length
    tj = treej.length
    t = newedge.length
    setLength!(newedge,ti+tj)
    if(approxEq(ti,0.0) && approxEq(tj,0.0))
        setLength!(treei,ti)
        setLength!(treej,tj)
    else
        setLength!(treei,(ti/(ti+tj))*t)
        setLength!(treej,(tj/(ti+tj))*t)
    end
    @debug begin printEdges(net); "printed edges" end
    @debug writeTopologyLevel1(net)
    if !undo
        if from_otheri
            # -- update partition
            indexPtreei = whichPartition(net,treei,node.number)
            ind = getIndex(newedge,net.partition[indexPtreei].edges)
            deleteat!(net.partition[indexPtreei].edges,ind)
            ind = getIndex(treei,net.partition[indexPtreei].edges)
            deleteat!(net.partition[indexPtreei].edges,ind)
            push!(net.partition[indexPtreei].edges,treej)
            edges = hybridEdges(otheri)
            for i in 1:3 #check of 3 edges inside hybridEdges
                if(!isEqual(edges[i],treei) && !isEqual(edges[i],newedge))
                    descendants = [edges[i]]
                    cycleNum = [node.inCycle]
                    getDescendants!(getOtherNode(edges[i],otheri),edges[i],descendants,cycleNum)
                    !isempty(descendants) || error("descendants is empty for node $(otheri.number)")
                    @debug "for node $(otheri.number), descendants are $([e.number for e in descendants]), and cycleNum is $(cycleNum)"
                    partition = Partition(unique(cycleNum),descendants) # create new partition
                    push!(net.partition, partition)
                    for e in descendants #delete edges from partition with tree originally
                        ind = getIndex(e,net.partition[indexPtreei].edges)
                        deleteat!(net.partition[indexPtreei].edges,ind)
                    end
                    net.partition[indexPtreei].cycle = union([node.inCycle],symdiff(net.partition[indexPtreei].cycle,cycleNum))
                    break
                end
            end
            # -- update in cycle
            newedge.inCycle = node.number
            switchCycleTree!(tree1,tree2,node)
            node.k += 1
            otheri.inCycle = node.number
        elseif(from_otherj)
            if(newedge.inCycle == node.number)
                # -- update partition
                indexPtreei = whichPartition(net,treei,node.number)
                push!(net.partition[indexPtreei].edges,treej)
                push!(net.partition[indexPtreei].edges,newedge)
                ind = getIndex(treei,net.partition[indexPtreei].edges)
                deleteat!(net.partition[indexPtreei].edges,ind)
                edges = hybridEdges(otherj)
                for i in 1:3 #check of 3 edges inside hybridEdges
                    if(!isEqual(edges[i],treej) && !isEqual(edges[i],newedge))
                        indexP = whichPartition(net,edges[i],node.number)
                        for e in net.partition[indexP].edges
                            push!(net.partition[indexPtreei].edges,e)
                        end
                        net.partition[indexPtreei].cycle = union(net.partition[indexPtreei].cycle,net.partition[indexP].cycle)
                        part = splice!(net.partition,indexP) #splice at the end to not change indexPtreei
                        break
                    end
                end
                # -- update in cycle
                newedge.inCycle = -1
                switchCycleTree!(tree1,tree2,node)
                otherj.inCycle = -1
                node.k -= 1
                return true
            else
                # -- update partition
                indexPnew = whichPartition(net,newedge,node.number)
                indexPtreei = whichPartition(net,treei,node.number)
                ind = getIndex(newedge,net.partition[indexPnew].edges)
                deleteat!(net.partition[indexPnew].edges,ind)
                ind = getIndex(treei,net.partition[indexPtreei].edges)
                deleteat!(net.partition[indexPtreei].edges,ind)
                push!(net.partition[indexPnew].edges,treei)
                push!(net.partition[indexPtreei].edges,newedge)
            end
        end
    else #yes undo
        if(newedge.inCycle == node.number)
            # -- update partition
            indexPtreei = whichPartition(net,treei,node.number)
            push!(net.partition[indexPtreei].edges,treej)
            push!(net.partition[indexPtreei].edges,newedge)
            ind = getIndex(treei,net.partition[indexPtreei].edges)
            deleteat!(net.partition[indexPtreei].edges,ind)
            edges = hybridEdges(otherj)
            for i in 1:3 #check of 3 edges inside hybridEdges
                if(!isEqual(edges[i],treej) && !isEqual(edges[i],newedge))
                    indexP = whichPartition(net,edges[i],node.number)
                    for e in net.partition[indexP].edges
                        push!(net.partition[indexPtreei].edges,e)
                    end
                    net.partition[indexPtreei].cycle = union(net.partition[indexPtreei].cycle,net.partition[indexP].cycle)
                    part = splice!(net.partition,indexP) #splice at the end
                    break
                end
            end
            # -- update in cycle
            newedge.inCycle = -1
            switchCycleTree!(tree1,tree2,node)
            node.k -= 1
            otherj.inCycle = -1
        elseif(newedge.inCycle == -1)
            if(newedgeincycle)
                # -- update partition
                indexPtreei = whichPartition(net,treei,node.number)
                ind = getIndex(newedge,net.partition[indexPtreei].edges)
                deleteat!(net.partition[indexPtreei].edges,ind)
                ind = getIndex(treei,net.partition[indexPtreei].edges)
                deleteat!(net.partition[indexPtreei].edges,ind)
                push!(net.partition[indexPtreei].edges,treej)
                edges = hybridEdges(otheri)
                for i in 1:3 #check of 3 edges inside hybridEdges
                    if(!isEqual(edges[i],treei) && !isEqual(edges[i],newedge))
                        descendants = [edges[i]]
                        cycleNum = [node.inCycle]
                        getDescendants!(getOtherNode(edges[i],otheri),edges[i],descendants,cycleNum)
                        !isempty(descendants) || error("descendants is empty for node $(otheri.number)")
                        @debug "for node $(otheri.number), descendants are $([e.number for e in descendants]), and cycleNum is $(cycleNum)"
                        partition = Partition(unique(cycleNum),descendants) # create new partition
                        push!(net.partition, partition)
                        for e in descendants #delete edges from partition with tree originally
                            ind = getIndex(e,net.partition[indexPtreei].edges)
                            deleteat!(net.partition[indexPtreei].edges,ind)
                        end
                        net.partition[indexPtreei].cycle = union([node.inCycle],symdiff(net.partition[indexPtreei].cycle,cycleNum))
                        break
                    end
                end
                # -- update in cycle
                newedge.inCycle = node.number
                switchCycleTree!(tree1,tree2,node)
                otheri.inCycle = node.number
                node.k += 1
            else
                # -- update partition
                indexPnew = whichPartition(net,newedge,node.number)
                indexPtreei = whichPartition(net,treei,node.number)
                ind = getIndex(newedge,net.partition[indexPnew].edges)
                deleteat!(net.partition[indexPnew].edges,ind)
                ind = getIndex(treei,net.partition[indexPtreei].edges)
                deleteat!(net.partition[indexPtreei].edges,ind)
                push!(net.partition[indexPnew].edges,treei)
                push!(net.partition[indexPtreei].edges,newedge)
            end
        end
    end
    return false
end

moveOrigin!(net::HybridNetwork,node::Node,othermin::Node,tree1::Edge, tree2::Edge,newedge::Edge) = moveOrigin!(net,node,othermin,tree1, tree2,newedge, false,false)

# function to switch the cycle of 2 tree edges in moveTarget
function switchCycleTree!(tree1::Edge, tree2::Edge, node::Node)
    (!tree1.hybrid && !tree2.hybrid) || error("tree1 edge $(tree1.number) and tree2 $(tree2.number) cannot be hybrid to switch incycle attribute")
    (tree1.inCycle == node.number && tree2.inCycle == -1) || (tree2.inCycle == node.number && tree1.inCycle == -1) || error("tree1 edge or tree2 edge must by in cycle by node $(node.number) and the other -1, but: $([tree1.inCycle, tree2.inCycle])")
    k = tree1.inCycle
    m = tree2.inCycle
    tree2.inCycle = k
    tree1.inCycle = m
end

# function to choose minor/major hybrid edge
# to move the origin or target
# random=false always chooses the minor edge
# returns the othermin node and majoredge for target
function chooseMinorMajor(node::Node, random::Bool, target::Bool)
    node.hybrid || error("node $(node.number) is not hybrid, so we cannot delete hybridization event around it")
    major,minor,tree = hybridEdges(node);
    #println("hybrid node $(node.number) has major $(major.number), minor $(minor.number) and tree $(tree.number)")
    if random
        r = rand()
        major.gamma >= 0.5 || error("strange major hybrid edge $(major.number) with gamma $(major.gamma) less than 0.5")
        @debug (major.gamma != 1.0 ? "" :
            "strange major hybrid edge $(major.number) with gamma $(major.gamma) equal to 1.0")
        othermin = r < major.gamma ? getOtherNode(minor,node) : getOtherNode(major,node)
        majoredge = r < major.gamma ? major : minor
        @debug (r < major.gamma ? "MOVE: will do move on minor hybrid edge" :
                                  "MOVE: will do move on major hybrid edge")
    else
        @debug "MOVE: will do move on minor hybrid edge"
        majoredge = major
        othermin = getOtherNode(minor,node);
    end
    if target
        return othermin, majoredge
    else
        return othermin
    end
end

chooseMinorMajor(node::Node, random::Bool) = chooseMinorMajor(node, random, false)

# function to give the vector of edges neighbor
# that includes all the suitable neighbors of othermin
# to move the origin
function getNeighborsOrigin(othernode::Node)
    othernode.hasHybEdge || error("other node $(othernode.number) should the attach to the hybrid edge whose origin will be moved")
    !othernode.hybrid || error("other node $(othernode.number) must not be the hybrid node to move origin")
    edges = hybridEdges(othernode)
    length(edges) == 3 || error("length of node.edge for node $(othernode.number) should be 3, not $(length(edges))")
    edges[1].hybrid || error("edge $(edges[1].number) should be hybrid because it is the one to move the origin from")
    neighbor = Edge[]
    n1 = getOtherNode(edges[2],othernode)
    n2 = getOtherNode(edges[3],othernode)
    if(!n1.leaf)
        e = hybridEdges(n1,edges[2])
        length(e) == 2 || error("node $(n1.number) is not a leaf but has $(length(e)+1) edges")
        push!(neighbor, e[1])
        push!(neighbor, e[2])
    end
    if(!n2.leaf)
        e = hybridEdges(n2,edges[3])
        length(e) == 2 || error("node $(n1.number) is not a leaf but has $(length(e)+1) edges")
        push!(neighbor, e[1])
        push!(neighbor, e[2])
    end
    return neighbor
end

# function to move the origin of a hybrid edge
# and update everything that needs update: gammaz
# input: network, hybrid node, othermin already chosen with chooseMinorMajor
# returns: success (bool), flag, nocycle, flag2
function moveOriginUpdate!(net::HybridNetwork, node::Node, othermin::Node, newedge::Edge)
    global CHECKNET
    node.hybrid || error("node $(node.number) is not hybrid, so we cannot delete hybridization event around it")
    @debug "MOVE: move Origin for hybrid node $(node.number)"
    in(newedge,net.edge) || error("newedge $(newedge.number) not in net.edge")
    edgebla, tree1, tree2 = hybridEdges(othermin);
    undoGammaz!(node,net);
    #println("othermin is $(othermin.number)")
    newedgeincycle = moveOrigin!(net,node,othermin,tree1,tree2,newedge)
    flag2, edgesGammaz = updateGammaz!(net,node)
    if(flag2)
        parameters!(net);
        @debug "MOVE: move Origin for hybrid node $(node.number) SUCCESSFUL"
        return true,flag2
    else
        @debug "MOVE: move Origin for hybrid node $(node.number) FAILED"
        CHECKNET && checkNet(net)
        isempty(edgesGammaz) || undoistIdentifiable!(edgesGammaz)
        undoGammaz!(node,net);
        @debug "MOVE: undoing move origin for conflict: gammaz"
        moveOrigin!(net,node,othermin,tree1,tree2,newedge,true,newedgeincycle);
        flag2, edgesGammaz = updateGammaz!(net,node)
        (flag2 || node.isVeryBadTriangle || node.isExtBadTriangle) || error("updating gammaz for undone moveOrigin, should not be any problem")
        CHECKNET && checkNet(net)
        return false, flag2
    end
end

# function to repeat moveOriginUpdate with all the neighbors
# until success or failure of all for a given hybrid node
# returns success (bool): failure from: neighbors not suitable to begin
# with, or conflicts incycle
function moveOriginUpdateRepeat!(net::HybridNetwork, node::Node, random::Bool)
    node.hybrid || error("cannot move origin because node $(node.number) is not hybrid")
    othermin = chooseMinorMajor(node,random)
    #println("othermin is $(othermin.number) with edges $([e.number for e in othermin.edge])")
    neighbor = getNeighborsOrigin(othermin)
    #println("neighbors list is $([n.number for n in neighbor])")
    success = false
    while(!isempty(neighbor) && !success)
        success1,newedge,ind = chooseEdgeOriginTarget!(net, neighbor,node)
        !isa(newedge,Nothing) || return false
        success1 || return false
        #println("newedge is $(newedge.number), success1 is $(success1)")
        in(newedge,net.edge) || error("newedge $(newedge.number) is not in net.edge")
        success,flag2 = moveOriginUpdate!(net, node, othermin, newedge)
        #println("after update, success is $(success)")
        if(!success)
            @debug "move origin failed, will delete that neighbor and try new one"
            deleteat!(neighbor,ind)
        end
        #println("neighbor list is $([n.number for n in neighbor])")
    end
    success || return false
    return true
end

# ------------------------------ move target of hybridization -------------------------------


# function to give the vector of edges neighbor
# that includes all the suitable neighbors of othermin
# to move the target
function getNeighborsTarget(node::Node,majoredge::Edge)
    node.hybrid || error("node $(node.number) must be the hybrid node to move target")
    major,minor,tree = hybridEdges(node)
    othermajor = getOtherNode(majoredge,node)
    #println("othermajor is $(othermajor.number)")
    neighbor = Edge[]
    n1 = getOtherNode(majoredge,node)
    n2 = getOtherNode(tree,node)
    if(!n1.leaf)
        e = hybridEdges(n1,majoredge)
        length(e) == 2 || error("node $(n1.number) is not a leaf but has $(length(e)+1) edges")
        push!(neighbor, e[1])
        push!(neighbor, e[2])
    end
    if(!n2.leaf)
        e = hybridEdges(n2,tree)
        length(e) == 2 || error("node $(n1.number) is not a leaf but has $(length(e)+1) edges")
        push!(neighbor, e[1])
        push!(neighbor, e[2])
    end
    return neighbor
end

# function to switch major and tree node in moveTarget
function switchMajorTree!(major::Edge, tree::Edge, node::Node)
    !tree.hybrid || error("tree edge $(tree.number) cannot be hybrid to switch to major")
    major.hybrid || error("major edge $(major.number) has to be hybrid to switch to tree")
    @debug "switch major $(major.number) tree $(tree.number), node $(node.number)"
    g = major.gamma #needed because changed inside makeEdgeTree
    cycle = major.inCycle #needed because changed inside makeEdgeTree
    makeEdgeTree!(major,node)
    makeEdgeHybrid!(tree,node,g)
    tree.inCycle = cycle
    major.inCycle = -1
end

# function to move the target of a hybrid edge
# warning: it does not do any update
# input: necessary nodes/edges, see moveOriginUpdate and
#        ipad notes, undo=true if it is undoing a previous moveTarget
# warning: it changes the branch lengths of newedge, tree1, tree2 to match the
#          optimum branch lengths in the corresponding other edge (see ipad notes)
# returns flag=true if newedge was incycle before, to be able to undo if needed (newedgeincycle)
# also returns alreadyNoRoot
function moveTarget!(net::HybridNetwork,node::Node, major::Edge, tree::Edge, newedge::Edge, undo::Bool, newedgeincycle::Bool)
    node.hybrid || error("cannot move origin of hybridization because node $(node.number) is not hybrid")
    length(newedge.node) == 2 || error("strange edge $(newedge.number) that has $(size(newedge.node,1)) nodes instead of 2")
    newedge.inCycle == node.number || newedge.inCycle == -1 || error("newedge in cycle must be -1 or $(node.number), not $(newedge.inCycle)")
    in(newedge,net.edge) || error("newedge $(newedge.number) not in net.edge")
    othermajor = getOtherNode(major,node);
    treenode = getOtherNode(tree,node);
    node1 = newedge.node[1]; # not waste of memory, needed step
    node2 = newedge.node[2];
    @debug "MOVE: hybrid node is $(node.number), major edge $(major.number), tree edge $(tree.number)"
    @debug "MOVE: newedge is $(newedge.number), othermajor $(othermajor.number), treenode $(treenode.number), node1 $(node1.number), node2 $(node2.number)"
    neighbor = false
    from_othermajor = false
    from_treenode = false
    if(isEqual(othermajor,node1))
        n1 = node1
        n2 = node2
        neighbor = true
        from_othermajor = true
    elseif(isEqual(othermajor,node2))
        n1 = node2
        n2 = node1
        neighbor = true
        from_othermajor = true
    elseif(isEqual(treenode,node1))
        n1 = node1
        n2 = node2
        neighbor = true
        from_treenode = true
    elseif(isEqual(treenode,node2))
        n1 = node2
        n2 = node1
        neighbor = true
        from_treenode = true
    end
    # -- alreadyNoRoot?
    if(newedge.inCycle != -1)
        alreadyNoRoot = false
    else
        if(neighbor && from_othermajor)
            alreadyNoRoot = !newedge.containRoot
        elseif(neighbor && from_treenode)
            for e in node.edge
                if(!isEqual(e,major) && !isEqual(e,tree)) #looking for minor edge
                    o = getOtherNode(e,node)
                    for e2 in o.edge
                        if(!e2.hybrid && e2.inCycle != -1)
                            alreadyNoRoot = !e2.containRoot
                            break
                        end
                    end
                    break
                end
            end
        elseif(!neighbor)
            alreadyNoRoot = false
        end
    end
    @debug "neighbor $(neighbor), from othermajor $(from_othermajor), from treenode $(from_treenode), n1 $(n1.number), n2 $(n2.number), alreadyNoRoot $(alreadyNoRoot)"
    if(neighbor && from_othermajor)
        #println("leaving n1 $(n1.number) as it is")
        #println("removing n2 $(n2.number) from newedge $(newedge.number) and viceversa")
        removeEdge!(n2,newedge)
        removeNode!(n2,newedge)
        #println("removing treenode $(treenode.number) from tree edge $(tree.number) and viceversa")
        removeEdge!(treenode,tree)
        removeNode!(treenode,tree)
        #println("setting newedge $(newedge.number) to treenode $(treenode.number) and viceversa")
        setNode!(newedge,treenode)
        setEdge!(treenode,newedge)
        #println("setting tree edge $(tree.number) to n2 $(n2.number) and viceversa")
        setNode!(tree,n2)
        setEdge!(n2,tree)
    elseif(neighbor && from_treenode)
        #println("leaving n1 $(n1.number) as it is")
        #println("removing n2 $(n2.number) from newedge $(newedge.number) and viceversa")
        removeEdge!(n2,newedge)
        removeNode!(n2,newedge)
        #println("removing othermajor $(othermajor.number) from major edge $(major.number) and viceversa")
        removeEdge!(othermajor,major)
        removeNode!(othermajor,major)
        #println("setting newedge $(newedge.number) to othermajor $(othermajor.number) and viceversa")
        setNode!(newedge,othermajor)
        setEdge!(othermajor,newedge)
        #println("setting major edge $(major.number) to n2 $(n2.number) and viceversa")
        setNode!(major,n2)
        setEdge!(n2,major)
    else
        error("move target to edge $(newedge.number) not neighbor to major edge $(major.number) nor tree edge $(tree.number), function not debugged!")
        #println("removing major edge $(major.number) from othermajor node $(othermajor.number) and viceverse")
        removeEdge!(othermajor,major)
        removeNode!(othermajor,major)
        #println("removing treenode $(treenode.number) from tree edge $(tree.number) and viceversa")
        removeEdge!(treenode,tree)
        removeNode!(treenode,tree)
        #println("removing node1 $(node1.number) from newedge $(newedge.number) and viceversa")
        removeEdge!(node1,newedge)
        removeNode!(node1,newedge)
        #println("removing node2 $(node2.number) from newedge $(newedge.number) and viceversa")
        removeEdge!(node2,newedge)
        removeNode!(node2,newedge)
        #println("setting newedge $(newedge.number) to othermajor $(othermajor.number) and viceversa")
        setNode!(newedge,othermajor)
        setEdge!(othermajor,newedge)
        #println("setting newedge $(newedge.number) to treenode $(treenode.number) and viceversa")
        setNode!(newedge,treenode)
        setEdge!(treenode,newedge)
        #println("setting major edge $(major.number) to node1 $(node1.number) and viceversa")
        setNode!(major,node1)
        setEdge!(node1,major)
        #println("setting tree edge $(tree.number) to node2 $(node2.number) and viceversa")
        setNode!(tree,node2)
        setEdge!(node2,tree)
    end
    t1 = major.length
    t2 = tree.length
    t = newedge.length
    setLength!(newedge,t1+t2)
    if(approxEq(t1,0.0) && approxEq(t2,0.0))
        setLength!(major,t1)
        setLength!(tree,t2)
    else
        setLength!(major,t1/(t1+t2)*t)
        setLength!(tree,t2/(t1+t2)*t)
    end
    if(!undo)
        if(from_treenode)
            # -- update partition
            indexPtree = whichPartition(net,tree,node.number)
            ind = getIndex(newedge,net.partition[indexPtree].edges)
            deleteat!(net.partition[indexPtree].edges,ind)
            ind = getIndex(tree,net.partition[indexPtree].edges)
            deleteat!(net.partition[indexPtree].edges,ind)
            push!(net.partition[indexPtree].edges,major)
            edges = hybridEdges(treenode)
            for i in 1:3 #check of 3 edges inside hybridEdges
                if(!isEqual(edges[i],tree) && !isEqual(edges[i],newedge))
                    descendants = [edges[i]]
                    cycleNum = [node.inCycle]
                    getDescendants!(getOtherNode(edges[i],treenode),edges[i],descendants,cycleNum)
                    !isempty(descendants) || error("descendants is empty for node $(treenode.number)")
                    @debug "for node $(treenode.number), descendants are $([e.number for e in descendants]), and cycleNum is $(cycleNum)"
                    partition = Partition(unique(cycleNum),descendants) # create new partition
                    push!(net.partition, partition)
                    for e in descendants #delete edges from partition with tree originally
                        ind = getIndex(e,net.partition[indexPtree].edges)
                        deleteat!(net.partition[indexPtree].edges,ind)
                    end
                    net.partition[indexPtree].cycle = union([node.inCycle],symdiff(net.partition[indexPtree].cycle,cycleNum))
                    break
                end
            end
            # -- update in cycle
            @debug "from treenode treatment, switch major $(major.number) to tree $(tree.number)"
            switchMajorTree!(major,tree,node)
            node.k += 1
            newedge.inCycle = node.number
            treenode.inCycle = node.number
        elseif(from_othermajor)
            if(newedge.inCycle == node.number)
                @debug "from othermajor and newedge incycle treatment, switch major $(major.number) to tree $(tree.number)"
                # -- update partition
                indexPtree = whichPartition(net,tree,node.number)
                push!(net.partition[indexPtree].edges,major)
                push!(net.partition[indexPtree].edges,newedge)
                ind = getIndex(tree,net.partition[indexPtree].edges)
                deleteat!(net.partition[indexPtree].edges,ind)
                edges = hybridEdges(othermajor)
                for i in 1:3 #check of 3 edges inside hybridEdges
                    if(!isEqual(edges[i],major) && !isEqual(edges[i],newedge))
                        indexP = whichPartition(net,edges[i],node.number)
                        for e in net.partition[indexP].edges
                            push!(net.partition[indexPtree].edges,e)
                        end
                        net.partition[indexPtree].cycle = union(net.partition[indexPtree].cycle,net.partition[indexP].cycle)
                        part = splice!(net.partition,indexP) #splice at the end
                        break
                    end
                end
                # -- update inCycle
                switchMajorTree!(major,tree,node)
                node.k -= 1
                newedge.inCycle = -1
                othermajor.inCycle = -1
                return true, alreadyNoRoot
            else
                # -- update partition
                indexPnew = whichPartition(net,newedge,node.number)
                indexPtree = whichPartition(net,tree,node.number)
                ind = getIndex(newedge,net.partition[indexPnew].edges)
                deleteat!(net.partition[indexPnew].edges,ind)
                ind = getIndex(tree,net.partition[indexPtree].edges)
                deleteat!(net.partition[indexPtree].edges,ind)
                push!(net.partition[indexPnew].edges,tree)
                push!(net.partition[indexPtree].edges,newedge)
                # -- update in cycle
                major.istIdentifiable = isEdgeIdentifiable(major)
            end
        end
    else
        if(from_treenode)
            # -- update partition
            indexPmajor = whichPartition(net,major,node.number)
            push!(net.partition[indexPmajor].edges,tree)
            push!(net.partition[indexPmajor].edges,newedge)
            ind = getIndex(major,net.partition[indexPmajor].edges)
            deleteat!(net.partition[indexPmajor].edges,ind)
            edges = hybridEdges(treenode)
            for i in 1:3 #check of 3 edges inside hybridEdges
                if(!isEqual(edges[i],tree) && !isEqual(edges[i],newedge))
                    indexP = whichPartition(net,edges[i],node.number)
                    for e in net.partition[indexP].edges
                        push!(net.partition[indexPmajor].edges,e)
                    end
                    net.partition[indexPmajor].cycle = union(net.partition[indexPmajor].cycle,net.partition[indexP].cycle)
                    part = splice!(net.partition,indexP) #splice at the end
                    break
                end
            end
            # -- update in cycle
            switchMajorTree!(tree,major,node)
            node.k -= 1
            newedge.inCycle = -1
            treenode.inCycle = -1
        elseif(from_othermajor)
            if(newedgeincycle)
                # -- update partition
                indexPmajor = whichPartition(net,major,node.number)
                ind = getIndex(major,net.partition[indexPmajor].edges)
                deleteat!(net.partition[indexPmajor].edges,ind) #delete major
                ind = getIndex(newedge,net.partition[indexPmajor].edges)
                deleteat!(net.partition[indexPmajor].edges,ind) # delete newedge
                push!(net.partition[indexPmajor].edges,tree) # add tree
                edges = hybridEdges(othermajor)
                for i in 1:3 #check of 3 edges inside hybridEdges
                    if(!isEqual(edges[i],major) && !isEqual(edges[i],newedge))
                        descendants = [edges[i]]
                        cycleNum = [node.inCycle]
                        getDescendants!(getOtherNode(edges[i],othermajor),edges[i],descendants,cycleNum)
                        !isempty(descendants) || error("descendants is empty for node $(othermajor.number)")
                        @debug "for node $(othermajor.number), descendants are $([e.number for e in descendants]), and cycleNum is $(cycleNum)"
                        partition = Partition(unique(cycleNum),descendants) # create new partition
                        push!(net.partition, partition)
                        for e in descendants #delete edges from partition with major originally
                            ind = getIndex(e,net.partition[indexPmajor].edges)
                            deleteat!(net.partition[indexPmajor].edges,ind)
                        end
                        net.partition[indexPmajor].cycle = union([node.inCycle],symdiff(net.partition[indexPmajor].cycle,cycleNum))
                        break
                    end
                end
                # -- update inCycle
                switchMajorTree!(tree,major,node)
                node.k += 1
                newedge.inCycle = node.number
                othermajor.inCycle = node.number
            else
                # -- update partition
                indexPnew = whichPartition(net,newedge,node.number)
                indexPtree = whichPartition(net,tree,node.number)
                ind = getIndex(newedge,net.partition[indexPnew].edges)
                deleteat!(net.partition[indexPnew].edges,ind)
                ind = getIndex(tree,net.partition[indexPtree].edges)
                deleteat!(net.partition[indexPtree].edges,ind)
                push!(net.partition[indexPnew].edges,tree)
                push!(net.partition[indexPtree].edges,newedge)
                # -- update in cycle
                major.istIdentifiable = isEdgeIdentifiable(major)
            end
        end
    end
    return false, alreadyNoRoot
end

moveTarget!(net::HybridNetwork,node::Node, major::Edge, tree::Edge, newedge::Edge) = moveTarget!(net,node, major, tree, newedge, false, false)

# function to move the target of a hybrid edge
# and update everything that needs update: gammaz, root
# input: network, hybrid node, othermin, majoredge (chosen with chooseMinorMajor), newedge
# returns: success (bool), flag2, flag3
function moveTargetUpdate!(net::HybridNetwork, node::Node, othermin::Node, majoredge::Edge, newedge::Edge)
    global CHECKNET
    node.hybrid || error("node $(node.number) is not hybrid, so we cannot delete hybridization event around it")
    @debug "MOVE: move Target for hybrid node $(node.number)"
    in(newedge,net.edge) || error("newedge $(newedge.number) not in net.edge")
    major,minor,tree = hybridEdges(node)
    undoGammaz!(node,net);
    edgesRoot = identifyContainRoot(net,node);
    edgesRoot = setdiff(edgesRoot,[tree]) #need to remove tree from edgesRoot or it will be undone/updated with containRoot
    #undoContainRoot!(edgesRoot); now in updateContainRootChangeDir
    #@debug "undoContainRoot for edges $([e.number for e in edgesRoot])"
    newedgeincycle, alreadyNoRoot = moveTarget!(net,node,majoredge,tree,newedge)
    flag2, edgesGammaz = updateGammaz!(net,node)
    if(flag2)
#        flag3,edgesroot = updateContainRoot!(net,node) now in updateContainRootChangeDir
        flag3,edgesroot = updateContainRootChangeDir!(net,node,edgesRoot, alreadyNoRoot)
        if(flag3)
            parameters!(net);
            @debug "MOVE: move Target for hybrid node $(node.number) SUCCESSFUL"
            return true,flag2, flag3
        else
            @debug "MOVE: move Target for hybrid node $(node.number) FAILED"
            @debug begin printEverything(net); "printed everything" end
            alreadyNoRoot || undoContainRoot!(edgesroot) #only undo edgesroot if there were changed: alreadyNoRoot=true => no change
            alreadyNoRoot || undoContainRoot!(edgesRoot,false);
            isempty(edgesGammaz) || undoistIdentifiable!(edgesGammaz)
            undoGammaz!(node,net);
            @debug "MOVE: undoing move target for conflict: updategammaz"
            moveTarget!(net,node,majoredge,tree,newedge,true,newedgeincycle);
            flag2, edgesGammaz = updateGammaz!(net,node)
            #@debug "undoContainRoot for edges $([e.number for e in edgesRoot])"
            (flag2 || node.isVeryBadTriangle || node.isExtBadTriangle) || error("updating gammaz/root for undone moveTarget, should not be any problem, but flag2 $(flag2) and node not very/ext bad triangle")
            CHECKNET && checkNet(net)
            return false, flag2,flag3
        end
    else
        @debug "MOVE: move Target for hybrid node $(node.number) FAILED"
        @debug begin printEverything(net); "printed everything" end
        if CHECKNET
            flag3,edgesroot = updateContainRootChangeDir!(net,node,edgesRoot, alreadyNoRoot) #only to be sure there are no errors in the modified net
            checkNet(net)
            alreadyNoRoot || undoContainRoot!(edgesroot) #only undo edgesroot if there were changed: alreadyNoRoot=true => no change
            alreadyNoRoot || undoContainRoot!(edgesRoot,false);
        end
        isempty(edgesGammaz) || undoistIdentifiable!(edgesGammaz)
        undoGammaz!(node,net);
        @debug "MOVE: undoing move target for conflict: updategammaz"
        moveTarget!(net,node,majoredge,tree,newedge,true,newedgeincycle);
        flag2, edgesGammaz = updateGammaz!(net,node)
        #@debug "undoContainRoot for edges $([e.number for e in edgesRoot])"
        (flag2 || node.isVeryBadTriangle || node.isExtBadTriangle) || error("updating gammaz/root for undone moveTarget, should not be any problem, but flag2 $(flag2) and node not very/ext bad triangle")
        CHECKNET && checkNet(net)
        return false, flag2, false
    end
end


# function to repeat moveTargetUpdate with all the neighbors
# until success or failure of all for a given hybrid node
# returns success (bool): failure from: neighbors not suitable to begin
# with, or conflicts incycle
function moveTargetUpdateRepeat!(net::HybridNetwork, node::Node, random::Bool)
    node.hybrid || error("cannot move origin because node $(node.number) is not hybrid")
    othermin,majoredge = chooseMinorMajor(node,random, true)
    #println("othermin is $(othermin.number), will move edge: majoredge $(majoredge.number)")
    neighbor = getNeighborsTarget(node,majoredge)
    #println("neighbors list is $([n.number for n in neighbor])")
    success = false
    while(!isempty(neighbor) && !success)
        success1,newedge,ind = chooseEdgeOriginTarget!(net, neighbor,node);
        success1 || return false
        #println("newedge is $(newedge.number), success1 is $(success1)")
        in(newedge,net.edge) || error("newedge $(newedge.number) not in net.edge")
        success,flag2 = moveTargetUpdate!(net, node, othermin, majoredge,newedge)
        #println("after update, success is $(success)")
        if(!success)
            @debug "move target failed, will delete neighbor and try new one"
            deleteat!(neighbor,ind)
        end
        #println("neighbor list is $([n.number for n in neighbor])")
    end
    success || return false
    return true
end


# ------------------------------------- tree move NNI -------------------------------------------

# function to check if an edge is internal edge
function isInternalEdge(edge::Edge)
    length(edge.node) == 2 || error("edge $(edge.number) has $(length(edge.node)) nodes, should be 2")
    return !edge.node[1].leaf && !edge.node[2].leaf
end

# function to check if the two node of an edge
# do not have hybrid edges
function hasNeighborHybrid(edge::Edge)
    length(edge.node) == 2 || error("edge $(edge.number) has $(length(edge.node)) nodes, should be 2")
    return edge.node[1].hasHybEdge || edge.node[2].hasHybEdge
end

# function to choose a tree edge for NNI
# its two nodes must have hasHybEdge=false
function chooseEdgeNNI(net::Network,N::Integer)
    N > 0 || error("N must be positive: $(N)")
    index1 = 0
    i = 0
    while((index1 == 0 || index1 > size(net.edge,1) || net.edge[index1].hybrid || hasNeighborHybrid(net.edge[index1]) || !isInternalEdge(net.edge[index1])) && i < N)
        index1 = round(Integer,rand()*size(net.edge,1));
        i += 1
    end
    if(i < N)
        return true,net.edge[index1]
    else
        @debug "cannot find suitable tree edge for NNI after $(N) attempts"
        return false,nothing
    end
end


# tree move NNI
# it does not consider keeping the same setup
# as a possibility
# input: edge from chooseEdgeNNI
# returns success; failure if cannot do nni without creating intersecting cycles
function NNI!(net::Network,edge::Edge)
    !edge.hybrid || error("cannot do NNI on hybrid edge: $(edge.number)")
    @debug "MOVE: NNI on tree edge $(edge.number)"
    edges1 = hybridEdges(edge.node[1],edge)
    edges2 = hybridEdges(edge.node[2],edge)
     (!edges1[1].hybrid && !edges1[2].hybrid) || error("cannot do tree move NNI if one of the edges is hybrid: check neighbors of edge $(edge.number)")
     (!edges2[1].hybrid && !edges2[2].hybrid) || error("cannot do tree move NNI if one of the edges is hybrid: check neighbors of edge $(edge.number)")
    if(rand() < 0.5)
        e1 = edges1[1]
        e2 = edges1[2]
    else
        e1 = edges1[2]
        e2 = edges1[1]
    end
     if(rand() < 0.5)
        e3 = edges2[1]
        e4 = edges2[2]
    else
        e3 = edges2[2]
        e4 = edges2[1]
    end
    t1 = e1.length
    t = edge.length
    t4 = e4.length
    n1 = edge.node[1]
    n2 = edge.node[2]
    nohybrid = false
    @debug "e1 is $(e1.number), e2 is $(e2.number), e3 is $(e3.number), e4 is $(e4.number), n1 is $(n1.number), n2 is $(n2.number)"
    if(edge.inCycle != -1)
        node = net.node[getIndexNode(edge.inCycle,net)]
        node.hybrid || error("edge $(edge.number) has incycle $(edge.inCycle) but node $(node.number) is not hybrid")
        (n1.inCycle == edge.inCycle && n2.inCycle == edge.inCycle) || error("edge $(edge.number) is in cycle $(edge.inCycle) but its nodes are not: $(n1.number), $(n2.number)")
        if((e2.inCycle == e3.inCycle == edge.inCycle && e1.inCycle == e4.inCycle == -1) || (e1.inCycle == e4.inCycle == edge.inCycle && e2.inCycle == e3.inCycle == -1))
            nothing
        elseif((e2.inCycle == e4.inCycle == edge.inCycle && e1.inCycle == e3.inCycle == -1) || (e1.inCycle == e3.inCycle == edge.inCycle && e2.inCycle == e4.inCycle == -1))
            # -- update partition
            if(e2.inCycle != -1)
                indexPe3 = whichPartition(net,e3,node.number)
                part = splice!(net.partition,indexPe3) # delete partition e3
                indexPe1 = whichPartition(net,e1,node.number)
                for e in part.edges
                    push!(net.partition[indexPe1].edges,e) #put into partition e1
                end
                push!(net.partition[indexPe1].edges,edge) #put edge into partition e1
                net.partition[indexPe1].cycle = union(net.partition[indexPe1].cycle,part.cycle)
            else
                indexPe2 = whichPartition(net,e2,node.number)
                part = splice!(net.partition,indexPe2) # delete partition e2
                indexPe4 = whichPartition(net,e4,node.number)
                @debug "deleted partition $([e.number for e in part.edges]) from net.partition"
                for e in part.edges
                    @debug "partition for e4 is $([e.number for e in net.partition[indexPe4].edges]), pushing edge $(e.number)"
                    push!(net.partition[indexPe4].edges,e) #put into partition e4
                end
                @debug "partition for e4 is $([e.number for e in net.partition[indexPe4].edges]), pushing edge $(edge.number)"
                push!(net.partition[indexPe4].edges,edge) #put edge into partition e4
                @debug "added the edges of such partition to partition of e4, so new partition is $([e.number for e in net.partition[indexPe4].edges])"
                net.partition[indexPe4].cycle = union(net.partition[indexPe4].cycle,part.cycle)
                @debug "new partition has cycle $(net.partition[indexPe4].cycle)"
            end
            # -- update inCycle
            @debug "incycle e1 $(e1.inCycle), e2 $(e2.inCycle), e3 $(e3.inCycle), e4 $(e4.inCycle), edge $(edge.inCycle)"
            if(e2.inCycle == e4.inCycle == edge.inCycle)
                edge.inCycle = -1
                n1.inCycle = -1
                node.k -= 1;
            elseif(e1.inCycle == e3.inCycle == edge.inCycle)
                edge.inCycle = -1
                n2.inCycle = -1
                node.k -= 1;
            end
        else
            error("internal edge $(edge.number) is in cycle $(edge.inCycle), but it is not consistent with other edges")
        end
    else #edge not in cycle
        (e1.inCycle == e2.inCycle && e3.inCycle == e4.inCycle) || error("both edges in edges1 $([e.number for e in edges1]) (or edges2 $([e.number for e in edges2])) must have the same inCycle attribute")
        if(e1.inCycle != -1 && e3.inCycle != -1)
            @debug "cannot do tree NNI because it will create intersecting cycles, nothing done. e1,e2,e3,e4: $([e1.number,e2.number,e3.number,e4.number])"
            return false
        elseif(e1.inCycle != -1 && e3.inCycle == -1)
            node = net.node[getIndexNode(e1.inCycle,net)]
            node.hybrid || error("edge $(ed1.number) has incycle $(ed1.inCycle) but node $(node.number) is not hybrid")
            # -- update partition
            indexP = whichPartition(net,edge,node.number) # find partition where edge is
            ind = getIndex(edge,net.partition[indexP].edges)
            deleteat!(net.partition[indexP].edges,ind)
            edges = hybridEdges(n2)
            for i in 1:3 #check of 3 edges inside hybridEdges
                if(!isEqual(edges[i],e4) && !isEqual(edges[i],edge))
                    descendants = [edges[i]]
                    cycleNum = [node.inCycle]
                    getDescendants!(getOtherNode(edges[i],n2),edges[i],descendants,cycleNum)
                    !isempty(descendants) || error("descendants is empty for node $(n2.number)")
                    @debug "for node $(n2.number), descendants are $([e.number for e in descendants]), and cycleNum is $(cycleNum)"
                    partition = Partition(unique(cycleNum),descendants) # create new partition
                    push!(net.partition, partition)
                    for e in descendants #delete edges from partition with tree originally
                        ind = getIndex(e,net.partition[indexP].edges)
                        deleteat!(net.partition[indexP].edges,ind)
                    end
                    net.partition[indexP].cycle = union([node.inCycle],symdiff(net.partition[indexP].cycle,cycleNum))
                    break
                end
            end
            # -- update in cycle
            edge.inCycle = e1.inCycle
            n2.inCycle = e1.inCycle
            node.k += 1
        elseif(e3.inCycle != -1 && e1.inCycle == -1)
            node = net.node[getIndexNode(e3.inCycle,net)]
            node.hybrid || error("edge $(ed3.number) has incycle $(ed3.inCycle) but node $(node.number) is not hybrid")
            # -- update partition
            indexP = whichPartition(net,edge,node.number) # find partition where edge is
            ind = getIndex(edge,net.partition[indexP].edges)
            deleteat!(net.partition[indexP].edges,ind)
            edges = hybridEdges(n1)
            for i in 1:3 #check of 3 edges inside hybridEdges
                if(!isEqual(edges[i],e1) && !isEqual(edges[i],edge))
                    descendants = [edges[i]]
                    cycleNum = [node.inCycle]
                    getDescendants!(getOtherNode(edges[i],n1),edges[i],descendants,cycleNum)
                    !isempty(descendants) || error("descendants is empty for node $(n1.number)")
                    @debug "for node $(n1.number), descendants are $([e.number for e in descendants]), and cycleNum is $(cycleNum)"
                    partition = Partition(unique(cycleNum),descendants) # create new partition
                    push!(net.partition, partition)
                    for e in descendants #delete edges from partition with tree originally
                        ind = getIndex(e,net.partition[indexP].edges)
                        deleteat!(net.partition[indexP].edges,ind)
                    end
                    net.partition[indexP].cycle = union([node.inCycle],symdiff(net.partition[indexP].cycle,cycleNum))
                    break
                end
            end
            # -- update in cycle
            edge.inCycle = e3.inCycle
            n1.inCycle = e3.inCycle
            node.k += 1
        else
            nohybrid = true
        end
    end
    removeNode!(n1,e2)
    removeEdge!(n1,e2)
    removeNode!(n2,e3)
    removeEdge!(n2,e3)
    setNode!(e3,n1)
    setEdge!(n1,e3)
    setNode!(e2,n2)
    setEdge!(n2,e2)
    if(rand() < 0.5) # update lengths
        r = rand()
        setLength!(e1,r*t1)
        setLength!(edge,(1-r)*t1)
        setLength!(e4,t4+t)
    else
        r = rand()
        setLength!(e1,t1+t)
        setLength!(edge,(1-r)*t4)
        setLength!(e4,r*t4)
    end
    if(!isTree(net) && !nohybrid)
        undoGammaz!(node,net);
        flag,edges = updateGammaz!(net,node);
        flag || error("cannot encounter bad triangles with NNI move")
    end
    parameters!(net);
    @debug "MOVE: NNI around edge $(edge.number) SUCCESSFUL"
    return true
end

# function to repeat NNI until success
function NNIRepeat!(net::HybridNetwork,N::Integer)
    N > 0 || error("N must be positive: $(N)")
    flag,edge = chooseEdgeNNI(net,N)
    flag || return false
    i = 0
    success = false
    while(!success && i < N)
        success = NNI!(net,edge)
        i += 1
    end
    success || return false
    return true
end

