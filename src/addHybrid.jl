# functions to add a hybridization
# originally in functions.jl
# Claudia March 2015


# ------------------------------------------------ add new hybridization------------------------------------

# function to add hybridization event
# input: edge1, edge2 are the edges to remove
#        edge3, edge4, edge5 are the new tree edges to add
#        net is the network
#        gamma is the gamma for the hybridization
# warning: assumes that edge1, edge2 are tree edges with inCycle=-1
#          assumes the hybrid edge goes from edge1 to edge2
#          sets minor hybrid edge length to zero
# this function create the hybrid node/edges and connect everything
# and deletes edge1,2 from the nodes, and removes the nodes from edge1,2
# returns the hybrid node to start future updates there
function createHybrid!(edge1::Edge, edge2::Edge, edge3::Edge, edge4::Edge, net::HybridNetwork, gamma::Float64)
    0 < gamma < 1 || error("gamma must be between 0 and 1: $(gamma)")
    (edge1.hybrid || edge2.hybrid) ? error("edges to delete must be tree edges") : nothing
    (edge3.hybrid || edge4.hybrid) ? error("edges to add must be tree edges") : nothing
    pushEdge!(net,edge3);
    pushEdge!(net,edge4);
    # create hybridization
    max_node = maximum([e.number for e in net.node]);
    max_edge = maximum([e.number for e in net.edge]);
    gamma < 0.5 || @warn "adding a major hybrid edge with gamma $(gamma), this can cause problems when updating incycle"
    hybrid_edge = Edge(max_edge+1,0.0,true,gamma,gamma>=0.5);
    pushEdge!(net,hybrid_edge);
    hybrid_node = Node(max_node+1,false,true,[edge2,hybrid_edge,edge4]);
    tree_node = Node(max_node+2,false,false,[edge1,edge3,hybrid_edge]);
    setNode!(hybrid_edge,[tree_node,hybrid_node]);
    setNode!(edge3,[tree_node,edge1.node[2]]);
    setNode!(edge4,[hybrid_node,edge2.node[2]]);
    setEdge!(edge1.node[2],edge3);
    setEdge!(edge2.node[2],edge4);
    removeEdge!(edge2.node[2],edge2);
    removeEdge!(edge1.node[2],edge1);
    removeNode!(edge1.node[2],edge1);
    setNode!(edge1,tree_node);
    removeNode!(edge2.node[2],edge2);
    #[n.number for n in edge2.node]
    setNode!(edge2,hybrid_node)
    pushNode!(net,hybrid_node);
    pushNode!(net,tree_node);
    #pushHybrid!(net,hybrid_node);
    return hybrid_node
end

# aux function for chooseEdgesGamma to identify
# if two edges are sisters and if they are cherry
# (have leaves)
# returns true/false for sisters, true/false for cherry
#         true/false for nonidentifiable (two leaves, k=1 node crossed by hybridization)
function sisterOrCherry(edge1::Edge,edge2::Edge)
    sisters = false
    cherry = false
    nonidentifiable = false
    node = nothing;
    if isEqual(edge1.node[1],edge2.node[1]) || isEqual(edge1.node[1],edge2.node[2])
        node = edge1.node[1];
    elseif isEqual(edge1.node[2],edge2.node[1]) || isEqual(edge1.node[2],edge2.node[2])
        node = edge1.node[2];
    end
    if node !== nothing
        size(node.edge,1) == 3 || error("node found $(node.number) that does not have exactly 3 edges, it has $(size(node.edge,1)) edges instead.")
        sisters = true
        if getOtherNode(edge1,node).leaf && getOtherNode(edge2,node).leaf
            cherry = true
        elseif getOtherNode(edge1,node).leaf || getOtherNode(edge2,node).leaf
            edge = nothing
            for e in node.edge
                if(!isEqual(e,edge1) && !isEqual(e,edge2))
                    edge = e
                end
            end
            if getOtherNode(edge,node).leaf
                nonidentifiable = true
            end
        end
    end
    return sisters, cherry, nonidentifiable
end

# aux function to addHybridization
# it chooses the edges in the network and the gamma value
# warning: chooses edge1, edge2, gamma randomly, but
#          we could do better later
# check: gamma is uniform(0,1/2) to avoid big gammas
# fixit: add different function to choose gamma
# fixit: how to stop from infinite loop if there are no options
# blacklist used for afterOptBLAll
# input: edges, list of edges from which to choose, default is net.edge
# warning: if edges is not net.edge, it still need to contain Edge objects from net (not deepcopies)
function chooseEdgesGamma(net::HybridNetwork, blacklist::Bool, edges::Vector{Edge})
    index1 = 1;
    index2 = 1;
    inlimits = false
    inblack = true
    cherry = false
    nonidentifiable = false
    while !inlimits || edges[index1].inCycle != -1 || edges[index2].inCycle != -1 || inblack || cherry || nonidentifiable
        index1 = round(Integer,rand()*size(edges,1));
        index2 = round(Integer,rand()*size(edges,1));
        if index1 != index2 && index1 != 0 && index2 != 0 && index1 <= size(edges,1) && index2 <= size(edges,1)
            inlimits = true
            sisters, cherry, nonidentifiable = sisterOrCherry(edges[index1],edges[index2]);
        else
            inlimits = false
        end
        if blacklist && !isempty(net.blacklist)
            length(net.blacklist) % 2 == 0 || error("net.blacklist should have even number of entries, not length: $(length(net.blacklist))")
            i = 1
            while i < length(net.blacklist)
                if edges[index1].number == net.blacklist[i]
                    if edges[index2].number == net.blacklist[i+1]
                        inblack = true
                    else
                        inblack = false
                    end
                elseif edges[index2].number == net.blacklist[i]
                    if edges[index1].number == net.blacklist[i+1]
                        inblack = true
                    else
                        inblack = false
                    end
                end
                i += 2
            end
        else
            inblack = false
        end
    end
    gamma = rand()*0.5;
    @debug "choose edges and gamma: from $(edges[index1].number) to $(edges[index2].number), $(gamma)"
    return edges[index1],edges[index2],gamma
end

chooseEdgesGamma(net::HybridNetwork) = chooseEdgesGamma(net, false, net.edge)
chooseEdgesGamma(net::HybridNetwork, blacklist::Bool) = chooseEdgesGamma(net, blacklist, net.edge)

# aux function for addHybridization
# that takes the output edge1, edge2, gamma from
# chooseEdgesGamma and created necessary edges
# returns edge3, edge4, and adjusts edge1, edge2 to shorter length
function parameters4createHybrid!(edge1::Edge, edge2::Edge,net::HybridNetwork)
    max_edge = maximum([e.number for e in net.edge]);
    t1 = rand()*edge1.length;
    t3 = edge1.length - t1;
    edge3 = Edge(max_edge+1,t3);
    edge1.length = t1;
    t1 = rand()*edge2.length;
    t3 = edge2.length - t1;
    edge4 = Edge(max_edge+2,t3);
    edge2.length = t1;
    edge3.containRoot = edge1.containRoot
    edge4.containRoot = edge2.containRoot
    return edge3, edge4
end

# aux function to add the hybridization
# without checking all the updates
# returns the hybrid node of the new hybridization
# calls chooseEdgesGamma, parameter4createHybrid and createHybrid
# blacklist used in afterOptBLAll
# usePartition=true if we use the information on net.partition, default true
function addHybridization!(net::HybridNetwork, blacklist::Bool, usePartition::Bool)
    if(net.numHybrids > 0 && usePartition)
        !isempty(net.partition) || error("net has $(net.numHybrids) but net.partition is empty")
        index = choosePartition(net)
        if(index == 0) #no place for new hybrid
            @debug "no partition suitable to place new hybridization"
            return nothing
        end
        partition = splice!(net.partition,index) #type partition
        @debug "add hybrid with partition $([n.number for n in partition.edges])"
        edge1, edge2, gamma = chooseEdgesGamma(net, blacklist,partition.edges);
    else
        edge1, edge2, gamma = chooseEdgesGamma(net, blacklist);
    end
    @debug "add hybridization between edge1, $(edge1.number) and edge2 $(edge2.number) with gamma $(gamma)"
    edge3, edge4 = parameters4createHybrid!(edge1,edge2,net);
    hybrid = createHybrid!(edge1, edge2, edge3, edge4, net, gamma);
    return hybrid
end

addHybridization!(net::HybridNetwork) = addHybridization!(net, false, true)

# function to update who is the major hybrid
# after a new hybridization is added and
# inCycle is updated
# warning: needs updateInCycle! run first
# can return the updated edge for when undoing network moves, not needed now
function updateMajorHybrid!(net::HybridNetwork, node::Node)
    node.hybrid || error("node $(node.number) is not hybrid, cannot update major hybrid after updateInCycle")
    length(node.edge) == 3 || error("hybrid node $(node.number) has $(length(node.edge)) edges, should have 3")
    hybedge = nothing
    edgecycle = nothing
    for e in node.edge
        if(e.hybrid)
            hybedge = e
        elseif(e.inCycle != -1 && !e.hybrid)
            edgecycle = e
        end
    end
    !isa(hybedge,Nothing) || error("hybrid node $(node.number) does not have hybrid edge")
    !isa(edgecycle,Nothing) || error("hybrid node $(node.number) does not have tree edge in cycle to update to hybrid edge after updateInCycle")
    #println("updating hybrid status to edgeincycle $(edgecycle.number) for hybedge $(hybedge.number)")
    makeEdgeHybrid!(edgecycle,node,1-hybedge.gamma)
end

# function to update everything of a new hybridization
# it follows the flow diagram in ipad
# input: new added hybrid, network,
#        updatemajor (bool) to decide if we need to update major edge
#        only need to update if new hybrid added, if read from file not needed
#        allow=true allows extreme/very bad triangles, needed when reading
#        updatePart = true will update PArtition at this moment, it makes sense with a newly added hybrid
#        but not if net just read (because in this case it needs all inCycle updated before)
# returns: success bool, hybrid, flag, nocycle, flag2, flag3
function updateAllNewHybrid!(hybrid::Node,net::HybridNetwork, updatemajor::Bool, allow::Bool, updatePart::Bool)
    flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,hybrid);
    if(nocycle)
        return false, hybrid, flag, nocycle, false, false
    else
        if(flag)
            if(updatemajor)
                updateMajorHybrid!(net,hybrid);
            end
            flag2, edgesGammaz = updateGammaz!(net,hybrid,allow);
            if(flag2)
                flag3, edgesRoot = updateContainRoot!(net,hybrid);
                if(updatePart)
                    updatePartition!(net,nodesInCycle)
                end
                if(flag3)
                    parameters!(net)
                    return true, hybrid, flag, nocycle, flag2, flag3
                else
                    #undoContainRoot!(edgesRoot);
                    #undoistIdentifiable!(edgesGammaz);
                    #undoGammaz!(hybrid,net);
                    #undoInCycle!(edgesInCycle, nodesInCycle);
                    return false, hybrid, flag, nocycle, flag2, flag3
                end
            else
                if(updatePart)
                    updatePartition!(net,nodesInCycle)
                end
                flag3, edgesRoot = updateContainRoot!(net,hybrid); #update contain root even if it is bad triangle to writeTopologyLevel1 correctly
                #undoistIdentifiable!(edgesGammaz);
                #undoGammaz!(hybrid,net);
                #undoInCycle!(edgesInCycle, nodesInCycle);
                return false, hybrid, flag, nocycle, flag2, flag3
            end
        else
            # no need to do updatePartition in this case, because we only call deleteHybrid after
            undoInCycle!(edgesInCycle, nodesInCycle);
            return false, hybrid, flag, nocycle, true, true
        end
    end
end

updateAllNewHybrid!(hybrid::Node,net::HybridNetwork, updatemajor::Bool) = updateAllNewHybrid!(hybrid,net, updatemajor, false, true)

# function to add a new hybridization event
# it calls chooseEdgesGamma and createHybrid!
# input: network
# check: assumes that one of the two possibilities for
#        major hybrid edge gives you a cycle, true?
# warning: "while" removed, it does not attempt to add until
#          success, it attempts to add once
# returns: success (bool), hybrid, flag, nocycle, flag2, flag3
# blacklist used in afterOptBLAll
function addHybridizationUpdate!(net::HybridNetwork, blacklist::Bool, usePartition::Bool)
    hybrid = addHybridization!(net,blacklist, usePartition);
    isa(hybrid,Nothing) && return false,nothing,false,false,false,false
    updateAllNewHybrid!(hybrid,net,true)
end

addHybridizationUpdate!(net::HybridNetwork) = addHybridizationUpdate!(net, false, true)
addHybridizationUpdate!(net::HybridNetwork, blacklist::Bool) = addHybridizationUpdate!(net, blacklist::Bool, true)




# function that will add a hybridization with addHybridizationUpdate,
# if success=false, it will try to move the hybridization before
# declaring failure
# blacklist used in afterOptBLAll
function addHybridizationUpdateSmart!(net::HybridNetwork, blacklist::Bool, N::Integer)
    global CHECKNET
    @debug "MOVE: addHybridizationUpdateSmart"
    success, hybrid, flag, nocycle, flag2, flag3 = addHybridizationUpdate!(net, blacklist)
    @debug begin
        printEverything(net)
        "success $(success), flag $(flag), flag2 $(flag2), flag3 $(flag3)"
    end
    i = 0
    if !success
        if isa(hybrid,Nothing)
            @debug "MOVE: could not add hybrid by any means"
        else
            while((nocycle || !flag) && i < N) #incycle failed
                @debug "MOVE: added hybrid causes conflict with previous cycle, need to delete and add another"
                deleteHybrid!(hybrid,net,true)
                success, hybrid, flag, nocycle, flag2, flag3 = addHybridizationUpdate!(net, blacklist)
            end
            if(nocycle || !flag)
                @debug "MOVE: added hybridization $(i) times trying to avoid incycle conflicts, but failed"
            else
                if(!flag3 && flag2) #containRoot failed
                    @debug "MOVE: added hybrid causes problems with containRoot, will change the direction to fix it"
                    success = changeDirectionUpdate!(net,hybrid) #change dir of minor
                elseif(!flag2 && flag3) #gammaz failed
                    @debug "MOVE: added hybrid has problem with gammaz (not identifiable bad triangle)"
                    if(flag3)
                        @debug "MOVE: we will move origin to fix the gammaz situation"
                        success = moveOriginUpdateRepeat!(net,hybrid,true)
                    else
                        @debug "MOVE: we will move target to fix the gammaz situation"
                        success = moveTargetUpdateRepeat!(net,hybrid,true)
                    end
                elseif(!flag2 && !flag3) #containRoot AND gammaz failed
                    @debug "MOVE: containRoot and gammaz both fail"
                end
            end
            if !success
                @debug "MOVE: could not fix the added hybrid by any means, we will delete it now"
                CHECKNET && checkNet(net)
                @debug begin printEverything(net); "printed everything" end
                deleteHybridizationUpdate!(net,hybrid)
                @debug begin printEverything(net); "printed everything" end
                CHECKNET && checkNet(net)
            end
        end
    end
    success && @debug "MOVE: added hybridization SUCCESSFUL: new hybrid $(hybrid.number)"
    return success
end

addHybridizationUpdateSmart!(net::HybridNetwork, N::Integer) = addHybridizationUpdateSmart!(net, false,N)


# ----------------------------------- add alternative hybridizations found in bootstrap ------------------------------------
## function to add the top hybridizations (not in best network)
## found in bootstrap networks
## it modifies the network and BSe to add the new edge number,
## it also adds a new column 0/1, with 1 if the edge is not present in the original estimated network
## cutoff=10, only show hybridizations with bootstrap support>10%
## top=3, only show the top three alternative hybridizations (sorted by bootstrap support)
## can we put in a different color in the plot?
## Cecile: you could add a new option to our RCall-based plot function to take a vector of colors.
"""
    addAlternativeHybridizations!(net::HybridNetwork, BSe::DataFrame;
                                  cutoff=10::Number, top=3::Int)

Modify the network `net` (the best network estimated with snaq) by adding other hybridizations
that are present in the bootstrap networks. By default, it will only consider hybrid edges with
more than 10% bootstrap support (`cutoff`) and it will only include the three top hybridizations
(`top`) sorted by bootstrap support.
The function also modifies the dataframe `BSe` obtained with `hybridBootstrapSupport`. In the original
`BSe` dataframe, hybrid edges that do not appear in the best network have a missing number.
After the hybrid edges are added with `addAlternativeHybridizations`, `BSe` is modified to include the
edge numbers of the newly added hybrid edges. Note that the function only adds the hybrid edges as minor to
keep the underlying tree topology.

# example

```julia
bootnet = readMultiTopology("bootstrap-networks.txt")
bestnet = readTopology("best.tre")
BSn, BSe, BSc, BSgam, BSedgenum = hybridBootstrapSupport(bootnet, bestnet);
addAlternativeHybridizations!(bestnet,BSe)
using PhyloPlots
plot(bestnet, edgeLabel=BSe[[:edge,:BS_hybrid_edge]])
```
"""
function addAlternativeHybridizations!(net::HybridNetwork,BSe::DataFrame; cutoff=10::Number,top=3::Int)
    top > 0 || error("top must be greater than 0")
    BSe[:alternative] = falses(length(BSe[:hybrid]))
    newBSe = BSe[BSe[:BS_hybrid_edge] .> cutoff,:]
    newBSe = newBSe[.!ismissing.(newBSe[:hybrid]) & .!ismissing.(newBSe[:sister]),:]
    newBSe = newBSe[ismissing.(newBSe[:edge]),:]
    newHyb = newBSe[1:top,:]

    if(size(newHyb,1) == 0)
        @warn "Did not find any alternative hybridizations with bootstrap support greater than the cutoff, so nothign added"
        return
    end

    for i in 1:top
        hybnum = newHyb[:hybrid][i]
        sisnum = newHyb[:sister][i]
        edgenum = addHybridBetweenClades!(hybnum,sisnum,net)
        ind1 = findall(x->!ismissing(x) && x==hybnum,BSe[:hybrid])
        ind2 = findall(x->!ismissing(x) && x==sisnum,BSe[:sister])
        ind = intersect(ind1,ind2)
        BSe[ind,:edge] = edgenum
        BSe[ind,:alternative] = true
    end
end


## function to a hybrid edge between hybrid clade (hybnum = node number) and sister clade (sisnum = node number)
## the function finds the parent edges to these nodes, and puts a hybrid edge between them
## the function modifies net, and it returns the number of the minor hybrid edge added
function addHybridBetweenClades!(hybnum::Number,sisnum::Number,net::HybridNetwork)
    hybind = getIndexNode(hybnum,net)
    sisind = getIndexNode(sisnum,net)

    ## hybridization ed1->ed2
    edge1 = getMajorParentEdge(net.node[sisind])
    edge2 = getMajorParentEdge(net.node[hybind])

    edge3,edge4 = parameters4createHybrid!(edge1, edge2,net)
    hybridnode = createHybrid!(edge1, edge2, edge3, edge4, net, 0.1) ## gamma=0.1, fixed later
    if(edge1.length < 0)
        setBranchLength!(edge1,-1.0)
        setBranchLength!(edge3,-1.0)
    end
    if(edge2.length < 0)
        setBranchLength!(edge2,-1.0)
        setBranchLength!(edge4,-1.0)
    end

    if(edge2.isChild1)
        edge4.hybrid = true
        setGamma!(edge4,0.9)
        edge4.isChild1 = true
    else
        edge2.hybrid = true
        setGamma!(edge2,0.9)
        edge2.isChild1 = false
    end
    ## used gamma=0.1 to make the new edge a minor edge, but we really do not have gamma value:
    emaj = getMajorParentEdge(hybridnode)
    emaj.gamma = -1
    e = getMinorParentEdge(hybridnode)
    e.gamma = -1
    return e.number
end
