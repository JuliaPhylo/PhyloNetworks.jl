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
    gamma < 0.5 || warn("adding a major hybrid edge with gamma $(gamma), this can cause problems when updating incycle")
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
    if(isEqual(edge1.node[1],edge2.node[1]) || isEqual(edge1.node[1],edge2.node[2]))
        node = edge1.node[1];
    elseif(isEqual(edge1.node[2],edge2.node[1]) || isEqual(edge1.node[2],edge2.node[2]))
        node = edge1.node[2];
    end
    if(!isa(node,Nothing))
        size(node.edge,1) == 3 || error("node found $(node.number) that does not have exactly 3 edges, it has $(size(node.edge,1)) edges instead.")
        sisters = true
        if(getOtherNode(edge1,node).leaf && getOtherNode(edge2,node).leaf)
            cherry = true
        elseif(getOtherNode(edge1,node).leaf || getOtherNode(edge2,node).leaf)
            edge = nothing
            for(e in node.edge)
                if(!isEqual(e,edge1) && !isEqual(e,edge2))
                    edge = e
                end
            end
            if(getOtherNode(edge,node).leaf)
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
    #println("net.edge length: $(length(net.edge))")
    index1 = 1;
    index2 = 1;
    inlimits = false
    inblack = true
    while(!inlimits || edges[index1].inCycle != -1 || edges[index2].inCycle != -1 || cherry || nonidentifiable || inblack)
        index1 = iround(rand()*size(edges,1));
        index2 = iround(rand()*size(edges,1));
        if(index1 != index2 && index1 != 0 && index2 != 0 && index1 <= size(edges,1) && index2 <= size(edges,1))
            inlimits = true
            sisters, cherry, nonidentifiable = sisterOrCherry(edges[index1],edges[index2]);
        else
            inlimits = false
        end
        if(blacklist && !isempty(net.blacklist))
            length(net.blacklist) % 2 == 0 || error("net.blacklist should have even number of entries, not length: $(length(net.blacklist))")
            i = 1
            while(i < length(net.blacklist))
                if(edges[index1].number == net.blacklist[i])
                    if(edges[index2].number == net.blacklist[i+1])
                        inblack = true
                    else
                        inblack = false
                    end
                elseif(edges[index2].number == net.blacklist[i])
                    if(edges[index1].number == net.blacklist[i+1])
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
    #println("from $(edge1.number) to $(edge2.number), $(gamma)");
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
    return edge3, edge4
end

## # aux function to add the hybridization
## # without checking all the updates
## # returns the hybrid node of the new hybridization
## # calls chooseEdgesGamma, parameter4createHybrid and createHybrid
## # blacklist used in afterOptBLAll
## function addHybridization!(net::HybridNetwork, blacklist::Bool)
##     edge1, edge2, gamma = chooseEdgesGamma(net, blacklist);
##     println("add hybridization between edge1, $(edge1.number) and edge2 $(edge2.number) with gamma $(gamma)")
##     edge3, edge4 = parameters4createHybrid!(edge1,edge2,net);
##     hybrid = createHybrid!(edge1, edge2, edge3, edge4, net, gamma);
##     return hybrid
## end

## addHybridization!(net::HybridNetwork) = addHybridization!(net, false)

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
        index == 0 && return nothing #no place for new hybrid
        partition = splice!(net.partition,index)
        edge1, edge2, gamma = chooseEdgesGamma(net, blacklist,partition);
    else
        edge1, edge2, gamma = chooseEdgesGamma(net, blacklist);
    end
    println("add hybridization between edge1, $(edge1.number) and edge2 $(edge2.number) with gamma $(gamma)")
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
    for(e in node.edge)
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
# returns: success bool, hybrid, flag, nocycle, flag2, flag3
function updateAllNewHybrid!(hybrid::Node,net::HybridNetwork, updatemajor::Bool, allow::Bool)
    flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,hybrid);
    if(nocycle)
        return false, hybrid, flag, nocycle, true, true
    else
        if(flag)
            if(updatemajor)
                updateMajorHybrid!(net,hybrid);
            end
            flag2, edgesGammaz = updateGammaz!(net,hybrid,allow);
            if(flag2)
                flag3, edgesRoot = updateContainRoot!(net,hybrid);
                if(flag3)
                    parameters!(net)
                    updatePartition!(net,nodesInCycle)
                    return true, hybrid, flag, nocycle, flag2, flag3
                else
                    #undoContainRoot!(edgesRoot);
                    #undoistIdentifiable!(edgesGammaz);
                    #undoGammaz!(hybrid,net);
                    #undoInCycle!(edgesInCycle, nodesInCycle);
                    return false, hybrid, flag, nocycle, flag2, flag3
                end
            else
                flag3, edgesRoot = updateContainRoot!(net,hybrid); #update contain root even if it is bad triangle to writeTopology correctly
                #undoistIdentifiable!(edgesGammaz);
                #undoGammaz!(hybrid,net);
                #undoInCycle!(edgesInCycle, nodesInCycle);
                return false, hybrid, flag, nocycle, flag2, true
            end
        else
            undoInCycle!(edgesInCycle, nodesInCycle);
            return false, hybrid, flag, nocycle, true, true
        end
    end
end

updateAllNewHybrid!(hybrid::Node,net::HybridNetwork, updatemajor::Bool) = updateAllNewHybrid!(hybrid,net, updatemajor, false)

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


# based on getDescendants on readData.jl but with vector of edges, instead of nodes
# finds the partition corresponding to the node and edge in the cycle
# used in chooseEdgesGamma and to set net.partition
function getDescendants!(node::Node, edge::Edge, descendants::Vector{Edge})
    ## if(node.leaf)
    ##     push!(descendants,edge)
    ## elseif(node.inCycle != -1)
    ##     push!(descendants,edge)
    ## else
    if(!node.leaf && node.inCycle == -1)
        for(e in node.edge)
            if(!isEqual(edge,e) && e.isMajor)
                push!(descendants,e)
                getDescendants!(getOtherNode(e,node),e,descendants)
            end
        end
    end
end

# function to update the net.partition attribute along a cycle formed
# by nodesChanged vector (obtained from updateInCycle)
function updatePartition!(net::HybridNetwork, nodesChanged::Vector{Node})
    length(nodesChanged) > 2 || error("incycle with only 2 nodes in it after updateGammaz")
    #println("nodesChanged are $([n.number for n in nodesChanged])")
    if(net.numHybrids == 0)
        net.partition = Vector{Edge}[Edge[]]
    end
    for(n in nodesChanged)
        edge = nothing
        for(e in n.edge)
            if(e.inCycle == -1)
                edge = e
            end
        end
        !isa(edge,Nothing) || error("one edge in n.edge for node $(n.number) should not be in cycle")
        descendants = [edge]
        getDescendants!(getOtherNode(edge,n),edge,descendants)
        !isempty(descendants) || error("descendants is empty for node $(n.number)")
        #println("for node $(n.number), descendants are $([e.number for e in descendants])")
        push!(net.partition, descendants)
    end
    if(isempty(net.partition[1]))
        net.partition = net.partition[2:end]
    end
    #println("length of net.partition at the end is $(length(net.partition))")
end

function choosePartition(net::HybridNetwork)
    all([length(n) == 1 for n in net.partition]) && return 0
    index1 = iround(rand()*size(net.partition,1));
    while(index1 == 0 || index1 > length(net.partition) || length(net.partition[index1]) < 2)
        index1 = iround(rand()*size(net.partition,1));
    end
    return index1
end

