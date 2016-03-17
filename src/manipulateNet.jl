# Rooting functions previously in descriptive.jl
# function to re root on a node
# resolve=true, a branch of length=0 is
# added if node is internal if node is leaf, root(net,outgroup) is
# called, so a new node is created on external edge
function root!(net::HybridNetwork, node::Node, resolve::Bool)
    node.hybrid && error("node $(node.number) is hybrid, cannot root network on a hybrid node")
    if(node.leaf)
        warn("node $(node.number) is a leaf, so we will root as an outgroup if possible")
        root!(net,node.name)
    else
        if(!isTree(net))
            if(!net.cleaned)
                DEBUG && println("net not cleaned inside root, need to run updateCR")
                for(n in net.hybrid)
                    flag,edges = updateContainRoot!(net,n)
                    flag || error("hybrid node $(n.hybrid) has conflicting containRoot")
                end
            end
        end
        if(canBeRoot(node))
            try
                ind = getIndex(node,net)
            catch
                error("cannot set node $(node.number) as root because it is not part of net")
            end
            ind = getIndex(node,net)
            net.root = ind
            if(resolve)
                resolve!(net,node)
            end
        else
            warn("node $(node.number) cannot be root, will leave root as is")
        end
    end
end

"""
`root!(net::HybridNetwork, nodeNumber::Int64, resolve::Bool)`

`root!(net::HybridNetwork, nodeNumber::Int64)`

Roots the network/tree object at the node with number 'nodeNumber'.
With resolve=true, the polytomy at the root is resolved arbitrarily (??) with a branch of length 0.
The second version uses resolve=false, that is, a polytomy is left at the root.
"""
function root!(net::HybridNetwork, nodeNum::Int64, resolve::Bool)
    try
        ind = getIndexNode(nodeNum,net)
    catch
        error("cannot set node $(nodeNum) as root because it is not part of net")
    end
    ind = getIndexNode(nodeNum,net)
    root!(net,net.node[ind],resolve)
end

root!(net::HybridNetwork, nodeNum::Int64) = root!(net, nodeNum, false)

# function to resolve an internal node for root function
function resolve!(net::HybridNetwork, node::Node)
    length(node.edge) == 3 || error("node $(node.number) has $(length(node.edge)) edges instead of 3")
    !node.hybrid || error("node $(node.number) is hybrid and cannot root/resolve on hybrid node")
    canBeRoot(node) || error("cannot resolve root in node $(node.number) because it cannot be root to begin with")
    if(net.cleaned)
        if(node.inCycle == -1) #node not in cycle
            done=false
            for(e in node.edge)
                if(e.containRoot)
                    done = true
                    newNodeResolve!(net,e,node)
                    net.root = length(net.node) #last node is root
                    break
                end
            end
            done || error("not found edge to contain root for node $(node.number)")
        else
            done = false
            for(e in node.edge)
                if(e.inCycle == -1 && e.containRoot)
                    done = true
                    newNodeResolve!(net,e,node)
                    net.root = length(net.node) #last node is root
                    break
                end
            end
            done || error("strange: not found edge to contain root not in cycle for node $(node.number) even when node can be root")
        end
    else # net not cleaned
        if(!node.hasHybEdge)
            done=false
            for(e in node.edge)
                if(!e.hybrid)
                    done = true
                    newNodeResolve!(net,e,node)
                    net.root = length(net.node) #last node is root
                    break
                end
            end
            done || error("not found edge to contain root for node $(node.number)")
        else
            done = false
            for(e in node.edge)
                if(e.inCycle == -1 && e.containRoot)
                    done = true
                    newNodeResolve!(net,e,node)
                    net.root = length(net.node) #last node is root
                    break
                end
            end
            done || error("strange: not found edge to contain root not in cycle for node $(node.number) even when node can be root")
        end
    end
end

# function to create new node in edge for resolve
function newNodeResolve!(net::HybridNetwork,e::Edge, node::Node)
    removeEdge!(node,e)
    removeNode!(node,e)
    max_edge = maximum([e.number for e in net.edge]);
    max_node = maximum([e.number for e in net.node]);
    newedge = Edge(max_edge+1)
    newnode = Node(max_node+1,false,false,[e,newedge])
    if(net.cleaned && !isTree(net))
        part = whichPartition(net,e)
        push!(net.partition[part].edges,newedge)
    end
    setNode!(e,newnode)
    setNode!(newedge,newnode)
    setEdge!(node,newedge)
    setNode!(newedge,node)
    pushEdge!(net,newedge)
    pushNode!(net,newnode)
    t = e.length
    setLength!(e,t/2)
    setLength!(newedge,t/2)
end

# function to root a network on an outgroup
# (single taxon)
"""
`root!(net::HybridNetwork, outgroup::AbstractString)`

Roots the network/tree object along the external edge leading to the taxon named 'outgroup'.
"""
function root!(net::HybridNetwork, outgroup::AbstractString)
    if(!isTree(net))
        if(!net.cleaned)
            DEBUG && println("net not cleaned inside root, need to run updateCR")
            for(n in net.hybrid)
                flag,edges = updateContainRoot!(net,n)
                flag || error("hybrid node $(n.hybrid) has conflicting containRoot")
            end
        end
    end
    updateRoot!(net,outgroup)
end

# function to root in an edge
function root!(net::HybridNetwork, edge::Edge)
    isEdgeNumIn(edge,net.edge) || error("edge $(edge.number) not in net")
    !edge.hybrid || error("cannot put root on hybrid edge at the moment")
    node1 = edge.node[1]
    node2 = edge.node[2]
    removeEdge!(node2,edge)
    removeNode!(node2,edge)
    max_edge = maximum([e.number for e in net.edge]);
    max_node = maximum([e.number for e in net.node]);
    newedge = Edge(max_edge+1)
    newnode = Node(max_node+1,false,false,[edge,newedge])
    setNode!(newedge,node2)
    setEdge!(node2,newedge)
    pushEdge!(net,newedge)
    pushNode!(net,newnode)
    if(edge.inCycle != -1)
        newedge.inCycle = edge.inCycle
        newnode.inCycle = edge.inCycle
    end
    root!(net,newnode,false)
end



# Claudia SL & Paul Bastide: November 2015, Cecile: Feb 2016

#################################################
# Direct Edges
#################################################

function traverseDirectEdges!(node::Node, edge::Edge)
    if(isEqual(node,edge.node[1]))
        DEBUG && println("changing edge $(edge.number)")
        edge.isChild1 = false
    else
        DEBUG && println("changing edge $(edge.number)")
        edge.isChild1 = true
    end
    node2 = getOtherNode(edge,node)
    if(!node2.leaf)
        for(e in node2.edge)
            if(!isEqual(e,edge) && !e.hybrid)
                traverseDirectEdges!(node2,e)
            end
        end
    end
end

function directEdges!(net::HybridNetwork)
    for(e in net.node[net.root].edge)
        traverseDirectEdges!(net.node[net.root],e)
    end
    net.isRooted = true
end

# fixit: add this to root! function

#################################################
## Topological sorting
#################################################

# function to get all parent nodes of a given node
# it assumes the isChild1 attributes are correct
function getParents(node::Node)
    parents = Node[]
    for(e in node.edge)
            if(isEqual(node,e.isChild1 ? e.node[1] : e.node[2])) #node is child of e
                push!(parents,getOtherNode(e,node))
            end
    end
    return parents
end

# function to get only one parent (the major if hybrid) of a given node
# assumes isChild1 and isMajor attributes are correct
function getMajorParent(n::Node)
    found = false
    for (e in node.edge)
        if (isEqual(n, e.isChild1 ? e.node[1] : e.node[2]) # n is child of e
            && e.isMajor) # in case n is a hybrid, e is major parent edge
            found = true
            break
        end
    end
    found || error("node $(n.number) has no major parent")
    return getOtherNode(e,n)
end


# function to order nodes in topological sorting
# saves vector of nodes in the right order in net.nodes_changed
"""
`preorder!(net::HybridNetwork)`

Updates attributes of the network to calculate a pre-ordering of the nodes
(also called topological sorting), such that each node is visited after its parent(s).
The edges' direction needs to be correct before calling preorder!, using directEdges!
"""
function preorder!(net::HybridNetwork)
    net.isRooted || error("net needs to be rooted for preorder, run root! or directEdges!")
    net.nodes_changed = Node[] # path of nodes in preorder.
    net.preorder_nodeIndex = Int64[] # corresponding order, but with node indices
    net.preorder_edgeIndex = Int64[] # Major edges only: major parent of corresponding node
    queue = Node[] # problem with PriorityQueue(): dequeue() takes a
                   # random member if all have the same priority 1.
    net.visited = [false for i = 1:size(net.node,1)];
    push!(net.nodes_changed,net.node[net.root]) #push root into path
    push!(net.preorder_nodeIndex,    net.root )
    push!(net.preorder_edgeIndex,    0) # root has no parent edge: giving bad index 0
    net.visited[net.root] = true   # visit root
    for(e in net.node[net.root].edge)
        if (!e.hybrid) # enqueue child of root, if child is not a hybrid.
            push!(queue,getOtherNode(e,net.node[net.root]))
        end # if child is hybrid, its second parent has not been visited yet.
    end
    # print("queued the root's children: "); @show queue
    while(!isempty(queue))
        #println("at this moment, queue is $([n.number for n in queue])")
        curr = pop!(queue); # deliberate choice over shift! for cladewise order
        # @show curr.number
        net.visited[getIndex(curr,net)] = true # visit curr node
        push!(net.nodes_changed,curr) #push curr into path
        push!(net.preorder_nodeIndex, getIndex(curr,net))
        for (e in curr.edge)
            if (e.isMajor && curr == e.node[e.isChild1 ? 1 : 2]) # e is major parent edge of curr node
                push!(net.preorder_edgeIndex, getIndex(e,net))
                break
            end
        end
        for(e in curr.edge)
            if(isEqual(curr,e.node[e.isChild1 ? 2 : 1])) # curr is the parent node if e
                other = getOtherNode(e,curr)
                if(!e.hybrid)
                    push!(queue,other)
                    # print("queuing: "); @show other.number
                else
                    for(e2 in other.edge) # find other hybrid parent edge for 'other'
                        if(e2.hybrid && !isEqual(e,e2))
                            parent = getOtherNode(e2,other)
                            if(net.visited[getIndex(parent,net)])
                                push!(queue,other)
                                # print("queuing: "); @show other.number
                            end
                            break
                        end
                    end
                end
            end
        end
    end
    # println("path of nodes is $([n.number for n in net.nodes_changed])")
end
