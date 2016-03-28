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
            directEdges!(net)
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
    directEdges!(net)
end

"""
`root!(net::HybridNetwork, edge::Edge)`

`root!(net::HybridNetwork, edgeNumber::Int64)`

Roots the network/tree object along an edge with number 'edgeNumber'.
This adds a new node (and a new edge) to the network.
Use plot(net, showEdgeNumber=true, showEdgeLength=false) to
visualize and identify an edge of interest.
"""
function root!(net::HybridNetwork, edgeNum::Int64)
    ind=0 # to declare outside of try/catch
    try
        ind = getIndexEdge(edgeNum,net)
    catch
        error("cannot set node $(nodeNum) as root because it is not part of net")
    end
    root!(net,net.edge[ind])
end

function root!(net::HybridNetwork, edge::Edge)
    isEdgeNumIn(edge,net.edge) || error("edge $(edge.number) not in net")
    !edge.hybrid || error("cannot put root on hybrid edge at the moment")
    node1 = edge.node[1]
    node2 = edge.node[2]
    removeEdge!(node2,edge)
    removeNode!(node2,edge)
    max_edge = maximum([e.number for e in net.edge]);
    max_node = maximum([e.number for e in net.node]);
    newedge = Edge(max_edge+1,0.0) # length 0 for new edge. half that of edge instead?
    newnode = Node(max_node+1,false,false,[edge,newedge])
    setNode!(newedge,newnode)
    setNode!(newedge,node2)
    setEdge!(node2,newedge)
    setNode!(edge, newnode)
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

"""
`directEdges!(net::HybridNetwork; checkMajor=true::Bool)`

Updates the edges' attribute isChild1, according to the root placement.
Relies on hybrid nodes having exactly 1 major hybrid parent edge,
but checks for that if checkMajor=true.

Warning: updates tree edges only. Assumes that isChild1 is correct on hybrid edges
(to avoid changing the identity of which nodes are hybrids and which are not).
"""
function directEdges!(net::HybridNetwork; checkMajor=true::Bool)
    if checkMajor # check each node has 2+ hybrid parent edges (if any), and exactly one major.
        for (n in net.node)
            nparents = 0 # 0 or 2 normally, but could be >2 if polytomy.
            nmajor = 0   # there should be exactly 1 major parent if nparents>0
            for (e in n.edge)
                if (e.hybrid && n == e.node[e.isChild1 ? 1 : 2])
                    nparents += 1
                    if (e.isMajor) nmajor +=1; end
                end
            end
            (nparents!=1) || error("node $(n.number) has exactly 1 hybrid parent edge")
            (nparents==0 || nmajor == 1) ||
              error("hybrid node $(n.number) has only 0 or 2+ major hybrid parents")
        end
    end
    for(e in net.node[net.root].edge)
        traverseDirectEdges!(net.node[net.root],e)
    end
    net.isRooted = true
end

function traverseDirectEdges!(node::Node, edge::Edge)
    if (edge.hybrid && node==edge.node[edge.isChild1 ? 1 : 2])
        error("the direction of hybrid edge $(edge.number) (isChild1) is incompatible with the root.")
    end
    if (node == edge.node[1])
        edge.isChild1 = false
        cn = edge.node[2] # cn = child node
    else
        edge.isChild1 = true
        cn = edge.node[1]
    end
    if (!cn.leaf && (!edge.hybrid || edge.isMajor))
        nchildren=0
        for (e in cn.edge)
            if e==edge continue; end
            if (e.hybrid && cn == e.node[e.isChild1 ? 1 : 2]) continue; end
            traverseDirectEdges!(cn,e)
            nchildren += 1
        end
        nchildren>0 || error("non-leaf node $(cn.number) had 0 children. It could be a hybrid
       whose parents' direction is in contradiction with the root.")
    end
    return nothing
end

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


"""
`preorder!(net::HybridNetwork)`

Updates attribute net.nodes_changed in which the nodes are pre-ordered
(also called topological sorting), such that each node is visited after its parent(s).
The edges' direction needs to be correct before calling preorder!, using directEdges!
"""
function preorder!(net::HybridNetwork)
    net.isRooted || error("net needs to be rooted for preorder!, run root! or directEdges!")
    net.nodes_changed = Node[] # path of nodes in preorder.
    queue = Node[] # problem with PriorityQueue(): dequeue() takes a
                   # random member if all have the same priority 1.
    net.visited = [false for i = 1:size(net.node,1)];
    push!(net.nodes_changed,net.node[net.root]) #push root into path
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


"""
`cladewiseorder!(net::HybridNetwork)`

Updates attribute net.cladewiseorder_nodeIndex. Used for plotting the network.
In the major tree, all nodes in a given clade are consecutive. On a tree, this function
also provides a pre-ordering of the nodes.
The edges' direction needs to be correct before calling cladewiseorder!, using directEdges!
"""
function cladewiseorder!(net::HybridNetwork)
    net.isRooted || error("net needs to be rooted for cladewiseorder!\n run root! or directEdges!")
    net.cladewiseorder_nodeIndex = Int64[]
    queue = Int64[] # index (in net) of nodes in the queue
    push!(net.cladewiseorder_nodeIndex, net.root)
    for (e in net.node[net.root].edge)
        if (e.isMajor) # follow the major tree only
            push!(queue, getIndex(getOtherNode(e,net.node[net.root]),net))
        end
    end
    # print("queued the root's children's indices: "); @show queue
    while (!isempty(queue))
        ni = pop!(queue); # deliberate choice over shift! for cladewise order
        # @show net.node[ni].number
        push!(net.cladewiseorder_nodeIndex, ni)
        for (e in net.node[ni].edge)
            if (isEqual(net.node[ni],e.node[e.isChild1 ? 2 : 1])) # net.node[ni] is parent node of e
                other = getOtherNode(e, net.node[ni])
                if (e.isMajor)
                    push!(queue, getIndex(other,net))
                    # print("queuing: "); @show other.number
                end
            end
        end
    end
end

"""
`rotate!(net::HybridNetwork, nodeNumber::Int64; edgeNumberOrder::Array{Int64,1})`

Rotates the order of the node's children edges. Useful for plotting,
to remove crossing edges.
If `node` is a tree node with no polytomy, the 2 children edges are switched
and the optional argument edgeNumberOrder is ignored.

Use plot(net, showNodeNumber=true, showEdgeNumber=false) to map node and edge numbers
on the network, as shown in the examples below.

Warning: assumes that edges are correctly directed (isChild1 updated). This is done
by plot(net). Otherwise run directEdges!(net).

# Example
```julia
net = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
plot(net, showNodeNumber=true)
rotate!(net, -4)
plot(net)

net=readTopology("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.136):10.0);");
plot(net, showNodeNumber=true, showEdgeNumber=true)
rotate!(net, -1, enumOrder=[1,12,9])
plot(net, showNodeNumber=true, showEdgeNumber=true)
rotate!(net, -3)
plot(net)
```
"""
function rotate!(net::HybridNetwork, nnum::Int64; enumOrder=Int64[]::Array{Int64,1})
    nind = 0
    try
        nind = getIndexNode(nnum,net)
    catch
        error("cannot find any node with number $nnum in network.")
    end
    n = net.node[nind]
    ci = Int[] # children edge indices
    for (i = 1:length(n.edge))
        if (n == n.edge[i].node[n.edge[i].isChild1? 2 : 1])
            push!(ci,i)
        end
    end
    if length(ci) < 2
        warn("no edge to rotate: node $nnum has $(length(ci)) children edge.")
    elseif length(ci)==2 || length(enumOrder)==0
        etmp          = n.edge[ci[1]]
        n.edge[ci[1]] = n.edge[ci[2]]
        n.edge[ci[2]] = etmp
    else # 3+ children edges and enumOrder provided
        length(enumOrder)==length(ci) || error("enumOrder $enumOrder should be of length $(length(ci))")
        length(unique(enumOrder))==length(ci) || error("enumOrder $enumOrder should not have duplicates")
        childrenedge = n.edge[ci]
        for i=1:length(ci)
            tmp = findin([e.number for e in childrenedge], enumOrder[i])
            length(tmp)==1 || error("edge number $(enumOrder[i]) not found as child of node $(n.number)")
            n.edge[ci[i]] = childrenedge[tmp[1]]
        end
    end
    return nothing
end
