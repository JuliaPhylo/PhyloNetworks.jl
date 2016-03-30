"""
`rootatnode!(HybridNetwork, nodeNumber::Int64; index=false::Bool)`

`rootatnode!(HybridNetwork, Node)`

`rootatnode!(HybridNetwork, nodeName::AbstractString)`

Roots the network/tree object at the node with name 'nodeName' or
number 'nodeNumber' (by default) or with index 'nodeNumber' if index=true.
Attributes isChild1 and containRoot are updated along the way.
Use `plot(net, showNodeNumber=true, showEdgeLength=false)` to
visualize and identify a node of interest.

Warnings:
- If the node is a leaf, the root will be placed along
  the edge adjacent to the leaf, with a message. This might add a new node.
- If the desired root placement is incompatible with one or more hybrids,
  then the network will still have some attributes modified.

Return value: like directEdges, returns `true` if successful,
`false` if the desired root placement was in
conflict with the direction of any hybrid edge.

See also: `rootonedge!`.
"""
function rootatnode!(net::HybridNetwork, node::Node)
    rootatnode!(net, node.number, index=false)
end

function rootatnode!(net::HybridNetwork, nodeName::AbstractString)
    tmp = findin([n.name for n in net.node], [nodeName])
    if length(tmp)==0
        error("node named $nodeName was not found in the network.")
    elseif length(tmp)>1
        error("several nodes were found with name $nodeName.")
    end
    rootatnode!(net, tmp[1], index=true)
end

function rootatnode!(net::HybridNetwork, nodeNumber::Int64; index=false::Bool)
    ind = nodeNumber # good if index=true
    if !index
      try
        ind = getIndexNode(nodeNumber,net)
      catch
        error("cannot set node $(nodeNumber) as root because it is not part of net")
      end
    elseif ind > length(net.node)
        error("node index $ind too large: the network only has $(length(net.node)) nodes.")
    end
    if net.node[ind].leaf
        println("node $(net.node[ind].number) is a leaf. Will create a new node if needed, to set taxon \"$(net.node[ind].name)\" as outgroup.")
        length(net.node[ind].edge)==1 || error("leaf has $(length(net.node[ind].edge)) edges!")
        pn = getOtherNode(net.node[ind].edge[1], net.node[ind])
        if length(pn.edge) <= 2 # if parent of leaf has degree 2, use it as new root
            rootatnode!(net, pn.number)
        else # otherwise, create a new node between leaf and its parent
            rootonedge!(net,net.node[ind].edge[1])
        end
    else
        rootsaved = net.root
        net.root = ind
        res = directEdges!(net)
        if !res # new root incompatible with hybrid directions: revert back.
            net.root = rootsaved
        elseif (net.root != rootsaved && length(net.node[rootsaved].edge)==2)
            fuseedgesat!(rootsaved,net) # remove old root node if degree 2
        end
        return(res)
    end
end


"""
`rootonedge!(HybridNetwork, edgeNumber::Int64; index=false::Bool)`

`rootonedge!(HybridNetwork, Edge)`

Roots the network/tree object along an edge with number 'edgeNumber' (by default)
or with index 'edgeNumber if index=true. Attributes isChild1 and containRoot are
updated along the way.

This adds a new node and a new edge to the network.
Use `plot(net, showEdgeNumber=true, showEdgeLength=false)` to
visualize and identify an edge of interest.

See also: `rootatnode!`.
"""
function rootonedge!(net::HybridNetwork, edge::Edge)
    rootonedge!(net, edge.number, index=false)
end

function rootonedge!(net::HybridNetwork, edgeNumber::Int64; index=false::Bool)
    ind = edgeNumber # good if index=true
    if !index
      try
        ind = getIndexEdge(edgeNumber,net)
      catch
        error("cannot set root along edge $(edgeNumber): no such edge in network")
      end
    elseif ind > length(net.edge)
        error("edge index $ind too large: the network only has $(length(net.edge)) edges.")
    end
    rootsaved = net.root
    net.root = breakedge!(net.edge[ind],net)
    res = directEdges!(net)
    if !res
        fuseedgesat!(net.root,net) # reverts breakedge!
        net.root = rootsaved
    elseif (net.root != rootsaved && length(net.node[rootsaved].edge)==2)
        fuseedgesat!(rootsaved,net) # remove old root node if degree 2
    end
    return(res)
end

"""
`breakedge!(Edge, HybridNetwork)`

breaks an edge into 2 edges (each of length half that of original edge).
creates new node of degree 2. Useful to root network along an edge.

warning: updates `isChild1` and `containRoot`, but
does NOT update attributes like: inCycle, partition, gammaz, etc.
"""
function breakedge!(edge::Edge, net::HybridNetwork)
    pn = edge.node[edge.isChild1 ? 2 : 1] # parent node
    # new child edge = old edge, same hybrid attribute
    removeEdge!(pn,edge)
    removeNode!(pn,edge)
    max_edge = maximum([e.number for e in net.edge]);
    max_node = maximum([e.number for e in net.node]);
    newedge = Edge(max_edge+1) # create new parent (tree) edge
    newnode = Node(max_node+1,false,false,[edge,newedge]) # tree node
    setNode!(edge,newnode) # newnode comes 2nd, and parent node along 'edge'
    edge.isChild1 = true
    setNode!(newedge,newnode) # newnode comes 1st in newedge, but child node
    newedge.isChild1 = true
    setEdge!(pn,newedge)
    setNode!(newedge,pn) # pn comes 2nd in newedge
    if edge.length == -1.0
        newedge.length = -1.0
    else
        edge.length /= 2
        newedge.length = edge.length
    end
    newedge.containRoot = edge.containRoot
    pushEdge!(net,newedge)
    pushNode!(net,newnode)
    return length(net.node) # index of newnode: last one pushed to net.node
end

"""
`fuseedgesat!(i::Int64,net::HybridNetwork)`

Removes `i`th node in net.node, if it is of degree 2.
The parent and child edges of this node are fused.
Reverts the action of breakedge!.
"""
function fuseedgesat!(i::Int64, net::HybridNetwork)
    i <= length(net.node) ||
      error("node index $i too large: only $(length(net.node)) nodes in the network.")
    length(net.node[i].edge) == 2 ||
      error("can't fuse edges at node number $(net.node[i].number): connected to $(length(net.node[i].edge)) edges.")
    !(net.node[i].edge[1].hybrid && net.node[i].edge[2].hybrid) ||
      error("can't fuse edges at node number $(net.node[i].number): connected to exactly 2 hybrid edges")
    j = indmax([e.number for e in net.node[i].edge])
    pe = net.node[i].edge[j] # edge to remove: pe.number > ce.number
    if pe.hybrid             #                 unless it's a hybrid
        ce = pe # keep ce, the hybrid edge -> keep its isMajor & gamma.
        pe = net.node[i].edge[j==1 ? 2 : 1] # remove non-hybrid edge
    else
        ce = net.node[i].edge[j==1 ? 2 : 1] # edge to keep
    end
    isnodeiparent = (net.node[i] == ce.node[ce.isChild1 ? 2 : 1])
    (!ce.hybrid || isnodeiparent) ||
      error("node $(net.node[i].number) has 1 tree edge ($(pe.number)) and 1 hybrid edge ($(ce.number)), but is child of the hybrid edge.")
    pn = getOtherNode(pe,net.node[i])
    removeEdge!(net.node[i],ce) # perhaps useless. in case gc() on ith node affects its edges.
    removeNode!(net.node[i],ce)
    removeEdge!(pn,pe)          # perhaps useless. in case gc() on pe affects its nodes.
    removeNode!(pn,pe)
    setEdge!(pn,ce)
    setNode!(ce,pn)             # pn comes 2nd in ce now.
    ce.isChild1 = isnodeiparent # to retain same direction as before.
    ce.length = addBL(ce.length, pe.length)
    deleteNode!(net,net.node[i])
    deleteEdge!(net,pe)
    return nothing
end

#-----------------------------------------------------------------#
#      older functions, mostly unused now.                        #
#                                                                 #
#      they heavily require branch lengths, a level-1 network,    #
#      and do a lot of work to update all associated attributes.  #
#      previously in descriptive.jl                               #
#-----------------------------------------------------------------#
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

#-----------------------------------------------------------------#
#        end of older (unused) functions                          #
#-----------------------------------------------------------------#


# Claudia SL & Paul Bastide: November 2015, Cecile: Feb 2016

#################################################
# Direct Edges
#################################################

"""
`directEdges!(net::HybridNetwork; checkMajor=true::Bool)`

Updates the edges' attribute `isChild1`, according to the root placement.
Also updates edges' attribute `containRoot`, for other possible root placements
compatible with the direction of existing hybrid edges.
Relies on hybrid nodes having exactly 1 major hybrid parent edge,
but checks for that if checkMajor=true.

Warning: Assumes that isChild1 is correct on hybrid edges
(to avoid changing the identity of which nodes are hybrids and which are not).

Return value: `true` if the process was successful,
`false` if the starting root placement (using net.root) was in
conflict with the direction of any hybrid edge.
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
              error("hybrid node $(n.number) has 0 or 2+ major hybrid parents")
            (nparents!=2 || n.hybrid) ||
              warn("node $(n.number) has 2 parents but its hybrid attribute is false.
It is not used in directEdges!, but might cause an error elsewhere.")
            # to fix this: change n.hybrid, net.hybrid, net.numHybrids etc.
            # none of those attributes are used here.
        end
    end
    for(e in net.node[net.root].edge)
        traverseDirectEdges!(net.node[net.root],e,true) ||
         return false # exits early with 'false' if rooting incompatibility
    end
    net.isRooted = true
    return true
end

# containroot = true until the path goes through a hybrid node, below which
# containroot is turned to false.
function traverseDirectEdges!(node::Node, edge::Edge, containroot::Bool)
    if (edge.hybrid && node==edge.node[edge.isChild1 ? 1 : 2])
        warn("the direction of hybrid edge $(edge.number) (isChild1) is incompatible with the root.
       iChild1 and containRoot were updated for a subset of edges in the network only.")
        return false
    end
    if (node == edge.node[1])
        edge.isChild1 = false
        cn = edge.node[2] # cn = child node
    else
        edge.isChild1 = true
        cn = edge.node[1]
    end
    edge.containRoot = containroot
    if (!cn.leaf && (!edge.hybrid || edge.isMajor)) # continue down recursion
        if edge.hybrid containroot=false; end # changes containroot locally, intentional.
        nchildren=0
        for (e in cn.edge)
            if e==edge continue; end
            if (e.hybrid && cn == e.node[e.isChild1 ? 1 : 2]) continue; end
            traverseDirectEdges!(cn,e,containroot) ||
              return false # exits early with 'false' if rooting incompatibility
            nchildren += 1
        end
        if nchildren==0
            warn("non-leaf node $(cn.number) had 0 children. It could be a hybrid
       whose parents' direction is in contradiction with the root.
       iChild1 and containRoot were updated for a subset of edges in the network only.")
            return false
        end
    end
    return true
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
`rotate!(net::HybridNetwork, nodeNumber::Int64; orderedEdgeNum::Array{Int64,1})`

Rotates the order of the node's children edges. Useful for plotting,
to remove crossing edges.
If `node` is a tree node with no polytomy, the 2 children edges are switched
and the optional argument `orderedEdgeNum` is ignored.

Use plot(net, showNodeNumber=true, showEdgeNumber=false) to map node and edge numbers
on the network, as shown in the examples below.

Warning: assumes that edges are correctly directed (isChild1 updated). This is done
by plot(net). Otherwise run directEdges!(net).

# Example
```julia
net = readTopology('(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;');
plot(net, showNodeNumber=true)
rotate!(net, -4)
plot(net)

net=readTopology('(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.136):10.0);');
plot(net, showNodeNumber=true, showEdgeNumber=true)
rotate!(net, -1, orderedEdgeNum=[1,12,9])
plot(net, showNodeNumber=true, showEdgeNumber=true)
rotate!(net, -3)
plot(net)
```
"""
function rotate!(net::HybridNetwork, nnum::Int64; orderedEdgeNum=Int64[]::Array{Int64,1})
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
    elseif length(ci)==2 || length(orderedEdgeNum)==0
        etmp          = n.edge[ci[1]]
        n.edge[ci[1]] = n.edge[ci[2]]
        n.edge[ci[2]] = etmp
    else # 3+ children edges and orderedEdgeNum provided
        length(orderedEdgeNum)==length(ci) || error("orderedEdgeNum $orderedEdgeNum should be of length $(length(ci))")
        length(unique(orderedEdgeNum))==length(ci) || error("orderedEdgeNum $orderedEdgeNum should not have duplicates")
        childrenedge = n.edge[ci]
        for i=1:length(ci)
            tmp = findin([e.number for e in childrenedge], orderedEdgeNum[i])
            length(tmp)==1 || error("edge number $(orderedEdgeNum[i]) not found as child of node $(n.number)")
            n.edge[ci[i]] = childrenedge[tmp[1]]
        end
    end
    return nothing
end
