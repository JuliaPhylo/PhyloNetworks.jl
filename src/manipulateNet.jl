"""
    undirectedOtherNetworks(net::HybridNetwork)

Return a vector of HybridNetwork objects, obtained by switching the placement
of each hybrid node to other nodes inside its cycle. This amounts to changing
the direction of a gene flow event (recursively to move around the whole cycle
of each reticulation).

Optional argument: `outgroup`, as a String. If an outgroup is specified,
then networks conflicting with the placement of the root are avoided.

Assumptions: `net` is assumed to be of level 1, that is, each blob has a
single cycle with a single reticulation.
All level-1 fields of `net` are assumed up-to-date.
# Example
```julia
julia> net = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
julia> vnet = undirectedOtherNetworks(net)
```
"""
function undirectedOtherNetworks(net0::HybridNetwork; outgroup="none"::AbstractString, insideSnaq=false::Bool)
# extra optional argument: "insideSnaq". When true, all level1-attributes are assumed up-to-date
# So far, undirectedOtherNetworks is called inside optTopRuns only
# Potential bug: if new node is -1, then inCycle will become meaningless: changed in readSubTree here
# WARNING: does not update partition, because only thing to change is hybrid node number
    if !insideSnaq
        net0 = readTopologyLevel1(writeTopologyLevel1(net0))
    end
    otherNet = HybridNetwork[]
    for i in 1:net0.numHybrids #need to do for by number, not node
        net = deepcopy(net0) # to avoid redoing attributes after each cycle is finished
        ## undo attributes at current hybrid node:
        hybrid = net.hybrid[i]
        nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,hybrid);
        @debug "nodesInCycle are: $([n.number for n in nodesInCycle])"
        !nocycle || error("the hybrid node $(hybrid.number) does not create a cycle")
        edgesRoot = identifyContainRoot(net,hybrid);
        edges = hybridEdges(hybrid);
        undoGammaz!(hybrid,net);
        othermaj = getOtherNode(edges[1],hybrid)
        edgesmaj = hybridEdges(othermaj)
        if edgesmaj[3].containRoot #if containRoot=true, then we need to undo
            undoContainRoot!(edgesRoot);
        end
        ## changes to new hybrid node:
        for newn in nodesInCycle
            if newn.number != hybrid.number # nodesInCycle contains the hybrid too
                newnet = deepcopy(net)
                newnocycle, newedgesInCycle, newnodesInCycle = identifyInCycle(newnet,newnet.hybrid[i]);
                !newnocycle || error("the hybrid node $(newnet.hybrid[i].number) does not create a cycle")
                ind = getIndexNode(newn.number,newnet) # find the newn node in the new network
                @debug "moving hybrid to node $(newnet.node[ind].number)"
                hybridatnode!(newnet, newnet.hybrid[i], newnet.node[ind])
                @debug begin printEdges(newnet); "printed edges" end
                @debug begin printNodes(newnet); "printed nodes" end
                undoInCycle!(newedgesInCycle, newnodesInCycle);
                @debug begin printEdges(newnet); "printed edges" end
                @debug begin printNodes(newnet); "printed nodes" end
                ##undoPartition!(net,hybrid, edgesInCycle)
                success, hybrid0, flag, nocycle, flag2, flag3 = updateAllNewHybrid!(newnet.node[ind], newnet, false,false,false)
                if success
                    @debug "successfully added new network: $(writeTopologyLevel1(newnet))"
                    push!(otherNet,newnet)
                else
                    println("the network obtained by putting the new hybrid in node $(newnet.node[ind].number) is not good, inCycle,gammaz,containRoot: $([flag,flag2,flag3]), we will skip it")
                end
            end
        end
    end
    # check root in good position
    if outgroup == "none"
        for n in otherNet
            !isTree(n) && checkRootPlace!(n, verbose=false)
        end
        return otherNet
    else ## root already in good place
        @debug "we will remove networks contradicting the outgroup in undirectedOtherNetworks"
        whichKeep = ones(Bool,length(otherNet)) # repeats 'true'
        i = 1
        for n in otherNet
            if !isTree(n)
                try
                    checkRootPlace!(n, verbose=true, outgroup=outgroup)
                catch
                    @debug "found one network incompatible with outgroup"
                    @debug "$(writeTopologyLevel1(n))"
                    whichKeep[i] = false
                end
            end
            i = i+1;
        end
        return otherNet[whichKeep]
    end
end

"""
    hybridatnode!(net::HybridNetwork, nodeNumber::Integer)

Change the status of edges in network `net`,
to move the hybrid node in a cycle to the node with number `nodeNumber`.
This node must be in one (and only one) cycle, otherwise an error will be thrown.


`net` is assumed to be of level 1, that is, each blob has a
single cycle with a single reticulation.
Check and update the nodes' field `inCycle`.

# Example
```julia
julia> net = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
julia> using PhyloPlots
julia> plot(net, showNodeNumber=true)
julia> hybridatnode!(net, -4)
julia> plot(net)
```
"""
function hybridatnode!(net::HybridNetwork, nodeNumber::Integer)
    undoInCycle!(net.edge, net.node)
    for n in net.hybrid
        flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,n);
        flag || error("not level1 network, hybrid $(n.number) cycle intersects another cycle")
        !nocycle || error("strange network without cycle for hybrid $(n.number)")
    end
    ind = 0
    try
        ind = getIndexNode(nodeNumber,net)
    catch
        error("cannot set node $(nodeNumber) as hybrid because it is not part of net")
    end
    net.node[ind].inCycle != -1 || error("node $(nodeNumber) is not part of any cycle, so we cannot make it hybrid")
    indhyb = 0
    try
        indhyb = getIndexNode(net.node[ind].inCycle,net)
    catch
        error("cannot find the hybrid node with number $(net.node[ind].inCycle)")
    end
    hybrid = net.node[indhyb]
    hybridatnode!(net,hybrid,net.node[ind])
    return net
end

"""
    hybridatnode!(net, hybrid::Node, newNode::Node)

Move the reticulation from `hybrid` to `newNode`,
which must in the same cycle. `net` is assumed to be of level 1,
but **no checks** are made and fields are supposed up-to-date.

Called by `hybridatnode!(net, node number)`, which is itself
called by [`undirectedOtherNetworks`](@ref).
"""
function hybridatnode!(net::HybridNetwork, hybrid::Node, newNode::Node)
    hybrid.hybrid || error("node $(hybrid.number) should be hybrid, but it is not")
    hybedges = hybridEdges(hybrid)
    makeEdgeTree!(hybedges[1],hybrid)
    makeEdgeTree!(hybedges[2],hybrid)
    hybedges[1].inCycle = hybrid.number #just to keep attributes ok
    hybedges[2].inCycle = hybrid.number
    switchHybridNode!(net,hybrid,newNode)
    found = false
    for e in newNode.edge
        if e.inCycle == hybrid.number
            if !found
                found = true
                makeEdgeHybrid!(e,newNode, 0.51, switchHyb=true) #first found, major edge, need to optimize gamma anyway
                ##e.gamma = -1
                ##e.containRoot = true ## need attributes like in snaq
            else
                makeEdgeHybrid!(e,newNode, 0.49, switchHyb=true) #second found, minor edge
                ##e.gamma = -1
                ##e.containRoot = true
            end
        end
    end
end

# function to change the hybrid node in a cycle
# does not assume that the network was read with readTopologyUpdate
# does not modify net0 because it needs to update all attributes
# so, it returns the new network
# Not used anywhere, but tested
@doc (@doc hybridatnode!) hybridatnode
function hybridatnode(net0::HybridNetwork, nodeNumber::Integer)
    net = readTopologyLevel1(writeTopologyLevel1(net0)) # we need inCycle attributes
    ind = 0
    try
        ind = getIndexNode(nodeNumber,net)
    catch
        error("cannot set node $(nodeNumber) as hybrid because it is not part of net")
    end
    net.node[ind].inCycle != -1 || error("node $(nodeNumber) is not part of any cycle, so we cannot make it hybrid")
    indhyb = 0
    try
        indhyb = getIndexNode(net.node[ind].inCycle,net)
    catch
        error("cannot find the hybrid node with number $(net.node[ind].inCycle)")
    end
    hybrid = net.node[indhyb]
    hybrid.hybrid || error("node $(hybrid.number) should be hybrid, but it is not")
    hybedges = hybridEdges(hybrid)
    makeEdgeTree!(hybedges[1],hybrid)
    makeEdgeTree!(hybedges[2],hybrid)
    hybedges[1].inCycle = hybrid.number #just to keep attributes ok
    hybedges[2].inCycle = hybrid.number
    switchHybridNode!(net,hybrid,net.node[ind])
    found = false
    for e in net.node[ind].edge
        if e.inCycle == hybrid.number
            if !found
                found = true
                makeEdgeHybrid!(e,net.node[ind], 0.51, switchHyb=true) #first found, major edge, need to optimize gamma anyway
                e.gamma = -1
                e.containRoot = true
            else
                makeEdgeHybrid!(e,net.node[ind], 0.49, switchHyb=true) #second found, minor edge
                e.gamma = -1
                e.containRoot = true
            end
        end
    end
    return net
end


"""
    rootatnode!(HybridNetwork, nodeNumber::Integer; index=false::Bool, verbose=true::Bool)
    rootatnode!(HybridNetwork, Node; verbose=true)
    rootatnode!(HybridNetwork, nodeName::AbstractString; verbose=true)

Root the network/tree object at the node with name 'nodeName' or
number 'nodeNumber' (by default) or with index 'nodeNumber' if index=true.
Attributes isChild1 and containRoot are updated along the way.
Use `plot(net, showNodeNumber=true, showEdgeLength=false)` to
visualize and identify a node of interest.
(see package [PhyloPlots](https://github.com/cecileane/PhyloPlots.jl))

Return the network.

Warnings:

- If the node is a leaf, the root will be placed along
  the edge adjacent to the leaf. This might add a new node.
- If the desired root placement is incompatible with one or more hybrids, then

  * a RootMismatch error is thrown; use `verbose=false` to silence
    the root mismatch info printed before the error is thrown.
  * the input network will still have some attributes modified.

See also: [`rootonedge!`](@ref).
"""
function rootatnode!(net::HybridNetwork, node::Node; kwargs...)
    rootatnode!(net, node.number; kwargs..., index=false)
end

function rootatnode!(net::HybridNetwork, nodeName::AbstractString; kwargs...)
    tmp = findall(n -> n.name == nodeName, net.node)
    if length(tmp)==0
        error("node named $nodeName was not found in the network.")
    elseif length(tmp)>1
        error("several nodes were found with name $nodeName.")
    end
    rootatnode!(net, tmp[1]; kwargs..., index=true)
end

function rootatnode!(net::HybridNetwork, nodeNumber::Integer; index=false::Bool, verbose=true::Bool)
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
        # @info "node $(net.node[ind].number) is a leaf. Will create a new node if needed, to set taxon \"$(net.node[ind].name)\" as outgroup."
        length(net.node[ind].edge)==1 || error("leaf has $(length(net.node[ind].edge)) edges!")
        pn = getOtherNode(net.node[ind].edge[1], net.node[ind])
        if length(pn.edge) <= 2 # if parent of leaf has degree 2, use it as new root
            rootatnode!(net, pn.number; verbose=verbose)
        else # otherwise, create a new node between leaf and its parent
            rootonedge!(net,net.node[ind].edge[1]; verbose=verbose)
        end
    else
        rootsaved = net.root
        net.root = ind
        try
          directEdges!(net)
        catch e
          if isa(e, RootMismatch) # new root incompatible with hybrid directions: revert back
            verbose && println("RootMismatch: reverting to old root position.")
            net.root = rootsaved
          end
          rethrow(e)
        end
        if (net.root != rootsaved && length(net.node[rootsaved].edge)==2)
            fuseedgesat!(rootsaved,net) # remove old root node if degree 2
        end
        return net
    end
end


"""
    rootonedge!(HybridNetwork, edgeNumber::Integer; index=false::Bool, verbose=true::Bool)
    rootonedge!(HybridNetwork, Edge; verbose=true::Bool)

Root the network/tree along an edge with number `edgeNumber` (by default)
or with index `edgeNumber` if `index=true`.
Attributes `isChild1` and `containRoot` are updated along the way.

This adds a new node and a new edge to the network.
Use `plot(net, showEdgeNumber=true, showEdgeLength=false)` to
visualize and identify an edge of interest.
(see package [PhyloPlots](https://github.com/cecileane/PhyloPlots.jl))

See also: [`rootatnode!`](@ref).
"""
function rootonedge!(net::HybridNetwork, edge::Edge; kwargs...)
    rootonedge!(net, edge.number, index=false; kwargs...)
end

function rootonedge!(net::HybridNetwork, edgeNumber::Integer; index=false::Bool, verbose=true::Bool)
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
    try
      directEdges!(net)
    catch e
      if isa(e, RootMismatch) # new root incompatible with hybrid directions: revert back
        verbose && println("RootMismatch: reverting to old root position.")
        fuseedgesat!(net.root,net) # reverts breakedge!
        net.root = rootsaved
      end
      rethrow(e)
    end
    if (net.root != rootsaved && length(net.node[rootsaved].edge)==2)
        fuseedgesat!(rootsaved,net) # remove old root node if degree 2
    end
    return net
end

"""
    breakedge!(Edge, HybridNetwork)

breaks an edge into 2 edges (each of length half that of original edge).
creates new node of degree 2. Useful to root network along an edge.

warning: updates `isChild1` and `containRoot`, but
does NOT update attributes like: inCycle, partition, gammaz, etc.

returns the index of the newly created node in the network
"""
function breakedge!(edge::Edge, net::HybridNetwork)
    pn = getParent(edge) # parent node
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
    fuseedgesat!(i::Integer,net::HybridNetwork)

Removes `i`th node in net.node, if it is of degree 2.
The parent and child edges of this node are fused.
Reverts the action of breakedge!.

returns the fused edge.
"""
function fuseedgesat!(i::Integer, net::HybridNetwork)
    i <= length(net.node) ||
      error("node index $i too large: only $(length(net.node)) nodes in the network.")
    length(net.node[i].edge) == 2 ||
      error("can't fuse edges at node number $(net.node[i].number): connected to $(length(net.node[i].edge)) edges.")
    !(net.node[i].edge[1].hybrid && net.node[i].edge[2].hybrid) ||
      error("can't fuse edges at node number $(net.node[i].number): connected to exactly 2 hybrid edges")
    j = argmax([e.number for e in net.node[i].edge])
    pe = net.node[i].edge[j] # edge to remove: pe.number > ce.number
    if pe.hybrid             #                 unless it's a hybrid
        ce = pe # keep ce, the hybrid edge -> keep its isMajor & gamma.
        pe = net.node[i].edge[j==1 ? 2 : 1] # remove non-hybrid edge
    else
        ce = net.node[i].edge[j==1 ? 2 : 1] # edge to keep
    end
    isnodeiparent = (net.node[i] ≡ getParent(ce))
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
    if net.root==i # if `node` was the root, new root = pn
        net.root = getIndex(pn,net)
    end
    deleteNode!(net,net.node[i])
    deleteEdge!(net,pe,part=false) # do not update partitions. irrelevant for networks of level>1.
    return ce
end

# Claudia SL & Paul Bastide: November 2015, Cecile: Feb 2016

#################################################
# Direct Edges
#################################################

"""
    directEdges!(net::HybridNetwork; checkMajor=true::Bool)

Updates the edges' attribute `isChild1`, according to the root placement.
Also updates edges' attribute `containRoot`, for other possible root placements
compatible with the direction of existing hybrid edges.
Relies on hybrid nodes having exactly 1 major hybrid parent edge,
but checks for that if checkMajor=true.

Warning: Assumes that isChild1 is correct on hybrid edges
(to avoid changing the identity of which nodes are hybrids and which are not).

Returns the network. Throws a 'RootMismatch' Exception if the root was found to
conflict with the direction of any hybrid edge.
"""
function directEdges!(net::HybridNetwork; checkMajor=true::Bool)
    if checkMajor # check each node has 2+ hybrid parent edges (if any), and exactly one major.
        for n in net.node
            nparents = 0 # 0 or 2 normally, but could be >2 if polytomy.
            nmajor = 0   # there should be exactly 1 major parent if nparents>0
            for e in n.edge
                if e.hybrid && n == getChild(e)
                    nparents += 1
                    if (e.isMajor) nmajor +=1; end
                end
            end
            (nparents!=1) || error("node $(n.number) has exactly 1 hybrid parent edge")
            (nparents==0 || nmajor == 1) ||
              error("hybrid node $(n.number) has 0 or 2+ major hybrid parents")
            (nparents!=2 || n.hybrid) ||
              @warn "node $(n.number) has 2 parents but its hybrid attribute is false.
It is not used in directEdges!, but might cause an error elsewhere."
            # to fix this: change n.hybrid, net.hybrid, net.numHybrids etc.
            # none of those attributes are used here.
        end
    end
    net.cleaned = false # attributed used by snaq! Will change isChild1 and containRoot
    for e in net.node[net.root].edge
        traverseDirectEdges!(net.node[net.root],e,true)
    end
    net.isRooted = true
    return net
end

# containroot = true until the path goes through a hybrid node, below which
# containroot is turned to false.
function traverseDirectEdges!(node::Node, edge::Edge, containroot::Bool)
    if edge.hybrid && node==getChild(edge)
        throw(RootMismatch(
"direction (isChild1) of hybrid edge $(edge.number) conflicts with the root.
isChild1 and containRoot were updated for a subset of edges in the network only."))
    end
    if node == edge.node[1]
        edge.isChild1 = false
        cn = edge.node[2] # cn = child node
    else
        edge.isChild1 = true
        cn = edge.node[1]
    end
    edge.containRoot = containroot
    if !cn.leaf && (!edge.hybrid || edge.isMajor) # continue down recursion
        if edge.hybrid containroot=false; end # changes containroot locally, intentional.
        nchildren=0
        for e in cn.edge
            if e==edge continue; end
            if (e.hybrid && cn == getChild(e)) continue; end
            traverseDirectEdges!(cn,e,containroot)
            nchildren += 1
        end
        if nchildren==0
            throw(RootMismatch("non-leaf node $(cn.number) had 0 children.
Could be a hybrid whose parents' direction conflicts with the root.
isChild1 and containRoot were updated for a subset of edges in the network only."))
        end
    end
    return nothing
end

#################################################
## Topological sorting
#################################################

"""
    getParents(node)

Get vector of all parent nodes of `n`, based on `isChild1` field (for edges).
To get the parent node of an edge: see [`getParent`](@ref).  
To get individual parent edges (rather than all parent *nodes*):
see [`getMajorParentEdge`](@ref) and `getMinorParentEdge`.
"""
@inline function getParents(node::Node)
    parents = Node[]
    for e in node.edge
            if node == getChild(e)
                push!(parents, getParent(e))
            end
    end
    return parents
end

# getParent, getMajorParent, getMinorParent: defined in auxiliary.jl

"""
    getMajorParentEdge(node)
    getMinorParentEdge(node)

return the parent edge of a given node: the major / minor if hybrid.  
**warning**: assume isChild1 and isMajor attributes are correct

To get all parent *nodes*: see [`getParents`](@ref).  
"""
@inline function getMajorParentEdge(n::Node)
    for ee in n.edge
        if n == ee.node[(ee.isChild1 ? 1 : 2)] && ee.isMajor
            return ee
        end
    end
    error("node $(n.number) has no major parent")
end
@doc (@doc getMajorParentEdge) getMinorParentEdge
@inline function getMinorParentEdge(n::Node)
    for ee in n.edge
        if !ee.isMajor && n == ee.node[(ee.isChild1 ? 1 : 2)]
            return ee
        end
    end
    error("node $(n.number) has no minor parent")
end

"""
    getChildren(node)

return a vector with all children *nodes* of `node`.  
**warning**: assume `isChild1` field (for edges) are correct

To get all parent *nodes*: see [`getParents`](@ref).  
"""
function getChildren(node::Node)
    children = Node[]
    for e in node.edge
        if node == getParent(e)
            push!(children, getChild(e))
        end
    end
    return children
end

"""
    preorder!(net::HybridNetwork)

Updates attribute net.nodes_changed in which the nodes are pre-ordered
(also called topological sorting), such that each node is visited after its parent(s).
The edges' direction needs to be correct before calling preorder!, using directEdges!
"""
function preorder!(net::HybridNetwork)
    net.isRooted || error("net needs to be rooted for preorder!, run root functions or directEdges!")
    net.nodes_changed = Node[] # path of nodes in preorder.
    queue = Node[] # problem with PriorityQueue(): dequeue() takes a
                   # random member if all have the same priority 1.
    net.visited = [false for i = 1:size(net.node,1)];
    push!(queue,net.node[net.root]) # push root into queue
    while !isempty(queue)
        #println("at this moment, queue is $([n.number for n in queue])")
        curr = pop!(queue); # deliberate choice over shift! for cladewise order
        currind = findfirst(x -> x===curr, net.node)
        # the "curr"ent node may have been already visited: because simple loop (2-cycle)
        !net.visited[currind] || continue
        net.visited[currind] = true # visit curr node
        push!(net.nodes_changed,curr) #push curr into path
        for e in curr.edge
            if curr == getParent(e)
                other = getChild(e)
                if !e.hybrid
                    push!(queue,other)
                    # print("queuing: "); @show other.number
                else
                    e2 = getPartner(e, other)
                    parent = getParent(e2)
                    if net.visited[findfirst(x -> x===parent, net.node)]
                      push!(queue,other)
                      # warning: if simple loop, the same node will be pushed twice: child of "curr" via 2 edges
                    end
                end
            end
        end
    end
    # println("path of nodes is $([n.number for n in net.nodes_changed])")
end


"""
    cladewiseorder!(net::HybridNetwork)

Updates attribute net.cladewiseorder_nodeIndex. Used for plotting the network.
In the major tree, all nodes in a given clade are consecutive. On a tree, this function
also provides a pre-ordering of the nodes.
The edges' direction needs to be correct before calling
[`cladewiseorder!`](@ref), using [`directEdges!`](@ref)
"""
function cladewiseorder!(net::HybridNetwork)
    net.isRooted || error("net needs to be rooted for cladewiseorder!\n run root functions or directEdges!")
    net.cladewiseorder_nodeIndex = Int[]
    queue = [net.root] # index (in net) of nodes in the queue
    # print("queued the root's children's indices: "); @show queue
    while !isempty(queue)
        ni = pop!(queue); # deliberate choice over shift! for cladewise order
        # @show net.node[ni].number
        push!(net.cladewiseorder_nodeIndex, ni)
        for e in net.node[ni].edge
            if net.node[ni] ≡ getParent(e) # net.node[ni] is parent node of e
                if e.isMajor
                    push!(queue, findfirst(isequal(getChild(e)), net.node))
                    # print("queuing: "); @show other.number
                end
            end
        end
    end
end

"""
    rotate!(net::HybridNetwork, nodeNumber::Integer; orderedEdgeNum::Array{Int,1})

Rotates the order of the node's children edges. Useful for plotting,
to remove crossing edges.
If `node` is a tree node with no polytomy, the 2 children edges are switched
and the optional argument `orderedEdgeNum` is ignored.

Use `plot(net, showNodeNumber=true, showEdgeNumber=false)` to map node and edge numbers
on the network, as shown in the examples below.
(see package [PhyloPlots](https://github.com/cecileane/PhyloPlots.jl))

Warning: assumes that edges are correctly directed (isChild1 updated). This is done
by `plot(net)`. Otherwise run `directEdges!(net)`.

# Example

```julia
julia> net = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
julia> using PhyloPlots
julia> plot(net, showNodeNumber=true)
julia> rotate!(net, -4)
julia> plot(net)
julia> net=readTopology("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.136):10.0);");
julia> plot(net, showNodeNumber=true, showEdgeNumber=true)
julia> rotate!(net, -1, orderedEdgeNum=[1,12,9])
julia> plot(net, showNodeNumber=true, showEdgeNumber=true)
julia> rotate!(net, -3)
julia> plot(net)
```
"""
function rotate!(net::HybridNetwork, nnum::Integer; orderedEdgeNum=Int[]::Array{Int,1})
    nind = 0
    nind = findfirst(n -> n.number == nnum, net.node)
    nind !== nothing || error("cannot find any node with number $nnum in network.")
    n = net.node[nind]
    ci = Int[] # children edge indices
    for i = 1:length(n.edge)
        if n == getParent(n.edge[i])
            push!(ci,i)
        end
    end
    if length(ci) < 2
        @warn "no edge to rotate: node $nnum has $(length(ci)) children edge."
    elseif length(ci)==2 || length(orderedEdgeNum)==0
        etmp          = n.edge[ci[1]]
        n.edge[ci[1]] = n.edge[ci[2]]
        n.edge[ci[2]] = etmp
    else # 3+ children edges and orderedEdgeNum provided
        length(orderedEdgeNum)==length(ci) || error("orderedEdgeNum $orderedEdgeNum should be of length $(length(ci))")
        length(unique(orderedEdgeNum))==length(ci) || error("orderedEdgeNum $orderedEdgeNum should not have duplicates")
        childrenedge = n.edge[ci] # makes a shallow copy, because of subsetting [ci]
        for i=1:length(ci)
            tmp = findall(x -> x.number == orderedEdgeNum[i], childrenedge)
            length(tmp)==1 || error("edge number $(orderedEdgeNum[i]) not found as child of node $(n.number)")
            n.edge[ci[i]] = childrenedge[tmp[1]]
        end
    end
    return nothing
end


### WARNING:
# deleteleaf! is similar but also very different from
# deleteLeaf! in pseudolik.jl, which
# - does not necessarily remove nodes of degree 2,
# - requires and updates all attributes for level-1 networks:
#   inCycle, partition, branch lengths, diamond/triangle types etc.
# - is used a lot within snaq! to extract quartets and retains info
#   on which parameters in the full network affect the quartet.
# deleteIntLeaf! somewhat similar to fuseedgesat!
"""
    deleteleaf!(HybridNetwork, leafName::AbstractString; simplify=true)
    deleteleaf!(HybridNetwork, Node; simplify=true)
    deleteleaf!(HybridNetwork, Integer; index=false, simplify=true)

Deletes a leaf node from the network, possibly from its name, number, or index
in the network's array of nodes.

simplify: if true and if deleting the node results in 2 hybrid edges
forming a cycle of k=2 nodes, then these hybrid edges are merged and
simplified as a single tree edge.

The first 2 versions require that `node` is a leaf.
The 3rd version does **not** require that `node` is a leaf.
If `node` has degree 3 or more, nothing happens. If it has degree 1 or 2, it is deleted.

Warning: does **not** update attributes related to level-1 networks,
such as inCycle, partition, gammaz, etc.
Does not require branch lengths, and designed to work on networks
of all levels.
"""
function deleteleaf!(net::HybridNetwork, node::Node; simplify=true::Bool)
    node.leaf || error("node number $(node.number) is not a leaf.")
    deleteleaf!(net, node.number, index=false, simplify=simplify)
end

function deleteleaf!(net::HybridNetwork, nodeName::AbstractString; simplify=true::Bool)
    tmp = findall(n -> n.name == nodeName, net.node)
    if length(tmp)==0
        error("node named $nodeName was not found in the network.")
    elseif length(tmp)>1
        error("several nodes were found with name $nodeName.")
    end
    deleteleaf!(net, tmp[1], index=true, simplify=simplify)
end

# recursive algorithm. nodes previously removed are all necessaily
# *younger* than current node to remove ("leaf"). Stated otherwise:
# edges previously removed all go "down" in time towards current node:
# - tree edge down to an original leaf,
# - 2 hybrid edges down to a hybrid node.
# hybrid edges from node to another node are not removed. fused instead.
# consequence: node having 2 hybrid edges away from node should not occur.
function deleteleaf!(net::HybridNetwork, nodeNumber::Integer;
                     index=false::Bool, simplify=true::Bool)
    i = nodeNumber # good if index=true
    if !index
      try
        i = getIndexNode(nodeNumber,net)
      catch
        error("cannot delete leaf number $(nodeNumber) because it is not part of net")
      end
    elseif i > length(net.node)
        error("node index $i too large: the network only has $(length(net.node)) nodes.")
    end
    if length(net.node[i].edge)==0
        length(net.node)==1 || error("leaf $(net.node[i].name) has no edge but network has $(length(net.node)) nodes (instead of 1).")
        println("Only 1 node. Removing it: the network will be empty")
        deleteNode!(net,net.node[i]) # empties the network
    elseif length(net.node[i].edge)==1
        pe = net.node[i].edge[1]
        pn = getOtherNode(pe, net.node[i]) # parent node of leaf
        # remove leaf and pe.
        removeNode!(pn,pe)  # perhaps useless. in case gc() on pe affects pn
        removeEdge!(pn,pe)
        deleteEdge!(net,pe,part=false)
        if net.root==i # if node was the root, new root = pn
            net.root = getIndex(pn,net)
        end
        deleteNode!(net,net.node[i])
        if pn.leaf # network had 2 nodes only: pn and the leaf
            length(net.edge)==0 || error("neighbor of leaf $(net.node[i].name) is another leaf, but network had $(length(net.edge)) edges (instead of 1).")
            length(pn.edge)==0 || error("neighbor of leaf $(net.node[i].name) is another leaf, which had $(length(pn.edge)) edges (instead of 1)")
            return nothing # all done: exit function
        end
        deleteleaf!(net, pn.number, simplify=simplify)
    elseif length(net.node[i].edge)==2
        e1 = net.node[i].edge[1]
        e2 = net.node[i].edge[2]
        if (e1.hybrid && e2.hybrid)
            (net.node[i] ≡ getChild(e1) && net.node[i] ≡ getChild(e2)) ||
              error("after removing leaf and descendants, node $(net.node[i].number) has 2 hybrid edges but is not the child of both.")
            p1 = getParent(e1) # find both parents of hybrid leaf
            p2 = getParent(e2)
            # remove the hybrid `node` and both e1, e2
            sameparent = (p1≡p2) # should not happen unless original network has k=2 cycles
            removeNode!(p1,e1);  removeNode!(p2,e2) # perhaps useless
            removeEdge!(p1,e1);  removeEdge!(p2,e2)
            deleteEdge!(net,e1,part=false); deleteEdge!(net,e2,part=false)
            if net.root==i net.root=getIndex(p1,net); end # should never occur though.
            deleteNode!(net,net.node[i])
            # recursive call on both p1 and p2.
            deleteleaf!(net, p1.number, simplify=simplify)
            sameparent || deleteleaf!(net, p2.number, simplify=simplify)
        else
            e1 = fuseedgesat!(i,net) # fused edge
            if simplify && e1.hybrid # check for cycle of k=2 nodes
                cn = getChild(e1)
                e2 = nothing
                for e in cn.edge # find companion hybrid edge
                    if e.hybrid && e ≢ e1 && cn ≡ getChild(e)
                        e2=e; break;
                    end
                end
                e2!=nothing || error("node $(cn.number) with a single parent hybrid edge")
                pn  = getParent(e1)
                if pn ≡ getParent(e2)
                    # e1 and e2 have same child and same parent. Remove e1.
                    e2.hybrid=false;
                    e2.isMajor=true; e2.gamma += e1.gamma
                    removeEdge!(pn,e1); removeEdge!(cn,e1)
                    deleteEdge!(net,e1,part=false)
                    # call recursion again because pn and/or cn might be of degree 2.
                    deleteleaf!(net, cn.number, simplify=simplify)
                end
            end
        end
    #else: do nothing. leaf has degree 3+. could occur through recursive calling reaching a polytomy.
    end
    return nothing
end

"""
    resetNodeNumbers!(net::HybridNetwork; checkPreorder=true, ape=true)

Change internal node numbers of `net` to consecutive numbers from 1 to the total
number of nodes.

keyword arguments:
- `ape`: if true, the new numbers satisfy the conditions assumed by the
  `ape` R package: leaves are 1 to n, the root is n+1, and internal nodes
  are higher consecutive integers. If false, nodes are numbered in post-order,
  with leaves from 1 to n (and the root last).
- `checkPreorder`: if false, the `isChild1` edge field and the `net.nodes_changed`
  network field are supposed to be correct (to get nodes in preorder)

# Examples

```jldoctest
julia> net = readTopology("(A,(B,(C,D)));");

julia> PhyloNetworks.resetNodeNumbers!(net)

julia> printNodes(net) # first column "node": root is 5
node leaf  hybrid hasHybEdge name inCycle edges'numbers
1    true  false  false      A    -1      1   
2    true  false  false      B    -1      2   
3    true  false  false      C    -1      3   
4    true  false  false      D    -1      4   
7    false false  false           -1      3    4    5   
6    false false  false           -1      2    5    6   
5    false false  false           -1      1    6   

julia> net = readTopology("(A,(B,(C,D)));");

julia> PhyloNetworks.resetNodeNumbers!(net; ape=false)

julia> printNodes(net) # first column "node": root is 7
node leaf  hybrid hasHybEdge name inCycle edges'numbers
1    true  false  false      A    -1      1   
2    true  false  false      B    -1      2   
3    true  false  false      C    -1      3   
4    true  false  false      D    -1      4   
5    false false  false           -1      3    4    5   
6    false false  false           -1      2    5    6   
7    false false  false           -1      1    6   
```
"""
function resetNodeNumbers!(net::HybridNetwork; checkPreorder=true::Bool, ape=true::Bool)
    if checkPreorder
      directEdges!(net)
      preorder!(net) # to create/update net.nodes_changed
    end
    lnum = 1 # first number
    for n in net.node
        n.leaf || continue
        n.number = lnum
        lnum += 1
    end
    if ape
        nodelist = net.nodes_changed # pre-order: root first
    else
        nodelist = reverse(net.nodes_changed) # post-order
    end
    for n in nodelist
        !n.leaf || continue
        n.number = lnum
        lnum += 1
    end
end

"""
    resetEdgeNumbers!(net::HybridNetwork)

Check that edge numbers of `net` are consecutive numbers from 1 to the total
number of edges. If not, reset the edge numbers to be so.
"""
function resetEdgeNumbers!(net::HybridNetwork)
    enum = [e.number for e in net.edge]
    ne = length(enum)
    unused = setdiff(1:ne, enum)
    if isempty(unused)
        return nothing # all good
    end
    @warn "resetting edge numbers to be from 1 to $ne"
    ind2change = findall(x -> x ∉ 1:ne, enum)
    length(ind2change) == length(unused) || error("can't reset edge numbers")
    for i in 1:length(unused)
        net.edge[ind2change[i]].number = unused[i]
    end
    return nothing;
end

"""
    norootbelow!(e::Edge)

Set `containRoot` to `false` for edge `e` and all edges below, recursively.
The traversal stops if `e.containRoot` is already `false`, assuming
that `containRoot` is already false all the way below down that edge.
"""
function norootbelow!(e::Edge)
    e.containRoot || return nothing # if already false: stop
    # if true: turn to false then move down to e's children
    e.containRoot = false
    cn = getChild(e) # cn = child node
    for ce in cn.edge
        ce !== e || continue # skip e
        getParent(ce) === cn || continue # skip edges that aren't children of cn
        norootbelow!(ce)
    end
    return nothing
end

"""
    allowrootbelow!(e::Edge)
    allowrootbelow!(n::Node, parent_edge_of_n::Edge)

Set `containRoot` to `true` for edge `e` and all edges below, recursively.
The traversal stops whenever a hybrid node is encountered:
if the child of `e` is a hybrid node (that is, if `e` is a hybrid edge)
or if `n` is a hybrid node, then the edges below `e` or `n` are *not* traversed.
"""
function allowrootbelow!(e::Edge)
    e.containRoot = true
    e.hybrid && return nothing # e hybrid edge <=> its child hybrid node: stop
    allowrootbelow!(getChild(e), e)
end
function allowrootbelow!(n::Node, pe::Edge)
    # pe assumed to be the parent of n
    for ce in n.edge
        ce !== pe || continue # skip parent edge of n
        getParent(ce) === n || continue # skip edges that aren't children of n
        allowrootbelow!(ce)
    end
    return nothing
end
