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

# example

```julia
net = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
using PhyloPlots
plot(net, :R, showNodeNumber=true); # to locate nodes and their numbers. D of hybrid origin
hybridatnode!(net, -4)
plot(net, :R, showNodeNumber=true); # hybrid direction reversed: now 2B of hybrid origin
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
    hybridatnode!(net::HybridNetwork, hybrid::Node, newNode::Node)

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
    breakedge!(net.edge[ind],net) # returns new node, new edge (last ones pushed)
    net.root = length(net.node)   # index of new node: was the last one pushed
    try
      directEdges!(net)
    catch e
      if isa(e, RootMismatch) # new root incompatible with hybrid directions: revert back
        verbose && println("RootMismatch: reverting to old root position.")
        fuseedgesat!(net.root,net) # reverts breakedge!
        net.root = rootsaved
        directEdges!(net)
      end
      rethrow(e)
    end
    if (net.root != rootsaved && length(net.node[rootsaved].edge)==2)
        fuseedgesat!(rootsaved,net) # remove old root node if degree 2
    end
    return net
end

"""
    breakedge!(edge::Edge, net::HybridNetwork)

Break an edge into 2 edges, each of length half that of original edge,
creating a new node of degree 2. Useful to root network along an edge.
Return the new node and the new edge, which is the "top" half of the
original starting edge.
These new node & edge are pushed last in `net.node` and `net.edge`.

If the starting edge was:
```
n1  --edge-->  n2
```
then we get this:
```
n1  --newedge-->  newnode  --edge-->  n2
```

`isChild1` and `containRoot` are updated, but not fields for level-1
networks like `inCycle`, `partition`, `gammaz`, etc.

# examples

```jldoctest
julia> net = readTopology("(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));");

julia> length(net.node)
19

julia> net.edge[4] # edge 4 goes from node -8 to 3
PhyloNetworks.Edge:
 number:4
 length:-1.0
 attached to 2 node(s) (parent first): -8 3


julia> newnode, newedge = PhyloNetworks.breakedge!(net.edge[4], net);

julia> length(net.node) # one more than before
20

julia> newedge # new edge 21 goes from node -8 and 11 (new)
PhyloNetworks.Edge:
 number:21
 length:-1.0
 attached to 2 node(s) (parent first): -8 11


julia> net.edge[4] # original edge 4 now goes from node 11 (new) to 3
PhyloNetworks.Edge:
 number:4
 length:-1.0
 attached to 2 node(s) (parent first): 11 3


julia> writeTopology(net) # note extra pair of parentheses around S1
"(((S8,S9),((((S4,(S1)),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
```
"""
function breakedge!(edge::Edge, net::HybridNetwork)
    pn = getParent(edge) # parent node
    # new child edge = old edge, same hybrid attribute
    removeEdge!(pn,edge)
    removeNode!(pn,edge)
    max_edge = maximum(e.number for e in net.edge)
    max_node = maximum(n.number for n in net.node)
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
    return newnode, newedge
end

"""
    fuseedgesat!(i::Integer,net::HybridNetwork, multgammas=false::Bool)

Removes `i`th node in net.node, if it is of degree 2.
The parent and child edges of this node are fused.
Reverts the action of breakedge!.

returns the fused edge.
"""
function fuseedgesat!(i::Integer, net::HybridNetwork, multgammas=false::Bool)
    i <= length(net.node) ||
      error("node index $i too large: only $(length(net.node)) nodes in the network.")
    nodei = net.node[i]
    length(nodei.edge) == 2 ||
      error("can't fuse edges at node number $(nodei.number): connected to $(length(nodei.edge)) edges.")
    !(nodei.edge[1].hybrid && nodei.edge[2].hybrid) ||
      error("can't fuse edges at node number $(nodei.number): connected to exactly 2 hybrid edges")
    j = argmax([e.number for e in nodei.edge])
    pe = nodei.edge[j] # edge to remove: pe.number > ce.number
    ce = nodei.edge[j==1 ? 2 : 1]
    if pe.hybrid       # unless it's a hybrid: should be --tree--> node i --hybrid-->
        (ce,pe) = (pe,ce) # keep the hybrid edge: keep its isMajor
    end
    isnodeiparent = (nodei ≡ getParent(ce))
    (!ce.hybrid || isnodeiparent) ||
      error("node $(nodei.number) has 1 tree edge ($(pe.number)) and 1 hybrid edge ($(ce.number)), but is child of the hybrid edge.")
    pn = getOtherNode(pe,nodei)
    removeEdge!(nodei,ce) # perhaps useless. in case gc() on ith node affects its edges.
    removeNode!(nodei,ce)
    removeEdge!(pn,pe)
    removeNode!(pn,pe)    # perhaps useless. in case gc() on pe affects its nodes.
    setEdge!(pn,ce)
    setNode!(ce,pn)       # pn comes 2nd in ce now: ce.node is: [original, pn]
    ce.isChild1 = isnodeiparent # to retain same direction as before.
    ce.length = addBL(ce.length, pe.length)
    if multgammas
        ce.gamma = multiplygammas(ce.gamma, pe.gamma)
    end
    if net.root==i # isnodeiparent should be true, unless the root and ce's direction were not in sync
        newroot = pn
        if newroot.leaf && !ce.hybrid # then reverse ce's direction. pn.leaf and ce.hybrid should never both occur!
            newroot = ce.node[1] # getOtherNode(ce, pn)
            ce.isChild1 = false
        end
        net.root = findfirst(isequal(newroot), net.node)
    end
    deleteNode!(net,nodei)
    deleteEdge!(net,pe,part=false) # do not update partitions. irrelevant for networks of level>1.
    return ce
end

"""
    removedegree2nodes!(net::HybridNetwork)

Delete *all* nodes of degree two in `net`, fusing the two adjacent edges
together each time, and return the network.
If the network has a degree-2 root, then the root is eliminated as well,
leaving the network unrooted.

See [`fuseedgesat!`](@ref).

```jldoctest
julia> net = readTopology("(((((S1,(S2)#H1),(#H1,S3)))#H2),(#H2,S4));");

julia> PhyloNetworks.breakedge!(net.edge[3], net); # create a degree-2 node along hybrid edge

julia> PhyloNetworks.breakedge!(net.edge[3], net); # another one: 2 in a row

julia> PhyloNetworks.breakedge!(net.edge[10], net); # another one, elsewhere

julia> writeTopology(net) # extra pairs of parentheses
"((#H2,S4),(((((S1,(((S2)#H1))),(#H1,S3)))#H2)));"

julia> PhyloNetworks.removedegree2nodes!(net);

julia> writeTopology(net) # even the root is gone
"(#H2,S4,(((S1,(S2)#H1),(#H1,S3)))#H2);"

```
"""
function removedegree2nodes!(net::HybridNetwork)
    ndegree2nodes = sum(length(n.edge) == 2 for n in net.node)
    # caution: nodes and their indices in the 'current' network may change some of them are removed
    for ni in 1:ndegree2nodes # empty if 0 degree-2 nodes
        i = findfirst(n -> length(n.edge) == 2, net.node)
        i !== nothing || error("incorrect predicted number of degree-2 nodes to remove...")
        fuseedgesat!(i, net)
    end
    return net
end

"""
    addleaf!(net::HybridNetwork, node::Node, leafname::String, edgelength::Float64=-1.0)
    addleaf!(net::HybridNetwork, edge::Edge, leafname::String, edgelength::Float64=-1.0)

Add a new external edge between `node` or between the "middle" of `edge`
and a newly-created leaf, of name `leafname`.
By default, the new edge length is missing (-1).

# examples

```jldoctest
julia> net = readTopology("((S1,(((S2,(S3)#H1),(#H1,S4)))#H2),(#H2,S5));");

julia> net.node[6].name # leaf S4
"S4"

julia> PhyloNetworks.addleaf!(net, net.node[6], "4a"); # adding leaf to a node

julia> writeTopology(net, internallabel=true)
"((S1,(((S2,(S3)#H1),(#H1,(4a)S4)))#H2),(#H2,S5));"

julia> PhyloNetworks.addleaf!(net, net.node[6], "4b");

julia> writeTopology(net, internallabel=true)
"((S1,(((S2,(S3)#H1),(#H1,(4a,4b)S4)))#H2),(#H2,S5));"
```

```jldoctest
julia> net = readTopology("((S1,(((S2,(S3)#H1),(#H1,S4)))#H2),(#H2,S5));");

julia> [n.name for n in net.edge[7].node] # external edge to S4
2-element Array{String,1}:
 "S4"
 ""  

julia> PhyloNetworks.addleaf!(net, net.edge[7], "4a"); # adding leaf to an edge

julia> writeTopology(net, internallabel=true)
"((S1,(((S2,(S3)#H1),(#H1,(S4,4a))))#H2),(#H2,S5));"
```
"""
function addleaf!(net::HybridNetwork, speciesnode::Node, leafname::String, edgelength::Float64=-1.0)
    exterioredge = Edge(maximum(e.number for e in net.edge) + 1, edgelength) # isChild1 = true by default in edge creation
    pushEdge!(net, exterioredge)
    setEdge!(speciesnode, exterioredge)
    if speciesnode.hybrid || (speciesnode != net.node[net.root] && !getMajorParentEdge(speciesnode).containRoot)
        exterioredge.containRoot = false
    end
    newleaf = Node(maximum(n.number for n in net.node) + 1, true, false, [exterioredge]) # Node(number, leaf, hybrid, edge array)
    newleaf.name = leafname
    setNode!(exterioredge, [newleaf, speciesnode]) # [child, parent] to match isChild1 = true by default
    if speciesnode.leaf
        deleteat!(net.leaf,findfirst(isequal(speciesnode), net.leaf))
        speciesnode.leaf = false
        net.numTaxa -= 1
    end
    pushNode!(net, newleaf) # push node into network (see auxillary.jl)
    return net
end

function addleaf!(net::HybridNetwork, startingedge::Edge, leafname::String, edgelength::Float64=-1.0)
    newnode, newedge = breakedge!(startingedge, net)
    addleaf!(net, newnode, leafname, edgelength)
    return net
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

Warnings:
1. Assumes that isChild1 is correct on hybrid edges
(to avoid changing the identity of which nodes are hybrids and which are not).
2. Does not check for cycles (to maintain a network's DAG status)

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
        if node === getParent(e)
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

Note that `LinearAlgebra` also exports a function named `rotate!` in Julia v1.5.
If both packages need to be used in Julia v1.5 or higher,
usage of `rotate!` needs to be qualified, such as with `PhyloNetworks.rotate!`.
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
    deleteleaf!(HybridNetwork, leafName::AbstractString; ...)
    deleteleaf!(HybridNetwork, Node; ...)
    deleteleaf!(HybridNetwork, Integer; index=false, ...)

Delete a node from the network, possibly from its name, number, or index
in the network's array of nodes.
The first two versions require that the node is a leaf.
The third version does **not** require that the node is a leaf:
If it has degree 3 or more, nothing happens.
If it has degree 1 or 2, then it is deleted.

## keyword arguments

`simplify`: if true and if deleting the node results in 2 hybrid edges
forming a cycle of k=2 nodes, then these hybrid edges are merged and
simplified as a single tree edge.

`unroot`: if true, a root of degree 1 or 2 is deleted. If false,
the root is deleted if it is of degree 1 (no root edge is left),
but is kept if it is of degree 2. Deleting all leaves in an outgroup
clade or grade will leave the ingroup rooted
(that is, the new root will be of degree 2).

`nofuse`: if true, keep nodes (and edges) provided that they have at least
one descendant leaf, even if they are of degree 2.
This will keep two-cycles (forcing `simplify` to false).
Nodes without any descendant leaves are deleted.
If `nofuse` is false, edges adjacent to degree-2 nodes are fused.

`multgammas`: if true, the fused edge has γ equal to the product of
the hybrid edges that have been fused together, which may result in
tree edges with γ<1, or with reticulations in which the two parent
γ don't add up to 1.

`keeporiginalroot`: if true, keep the root even if it is of degree one.

Warning: does **not** update attributes related to level-1 networks,
such as inCycle, partition, gammaz, etc.
Does not require branch lengths, and designed to work on networks
of all levels.
"""
function deleteleaf!(net::HybridNetwork, node::Node; kwargs...)
    node.leaf || error("node number $(node.number) is not a leaf.")
    deleteleaf!(net, node.number; kwargs..., index=false)
end

function deleteleaf!(net::HybridNetwork, nodeName::AbstractString; kwargs...)
    tmp = findall(n -> n.name == nodeName, net.node)
    if length(tmp)==0
        error("node named $nodeName was not found in the network.")
    elseif length(tmp)>1
        error("several nodes were found with name $nodeName.")
    end
    deleteleaf!(net, tmp[1]; kwargs..., index=true)
end

# recursive algorithm. nodes previously removed are all necessaily
# *younger* than the current node to remove. Stated otherwise:
# edges previously removed all go "down" in time towards current node:
# - tree edge down to an original leaf,
# - 2 hybrid edges down to a hybrid node.
# hybrid edges from node to another node are not removed. fused instead.
# consequence: node having 2 hybrid edges away from node should not occur.
function deleteleaf!(net::HybridNetwork, nodeNumber::Integer;
                     index=false::Bool, nofuse=false::Bool,
                     simplify=true::Bool, unroot=false::Bool,
                     multgammas=false::Bool, keeporiginalroot=false::Bool)
    i = nodeNumber # good if index=true
    if !index
        i = findfirst(n -> n.number == nodeNumber, net.node)
        i !== nothing ||
            error("cannot delete leaf number $(nodeNumber) because it is not part of net")
    elseif i > length(net.node)
        error("node index $i too large: the network only has $(length(net.node)) nodes.")
    end
    nodei = net.node[i]
    nodeidegree = length(nodei.edge)
    if nodeidegree == 0
        length(net.node)==1 || error("leaf $(nodei.name) has no edge but network has $(length(net.node)) nodes (instead of 1).")
        println("Only 1 node. Removing it: the network will be empty")
        deleteNode!(net,nodei) # empties the network
    elseif nodeidegree == 1
        pe = nodei.edge[1]
        pn = getOtherNode(pe, nodei) # parent node of leaf
        if net.root == i && keeporiginalroot
            return nothing
        end
        # remove leaf and pe.
        removeNode!(pn,pe)  # perhaps useless. in case gc() on pe affects pn
        removeEdge!(pn,pe)
        deleteEdge!(net,pe,part=false)
        if net.root==i # if node was the root, new root = pn
            net.root = findfirst(x -> x===pn, net.node)
        end
        deleteNode!(net,nodei) # this updates the index net.root
        if pn.leaf # network had 2 nodes only: pn and the leaf
            length(net.edge)==0 || error("neighbor of leaf $(nodei.name) is another leaf, but network had $(length(net.edge)) edges (instead of 1).")
            length(pn.edge)==0 || error("neighbor of leaf $(nodei.name) is another leaf, which had $(length(pn.edge)) edges (instead of 1)")
            return nothing # all done: exit function
        end
        deleteleaf!(net, pn.number; nofuse = nofuse, simplify=simplify, unroot=unroot, multgammas=multgammas,
                    keeporiginalroot=keeporiginalroot)
        return nothing
    elseif nodeidegree > 2
        # do nothing: nodei has degree 3+ (through recursive calling)
        return nothing
    end
    # if we get to here, nodei has degree 2 exactly: --e1-- nodei --e2--
    if !unroot && i==net.root
        return nothing # node = root of degree 2 and we don't want to unroot
    end
    e1 = nodei.edge[1]
    e2 = nodei.edge[2]
    if e1.hybrid && e2.hybrid
        (nodei ≡ getChild(e1) && nodei ≡ getChild(e2)) ||
            error("after removing descendants, node $(nodei.number) has 2 hybrid edges but is not the child of both.")
        p1 = getParent(e1) # find both parents of hybrid leaf
        p2 = getParent(e2)
        # remove node1 and both e1, e2
        sameparent = (p1≡p2) # 2-cycle
        removeNode!(p1,e1);  removeNode!(p2,e2) # perhaps useless
        removeEdge!(p1,e1);  removeEdge!(p2,e2)
        deleteEdge!(net,e1,part=false); deleteEdge!(net,e2,part=false)
        if net.root==i net.root=getIndex(p1,net); end # should never occur though.
        deleteNode!(net,nodei)
        # recursive call on both p1 and p2.
        deleteleaf!(net, p1.number; nofuse = nofuse, simplify=simplify, unroot=unroot,
                    multgammas=multgammas, keeporiginalroot=keeporiginalroot)
        sameparent ||
            deleteleaf!(net, p2.number; nofuse=nofuse, simplify=simplify,
                        unroot=unroot, multgammas=multgammas, keeporiginalroot=keeporiginalroot)
    elseif !nofuse #if keeepNodes, do not fuseedges. The recursion should stop.
        e1 = fuseedgesat!(i,net, multgammas) # fused edge
        if simplify && e1.hybrid # check for 2-cycle at new hybrid edge
            cn = getChild(e1)
            e2 = getPartner(e1, cn) # companion hybrid edge
            pn  = getParent(e1)
            if pn ≡ getParent(e2)
                # e1 and e2 have same child and same parent. Remove e1.
                e2.hybrid=false;
                e2.isMajor=true;
                e2.gamma = addBL(e1.gamma, e2.gamma)
                removeEdge!(pn,e1); removeEdge!(cn,e1)
                deleteEdge!(net,e1,part=false)
                # call recursion again because pn and/or cn might be of degree 2.
                deleteleaf!(net, cn.number; nofuse = nofuse, simplify=simplify, unroot=unroot,
                            multgammas=multgammas, keeporiginalroot=keeporiginalroot)
            end
        end
    end
    return nothing
end

"""
    resetNodeNumbers!(net::HybridNetwork; checkPreorder=true, type=:ape)

Change internal node numbers of `net` to consecutive numbers from 1 to the total
number of nodes.

keyword arguments:
- `type`: default is `:ape`, to get numbers that satisfy the conditions assumed by the
  `ape` R package: leaves are 1 to n, the root is n+1, and internal nodes
  are higher consecutive integers.
  If `:postorder`, nodes are numbered in post-order,
  with leaves from 1 to n (and the root last).
  If `:internalonly`, leaves are unchanged. Only internal nodes are modified,
  to take consecutive numbers from (max leaf number)+1 and up. With this
  last option, the post-ordering of nodes is by-passed.
- `checkPreorder`: if false, the `isChild1` edge field and the `net.nodes_changed`
  network field are supposed to be correct (to get nodes in preorder).
  This is not needed when `type=:internalonly`.

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

julia> PhyloNetworks.resetNodeNumbers!(net; type=:postorder)

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
function resetNodeNumbers!(net::HybridNetwork;
    checkPreorder=true::Bool,
    type::Symbol=:ape)
    if checkPreorder
      directEdges!(net)
      preorder!(net) # to create/update net.nodes_changed
    end
    # first: re-number the leaves
    if type == :internalonly
        lnum = maximum(n.number for n in net.node if n.leaf) + 1
    else
        lnum = 1 # first number
        for n in net.node
            n.leaf || continue
            n.number = lnum
            lnum += 1
        end
    end
    # second: re-number internal nodes
    if type == :ape
        nodelist = net.nodes_changed # pre-order: root first
    elseif type == :postorder
        nodelist = reverse(net.nodes_changed) # post-order
    elseif type == :internalonly
        nodelist = net.node
    end
    for n in nodelist
        !n.leaf || continue
        n.number = lnum
        lnum += 1
    end
end

"""
    resetEdgeNumbers!(net::HybridNetwork, verbose=true)

Check that edge numbers of `net` are consecutive numbers from 1 to the total
number of edges. If not, reset the edge numbers to be so.
"""
function resetEdgeNumbers!(net::HybridNetwork, verbose=true::Bool)
    enum = [e.number for e in net.edge]
    ne = length(enum)
    unused = setdiff(1:ne, enum)
    if isempty(unused)
        return nothing # all good
    end
    verbose && @warn "resetting edge numbers to be from 1 to $ne"
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

"""
    unzip_canonical!(net::HybridNetwork)

Unzip all reticulations: set the length of child edge to 0, and increase
the length of both parent edges by the original child edge's length,
to obtain the canonical version of the network according to
Pardi & Scornavacca (2015).

Output: vector of hybrid node in postorder, vector of child edges
whose length is constrained to be 0, and vector of their original
branch lengths to re-zip if needed using [`rezip_canonical!`](@ref).

Assumption: `net.hybrid` is correct, but a preordering of all nodes
is *not* assumed.

Note: This unzipping is not as straightforward as it might seem, because
of "nested" zippers: when the child of a hybrid node is itself a hybrid node.
The unzipping is propagated all the way through.
"""
function unzip_canonical!(net::HybridNetwork)
    hybchild = Dict{Int,Tuple{Edge,Float64}}() # keys: hybrid node number (immutable)
    hybladder = Dict{Int,Union{Nothing,Node}}()
    for h in net.hybrid
        ce = getChildEdge(h)
        hybchild[h.number] = (ce, ce.length) # original lengths before unzipping
        hybladder[h.number] = (ce.hybrid ? getChild(ce) : nothing)
    end
    hybrid_po = Node[] # will list hybrid nodes with partial post-order
    hybpriority = PriorityQueue(h => i for (i,h) in enumerate(net.hybrid))
    nextpriority = length(hybpriority)+1
    while !isempty(hybpriority)
        h,p = peek(hybpriority)
        hl = hybladder[h.number]
        if isnothing(hl) || !haskey(hybpriority,hl)
            # hl no longer key because had priority < p, so already dequeued
            push!(hybrid_po, h)
            delete!(hybpriority, h)
        else
            hybpriority[h] = nextpriority
            nextpriority += 1
        end
    end
    zeroedge = Edge[]
    originallength = Float64[]
    for h in hybrid_po # partial post-order: child < parent if hybrid ladder
        ce, ol = hybchild[h.number] # child edge, original length
        push!(zeroedge, ce)
        push!(originallength, ol)
        unzipat_canonical!(h,ce)
    end
    return hybrid_po, zeroedge, originallength
end

"""
    unzipat_canonical!(hybnode::Node, childedge::Edge)

Unzip the reticulation a node `hyb`. See [`unzip_canonical!`](ref).
Warning: no check that `hyb` has a single child.

Output: constrained edge (child of `hyb`) and its original length.
"""
function unzipat_canonical!(hybnode::Node, childedge::Edge)
    clen = childedge.length # might be different from original if hybrid ladder
    for e in hybnode.edge
        if e === childedge
            e.length = 0.0
        else
            e.length += clen
        end
    end
    return clen
end

"""
    rezip_canonical!(hybridnodes::Vector{Node}, childedges::Vector{Edge},
                     originallengths::Vector{Float64})

Undo [`unzip_canonical!`](@ref).
"""
function rezip_canonical!(hybridnode::Vector{Node}, childedge::Vector{Edge},
                          originallength::Vector{Float64})
    for (i,h) in enumerate(hybridnode) # assumed in post-order
        ce = childedge[i]
        ol = originallength[i]
        lendiff = ol - ce.length # ce.length might be temporarily < 0 if hybrid ladder
        for e in h.edge
            if e === ce
                e.length = ol
            else
                e.length -= lendiff
            end
        end
    end
    return nothing
end
