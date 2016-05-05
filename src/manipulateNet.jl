# function to give all the networks obtained from moving the hybrid node
# inside its cycle
# WARNING: assumes net has all the attributes. It is called inside optTopRuns only
function undirectedOtherNetworks(net0::HybridNetwork; outgroup="none"::AbstractString)
    otherNet = HybridNetwork[]
    for(i in 1:net0.numHybrids) #need to do for by number, not node
        net = deepcopy(net0) # to avoid redoing attributes after each cycle is finished
        ## undo attributes at current hybrid node:
        hybrid = net.hybrid[i]
        nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,hybrid);
        !nocycle || error("the hybrid node $(hybrid.number) does not create a cycle")
        edgesRoot = identifyContainRoot(net,hybrid);
        edges = hybridEdges(hybrid);
        undoGammaz!(hybrid,net);
        othermaj = getOtherNode(edges[1],hybrid)
        edgesmaj = hybridEdges(othermaj)
        DEBUG && println("edgesmaj[3] $(edgesmaj[3].number) is the one to check if containRoot=false already: $(edgesmaj[3].containRoot)")
        if(edgesmaj[3].containRoot) #if containRoot=true, then we need to undo
            undoContainRoot!(edgesRoot);
        end
        ## changes to new hybrid node:
        for(newn in nodesInCycle)
            if(newn.number != hybrid.number) # nodesInCycle contains the hybrid too
                newnet = deepcopy(net)
                newnocycle, newedgesInCycle, newnodesInCycle = identifyInCycle(newnet,newnet.hybrid[i]);
                ind = getIndexNode(newn.number,newnet) # find the newn node in the new network
                hybridatnode!(newnet, newnet.hybrid[i], newnet.node[ind])
                undoInCycle!(newedgesInCycle, newnodesInCycle);
                DEBUG && println("")
                ##undoPartition!(net,hybrid, edgesInCycle)
                success, hybrid0, flag, nocycle, flag2, flag3 = updateAllNewHybrid!(newnet.node[ind], newnet, false,false,false)
                if(success)
                    push!(otherNet,newnet)
                else
                    println("the network obtained by putting the new hybrid in node $(newnet.node[ind].number) is not good, inCycle,gammaz,containRoot: $([flag,flag2,flag3]), we will skip it")
                end
            end
        end
    end
    # check root in good position
    if(outgroup == "none")
        for(n in otherNet)
            !isTree(n) && checkRootPlace!(n, verbose=true)
        end
        return otherNet
    else ## root already in good place
        whichKeep = rep(true,length(otherNet))
        i = 1
        for(n in otherNet)
            if(!isTree(n))
                try
                    checkRootPlace!(n, verbose=true, outgroup=outgroup)
                catch
                    whichKeep[i] = false
                end
            end
            i = i+1;
        end
        return otherNet[whichKeep]
    end
end

# function to change the hybrid node in a cycle
# will try to update incycle inside
"""
`hybridatnode!(net::HybridNetwork, nodeNumber::Int64)`

Changes the hybrid in a cycle to the node defined in nodeNumber. The
node with nodeNumber must be in a cycle. If the node is not in a
cycle, this function will prompt an error.

# Example #"
```julia
julia> net = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
julia> plot(net, showNodeNumber=true)
julia> hybridatnode!(net, -4)
julia> plot(net)
```
""" #"
function hybridatnode!(net::HybridNetwork, nodeNumber::Int64)
    undoInCycle!(net.edge, net.node)
    for(n in net.hybrid)
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

# auxiliary function to change the hybrid node inside a cycle to another node
# WARNING: it assumes all the attributes are correct
# it is called by hybridatnode! with node number as input, and it is called
# by undirectedOtherNetworks
function hybridatnode!(net::HybridNetwork, hybrid::Node, newNode::Node)
    hybrid.hybrid || error("node $(hybrid.number) should be hybrid, but it is not")
    hybedges = hybridEdges(hybrid)
    makeEdgeTree!(hybedges[1],hybrid)
    makeEdgeTree!(hybedges[2],hybrid)
    hybedges[1].inCycle = hybrid.number #just to keep attributes ok
    hybedges[2].inCycle = hybrid.number
    switchHybridNode!(net,hybrid,newNode)
    found = false
    for(e in newNode.edge)
        if(e.inCycle == hybrid.number)
            if(!found)
                found = true
                makeEdgeHybrid!(e,newNode, 0.51, switchHyb=true) #first found, major edge, need to optimize gamma anyway
                ##e.gamma = -1
                e.containRoot = true
            else
                makeEdgeHybrid!(e,newNode, 0.49, switchHyb=true) #second found, minor edge
                ##e.gamma = -1
                e.containRoot = true
            end
        end
    end
end

# function to change the hybrid node in a cycle
# does not assume that the network was read with readTopologyUpdate
# does not modify net0 because it needs to update all attributes
# so, it returns the new network
function hybridatnode(net0::HybridNetwork, nodeNumber::Int64)
    net = readTopologyUpdate(writeTopologyLevel1(net0)) # we need inCycle attributes
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
    for(e in net.node[ind].edge)
        if(e.inCycle == hybrid.number)
            if(!found)
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

Returns the network. Gives a message if the desired root placement was in
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
        info("node $(net.node[ind].number) is a leaf. Will create a new node if needed, to set taxon \"$(net.node[ind].name)\" as outgroup.")
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
        try
          directEdges!(net)
        catch e
          if isa(e, RootMismatch) # new root incompatible with hybrid directions: revert back
            println("RootMismatch: ", e.msg, "\nReverting to old root position.")
            net.root = rootsaved
          else rethrow(e); end
        end
        if (net.root != rootsaved && length(net.node[rootsaved].edge)==2)
            fuseedgesat!(rootsaved,net) # remove old root node if degree 2
        end
        return net
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
    try
      directEdges!(net)
    catch e
      if isa(e, RootMismatch) # new root incompatible with hybrid directions: revert back
        println("RootMismatch: ", e.msg, "\nReverting to old root position.")
        fuseedgesat!(net.root,net) # reverts breakedge!
        net.root = rootsaved
      else rethrow(e); end
    end
    if (net.root != rootsaved && length(net.node[rootsaved].edge)==2)
        fuseedgesat!(rootsaved,net) # remove old root node if degree 2
    end
    return net
end

"""
`breakedge!(Edge, HybridNetwork)`

breaks an edge into 2 edges (each of length half that of original edge).
creates new node of degree 2. Useful to root network along an edge.

warning: updates `isChild1` and `containRoot`, but
does NOT update attributes like: inCycle, partition, gammaz, etc.

returns the index of the newly created node in the network
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

returns the fused edge.
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
    if net.root==i # if `node` was the root, new root = pn
        net.root = getIndex(pn,net)
    end
    deleteNode!(net,net.node[i])
    deleteEdge!(net,pe,part=false) # do not update partitions. irrelevant for networks of level>1.
    return ce
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
        info("node $(node.number) is a leaf, so we will root as an outgroup if possible")
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
# !!bug warning!! If 'outgroup' is already an outgroup, then an extra node
#    of degree 2 is created, in addition to the root (already of degree 2)
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

Returns the network. Throws a 'RootMismatch' Exception if the root was found to
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
    net.cleaned = false # attributed used by snaq! Will change isChild1 and containRoot
    for(e in net.node[net.root].edge)
        traverseDirectEdges!(net.node[net.root],e,true)
    end
    net.isRooted = true
    return net
end

# containroot = true until the path goes through a hybrid node, below which
# containroot is turned to false.
function traverseDirectEdges!(node::Node, edge::Edge, containroot::Bool)
    if (edge.hybrid && node==edge.node[edge.isChild1 ? 1 : 2])
        throw(RootMismatch(
"direction (isChild1) of hybrid edge $(edge.number) conflicts with the root.
isChild1 and containRoot were updated for a subset of edges in the network only."))
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
    push!(queue,net.node[net.root]) # push root into queue
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

# Example #"
```julia
julia> net = readTopology("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
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
""" #"
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


"""
`deleteleaf!(HybridNetwork,Int64; index=false, simplify=true)`
`deleteleaf!(HybridNetwork,leafName::AbstractString; simplify=true)`
`deleteleaf!(HybridNetwork,Node; simplify=true)`

Deletes a leaf node from the network, possibly from its name, number, or index
in the network's array of nodes.

simplify: if true and if deleting the node results in 2 hybrid edges
forming a cycle of k=2 nodes, then these hybrid edges are merged and
simplified as a single tree edge.

The first version does **not** require that `node` is a leaf,
so might be used to remove nodes of degree 2.
The other versions do, and use the default simplify=true.

Warning: does **not** update attributes related to level-1 networks,
such as inCycle, partition, gammaz, etc.
Does not require branch lengths, and designed to work on networks
of all levels.
"""
### WARNING: this is similar but also very different from
# deleteLeaf! in pseudolik.jl, which
# - does not necessarily remove nodes of degree 2,
# - requires and updates all attributes for level-1 networks:
#   inCycle, partition, branch lengths, diamond/triangle types etc.
# - is used a lot within snaq! to extract quartets and retains info
#   on which parameters in the full network affect the quartet.
# deleteIntLeaf! somewhat similar to fuseedgesat!
function deleteleaf!(net::HybridNetwork, node::Node; simplify=true::Bool)
    node.leaf || error("node number $(net.node[i].number) is not a leaf.")
    deleteleaf!(net, node.number, index=false, simplify=simplify)
end

function deleteleaf!(net::HybridNetwork, nodeName::AbstractString; simplify=true::Bool)
    tmp = findin([n.name for n in net.node], [nodeName])
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
function deleteleaf!(net::HybridNetwork, nodeNumber::Int64;
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
            (net.node[i]==e1.node[e1.isChild1?1:2] && net.node[i]==e2.node[e2.isChild1?1:2]) ||
              error("after removing leaf and descendants, node $(net.node[i].number) has 2 hybrid edges but is not the child of both.")
            p1 = getOtherNode(e1, net.node[i]) # find both parents of hybrid leaf
            p2 = getOtherNode(e2, net.node[i])
            # remove the hybrid `node` and both e1, e2
            removeNode!(p1,e1);  removeNode!(p2,e2) # perhaps useless
            removeEdge!(p1,e1);  removeEdge!(p2,e2)
            deleteEdge!(net,e1,part=false); deleteEdge!(net,e2,part=false)
            if net.root==i net.root=getIndex(p1,net); end # should never occur though.
            deleteNode!(net,net.node[i])
            # recursive call on both p1 and p2.
            deleteleaf!(net, p1.number, simplify=simplify)
            deleteleaf!(net, p2.number, simplify=simplify)
        else
            e1 = fuseedgesat!(i,net) # fused edge
            if simplify && e1.hybrid # check for cycle of k=2 nodes
                cn = e1.node[e1.isChild1?1:2]
                e2 = nothing
                for (e in cn.edge) # find companion hybrid edge
                    if (e.hybrid && e!=e1 && cn==e.node[e.isChild1?1:2])
                        e2=e; break;
                    end
                end
                e2!=nothing || error("node $(cn.number) with a single parent hybrid edge")
                pn  = e1.node[e1.isChild1?2:1]
                if pn == e2.node[e2.isChild1?2:1]
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
    #else: do nothing. leaf has degree 3. could occur
    #      through recursive calling reaching a polytomy.
    end
    return nothing
end


