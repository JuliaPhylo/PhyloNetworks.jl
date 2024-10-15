# auxiliary functions for all the other methods
# originally in functions.jl
# Claudia February 2015
#####################



# ----- aux general functions ---------------


function approxEq(a::Number,b::Number,absTol::Number,relTol::Number)
    if(a<eps() || b<eps())
        abs(a-b) < absTol
    else
        abs(a-b) < relTol*eps(abs(a)+abs(b))
    end
end

approxEq(a::Number,b::Number) = approxEq(a,b,1e-5,100)

# isEqual functions: to test if 2 edges (or 2 nodes etc.) "look" alike.
#                    Useful after a deepcopy of a network.
# For nodes (or edges etc.) in the same network, use instead n1 == n2 or n1 != n2.
function isEqual(n1::Node,n2::Node)
    return (n1.number == n2.number && approxEq(n1.gammaz,n2.gammaz) && n1.inCycle == n2.inCycle)
end

function isEqual(n1::Edge,n2::Edge)
    return (n1.number == n2.number && approxEq(n1.length,n2.length))
end

function isEqual(net1::HybridNetwork, net2::HybridNetwork)
    result = true
    result &= (net1.numTaxa == net2.numTaxa)
    result &= (net1.numNodes == net2.numNodes)
    result &= (net1.numEdges == net2.numEdges)
    ## result &= (net1.node == net2.node)
    ## result &= (net1.edge == net2.edge)
    result &= (net1.root == net2.root)
    result &= (net1.names == net2.names)
##    result &= (net1.hybrid == net2.hybrid)
    result &= (net1.numHybrids == net2.numHybrids)
##    result &= (net1.leaf == net2.leaf)
    result &= (net1.ht == net2.ht)
    result &= (net1.numht == net2.numht)
    result &= (net1.numBad == net2.numBad)
    result &= (net1.hasVeryBadTriangle == net2.hasVeryBadTriangle)
    result &= (net1.index == net2.index)
    result &= (net1.loglik == net2.loglik)
    return result
end


#------------- functions to allow for ------------#
#              missing values (lengths or gammas) #

# adds x+y but interprets -1.0 as missing: so -1.0 + x = -1.0 here.
function addBL(x::Number,y::Number)
    (x==-1.0 || y==-1.0) ? -1.0 : x+y
end
function multiplygammas(x::Number,y::Number)
    (x==-1.0 || y==-1.0) ? -1.0 : x * y
end

#------------- EDGE functions --------------------#

# warning: node needs to be defined as hybrid before adding to a
#          hybrid edge. First, an edge is defined as hybrid, and then
#          the nodes are added to it. If the node added is leaf, the
#          edge length is set unidentifiable (as it is external edge)
function setNode!(edge::Edge, node::Node)
    size(edge.node,1)  !=  2 || error("vector of nodes already has 2 values");
    push!(edge.node,node);
    if size(edge.node,1) == 1
        if edge.hybrid
            edge.isChild1 = node.hybrid
        end
        edge.istIdentifiable = !node.leaf
    else
        if node.leaf
            !edge.node[1].leaf || error("edge $(edge.number) has two leaves")
            edge.istIdentifiable = false;
        else
          if edge.hybrid
            if node.hybrid
                # @debug (edge.node[1].hybrid ? "hybrid edge $(edge.number) has two hybrid nodes" : "")
                edge.isChild1 = false;
            else
                edge.node[1].hybrid || error("hybrid edge $(edge.number) has no hybrid nodes");
                edge.isChild1 = true;
            end
          else #edge is tree
            if !edge.node[1].leaf
                if !node.hybrid && !edge.node[1].hybrid
                    edge.istIdentifiable = !edge.fromBadDiamondI
                else
                    if node.hybrid && (node.isBadDiamondI || node.isBadDiamondII || node.isBadTriangle)
                        edge.istIdentifiable = false
                    elseif edge.node[1].hybrid && (edge.node[1].isBadDiamondI ||edge.node[1].isBadDiamondII || edge.node[1].isBadTriangle)
                        edge.istIdentifiable = false
                    else
                        edge.istIdentifiable = true
                    end
                end
            else
                edge.istIdentifiable = false
            end
          end
        end
    end
end

# warning: node needs to be defined as hybrid before adding to a hybrid edge.
#          First, an edge is defined as hybrid, and then the nodes are added to it.
#          If there is a leaf in node, the edge.istIdentifiable=false
function setNode!(edge::Edge,node::Array{Node,1})
    size(node,1) ==  2 || error("vector of nodes must have exactly 2 values")
    edge.node = node;
    if(edge.hybrid)
      if(node[1].hybrid)
          edge.isChild1 = true;
      else
          node[2].hybrid || error("hybrid edge without hybrid node");
          edge.isChild1 = false;
      end
    end
    if(edge.node[1].leaf || edge.node[2].leaf)
        edge.istIdentifiable = false;
    else
        edge.istIdentifiable = true;
    end
end

"""
    getroot(net)

Node used to root `net`. If `net` is to be considered as semi-directed or
unrooted, this root node is used to write the networks' Newick parenthetical
description or for network traversals.

See also: [`isrootof`](@ref)
"""
getroot(net::HybridNetwork) = net.node[net.root]

"""
    isrootof(node, net)

`true` if `node` is the root of `net` (or used as such for network traversals
in case the network is considered as semi-directed); `false` otherwise.

    isleaf(node)
    isexternal(edge)

`true` if `node` is a leaf or `edge` is adjacent to a leaf, `false` otherwise.

See also: [`getroot`](@ref),
[`getparent`](@ref), [`getchild`](@ref)
"""
isrootof(node::Node, net::HybridNetwork) = node === getroot(net)

@doc (@doc isrootof) isleaf
isleaf(node::Node) = node.leaf

@doc (@doc isrootof) isexternal
isexternal(edge::Edge) = any(isleaf.(edge.node))

"""
    isparentof(node, edge)
    ischildof(node, edge)

`true` if `node` is the tail / head, or parent / child, of `edge`; `false` otherwise.
Assumes that the edge's direction is correct, meaning it's field `isChild1` is
reliable (in sync with the rooting).

See also: [`getparent`](@ref), [`getchild`](@ref), [`isrootof`](@ref)
"""
isparentof(node::Node, edge::Edge) = node === getparent(edge)
@doc (@doc isparentof) ischildof
ischildof( node::Node, edge::Edge) = node === getchild(edge)

"""
    hassinglechild(node)

`true` if `node` has a single child, based on the edges' `isChild1` field;
`false` otherwise.

See also: [`getchild`](@ref), [`getparent`](@ref)
"""
hassinglechild(node::Node) = sum(e -> getparent(e) === node, node.edge) == 1

"""
    getchild(edge)
    getchild(node)
    getchildren(node)

Get child(ren) **node(s)**.
- `getchild`: single child node of `edge`, or of `node` after checking that
  `node` has a single child.
- `getchildren`: vector of all children *nodes* of `node`.


    getchildedge(node)

Single child **edge** of `node`. Checks that it's a single child.

*Warning*: these functions rely on correct edge direction, via their `isChild1` field.

See also:
[`getparent`](@ref),
[`getpartneredge`](@ref),
[`isparentof`](@ref),
[`hassinglechild`](@ref).
"""
getchild(edge::Edge) = edge.node[edge.isChild1 ? 1 : 2]
getchild(node::Node) = getchild(getchildedge(node))

@doc (@doc getchild) getchildren
function getchildren(node::Node)
    children = Node[]
    for e in node.edge
        if isparentof(node, e)
            push!(children, getchild(e))
        end
    end
    return children
end

@doc (@doc getchild) getchildedge
function getchildedge(node::Node)
    ce_ind = findall(e -> isparentof(node, e), node.edge)
    length(ce_ind) == 1 || error("node number $(node.number) has $(length(ce_ind)) children instead of 1 child")
    return node.edge[ce_ind[1]]
end

"""
    getparent(edge)
    getparent(node)
    getparentminor(node)
    getparents(node)

Get parental **node(s)**.
- `getparent`: **major** (or only) parent node of `edge` or `node`
- `getparentminor`: minor parent node of `node`
- `getparents`: vector of all parent nodes of `node`.


    getparentedge(node)
    getparentedgeminor(node)

Get one parental **edge** of a `node`.
- `getparentedge`: major parent edge. For a tree node, it's its only parent edge.
- `getparentedgeminor`: minor parent edge, if `node` is hybrid
  (with an error if `node` has no minor parent).
If `node` has multiple major (resp. minor) parent edges, the first one would be
returned without any warning or error.

*Warning*: these functions use the field `isChild1` of edges.

See also: [`getchild`](@ref),
[`getpartneredge`](@ref).
"""
getparent(edge::Edge) = edge.node[edge.isChild1 ? 2 : 1]
@inline function getparent(node::Node)
    for e in node.edge
        if e.isMajor && ischildof(node, e)
            return getparent(e)
        end
    end
    error("could not find major parent of node $(node.number)")
end

@doc (@doc getparent) getparentminor
@inline function getparentminor(node::Node)
    for e in node.edge
        if !e.isMajor && node == getchild(e)
            return getparent(e)
        end
    end
    error("could not find minor parent of node $(node.number)")
end

@doc (@doc getparent) getparents
@inline function getparents(node::Node)
    parents = Node[]
    for e in node.edge
        if ischildof(node, e)
            push!(parents, getparent(e))
        end
    end
    return parents
end

@doc (@doc getparent) getparentedge
@inline function getparentedge(n::Node)
    for ee in n.edge
        if ee.isMajor && ischildof(n,ee)
            return ee
        end
    end
    error("node $(n.number) has no major parent")
end
@doc (@doc getparent) getparentedgeminor
@inline function getparentedgeminor(n::Node)
    for ee in n.edge
        if !ee.isMajor && n == ee.node[(ee.isChild1 ? 1 : 2)]
            return ee
        end
    end
    error("node $(n.number) has no minor parent")
end

"""
    getpartneredge(edge::Edge)
    getpartneredge(edge::Edge, node::Node)

Edge that is the hybrid partner of `edge`, meaning that is has the same child
`node` as `edge`. This child `node` is given as an argument in the second method.
Assumptions, not checked:

- no in-coming polytomy: a node has 0, 1 or 2 parents, no more
- when `node` is given, it is assumed to be the child of `edge`
  (the first method calls the second).

See also: [`getparent`](@ref), [`getchild`](@ref)
"""
@inline function getpartneredge(edge)
    node = getchild(edge)
    getpartneredge(edge, node)
end
@inline function getpartneredge(edge::Edge, node::Node)
    for e in node.edge
        if e.hybrid && e !== edge && node === getchild(e)
            return e
        end
    end
    error("did not find a partner for edge $(edge.number)")
end

"""
    edgerelation(e::Edge, node::Node, origin::Edge)

Return a symbol:

- `:origin` if `e` is equal to `origin`, and otherwise:
- `:parent` if `e` is a parent of `node`,
- `:child` if `e` is a child of `node`

using the `isChild1` attribute of edges.
Useful when `e` iterates over all edges adjacent to `node` and when
`origin` is one of the edges adjacent to `node`,
to known the order in which these edges come.

example:
```julia
labs = [edgerelation(e, u, uv) for e in u.edge] # assuming u is a node of edge uv
parentindex = findfirst(isequal(:parent), labs) # could be 'nothing' if no parent
childindices = findall( isequal(:child), labs)  # vector. could be empty
```
"""
function edgerelation(ee::Edge, n::Node, origin::Edge)
    (ee===origin ? :origin : (n===getchild(ee) ? :parent : :child))
end

# -------------- NODE -------------------------#

function setEdge!(node::Node,edge::Edge)
   push!(node.edge,edge);
   node.hasHybEdge = any(e -> e.hybrid, node.edge)
end

function getOtherNode(edge::Edge, node::Node)
  edge.node[1] === node ? edge.node[2] : edge.node[1]
end



# -------------- NETWORK ----------------------- #

function getIndex(node::Node, net::Network)
    i = 1;
    while(i<= size(net.node,1) && !isEqual(node,net.node[i]))
        i = i+1;
    end
    i <= size(net.node,1) || error("node $(node.number) not in network")
    return i
end

function getIndex(edge::Edge, net::Network)
    i = 1;
    while(i<= size(net.edge,1) && !isEqual(edge,net.edge[i]))
        i = i+1;
    end
    i <= size(net.edge,1) || error("edge $(edge.number) not in network")
    return i
end

function getIndex(edge::Edge, edges::Vector{Edge})
    i = 1;
    while(i<= size(edges,1) && !isEqual(edge,edges[i]))
        i = i+1;
    end
    i <= size(edges,1) || error("edge $(edge.number) not in array of edges")
    return i
end

# aux function to find the index of a node in a
# node array
function getIndex(name::Node, array::Array{Node,1})
    i = 1;
    while(i<= size(array,1) && !isequal(name,array[i]))
        i = i+1;
    end
    i <= size(array,1) || error("$(name.number) not in array")
    return i
end


function getIndexNode(number::Integer,net::Network)
    ind = findfirst(n -> n.number == number, net.node)
    if ind === nothing
        error("node number not in net.node")
    end
    return ind
end

function getIndexEdge(number::Integer,net::Network)
    ind = findfirst(x -> x.number == number, net.edge)
    if ind === nothing
        error("edge number not in net.edge")
    end
    return ind
end

# find the index of an edge in node.edge
function getIndexEdge(edge::Edge,node::Node)
    findfirst(e -> isequal(edge,e), node.edge)
end

# find the index of an edge with given number in node.edge
# bug found & fixed 2019-08-22. Unused function?
function getIndexEdge(number::Integer,node::Node)
    findfirst(e -> isequal(number,e.number), node.edge)
end

# find the index of a node in edge.node
function getIndexNode(edge::Edge,node::Node)
    size(edge.node,1) == 2 || @warn "this edge $(edge.number) has more or less than 2 nodes: $([n.number for n in edge.node])"
    if isequal(node,edge.node[1])
        return 1
    elseif isequal(node,edge.node[2])
        return 2
    else
        error("node not in edge.node")
    end
end

# function to find hybrid index in net.hybrid
function getIndexHybrid(node::Node, net::Network)
    node.hybrid || error("node $(node.number) is not hybrid so it cannot be in net.hybrid")
    i = 1;
    while(i<= size(net.hybrid,1) && !isEqual(node,net.hybrid[i]))
        i = i+1;
    end
    if i>size(net.hybrid,1) error("hybrid node not in network"); end
    return i
end

"""
    getconnectingedge(node1::Node, node2::Node)

Edge shared by (or connecting) `node1` and `node2`, that is: edge incident
to both nodes. An error is thrown if the 2 nodes are not connected.

See also [`isconnected`](@ref)
"""
function getconnectingedge(node1::Node, node2::Node)
    for e1 in node1.edge
        for e2 in node2.edge
            e1 === e2 && return e1
        end
    end
    error("nodes not connected")
    return nothing
end

"""
    isconnected(node1::Node, node2::Node)

Check if two nodes are connected by an edge. Return true if connected, false
if not connected.

See also [`getconnectingedge`](@ref)
"""
function isconnected(node1, node2)
    !isdisjoint(node1.edge, node2.edge) # requires Julia v1.5
end



function isTree(net::HybridNetwork)
    net.numHybrids == length(net.hybrid) || error("numHybrids does not match to length of net.hybrid")
    net.numHybrids != 0 || return true
    return false
end

# function to push a Node in net.node and
# update numNodes and numTaxa
function pushNode!(net::Network, n::Node)
    push!(net.node,n);
    net.numNodes += 1;
    if(n.leaf)
        net.numTaxa += 1
        push!(net.leaf,n);
    end
    if(n.hybrid)
        pushHybrid!(net,n)
    end
end

# function to push an Edge in net.edge and
# update numEdges
function pushEdge!(net::Network, e::Edge)
    push!(net.edge,e);
    net.numEdges += 1;
end


# function to push a hybrid Node in net.hybrid and
# update numHybrids
function pushHybrid!(net::Network, n::Node)
    if(n.hybrid)
        push!(net.hybrid,n);
        net.numHybrids += 1;
    else
        error("node $(n.number) is not hybrid, so cannot be pushed in net.hybrid")
    end
end

"""
    deleteNode!(net::HybridNetwork, n::Node)

Delete node `n` from a network, i.e. removes it from
net.node, and from net.hybrid or net.leaf as appropriate.
Update attributes `numNodes`, `numTaxa`, `numHybrids`.


Warning: if the root is deleted, the new root is arbitrarily set to the
first node in the list. This is intentional to save time because this function
is used frequently in snaq!, which handles semi-directed (unrooted) networks.
"""
function deleteNode!(net::HybridNetwork, n::Node)
    index = findfirst(no -> no===n, net.node)
    # warning: isequal does ===
    #          isEqual (from above) could match nodes across different networks
    index !== nothing || error("Node $(n.number) not in network");
    deleteat!(net.node,index);
    net.numNodes -= 1;
    if net.root == index  # do not check containRoot to save time in snaq!
        net.root = 1      # arbitrary
    elseif net.root > index
        net.root -= 1
    end
    if n.hybrid
       removeHybrid!(net,n)
    end
    if n.leaf
        removeLeaf!(net,n)
    end
end


"""
    deleteEdge!(net::HybridNetwork,  e::Edge; part=true)

Delete edge `e` from `net.edge` and update `net.numEdges`.
If `part` is true, update the network's partition field.
"""
function deleteEdge!(net::HybridNetwork, e::Edge; part=true::Bool)
    if part
        if e.inCycle == -1 && !e.hybrid && !isempty(net.partition) && !isTree(net)
            ind = whichPartition(net,e)
            indE = getIndex(e,net.partition[ind].edges)
            deleteat!(net.partition[ind].edges,indE)
        end
    end
    i = findfirst(x -> x===e, net.edge)
    i !== nothing || error("edge $(e.number) not in network: can't delete");
    deleteat!(net.edge, i);
    net.numEdges -= 1;
end




# function to delete a leaf node in net.leaf
# and update numTaxa
function removeLeaf!(net::Network,n::Node)
    n.leaf || error("cannot delete node $(n.number) from net.leaf because it is not leaf")
    index = findfirst(no -> no === n, net.leaf)
    index !== nothing || error("leaf node $(n.number) not in network")
    deleteat!(net.leaf,index)
    net.numTaxa -= 1
end

# function to delete an internal node with only 2 edges
function deleteIntNode!(net::Network, n::Node)
    size(n.edge,1) == 2 || error("node $(n.number) does not have only two edges")
    index = n.edge[1].number < n.edge[2].number ? 1 : 2;
    edge1 = n.edge[index]; # edge1 will be kept
    edge2 = n.edge[index==1 ? 2 : 1] # we will delete edge2 and n, except if edge2 is hybrid
    if edge2.hybrid
        (edge2, edge1) = (edge1, edge2)
        if getchild(edge1) === n || edge2.hybrid
            @error "node with incoming hybrid edge or incident to 2 hybrid edges: will not be removed"
            return nothing
        end
    end
    node2 = getOtherNode(edge2,n)
    removeEdge!(node2,edge2)
    removeNode!(n,edge1)
    setEdge!(node2,edge1)
    setNode!(edge1,node2)
    deleteNode!(net,n)
    deleteEdge!(net,edge2)
    return nothing
end


# search the hybrid node(s) in network: returns the hybrid node(s)
# in an array
# throws error if no hybrid in network
function searchHybridNode(net::Network)
    a = [n for n in net.node if n.hybrid]
    suma = length(a)
    suma != 0 || error("network has no hybrid node")
    return a
end

# search and returns the hybrid edges in network
# throws error if no hybrid in network
function searchHybridEdge(net::Network)
    a = [n for n in net.edge if n.hybrid]
    suma = length(a)
    suma != 0 || error("network has no hybrid edge")
    return a
end



# print for every node, inCycle and edges
"""
    printNodes(net)
    printNodes(io, net)

Print information on the nodes of a `HybridNetwork` net: node number,
whether it's a leaf, whether it's a hybrid node, whether it's connected to one
or more hybrid edges, it's name (label),
the cycle in which it is belong (-1 if no cycle; makes sense for level-1 networks),
and the list of edges attached to it, by their numbers.
"""
printNodes(x) = printNodes(stdout::IO, x)
function printNodes(io::IO, net::Network)
    namepad = max(4, maximum(length.([n.name for n in net.node])))
    println(io, "node leaf  hybrid hasHybEdge ", rpad("name", namepad), " inCycle edges'numbers")
    for n in net.node
        @printf(io, "%-4d %-5s %-6s %-10s ", n.number, n.leaf, n.hybrid, n.hasHybEdge)
        print(io, rpad(n.name,namepad))
        @printf(io, " %-7d", n.inCycle)
        for e in n.edge
            @printf(io, " %-4d", e.number)
        end
        print(io, "\n")
    end
end

"""
    hybridEdges(node::Node)

Return the 3 edges attached to `node` in a specific order [e1,e2,e3].
**Warning**: assume a level-1 network with node field `hasHybEdge`
and edge field `inCycle` up-to-date.

If `node` is a hybrid node:

- e1 is the major hybrid parent edge of `node`
- e2 is the minor hybrid parent edge
- e3 is the tree edge, child of `node`.

If `node` is a tree node parent of one child edge:

- e1 is the hybrid edge, child of `node`
- e2 is the tree edge that belongs to the cycle created by e1
- e3 is the other tree edge attached to `node` (not in a cycle)

Otherwise:

- e3 is an external edge from `node` to a leaf, if one exists.
"""
function hybridEdges(node::Node)
    size(node.edge,1) == 3 || error("node $(node.number) has $(size(node.edge,1)) edges instead of 3");
    if(node.hybrid)
        hybmajor = nothing;
        hybminor = nothing;
        tree = nothing;
        for e in node.edge
            (e.hybrid && e.isMajor) ? hybmajor = e : nothing
            (e.hybrid && !e.isMajor) ? hybminor = e : nothing
            !e.hybrid ? tree = e : nothing
        end
        return hybmajor, hybminor, tree
    elseif(node.hasHybEdge)
        hybrid = nothing;
        treecycle = nothing;
        tree = nothing;
        for e in node.edge
            (e.hybrid) ? hybrid = e : nothing
            (!e.hybrid && e.inCycle != -1) ? treecycle = e : nothing
            (!e.hybrid && e.inCycle == -1) ? tree = e : nothing
        end
        return hybrid, treecycle, tree
    else
        #@warn "node $(node.number) is not hybrid $(node.hybrid) nor tree with hybrid edges (hasHybEdge) $(node.hasHybEdge), return the node.edge in order, unless a leaf is attached, then the edge attached to leaf is last";
        edge1 = nothing
        edge2 = nothing
        edge3 = nothing
        leaffound = false
        ind = 1
        for i in 1:3
            if(getOtherNode(node.edge[i],node).leaf)
                leaffound = true
                edge3 = node.edge[i]
                ind = i
                break
            end
        end
        if(leaffound)
            if(ind == 1)
                return node.edge[2], node.edge[3], edge3
            elseif(ind == 2)
                return node.edge[1], node.edge[3], edge3
            elseif(ind == 3)
                return node.edge[1], node.edge[2], edge3
            end
        else
            return node.edge[1], node.edge[2], node.edge[3]
        end
    end
end

"""
    hybridEdges(node::Node, e::Edge)

Return the 2 edges connected to `node` other than `e`,
in the same order as `node.edge`,
except that `e` absent from the list.

Despite what the name suggest, `node` need not be a hybrid node!
`node` is assumed to have 3 edges, though.
"""
function hybridEdges(node::Node, edge::Edge)
    size(node.edge,1) == 3 || error("node $(node.number) has $(size(node.edge,1)) edges instead of 3")
    edge1 = nothing
    edge2 = nothing
    for e in node.edge
        if(!isequal(e,edge))
            isa(edge1,Nothing) ? edge1 = e : edge2 = e
        end
    end
    return edge1,edge2
end


# function to remove an edge from a node
# warning: deletion is final, you can only
#          have edge back by pushing it again
# warning: if the edge removed is hybrid and node is tree,
#          node.hasHybEdge is set to false
#          assuming any tree node can only have one
#          one hybrid edge
function removeEdge!(node::Node, edg::Edge)
    index = findfirst(x -> x === edg, node.edge)
    index !== nothing || error("edge $(edg.number) not in node $(node.number)")
    deleteat!(node.edge,index)
    node.hasHybEdge = any(e -> e.hybrid, node.edge)
end

# function to remove a node from a edge
# warning: deletion is final, you can only
#          have node back by pushing it again
# warning: only removes node from edge, edge might still
#          be in node.edge
function removeNode!(nod::Node, edge::Edge)
    index = findfirst(x -> x === nod, edge.node)
    index !== nothing || error("node $(nod.number) not in edge")
    deleteat!(edge.node,index);
end


# ----------------------------------------------------------------------------------------



"""
    setGamma!(Edge, new γ)
    setGamma!(Edge, new γ, change_other=true::Bool)

Set inheritance probability γ for an edge, which must be a hybrid edge.
The new γ needs to be in [0,1]. The γ of the "partner" hybrid edge is changed
accordingly, to 1-γ. The field `isMajor` is also changed accordingly.
If the new γ is approximately 0.5, `Edge` is set to the major parent,
its partner is set to the minor parent.

If `net` is a HybridNetwork object, `printEdges(net)` will show the list of edges
and their γ's. The γ of the third hybrid edge (say) can be changed to 0.2 with
`setGamma!(net.edge[3],0.2)`.
This will automatically set γ of the partner hybrid edge to 0.8.

The last argument is true by default. If false: the partner edge is not updated.
This is useful if the new γ is 0.5, and the partner's γ is already 0.5,
in which case the `isMajor` attributes can remain unchanged.
"""
setGamma!(edge::Edge, new_gamma::Float64) = setGamma!(edge, new_gamma, true)

# warning in the bad diamond/triangle cases because gamma is not identifiable
# changeOther = true, looks for the other hybrid edge and changes gamma too

function setGamma!(edge::Edge, new_gamma::Float64, changeOther::Bool)
    new_gamma >= 0.0 || error("gamma has to be positive: $(new_gamma)")
    new_gamma <= 1.0 || error("gamma has to be less than 1: $(new_gamma)")
    edge.hybrid || error("cannot change gamma in a tree edge");
    node = getchild(edge) # child of hybrid edge
    node.hybrid || @warn "hybrid edge $(edge.number) not pointing at hybrid node"
    # @debug (node.isBadDiamondI ? "bad diamond situation: gamma not identifiable" : "")
    partner = Edge[] # list of other hybrid parents of node, other than edge
    for e in node.edge
        if e.hybrid && e != edge && node == getchild(e)
            push!(partner, e)
        end
    end
    length(partner) == 1 ||
      error("strange hybrid node $(node.number) with $(length(partner)+1) hybrid parents")
    e2 = partner[1]
    onehalf = isapprox(new_gamma,0.5)
    if onehalf new_gamma=0.5; end
    new_ismajor = new_gamma >= 0.5
    edge.gamma = new_gamma
    if changeOther
        edge.isMajor = new_ismajor
        e2.gamma = 1.0 - new_gamma
        e2.isMajor = !new_ismajor
    else
        if onehalf # who is major is arbitrary: so we pick what's consistent with the partner
            edge.isMajor = !e2.isMajor
        else
            edge.isMajor = new_ismajor
        end
    end
    return nothing
end

@inline function setmultiplegammas!(edges::Vector{Edge}, gammas::Vector{Float64})
    for (e,g) in zip(edges, gammas)
        setGamma!(e, g)
    end
end

"""
    remove_edgeLengthsGammas!(net::HybridNetwork)

Reset all edge lengths and all hybrid edge γs to be missing (coded as -1.0).
"""
function remove_edgeLengthsGammas!(net::HybridNetwork)
    for e in net.edge
        e.length = -1.0
        if e.hybrid
            e.gamma = -1.0
        end
    end
end

"""
    check_nonmissing_nonnegative_edgelengths(net, str="")

Throw an Exception if `net` has undefined edge lengths (coded as -1.0) or
negative edge lengths. The error message indicates the number of the offending
edge(s), followed by `str`.
"""
function check_nonmissing_nonnegative_edgelengths(net::HybridNetwork, str="")
    if any(e.length == -1.0 for e in net.edge)
        undefined = [e.number for e in net.edge if e.length == -1.0]
        error(string("Branch(es) number ", join(undefined,","), " have no length.\n", str))
    end
    if any(e.length < 0 for e in net.edge)
        negatives = [e.number for e in net.edge if e.length < 0.0]
        error(string("Branch(es) number ", join(negatives,","), " have negative length.\n", str))
    end
end




"""
    setGammas!(net, γ vector)

Set inheritance γ's of hybrid edges, using input vector for *major* edges.
Assume pre-order calculated already, with up-to-date field `nodes_changed`.
See [`getGammas`](@ref).

**Warning**: very different from [`setGamma!`](@ref), which focuses on a
single hybrid event,
updates the field `isMajor` according to the new γ, and is not used here.

**Assumption**: each hybrid node has only 2 parents, a major and a minor parent
(according to the edges' field `isMajor`).
"""
function setGammas!(net::HybridNetwork, gammas::Vector)
    for (i,nod) in enumerate(net.nodes_changed)
        nod.hybrid || continue # skip tree nodes: nothing to do
        majorhyb = getparentedge(nod) # major
        minorhyb = getparentedgeminor(nod) # error if doesn't exit
        majorhyb.gamma = gammas[i]
        minorhyb.gamma = 1 - gammas[i]
    end
    return nothing
end

"""
    getGammas(net)

Vector of inheritance γ's of all major edges (tree edges and major hybrid edges),
ordered according to the pre-order index of their child node,
assuming this pre-order is already calculated
(with up-to-date field `nodes_changed`).
Here, a "major" edge is an edge with field `isMajor` set to true,
regardless of its actual γ (below, at or above 0.5).

See [`setGammas!`](@ref)
"""
function getGammas(net::HybridNetwork)
    gammas = ones(length(net.nodes_changed))
    for (i,node) in enumerate(net.nodes_changed)
        node.hybrid || continue # skip tree nodes: their gamma is already set to 1
        majorhybedge = getparentedge(node) # major
        gammas[i] = majorhybedge.gamma
    end
    return gammas
end



"""
    getHeights(net, checkpreorder::Bool=true)
    getHeights!(net, checkpreorder::Bool=true)

Return the height (distance to the root) of all nodes. The function will return an error if the network is not time-consistent.
Ultrametricity is not assumed: tips need not all be at the same distance from the root.
If `checkpreorder=false`, assumes the network has already been preordered
with [`preorder!`](@ref).

[`getHeights`](@ref) returns an error upon finding a missing edge (i.e., an edge with length `-1`)
The second form [`getHeights!`](@ref) will attempt to fill in missing edge values. 

If a tree edge has a missing length, then the function throws an error. In general, resolving tree edges is nontrivial.

It is assumed that hybrid nodes are not leaves, such that external edges
are necessarily tree edges.
If a hybrid edge has a missing length, this length is changed as follows:
- If both partner hybrid edges lack a length: the shortest lengths are assigned
  to make the network time-consistent at the hybrid node. In particular,
  either the major edge or the minor edge is assigned length 0.0.
- Otherwise: the value needed to make the network time-consistent considering
  based on the partner edge's length. If this value is negative, then an error is thrown.

Output: vector of node heights, one per node, in the same order as in
`net.nodes_changed`. Examples:

```jldoctest
julia> net = readTopology("(((C:1,(A:1)#H1:1.5::0.7):1,(#H1:0.3::0.3,E:2.0):2.2):1.0,O:5.2);");

julia> # using PhyloPlots; plot(net, useedgelength=true, showedgelength=true, shownodenumber=true); # to see

julia> nodeheight = PhyloTraits.getHeights(net)
9-element Vector{Float64}:
 0.0
 5.2
 1.0
 3.2
 5.2
 2.0
 3.5
 4.5
 3.0

julia> [node.number => (nodeheight[i], node.name) for (i,node) in enumerate(net.nodes_changed)]
9-element Vector{Pair{Int64, Tuple{Float64, String}}}:
 -2 => (0.0, "")
  5 => (5.2, "O")
 -3 => (1.0, "")
 -6 => (3.2, "")
  4 => (5.2, "E")
 -4 => (2.0, "")
  3 => (3.5, "H1")
  2 => (4.5, "A")
  1 => (3.0, "C")

```
"""
function getHeights(net::HybridNetwork, checkpreorder::Bool=true)
    getHeightshelper(net,false,inconsistencyerror,false,checkpreorder)
end

function getHeights!(net::HybridNetwork, checkpreorder::Bool=true)
    getHeightshelper(net,true,inconsistencyerror,false,checkpreorder)
end

function getHeightshelper(
    net::HybridNetwork,
    fixmissing::Bool,
    inconsistencyhandler::Function,
    stopearly::Bool,
    checkpreorder::Bool=true,)

    checkpreorder && preorder!(net)

    missing_e = Tuple{Edge,Float64}[] ##Vector of tuples that has missing edges and their 'fixed' values as tuples
    isconsistent = Ref(true)

    rootdistance=recursion_preorder(
        net.nodes_changed,
        get_heights_init,
        updateRecursion_default!,
        update_getheights_tree!,
        update_getheights_hybrid!,
        fixmissing,
        missing_e,
        isconsistent,
        inconsistencyhandler,
        stopearly)

    if fixmissing ##assign values to fixed edges
        (x-> x[1].length=x[2]).(missing_e)
        #for e in missing_e
        #    e[1].length = e[2]
        #end
    end

    rootdistance
end

inconsistencyerror = function(
    paredges::Vector{Edge},
    distanceroot_i_major::Float64,
    distanceroot_i_minor::Float64,
    i::Int,
    distanceroot::Vector{Float64} )
    if approxEq(distanceroot_i_major,distanceroot_i_minor) # They should only be aporximately equal if one of the edges has negative branch length
        error("Edges $((x ->x.length).(paredges)) cannot be fixed without creating a negative branch length")
    else # we errored because they are not equal
        error("The paths that lead to node number $(getchild(paredges[1]).number) both have different heights")
    end
end

function get_heights_init(nodes::Vector{Node},params...)
    n = length(nodes)
    return zeros(Float64,n) 
end


function update_getheights_tree!(rootdistance::Vector{Float64},
    i::Int,
    parind::Int,
    paredge::Edge,
    fixmissing::Bool,
    params...
    ) 
    if paredge.length == -1.0 ##The branch is missing a value
            !fixmissing && error("Edge $(paredge.number) has a missing edge length")
            error("Edge $(paredge.number) has a missing edge length and tree edges are nontrivial to fix, if fixable")
    end
    rootdistance[i] = rootdistance[parind] + paredge.length 
    return true
end

function update_getheights_hybrid!(
    rootdistance::Vector{Float64},
    i,
    parinds::Vector{Int},
    paredges::Vector{Edge},
    fixmissing::Bool,
    missing_e::Vector{Tuple{Edge,Float64}},
    isconsistent::Ref{Bool},
    inconsistencyhandler::Function,
    stopearly::Bool,)

    maj_ind = !paredges[1].isMajor + 1 ## If the first is not major then the 2nd is
    min_ind = paredges[1].isMajor + 1
    
    maj_e = paredges[maj_ind]
    maj_rootdistance = rootdistance[parinds[maj_ind]] 

    min_e = paredges[min_ind]
    min_rootdistance = rootdistance[parinds[min_ind]]


    majmissing = maj_e.length == -1.0
    if majmissing
        !fixmissing &&  error("Edge $(maj_e.number) has a missing edge length.")
    end
    minmissing = min_e.length == -1.0
    if minmissing
        !fixmissing && error("Edge $(min_e.number) has a missing edge length.")
    end

    ##Attempt to fix missing edges
    if majmissing && minmissing # both parent edges lack a length
        if maj_rootdistance < min_rootdistance
            min_e_len = 0.0
            maj_e_len = min_rootdistance - maj_rootdistance
            rootdistance[i] = min_rootdistance
        else
            maj_e_len = 0.0
            min_e_len = maj_rootdistance - min_rootdistance
            rootdistance[i] = maj_rootdistance
        end
        push!(missing_e,(min_e,min_e_len))
        push!(missing_e,(maj_e,maj_e_len))
    elseif majmissing || minmissing # one of the two branches is missing 
        ##Attempt to fix the missing edge length. compute rootdistance[i] given the nonmissing edge and see if we can specify the missing edge s.t. it is time-consistent
        maj_e_len=Inf
        min_e_len=Inf
        if majmissing
            dr_i_min = min_rootdistance + min_e.length #compute the root distance to i from min side
            maj_e_len= dr_i_min - maj_rootdistance 
            dr_i_maj = maj_rootdistance + maj_e_len
            push!(missing_e,(maj_e,maj_e_len))
        else ##minmissing
            dr_i_maj = maj_rootdistance + maj_e.length
            min_e_len = dr_i_maj - min_rootdistance
            dr_i_min = min_rootdistance + min_e_len
            push!(missing_e,(min_e,min_e_len))
        end
        if maj_e_len < 0.0 || min_e_len < 0.0 # one of the branches could not be resolved 
            isconsistent[]=false
            inconsistencyhandler(paredges,dr_i_maj,dr_i_min,i,rootdistance)    
        end
        rootdistance[i] = dr_i_maj # should be the same as dr_i_min
    else # both parent edges have a length
        dr_i_maj = maj_rootdistance+maj_e.length
        dr_i_min = min_rootdistance+min_e.length
        if !isapprox(dr_i_maj, dr_i_min) # root distances not compatible
            isconsistent[]=false
            inconsistencyhandler(paredges,dr_i_maj,dr_i_min,i,rootdistance)  
        else # root distances are compatible 
            rootdistance[i]=dr_i_maj 
        end
    end
    return !stopearly && isconsistent[] ## tell recursion to end if we want to stop early and the phylogeny is not time-consistent
end




function numTreeEdges(net::HybridNetwork)
    2*net.numTaxa - 3 + net.numHybrids
end

function numIntTreeEdges(net::HybridNetwork)
    2*net.numTaxa - 3 + net.numHybrids - net.numTaxa
end


# function to get the partition where an edge is
# returns the index of the partition, or error if not found
# better to return the index than the partition itself, because we need the index
# to use splice and delete it from net.partition later on
# cycle: is the number to look for partition on that cycle only
function whichPartition(net::HybridNetwork,edge::Edge,cycle::Integer)
    !edge.hybrid || error("edge $(edge.number) is hybrid so it cannot be in any partition")
    edge.inCycle == -1 || error("edge $(edge.number) is in cycle $(edge.inCycle) so it cannot be in any partition")
    @debug "search partition for edge $(edge.number) in cycle $(cycle)"
    in(edge,net.edge) || error("edge $(edge.number) is not in net.edge")
    for i in 1:length(net.partition)
        @debug "looking for edge $(edge.number) in partition $(i): $([e.number for e in net.partition[i].edges])"
        if(in(cycle,net.partition[i].cycle))
            @debug "looking for edge $(edge.number) in partition $(i), with cycle $(cycle): $([e.number for e in net.partition[i].edges])"
            if in(edge,net.partition[i].edges)
                @debug "partition for edge $(edge.number) is $([e.number for e in net.partition[i].edges])"
                return i
            end
        end
    end
    @debug begin; printPartitions(net); "" end
    error("edge $(edge.number) is not hybrid, nor part of any cycle, and it is not in any partition")
end

# function to get the partition where an edge is
# returns the index of the partition, or error if not found
# better to return the index than the partition itself, because we need the index
# to use splice and delete it from net.partition later on
function whichPartition(net::HybridNetwork,edge::Edge)
    !edge.hybrid || error("edge $(edge.number) is hybrid so it cannot be in any partition")
    edge.inCycle == -1 || error("edge $(edge.number) is in cycle $(edge.inCycle) so it cannot be in any partition")
    @debug "search partition for edge $(edge.number) without knowing its cycle"
    in(edge,net.edge) || error("edge $(edge.number) is not in net.edge")
    for i in 1:length(net.partition)
        @debug "looking for edge $(edge.number) in partition $(i): $([e.number for e in net.partition[i].edges])"
        if(in(edge,net.partition[i].edges))
            @debug "partition for edge $(edge.number) is $([e.number for e in net.partition[i].edges])"
            return i
        end
    end
    @debug begin printPartitions(net); "printed partitions" end
    error("edge $(edge.number) is not hybrid, nor part of any cycle, and it is not in any partition")
end

# function that will print the partition of net
function printPartitions(net::HybridNetwork)
    println("partition.cycle\t partition.edges")
    for p in net.partition
        println("$(p.cycle)\t\t $([e.number for e in p.edges])")
    end
end

# function to find if a given partition is in net.partition
function isPartitionInNet(net::HybridNetwork,desc::Vector{Edge},cycle::Vector{Int})
    for p in net.partition
        if(sort(cycle) == sort(p.cycle))
            if(sort([e.number for e in desc]) == sort([e.number for e in p.edges]))
                return true
            end
        end
    end
    return false
end




X
"""
    assignhybridnames!(net)

Assign names to hybrid nodes in the network `net`.
Hybrid nodes with an empty `name` field ("") are modified with a name that
does not conflict with other hybrid names in the network. The preferred name
is "H3" if the node number is 3 or -3, but an index other than 3 would be used
if "H3" were the name of another node already.

If two hybrid nodes have non-empty and equal names, the name of one of them is changed and
re-assigned as described above (with a warning).
"""
function assignhybridnames!(net::HybridNetwork)
    rx = r"^H(\d+)$"
    # prep: collect indices 'i' of any tree nodes named like Hi
    trenum = Int[]  # indices 'i' in tree node name, in case some are named Hi
    for n in net.node
        !n.hybrid || continue # do nothing if node n is hybrid
        m = match(rx, n.name)
        m === nothing || push!(trenum, parse(Int, m[1]))
    end
    # first: go through *all* existing non-empty names
    hybnum = Int[]  # indices 'i' in hybrid names: Hi
    for ih in 1:length(net.hybrid)
        hnode = net.hybrid[ih]
        lab = hnode.name
        lab != "" || continue # do nothing if label is missing
        jh = findfirst(isequal(lab), [net.hybrid[j].name for j in 1:ih-1])
        if jh !== nothing # set repeated names to ""
            @warn "hybrid nodes $(hnode.number) and $(net.hybrid[jh].number) have the same label: $lab. Will change the name of the former."
            hnode.name = ""
        else # fill in list of existing indices "i" in Hi
            m = match(rx, lab)
            m !== nothing || continue # skip the rest if name is not of the form Hi
            ind = parse(Int, m[1])
            if ind in trenum
                @warn "hybrid node $(hnode.number) had same label as a tree node: H$ind. Will change hybrid name."
                hnode.name = ""
            else
                push!(hybnum, ind)
            end
        end
    end
    # second: assign empty names to "Hi" for some i
    hnext = 1
    for ih in 1:length(net.hybrid)
        net.hybrid[ih].name == "" || continue # do nothing if non-empty label
        hnum = abs(net.hybrid[ih].number)
        while hnum in hybnum || hnum in trenum
            hnum = hnext  # not efficient, but rare
            hnext += 1    # and okay on small networks
        end
        push!(hybnum, hnum)
        net.hybrid[ih].name = "H$hnum"
    end
end



"""
    setlengths!(edges::Vector{Edge}, lengths::AbstractVector)

Assign new lengths to a vector of `edges`.
Warning: does *not* make any checks that the new edge lengths are non-negative
(except for -1 values to be interpreted as missing).
"""
@inline function setlengths!(edges::Vector{Edge}, lengths::AbstractVector)
    for (e,l) in zip(edges, lengths)
        e.length = l
    end
end

"""
    getlengths(edges::Vector{Edge})

Vector of edge lengths for a vector of `edges`.
"""
getlengths(edges::Vector{Edge}) = [e.length for e in edges]

"""
    hashybridladder(net::HybridNetwork)

Return true if `net` contains a hybrid ladder: where a hybrid node's
child is itself a hybrid node.
This makes the network not treechild, assuming it is fully resolved.
(One of the nodes does not have any tree-node child).
"""
function hashybridladder(net::HybridNetwork)
    for h in net.hybrid
        if any(n.hybrid for n in getparents(h))
            return true
        end
    end
    return false
end

"""
    shrinkedge!(net::HybridNetwork, edge::Edge)

Delete `edge` from net, provided that it is a non-external tree edge.
Specifically: delete its child node (as determined by `isChild1`) and connect
all edges formerly incident to this child node to the parent node of `edge`,
thus creating a new polytomy, unless the child was of degree 2.

Warning: it's best for `isChild1` to be in sync with the root for this. If not,
the shrinking may fail (if `edge` is a tree edge but its "child" is a hybrid)
or the root may change arbitrarily (if the child of `edge` is the root).

Output: true if the remaining node (parent of `edge`) becomes a hybrid node with
more than 1 child after the shrinking; false otherwise (e.g. no polytomy was
created, or the new polytomy is below a tree node)
"""
function shrinkedge!(net::HybridNetwork, edge2shrink::Edge)
    edge2shrink.hybrid && error("cannot shrink hybrid edge number $(edge2shrink.number)")
    cn = getchild(edge2shrink)
    cn.hybrid && error("cannot shrink tree edge number $(edge2shrink.number): its child node is a hybrid. run directEdges! ?")
    pn = getparent(edge2shrink)
    isexternal(edge2shrink) &&  # (isleaf(cn) || isleaf(pn)) &&
      error("won't shrink edge number $(edge2shrink.number): it is incident to a leaf")
    removeEdge!(pn,edge2shrink)
    empty!(edge2shrink.node) # should help gc
    for ee in cn.edge
        ee !== edge2shrink || continue
        cn_index = findfirst(x -> x === cn, ee.node)
        ee.node[cn_index] = pn # ee.isChild1 remains synchronized
        push!(pn.edge, ee)
    end
    pn.hasHybEdge = any(e -> e.hybrid, pn.edge)
    empty!(cn.edge) # should help to garbage-collect cn
    deleteEdge!(net, edge2shrink; part=false)
    deleteNode!(net, cn)
    badpolytomy = false
    if pn.hybrid # count the number of pn's children, without relying on isChild1 of tree edges
        nc = sum((!e.hybrid || getchild(e) !== pn) for e in pn.edge)
        badpolytomy = (nc > 1)
    end
    return badpolytomy
end

@doc raw"""
    shrink2cycles!(net::HybridNetwork, unroot=false::Bool)

If `net` contains a 2-cycle, collapse the cycle into one edge of length
tA + γt1+(1-γ)t2 + tB (see below), and return true.
Return false otherwise.
A 2-cycle is a set of 2 parallel hybrid edges, from the same parent node to the
same hybrid child node.

           A                A
           | tA             |
         parent             |
           | \              |
    t2,1-γ |  | t1,γ        | tA + γ*t1 + (1-γ)*t2 + tB
           | /              |
         hybrid             |
           | tB             |
           B                B

If any of the lengths or gammas associated with a 2-cycle are missing,
the combined length is missing. If γ is missing, branch lengths
are calculated using γ=0.5.

If `unroot` is false and the root is up for deletion, it will be kept only if it
is has degree 2 or more. If `unroot` is true and the root is up for deletion, it
will be kept only if it has degree 3 or more. A root node with degree 1 will be
deleted in both cases.
"""
function shrink2cycles!(net::HybridNetwork, unroot=false::Bool)
    foundcycle = false
    nh = length(net.hybrid)
    ih = nh # hybrids deleted from the end
    while ih > 0
        h = net.hybrid[ih]
        minor = getparentedgeminor(h)
        major = getparentedge(h)
        pmin = getparent(minor) # minor parent node
        pmaj = getparent(major) # major parent node
        if pmin !== pmaj # no 2-cycle
            ih -= 1
            continue
        end
        # 2-cycle
        foundcycle = true
        shrink2cycleat!(net, minor, major, unroot)
        nh = length(net.hybrid)
        ih = nh
        # we re-do if a cycle was removed: a new cycle might have appeared
    end
    return foundcycle
end

"""
    shrink2cycleat!(net::HybridNetwork, minor::Edge, major::Edge, unroot::Bool)

Remove `minor` edge then update the branch length of the remaining `major` edge.
Called by [`shrink2cycles!`](@ref)

Assumption: `minor` and `major` do form a 2-cycle. That is, they start and end
at the same node.
"""
function shrink2cycleat!(net::HybridNetwork, minor::Edge, major::Edge,
                         unroot::Bool)
    g = minor.gamma
    if g == -1.0 g=.5; end
    major.length = addBL(multiplygammas(    g, minor.length),
                         multiplygammas(1.0-g, major.length))
    deletehybridedge!(net, minor, false,unroot,false,false,false) # nofuse,unroot,multgammas,simplify
    return nothing
end

"""
    shrink3cycles!(net::HybridNetwork, unroot=false::Bool)

Remove all 2- and 3-cycles from a network.

Return true if `net` contains a 2-cycle or a 3-cycle; false otherwise.
A 3-cycle (2-cycle) is a set of 3 (2) nodes that are all connected.
One of them must be a hybrid node, since `net` is a DAG.

If `unroot` is false and the root is up for deletion, it will be kept only if it
is has degree 2 or more. If `unroot` is true and the root is up for deletion, it
will be kept only if it has degree 3 or more. A root node with degree 1 will be
deleted in both cases.

See [`shrink3cycleat!`](@ref) for details on branch lengths and
inheritance values.
"""
function shrink3cycles!(net::HybridNetwork, unroot=false::Bool)
    foundcycle = false
    nh = length(net.hybrid)
    ih = nh # hybrids deleted from the end
    while ih > 0
        h = net.hybrid[ih]
        minor = getparentedgeminor(h)
        major = getparentedge(h)
        pmin = getparent(minor) # minor parent node
        pmaj = getparent(major) # major parent node
        if pmin === pmaj # 2-cycle
            foundcycle = true
            shrink2cycleat!(net, minor, major, unroot)
            nh = length(net.hybrid)
            ih = nh + 1 # start over if a cycle was removed by setting ih = nh + 1.
                        # Shrinking could have created a new cycle.
        else # 3-cycle
            result = shrink3cycleat!(net, h, minor, major, pmin, pmaj, unroot)
            if result
                foundcycle = true
                nh = length(net.hybrid)
                ih = nh + 1 # start over as above
            end
        end
        ih -= 1
    end
    return foundcycle
end

@doc raw"""
    shrink3cycleat!(net::HybridNetwork, hybrid::Node, edge1::Edge, edge2::Edge,
                    node1::Node, node2::Node, unroot::Bool)

Replace a 3-cycle at a given `hybrid` node by a single node, if any.
Assumption: `edge1` (`node1`) and `edge2` (`node2`) are the parent edges (nodes)
of `hybrid`. Return true if a 3-cycle is found and removed, false otherwise.
There is a 3-cycle if nodes 1 & 2 are connected, by an edge called `e3` below.

There are two cases, with differing effects on the γ inheritance
values and branch lengths.

**Hybrid case**: the 3-cycle is shrunk to a hybrid node, which occurs if
either node 1 or 2 is a hybrid node (that is, e3 is hybrid). If e3 goes
from node 1 to node 2, the 3-cycle (left) is shrunk as on the right:

    \eA      /eB           \eA  /eB
     1--e3->2       γ1+γ2γ3 \  / γ2(1-γ3)
      \    /               hybrid
     γ1\  /γ2
      hybrid

with new branch lengths:
new tA = tA + (γ1.t1 + γ2γ3.(t2+t3))/(γ1+γ2γ3),
new tB = tB + t2,
provided that γ1, γ2=1-γ1, and γ3 are not missing. If one of them is missing
then γ1 and γ2 remain as is, and e3 is deleted naively,
such that new tA = tA + t1 and new tB = tB + t2.
If γ's are not missing but one of t1,t2,t3 is missing, then the γ's are
updated to γ1+γ2γ3 and γ2(1-γ3), but t's are update naively.

**Tree case**: the 3-cycle is shrunk to a tree node, which occurs if node 1 & 2
are both tree nodes (that is, e3 is a tree edge). If eC is the child edge of
`hybrid`, the 3-cycle (left) is shrunk as on the right:

    \eA                  \eA
     1--e3--2--eB--       \
      \    /               n--eB--
     γ1\  /γ2              |
      hybrid               |eC
        |
        |eC

with new branch lengths:
new tA = tA + γ2.t3,
new tB = tB + γ1.t3,
new tC = tC + γ1.t1 + γ2.t2,
provided that γ1, γ2=1-γ1, t1, t2 and t3 are not missing.
If one is missing, then e1 is deleted naively such that
tB is unchanged, new tC = tC + t2 and new tA = tA + t3.
"""
function shrink3cycleat!(net::HybridNetwork, hybrid::Node, edge1::Edge,
                         edge2::Edge, node1::Node, node2::Node, unroot::Bool)
    # check for presence of 3 cycle
    edge3 = nothing
    for e in node1.edge # find edge connecting node1 and node2
        e !== edge1 || continue
        n = getOtherNode(e, node1)
        if n === node2
            edge3 = e
            break
        end
    end
    !isnothing(edge3) || return false # no 3-cycle at node h
    # identify case type
    if edge3.hybrid # one of the parent nodes is a hybrid
        # to shrink this, delete edge connecting these two nodes (edge3 here)
        if getchild(edge3) === node1
            node1, node2 = node2, node1
            edge1, edge2 = edge2, edge1
        end # now: node1 --edge3--> node2
        edgeA = nothing
        for e in node1.edge
            if e!== edge1 && e !== edge2
                edgeA = e
                break
            end
        end
        g1 = edge1.gamma
        g2g3 = multiplygammas(edge2.gamma, edge3.gamma)
        g1tilde = addBL(g1, g2g3)
        if g1tilde != -1.0 # else one of the γ is missing: do nothing with γs and ts
            edge1.gamma = g1tilde
            edge2.gamma = 1.0-g1tilde
            edge1.isMajor = g1tilde >= 0.5
            edge2.isMajor = !edge1.isMajor
            if edge1.length != -1.0 && edge2.length != -1.0 && edge3.length != -1.0
                edge1.length = (edge1.length *g1 + (edge3.length + edge2.length)*g2g3)/g1tilde
            end
        end
        deletehybridedge!(net, edge3, false,unroot,false,false) # nofuse,unroot,multgammas,simplify
    else # parent nodes 1 and 2 are both tree nodes
        edgeB = nothing
        for e in node2.edge
            if e !== edge1 && e !==edge3
                edgeB = e
                break
            end
        end
        g1t1 = multiplygammas(edge1.gamma, edge1.length)
        t3 = edge3.length
        if g1t1 != -1.0 && edge2.length != -1.0 && t3 != -1.0 # else do nothing: keep tA, tB, tC as is
            edgeB.length = addBL(edgeB.length, edge1.gamma * t3)
            edge3.length = t3 * edge2.gamma
            edge2.length = g1t1 + edge2.gamma * edge2.length
        end
        deletehybridedge!(net, edge1, false,unroot,false,false) # nofuse,unroot,multgammas,simplify
    end
    return true
end

"""
    adjacentedges(centeredge::Edge)

Vector of all edges that share a node with `centeredge`.
"""
function adjacentedges(centeredge::Edge)
    n = centeredge.node
    length(n) == 2 || error("center edge is connected to $(length(n)) nodes")
    @inbounds edges = copy(n[1].edge) # shallow copy, to avoid modifying the first node
    @inbounds for ei in n[2].edge
        ei === centeredge && continue # don't add the center edge again
        getOtherNode(ei, n[2]) === n[1] && continue # any parallel edge would have been in `edges` already
        push!(edges, ei)
    end
    return edges
end

#------------------------------------
function citation()
    bibfile = joinpath(@__DIR__, "..", "CITATION.bib")
    out = readlines(bibfile)
    println("Bibliography in bibtex format also in CITATION.bib")
    println(join(out,'\n'))
end
