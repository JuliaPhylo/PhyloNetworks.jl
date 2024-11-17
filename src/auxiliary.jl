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
    return (n1.number == n2.number && approxEq(n1.fvalue,n2.fvalue) && n1.intn1 == n2.intn1)
end

function isEqual(n1::Edge,n2::Edge)
    return (n1.number == n2.number && approxEq(n1.length,n2.length))
end

function isEqual(net1::HybridNetwork, net2::HybridNetwork)
    result = true
    result &= (net1.numtaxa == net2.numtaxa)
    result &= (net1.numnodes == net2.numnodes)
    result &= (net1.numedges == net2.numedges)
    ## result &= (net1.node == net2.node)
    ## result &= (net1.edge == net2.edge)
    result &= (net1.rooti == net2.rooti)
    result &= (net1.names == net2.names)
##    result &= (net1.hybrid == net2.hybrid)
    result &= (net1.numhybrids == net2.numhybrids)
##    result &= (net1.leaf == net2.leaf)
    result &= (net1.vec_float == net2.vec_float)
    result &= (net1.vec_int2 == net2.vec_int2)
    result &= (net1.intg1 == net2.intg1)
    result &= (net1.boolg1 == net2.boolg1)
    result &= (net1.vec_int3 == net2.vec_int3)
    result &= (net1.fscore == net2.fscore)
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

#= warning: node needs to be defined as hybrid before adding to a hybrid edge.
First, an edge is defined as hybrid, and then the nodes are added to it.
The edge length is *not* identifiable via edge.boole1, e.g. when the edge is
- external, or
- tree edge in "fromBadDiamondI"  via edge.boole2
- possibly, if in a bad 3- or 4-cycle via .booln2, .booln3, booln6
=#
function setNode!(edge::Edge, node::Node)
    size(edge.node,1)  !=  2 || error("vector of nodes already has 2 values");
    push!(edge.node,node);
    if size(edge.node,1) == 1
        if edge.hybrid
            edge.ischild1 = node.hybrid
        end
        edge.boole1 = !node.leaf # is edge identifiable?
    else
        if node.leaf
            !edge.node[1].leaf || error("edge $(edge.number) has two leaves")
            edge.boole1 = false;
        else
          if edge.hybrid
            if node.hybrid
                # @debug (edge.node[1].hybrid ? "hybrid edge $(edge.number) has two hybrid nodes" : "")
                edge.ischild1 = false;
            else
                edge.node[1].hybrid || error("hybrid edge $(edge.number) has no hybrid nodes");
                edge.ischild1 = true;
            end
          else #edge is tree
            if !edge.node[1].leaf
                if !node.hybrid && !edge.node[1].hybrid
                    edge.boole1 = !edge.boole2
                else
                    if node.hybrid && (node.booln2 || node.booln3 || node.booln6)
                        edge.boole1 = false
                    elseif edge.node[1].hybrid && (edge.node[1].booln2 ||edge.node[1].booln3 || edge.node[1].booln6)
                        edge.boole1 = false
                    else
                        edge.boole1 = true
                    end
                end
            else
                edge.boole1 = false
            end
          end
        end
    end
end

# warning: node needs to be defined as hybrid before adding to a hybrid edge.
#          First, an edge is defined as hybrid, and then the nodes are added to it.
#          If there is a leaf in node, then edge.boole1=false
#          to mean that its length is not identifiable.
function setNode!(edge::Edge,node::Array{Node,1})
    size(node,1) ==  2 || error("vector of nodes must have exactly 2 values")
    edge.node = node;
    if(edge.hybrid)
      if(node[1].hybrid)
          edge.ischild1 = true;
      else
          node[2].hybrid || error("hybrid edge without hybrid node");
          edge.ischild1 = false;
      end
    end
    if(edge.node[1].leaf || edge.node[2].leaf)
        edge.boole1 = false;
    else
        edge.boole1 = true;
    end
end

"""
    getroot(net)

Node used to root `net`. If `net` is to be considered as semi-directed or
unrooted, this root node is used to write the networks' Newick parenthetical
description or for network traversals.

See also: [`isrootof`](@ref)
"""
getroot(net::HybridNetwork) = net.node[net.rooti]

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
Assumes that the edge's direction is correct, meaning its field `ischild1` is
reliable (in sync with the rooting).

See also: [`getparent`](@ref), [`getchild`](@ref), [`isrootof`](@ref)
"""
isparentof(node::Node, edge::Edge) = node === getparent(edge)
@doc (@doc isparentof) ischildof
ischildof( node::Node, edge::Edge) = node === getchild(edge)

"""
    hassinglechild(node)

`true` if `node` has a single child, based on the edges' `ischild1` field;
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

*Warning*: these functions rely on correct edge direction, via their `ischild1` field.

See also:
[`getparent`](@ref),
[`getpartneredge`](@ref),
[`isparentof`](@ref),
[`hassinglechild`](@ref).
"""
getchild(edge::Edge) = edge.node[edge.ischild1 ? 1 : 2]
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

*Warning*: these functions use the field `ischild1` of edges.

See also: [`getchild`](@ref),
[`getpartneredge`](@ref).
"""
getparent(edge::Edge) = edge.node[edge.ischild1 ? 2 : 1]
@inline function getparent(node::Node)
    for e in node.edge
        if e.ismajor && ischildof(node, e)
            return getparent(e)
        end
    end
    error("could not find major parent of node $(node.number)")
end

@doc (@doc getparent) getparentminor
@inline function getparentminor(node::Node)
    for e in node.edge
        if !e.ismajor && node == getchild(e)
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
        if ee.ismajor && ischildof(n,ee)
            return ee
        end
    end
    error("node $(n.number) has no major parent")
end
@doc (@doc getparent) getparentedgeminor
@inline function getparentedgeminor(n::Node)
    for ee in n.edge
        if !ee.ismajor && n == ee.node[(ee.ischild1 ? 1 : 2)]
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

using the `ischild1` attribute of edges.
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
   node.booln1 = any(e -> e.hybrid, node.edge)
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
    net.numhybrids == length(net.hybrid) || error("numhybrids does not match to length of net.hybrid")
    net.numhybrids != 0 || return true
    return false
end

# function to push a Node in net.node and
# update numnodes and numtaxa
function pushNode!(net::Network, n::Node)
    push!(net.node,n);
    net.numnodes += 1;
    if(n.leaf)
        net.numtaxa += 1
        push!(net.leaf,n);
    end
    if(n.hybrid)
        pushHybrid!(net,n)
    end
end

# function to push an Edge in net.edge and
# update numedges
function pushEdge!(net::Network, e::Edge)
    push!(net.edge,e);
    net.numedges += 1;
end


# function to push a hybrid Node in net.hybrid and
# update numhybrids
function pushHybrid!(net::Network, n::Node)
    if(n.hybrid)
        push!(net.hybrid,n);
        net.numhybrids += 1;
    else
        error("node $(n.number) is not hybrid, so cannot be pushed in net.hybrid")
    end
end

"""
    deleteNode!(net::HybridNetwork, n::Node)

Delete node `n` from a network, i.e. removes it from
net.node, and from net.hybrid or net.leaf as appropriate.
Update attributes `numnodes`, `numtaxa`, `numhybrids`.


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
    net.numnodes -= 1;
    if net.rooti == index  # do not check containroot to save time in snaq!
        net.rooti = 1      # arbitrary
    elseif net.rooti > index
        net.rooti -= 1
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

Delete edge `e` from `net.edge` and update `net.numedges`.
If `part` is true, update the network's partition field.

*Warning*: if `part` is true (the default), then `net` is assumed to be of
level-1, with valid internal fields `e.inte1` (to track which cycle
`e` may be in) and valid `net.partition`, which then gets updated.
"""
function deleteEdge!(net::HybridNetwork, e::Edge; part::Bool=true)
    if part
        if e.inte1 == -1 && !e.hybrid && !isempty(net.partition) && !isTree(net)
            ind = whichPartition(net,e)
            indE = getIndex(e,net.partition[ind].edges)
            deleteat!(net.partition[ind].edges,indE)
        end
    end
    i = findfirst(x -> x===e, net.edge)
    i !== nothing || error("edge $(e.number) not in network: can't delete");
    deleteat!(net.edge, i);
    net.numedges -= 1;
end

"""
    removeHybrid!(net::Network, n::Node)

Delete a hybrid node `n` from `net.hybrid`, and update `net.numHybrid`.
The actual node `n` is not deleted. It is kept in the full list `net.node`.
Very internal function, used by [`deletehybridedge!`](@ref) and others.
"""
function removeHybrid!(net::Network, n::Node)
    n.hybrid || error("cannot delete node $(n.number) from net.hybrid because it is not hybrid")
    i = findfirst(x -> x===n, net.hybrid)
    i !== nothing || error("hybrid node $(n.number) not in the network's list of hybrids");
    deleteat!(net.hybrid, i);
    net.numhybrids -= 1;
end

# function to delete a leaf node in net.leaf
# and update numtaxa
function removeLeaf!(net::Network,n::Node)
    n.leaf || error("cannot delete node $(n.number) from net.leaf because it is not leaf")
    index = findfirst(no -> no === n, net.leaf)
    index !== nothing || error("leaf node $(n.number) not in network")
    deleteat!(net.leaf,index)
    net.numtaxa -= 1
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

"""
    printedges(net)
    printedges(io::IO, net)

Print information on the edges of a `HybridNetwork`
`net`: edge number, numbers of nodes attached to it, edge length, whether it's
a hybrid edge, its γ inheritance value, whether it's a major edge,
if it could contain the root (this field is not always updated, though)
and one more attribute pertaining to level-1 networks used in SNaQ:
in which cycle it is contained (-1 if no cycle).
"""
printedges(x) = printedges(stdout::IO, x)
function printedges(io::IO, net::HybridNetwork)
    if net.intg1 > 0
        println(io, "net has $(net.intg1) bad diamond I. Some γ and edge lengths t are not identifiable, although their γ * (1-exp(-t)) are.")
    end
    miss = ""
    println(io, "edge parent child  length  hybrid ismajor gamma   containroot i_cycle")
    for e in net.edge
        @printf(io, "%-4d %-6d %-6d ", e.number, getparent(e).number, getchild(e).number)
        if e.length==-1.0 @printf(io, "%-7s ", miss); else @printf(io, "%-7.3f ", e.length); end
        @printf(io, "%-6s %-7s ", e.hybrid, e.ismajor)
        if e.gamma==-1.0  @printf(io, "%-7s ", miss); else @printf(io, "%-7.4g ", e.gamma); end
        @printf(io, "%-11s %-7d\n", e.containroot, e.inte1)
    end
end


"""
    printnodes(net)
    printnodes(io, net)

Print information on the nodes of a `HybridNetwork` net: node number,
whether it's a leaf, whether it's a hybrid node, it's name (label),
its `intn1` field (for level-1 networks in SNaQ: number given to the cycle in
which the node might be, -1 if the node it *not* in a cycle cycle),
and the list of edges attached to it, by their numbers.
"""
printnodes(x) = printnodes(stdout::IO, x)
function printnodes(io::IO, net::Network)
    namepad = max(4, maximum(length.([n.name for n in net.node])))
    println(io, "node leaf  hybrid ", rpad("name", namepad), " i_cycle edges'numbers")
    for n in net.node
        @printf(io, "%-4d %-5s %-6s ", n.number, n.leaf, n.hybrid)
        print(io, rpad(n.name,namepad))
        @printf(io, " %-7d", n.intn1)
        for e in n.edge
            @printf(io, " %-4d", e.number)
        end
        print(io, "\n")
    end
end

"""
    hybridEdges(node::Node)

Return the 3 edges attached to `node` in a specific order [e1,e2,e3].
**Warning**: assume a level-1 network with up-to-date node fields
`booln1` (tracking whether the node is incident to a hybrid edge
and edge) and field `intn1` (tracking the number given to the cycle in which
the node might be).

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
    if node.hybrid
        hybmajor = nothing;
        hybminor = nothing;
        tree = nothing;
        for e in node.edge
            (e.hybrid && e.ismajor) ? hybmajor = e : nothing
            (e.hybrid && !e.ismajor) ? hybminor = e : nothing
            !e.hybrid ? tree = e : nothing
        end
        return hybmajor, hybminor, tree
    elseif node.booln1
        hybrid = nothing;
        treecycle = nothing;
        tree = nothing;
        for e in node.edge
            (e.hybrid) ? hybrid = e : nothing
            (!e.hybrid && e.inte1 != -1) ? treecycle = e : nothing
            (!e.hybrid && e.inte1 == -1) ? tree = e : nothing
        end
        return hybrid, treecycle, tree
    else
        #@warn "node $(node.number) is not hybrid $(node.hybrid) nor tree with hybrid edges (booln1) $(node.booln1), return the node.edge in order, unless a leaf is attached, then the edge attached to leaf is last";
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
#          node.booln1 is set to false
#          assuming any tree node can only have one
#          one hybrid edge
function removeEdge!(node::Node, edg::Edge)
    index = findfirst(x -> x === edg, node.edge)
    index !== nothing || error("edge $(edg.number) not in node $(node.number)")
    deleteat!(node.edge,index)
    node.booln1 = any(e -> e.hybrid, node.edge)
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
    setgamma!(Edge, new γ, change_other=true)

Set inheritance probability γ for an edge, which must be a hybrid edge.
The new γ needs to be in [0,1]. The γ of the "partner" hybrid edge is changed
accordingly, to 1-γ. The field `ismajor` is also changed accordingly.
If the new γ is approximately 0.5, `Edge` is set to the major parent,
its partner is set to the minor parent.

If `net` is a HybridNetwork object, `printedges(net)` will show the list of edges
and their γ's. The γ of the third hybrid edge (say) can be changed to 0.2 with
`setgamma!(net.edge[3],0.2)`.
This will automatically set γ of the partner hybrid edge to 0.8.

The last argument is true by default. If false: the partner edge is not updated.
This is useful if the new γ is 0.5, and the partner's γ is already 0.5,
in which case the `ismajor` attributes can remain unchanged.

See also [`PhyloNetworks.setmultiplegammas!`](@ref)
"""
function setgamma!(edge::Edge, new_gamma::Float64, changeOther::Bool=true)
    new_gamma >= 0.0 || error("gamma has to be positive: $(new_gamma)")
    new_gamma <= 1.0 || error("gamma has to be less than 1: $(new_gamma)")
    edge.hybrid || error("cannot change gamma in a tree edge");
    node = getchild(edge) # child of hybrid edge
    node.hybrid || @warn "hybrid edge $(edge.number) not pointing at hybrid node"
    # @debug (node.booln2 ? "bad diamond I situation: gamma not identifiable" : "")
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
        edge.ismajor = new_ismajor
        e2.gamma = 1.0 - new_gamma
        e2.ismajor = !new_ismajor
    else
        if onehalf # who is major is arbitrary: so we pick what's consistent with the partner
            edge.ismajor = !e2.ismajor
        else
            edge.ismajor = new_ismajor
        end
    end
    return nothing
end

"""
    setmultiplegammas!(edges::Vector{Edge}, γs::Vector{Float64})

Set the inheritance of the ith edge to the ith γ value,
calling [`setgamma!`](@ref).
"""
@inline function setmultiplegammas!(edges::Vector{Edge}, gammas::Vector{Float64})
    for (e,g) in zip(edges, gammas)
        setgamma!(e, g)
    end
end

"""
    remove_edgelengthsgammas!(net::HybridNetwork)

Reset all edge lengths and all hybrid edge γs to be missing (coded as -1.0).
"""
function remove_edgelengthsgammas!(net::HybridNetwork)
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
    getnodeheights(net, checkpreorder::Bool=true)
    getnodeheights!(net, checkpreorder::Bool=true)

Vector of node heights, that is: the distance of each node to the root.
An error is thrown if the network is not time-consistent.
A network is time-consistent if, for any node `v`, all paths from the root to `v`
have the same length. (It is sufficient to check this condition at hybrid nodes).
Ultrametricity is not assumed: tips need not all be at the same distance from the root.
If `checkpreorder=false`, assumes the network has already been preordered
with [`preorder!`](@ref).

If a tree edge has a missing length (coded as -1), both functions throw an error.
In general, there may be an exponential number of ways to assign tree edge
lengths that make the network time-consistent.

`getnodeheights` sends a warning upon finding a missing hybrid edge length,
otherwises proceeds as `getnodeheights!` but without modifying the network.
`getnodeheights!` will attempt to assign values to missing lengths, for hybrid
edges only, so as to make the network time-consistent.

If a hybrid edge `e` has a missing length, `getnodeheights!` proceeds as follows
at its child hybrid node `h`:
- If all of `h`'s parent edges lack a length: the shortest non-negative lengths
  are assigned to make the network time-consistent at `h`. In particular, one of
  the partner edges is assigned length 0, and `h` is made as old as possible,
  that is, as close to the root as possible: the reticulation is "zipped-up".
- Otherwise: the length of `e` is set to the unique value that makes the network
  time-consistent at `h`, based on the partner edge's length.
  If this value is negative, then an error is thrown.

Output: vector of node heights, one per node, in the same order as in
`net.vec_node`.

See also: [`istimeconsistent`](@ref) and [`getnodeheights_average`](@ref).

Examples:

```jldoctest
julia> net = readnewick("(((C:1,(A:1)#H1:1.5::0.7):1,(#H1:0.3::0.3,E:2.0):2.2):1.0,O:5.2)root;");

julia> # using PhyloPlots; plot(net, useedgelength=true, showedgelength=true, shownodenumber=true); # to see

julia> nodeheight = getnodeheights(net)
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

julia> [node.number => (height, node.name) for (height,node) in zip(nodeheight, net.vec_node)]
9-element Vector{Pair{Int64, Tuple{Float64, String}}}:
 -2 => (0.0, "root")
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
getnodeheights(net::HybridNetwork, checkpreorder::Bool=true) =
    _getnodeheights(net, false, timeinconsistency_error, checkpreorder)[2]
@doc (@doc getnodeheights) getnodeheights!
getnodeheights!(net::HybridNetwork, checkpreorder::Bool=true) =
    _getnodeheights(net, true, timeinconsistency_error, checkpreorder)[2]

"""
    getnodeheights_average(net, checkpreorder::Bool=true; warn=true)

Vector of average node heights, that is: the average distance from the root to
each node. The average is a weighted average with weights taken to be the
hybrid edges' inheritance values γ, if available. Equal weights are used at
hybrid nodes with some parents lacking a γ inheritance value (with a warning).

missing edge lengths:
- An error is thrown if a tree edge has a missing edge length.
- If all parent hybrid edges have missing lengths at a given hybrid node, then
  the hybrid node is assumed to be as close to the root as possible, that is,
  the reticulation is assumed "zipped-up" with one of its hybrid edges of length 0.
- If some but not all parent hybrid edges have a missing length, then the
  average node height is calculated based on the non-missing parents only.
  If the hybrid node height turns out to be lower than one of the parent's height
  (such that some missing length would need to be negative) then a warning is
  issued.

A warning is issued, unless `warn=false`, if the network is not time-consistent.

See also: [`istimeconsistent`](@ref), [`getnodeheights`](@ref), and `getnodeheights_majortree`](@ref).
"""
function getnodeheights_average(
    net::HybridNetwork,
    checkpreorder::Bool=true;
    warn::Bool=true
)
    (isTC, nh) = _getnodeheights(net, false, timeinconsistency_average, checkpreorder)
    warn && !isTC && @warn "the network is not time consistent"
    return nh
end

"""
    getnodeheights_majortree(net, checkpreorder::Bool=true; warn=true)

Vector of node heights from the major tree, that is: the distance from the root to
each node when considering the major tree for node heights. 

missing edge lengths:
- An error is thrown if a tree edge has a missing edge length.
- If all parent hybrid edges have missing lengths at a given hybrid node, then
  the hybrid node is assumed to be as close to the root as possible, that is,
  the reticulation is assumed "zipped-up" with one of its hybrid edges of length 0.
- If a major hybrid edge has a missing length, then the hybrid node height will
  be calculated using the node height and edge length of the minor parent with
  the largest inheritance γ (with a warning). If the major hybrid edge lacks a length and
  all non-missing minor edges lack an inheritance γ or have the same value,
  then an error is thrown.

A warning is issued, unless `warn=false`, if the network is not time-consistent.

See also: [`istimeconsistent`](@ref), [`getnodeheights`](@ref) and
[`getnodeheights_average`](@ref).

```jldoctest
#node heights of time-consistent networks are the same 
julia> consistent_net = readnewick("((A:2.5,#H1:1.5::0.4):0.25,(C:1.5,(B:1)#H1:0.5::0.6):1.25);");

julia> heights = getnodeheights(consistent_net)
7-element Vector{Float64}:
 0.0
 1.25
 2.75
 0.25
 1.75
 2.75
 2.75

julia> heights_average = getnodeheights_average(consistent_net);

julia> heights_major = getnodeheights_majortree(consistent_net);

julia> heights == heights_average == heights_major  
true
```

```jldoctest
#inconsistent networks give different results
julia> inconsistent_net = readnewick("((A:2.5,#H1:1.5::0.4):0.25,(C:1.5,(B:1)#H1:2.5::0.6):1.25);");

julia> getnodeheights_average(inconsistent_net;warn=false)
7-element Vector{Float64}:
 0.0
 1.25
 2.75
 0.25
 2.95
 3.95
 2.75

julia> getnodeheights_majortree(inconsistent_net;warn=false) 
7-element Vector{Float64}:
 0.0
 1.25
 2.75
 0.25
 3.75
 4.75
 2.75

```
"""
function getnodeheights_majortree(net::HybridNetwork, checkpreorder::Bool=true; warn::Bool=true)
    (isTC, nh) = _getnodeheights(net, false, timeinconsistency_majortree, checkpreorder)
    warn && !isTC && @warn "the network is not time consistent"
    return nh
end

"""
    istimeconsistent(net, checkpreorder::Bool=true)

True (resp. false) if `net` network is (resp. is not) time-consistent.
A network is time-consistent if for any node `v`, all paths from the root to `v`
have the same length.
It is sufficient to check this condition at nodes `v` that are hybrid nodes.

See also [`getnodeheights`](@ref) and [`getnodeheights_average`](@ref).
"""
istimeconsistent(net::HybridNetwork, checkpreorder::Bool=true) =
    _getnodeheights(net, false, timeinconsistency_check, checkpreorder)[1]


"""
    _getnodeheights(net::HybridNetwork, fixmissing::Bool,
                    inconsistencyhandler::Function, checkpreorder::Bool=true)

Helper to determine time-consistency and calculate node heights
(distance from the root), used by [`getnodeheights`](@ref) for example.

output: `(isconsistent, nodes_distance_from_root)`

Arguments:

- `fixmissing`:
  * if `false`, any missing hybrid edge length will cause a warning, the network
    is *not* modified, and the best-case is assumed to determine time-consistency
    (as explained below)
  * if `true`, will attempt to find values for missing lengths of *hybrid* edges,
    if any, to make the network time-consistent, placing hybrid nodes as close
    to the root as possible (as this gives most chances to find a time-consistent
    assignment of all missing hybrid edge lengths).
    If there is a time-consistent assignment, then the network is modified
    (with missing hybrid edge lengths set to time-consistent values).
    Otherwise, the network is not modified.

  This option is passed to `update_getnodeheights_hybrid!` that handles
  1 hybrid node at a time.

- `inconsistencyhandler`: function to check & handle time-consistency as desired
  as a given hybrid node `h`, and to decide if the traversal should continue.
  It should take as input:
  * a vector of ≥1 candicate heights for `h` from parent edges with non-missing
    length, and
  * a vector of ≥0 heights of parent nodes whose child edge to `h` has no length
  * `isconsistent`: a boolean that is modified to `false` if the network is
    not time-consistent at `h` (unless an error is thrown anyway!).
  
  This handler function decides what to do if the candidate heights are not
  all equal (the network is time-inconsistent), and if the values to be
  assigned to missing edge lengths would be negative. Its output should be:
  `(keepgoing_boolean, hybrid_node_height)`.

  Examples:
  [`timeinconsistency_error`](@ref) is conservative and throws an
  error in both cases.
  [`timeinconsistency_average`](@ref) is lenient: only throws warnings,
  but keeps going and returns γ-weighted average node heights
"""
function _getnodeheights(
    net::HybridNetwork,
    fixmissing::Bool,
    inconsistencyhandler::Function,
    checkpreorder::Bool=true,
)
    checkpreorder && preorder!(net)
    missing_e = Tuple{Edge,Float64}[] # vector of (edge, fixed_length) for edges with missing length
    isconsistent = Ref(true)
    rootdistance = traversal_preorder(
        net.vec_node,
        getnodeheights_init,
        traversalupdate_default!, # nothing to do at the root
        update_getnodeheights_tree!,
        update_getnodeheights_hybrid!,
        missing_e,
        isconsistent,
        inconsistencyhandler)
    # assign values to fix missing edge lengths: ! modifies net !
    if fixmissing
        (x-> x[1].length=x[2]).(missing_e)
    else
        isempty(missing_e) || @warn "some hybrid edge length is missing"
    end
    return (isconsistent[], rootdistance)
end

function getnodeheights_init(nodes::Vector{Node}, params...)
    n = length(nodes)
    return zeros(Float64,n)
end

function update_getnodeheights_tree!(
    rootdistance::Vector{Float64},
    i::Int,
    parind::Int,
    paredge::Edge,
    params...
) 
    if paredge.length == -1.0 # interpreted as missing
        paredge.hybrid && @error("weird, hybrid edge at tree node")
        error("Edge $(paredge.number) has a missing edge length: node height cannot be determined")
    end
    rootdistance[i] = rootdistance[parind] + paredge.length 
    return true
end

function update_getnodeheights_hybrid!(
    rootdistance::Vector{Float64},
    i::Int,
    parinds::Vector{Int},
    paredges::Vector{Edge},
    missing_e::Vector{Tuple{Edge,Float64}},
    isconsistent::Ref{Bool},
    inconsistencyhandler::Function,
)
    keepgoing = true
    candidate_nodeheight = Float64[]
    missingparent_height = Float64[]
    missingparent_j = Int[]
    nonmissingparent_j = Int[]
    for (pj,pj_ind) in enumerate(parinds)
        if paredges[pj].length == -1 # missing parent edge length
            push!(missingparent_j, pj)
            push!(missingparent_height, rootdistance[pj_ind])
        else
            push!(nonmissingparent_j, pj)
            push!(candidate_nodeheight, rootdistance[pj_ind] + paredges[pj].length)
        end
    end
    if isempty(candidate_nodeheight) # all edge lengths are missing: then we can
        # set them in a time-consistent way, zipped up: hybrid as close to the root as possible
        nodehght = maximum(missingparent_height)
    else # ≥1 length: check & handle time-consistency as desired (e.g. take average)
        keepgoing, nodehght = inconsistencyhandler(
            candidate_nodeheight, missingparent_height,
            isconsistent, paredges, nonmissingparent_j)
    end
    rootdistance[i] = nodehght
    for (pj, ph) in zip(missingparent_j, missingparent_height)
        push!(missing_e, (paredges[pj], nodehght - ph))
    end
    return keepgoing # stops the traversal early is not time-consistent
end

"""
    timeinconsistency_error(
        candidate_nodeheights,
        missingparent_heights,
        args...;
        atol::Real=1e-8, rtol::Real=√eps(Float64))
    timeinconsistency_check

Check that all candidate node heights are approximately equal to one another,
and that this shared value `nodeheight` is higher (farther from the root) than the
height of parents connected by edge of missing length `missingparent_heights`
to ensure that these edge lengths would be assigned non-negative values.

If any of these conditions is not met, `timeinconsistency_error` throws an error.
Otherwise, it returns `(true, nodeheight)` where `true` means that
the network is (or could be) time-consistent at the node being considered.
`timeinconsistency_check` returns `(is_timeconsistent, nodeheight)` but does
*not* throw an error the `is_timeconsistent` if false (for either reason).

Assumption: `candidate_nodeheight` is not empty, that is, the node has at least
one parent edge with a non-missing length.

See also [`_getnodeheights`](@ref PhyloNetworks._getnodeheights)
"""
function timeinconsistency_error(
    candidate_nodeheight::AbstractVector{T},
    missingparent_height::AbstractVector{T},
    args...;
    atol::Real=1e-8, # more lenient than default 0 in isapprox
    rtol::Real=√eps(T),
) where T<:Real
    min_nh, max_nh = extrema(candidate_nodeheight)
    length(candidate_nodeheight) == 1 ||
        isapprox(min_nh, max_nh; atol=atol, rtol=rtol) ||
        error("the network is not time consistent. paths of different lengths: $candidate_nodeheight")
    if any(missingparent_height .> max_nh) # may be empty vector
        error("""a missing edge length would be need to be set to a negative value
        to make the network time-consistent""")
    end
    return (true, max_nh)
end

@doc (@doc timeinconsistency_error) timeinconsistency_check
function timeinconsistency_check(
    candidate_nodeheight::AbstractVector{T},
    missingparent_height::AbstractVector{T},
    isconsistent::Ref{Bool},
    args...;
    atol::Real=1e-8, # more lenient than default 0 in isapprox
    rtol::Real=√eps(T),
) where T<:Real
    min_nh, max_nh = extrema(candidate_nodeheight)
    timecons = length(candidate_nodeheight) == 1 ||
        isapprox(min_nh, max_nh; atol=atol, rtol=rtol)
    if any(missingparent_height .> max_nh) # may be empty vector
        timecons = false
    end
    isconsistent[] &= timecons
    return (timecons, max_nh) # stop the traversal if not time consistent
end

"""
    timeinconsistency_average(
        candidate_nodeheights,
        missingparent_heights,
        isconsistent::Ref{Bool},
        parent_edges,
        nonmissingparent_j;
        atol::Real=1e-8, rtol::Real=√eps(Float64))

Calculate the γ-weighted average node height of a given hybrid node `h`, based
on the candidate node heights from its parents with non-missing edge lengths.
- If some of these parent edges have a missing γ, then equal weights are used
  and a warning is issued.
- If the hybrid node's average height (calculate from non-missing lengths)
  turns out to be lower than one of the parent's height with a missing length
  (such that this parent edge would need to be assigned a negative value)
  then a warning is issued.

Outcome:
- update `isconsistent` to false if the candidate node heights are not all equal
  or if some missing edge length would have to be assigned a negative value
  to make the network time-consistent
- returns `(true, nodeheight)`

Assumption: `candidate_nodeheight` is not empty, that is, the node has at least
one parent edge with a non-missing length.

See also [`_getnodeheights`](@ref) and [`getnodeheights_average`](@ref)
"""
function timeinconsistency_average(
    candidate_nodeheight::AbstractVector{T},
    missingparent_height::AbstractVector{T},
    isconsistent::Ref{Bool},
    paredges::Vector{Edge},
    nm_ind::AbstractVector;
    atol::Real=1e-8, # more lenient than default 0
    rtol::Real=√eps(T),
) where T<:Real
    min_nh, max_nh = extrema(candidate_nodeheight)
    timecons = length(candidate_nodeheight) == 1 ||
        isapprox(min_nh, max_nh; atol=atol, rtol=rtol)
    if timecons # if time-consistent: don't calculate the weighted average
        nh = max_nh
    else # weighted average, with equal weights if some γ's are missing
        if any(j -> paredges[j].gamma == -1, nm_ind) # ≥ 1 missing γ
            @warn "missing γ: will use equal weights"
            nh = sum(candidate_nodeheight) / length(candidate_nodeheight)
        else
            nh = zero(T); gamma_sum = zero(T)
            for (cnh, pj) in zip(candidate_nodeheight, nm_ind)
                nh += cnh * paredges[pj].gamma
                gamma_sum += paredges[pj].gamma
            end
            nh /= gamma_sum
        end
    end
    if any(missingparent_height .> nh) # may be empty vector
        @warn """a missing edge length would be need to be set to a negative value
        for these average node heights"""
        timecons = false
    end
    isconsistent[] &= timecons
    return (true, nh) # keep going, even if the network is inconsistent
end

"""
    timeinconsistency_majortree(
        candidate_nodeheights,
        missingparent_heights,
        isconsistent::Ref{Bool},
        parent_edges,
        nonmissingparent_j;
        atol::Real=1e-8, rtol::Real=√eps(Float64))

Calculate node height of a given hybrid node `h`, based on the its major parent
node height, if its major parent edge has a non-missing length. If missing, then
the non-missing edge length with the largest γ is used and a warning is issued.
If all parent edges have missing γ values then an error is thrown.

Outcome:
- update `isconsistent` to false if the candidate node heights are not all equal
  or if some missing edge length would have to be assigned a negative value
  to make the network time-consistent
- returns `(true, nodeheight)`

Assumption: `candidate_nodeheight` is not empty, that is, the node has at least
one parent edge with a non-missing length.

See also [`_getnodeheights`](@ref) and [`getnodeheights_average`](@ref)
"""
function timeinconsistency_majortree(
    candidate_nodeheight::AbstractVector{T},
    missingparent_height::AbstractVector{T},
    isconsistent::Ref{Bool},
    paredges::Vector{Edge},
    nm_ind::AbstractVector;
    atol::Real=1e-8, # more lenient than default 0
    rtol::Real=√eps(T),
) where T<:Real
    min_nh, max_nh = extrema(candidate_nodeheight)
    timecons = length(candidate_nodeheight) == 1 ||
        isapprox(min_nh, max_nh; atol=atol, rtol=rtol)
    maj_cantidate = findfirst(x->x.ismajor, paredges[nm_ind])
    if !isnothing(maj_cantidate) #  the major edge is among candidates
        nh = candidate_nodeheight[maj_cantidate]
    else # find candidate with largest γ
        @warn "major hybrid edge missing a length. Using non-missing minor edge with largest γ"
        gammas = (x-> x.gamma).(paredges[nm_ind])
        max_gamma = maximum(gammas)
        max_gamma == -1 && error("major edge lacks a length and all non-missing edges lack a γ")
        cantidate_ind = findall(gammas .== max_gamma)
        length(cantidate_ind) > 1 && error("major edge lacks a length and two alternative edges have the same γ")
        nh = candidate_nodeheight[cantidate_ind[1]]
    end
    if any(missingparent_height .> nh) # may be empty vector
        @warn """a missing edge length would be need to be set to a negative value
        to use the cantidate node height of the major edge"""
        timecons = false
    end
    isconsistent[] &= timecons
    return (true, nh) # keep going, even if the network is inconsistent
end


function numTreeEdges(net::HybridNetwork)
    2*net.numtaxa - 3 + net.numhybrids
end

function numIntTreeEdges(net::HybridNetwork)
    2*net.numtaxa - 3 + net.numhybrids - net.numtaxa
end


# function to get the partition where an edge is
# returns the index of the partition, or error if not found
# better to return the index than the partition itself, because we need the index
# to use splice and delete it from net.partition later on
# cycle: is the number to look for partition on that cycle only
function whichPartition(net::HybridNetwork,edge::Edge,cycle::Integer)
    !edge.hybrid || error("edge $(edge.number) is hybrid so it cannot be in any partition")
    edge.inte1 == -1 || error("edge $(edge.number) is in cycle $(edge.inte1) so it cannot be in any partition")
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
    edge.inte1 == -1 || error("edge $(edge.number) is in cycle $(edge.inte1) so it cannot be in any partition")
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
    setlength!(edge::Edge, new_length)

Assign new length to `edge`. `new_length` should be non-negative,
or `missing` (or -1, interpreted as missing).
"""
@inline function setlength!(edge::Edge, new_length)
    if ismissing(new_length) || new_length == -1
        edge.length = -1
    else
        new_length >= 0.0 || error("edge length must be non negative: $(new_length)")
        edge.length = new_length
    end
    return nothing
end

"""
    setlengths!(edges::Vector{Edge}, lengths::AbstractVector)

Assign new lengths to a vector of `edges`.
Checks that the new edge lengths are non-negative or `missing` (or -1 to be
interpreted as missing).
"""
@inline function setlengths!(edges::Vector{Edge}, lengths::AbstractVector)
    for (e,l) in zip(edges, lengths)
        setlength!(e, l)
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
Specifically: delete its child node (as determined by `ischild1`) and connect
all edges formerly incident to this child node to the parent node of `edge`,
thus creating a new polytomy, unless the child was of degree 2.

Warning: it's best for `ischild1` to be in sync with the root for this. If not,
the shrinking may fail (if `edge` is a tree edge but its "child" is a hybrid)
or the root may change arbitrarily (if the child of `edge` is the root).

Output: true if the remaining node (parent of `edge`) becomes a hybrid node with
more than 1 child after the shrinking; false otherwise (e.g. no polytomy was
created, or the new polytomy is below a tree node)
"""
function shrinkedge!(net::HybridNetwork, edge2shrink::Edge)
    edge2shrink.hybrid && error("cannot shrink hybrid edge number $(edge2shrink.number)")
    cn = getchild(edge2shrink)
    cn.hybrid && error("cannot shrink tree edge number $(edge2shrink.number): its child node is a hybrid. run directedges! ?")
    pn = getparent(edge2shrink)
    isexternal(edge2shrink) &&  # (isleaf(cn) || isleaf(pn)) &&
      error("won't shrink edge number $(edge2shrink.number): it is incident to a leaf")
    removeEdge!(pn,edge2shrink)
    empty!(edge2shrink.node) # should help gc
    for ee in cn.edge
        ee !== edge2shrink || continue
        cn_index = findfirst(x -> x === cn, ee.node)
        ee.node[cn_index] = pn # ee.ischild1 remains synchronized
        push!(pn.edge, ee)
    end
    pn.booln1 = any(e -> e.hybrid, pn.edge)
    empty!(cn.edge) # should help to garbage-collect cn
    deleteEdge!(net, edge2shrink; part=false)
    deleteNode!(net, cn)
    badpolytomy = false
    if pn.hybrid # count the number of pn's children, without relying on ischild1 of tree edges
        nc = sum((!e.hybrid || getchild(e) !== pn) for e in pn.edge)
        badpolytomy = (nc > 1)
    end
    return badpolytomy
end

@doc raw"""
    shrink2cycles!(net::HybridNetwork, unroot::Bool=false)

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
function shrink2cycles!(net::HybridNetwork, unroot::Bool=false)
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
    shrink3cycles!(net::HybridNetwork, unroot::Bool=false)

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
function shrink3cycles!(net::HybridNetwork, unroot::Bool=false)
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
            edge1.ismajor = g1tilde >= 0.5
            edge2.ismajor = !edge1.ismajor
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
