


"""
    rootatnode!(HybridNetwork, nodeNumber::Integer; index::Bool=false)
    rootatnode!(HybridNetwork, Node)
    rootatnode!(HybridNetwork, nodeName::AbstractString)

Root the network/tree object at the node with name 'nodeName' or
number 'nodeNumber' (by default) or with index 'nodeNumber' if index=true.
Attributes ischild1 and containroot are updated along the way.
Use `plot(net, shownodenumber=true, showedgelength=false)` to
visualize and identify a node of interest.
(see package [PhyloPlots](https://github.com/juliaphylo/PhyloPlots.jl))

Return the network.

Warnings:

- If the node is a leaf, the root will be placed along
  the edge adjacent to the leaf. This might add a new node.
- If the desired root placement is incompatible with one or more hybrids, then

  * the original network is restored with its old root and edges' direction.
  * a RootMismatch error is thrown.

See also: [`rootonedge!`](@ref), and [`suppressroot!`](@ref) to undo.
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

function rootatnode!(net::HybridNetwork, nodeNumber::Integer; index::Bool=false)
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
            rootatnode!(net, pn.number)
        else # otherwise, create a new node between leaf and its parent
            rootonedge!(net,net.node[ind].edge[1])
        end
    else
        rootsaved = net.rooti
        net.rooti = ind
        try
          directedges!(net)
        catch e
          if isa(e, RootMismatch) # new root incompatible with hybrid directions: revert back
            net.rooti = rootsaved
            directedges!(net)
          end
          throw(RootMismatch("""the desired root is below a reticulation,
                                reverting to old root position."""))
        end
        if (net.rooti != rootsaved && length(net.node[rootsaved].edge)==2)
            fuseedgesat!(rootsaved,net) # remove old root node if degree 2
        end
        return net
    end
end


"""
    rootonedge!(HybridNetwork, edgeNumber::Integer; index::Bool=false)
    rootonedge!(HybridNetwork, Edge)

Root the network/tree along an edge with number `edgeNumber` (by default)
or with index `edgeNumber` if `index=true`.
Attributes `ischild1` and `containroot` are updated along the way.

This adds a new node and a new edge to the network.
Use `plot(net, showedgenumber=true, showedgelength=false)` to
visualize and identify an edge of interest.
(see package [PhyloPlots](https://github.com/juliaphylo/PhyloPlots.jl))

See also: [`rootatnode!`](@ref), and [`suppressroot!`](@ref) to undo.
"""
function rootonedge!(net::HybridNetwork, edge::Edge; kwargs...)
    rootonedge!(net, edge.number, index=false; kwargs...)
end

function rootonedge!(net::HybridNetwork, edgeNumber::Integer; index::Bool=false)
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
    rootsaved = net.rooti
    breakedge!(net.edge[ind],net) # returns new node, new edge (last ones pushed)
    net.rooti = length(net.node)   # index of new node: was the last one pushed
    try
      directedges!(net)
    catch e
      if isa(e, RootMismatch) # new root incompatible with hybrid directions: revert back
        fuseedgesat!(net.rooti,net) # reverts breakedge!
        net.rooti = rootsaved
        directedges!(net)
      end
      throw(RootMismatch("""the desired root is below a reticulation,
                                reverting to old root position."""))
    end
    if (net.rooti != rootsaved && length(net.node[rootsaved].edge)==2)
        fuseedgesat!(rootsaved,net) # remove old root node if degree 2
    end
    return net
end

"""
    breakedge!(edge::Edge, net::HybridNetwork;
               lengthratio::Union{Nothing,Real}=nothing,
               lengths=nothing)

Break an edge into 2 edges, creating a new node of degree 2.
If the starting edge was:
```
n0  --edge-->  n1
```
then we get this, where the newly-created edge `newedge` is the "top" portion
of the original starting edge:
```
n0  --newedge-->  newnode  --edge-->  n1
```
These new node & edge are pushed last in `net.node` and `net.edge`.

Output: `(newnode, newedge)`.

`ischild1` and `containroot` are updated, but not internal fields
(e.g. not those used by SNaQ for level-1 networks).

Branch lengths are assigned with either one of 2 options. The default is
`lengthratio=0.5`, which means that `edge` is broken into 2 equal halves (unless
its original length is missing, in which case the new edge length is set to be
missing also).

1. `lengthratio`: number in [0,1]. If the original length of `edge` was `l`,
   the `newedge` is assigned length `lengthratio × l` and
   `edge` has its length updated to be `(1-lengthratio) × l`. In other words:
   - the original distance between `n0` and `n1` is preserved,
   - a length ratio of 0 corresponds to placing `newnode` close to `n0`, and
   - a length ratio of 1 corresponds to placing `newnode` close to `n1`.
2. `lengths`: 1 length (new) or 2 lengths (new,old), to be assigned to
   `newedge` (first one) and `edge` (second one, if any). If only the new edge
   length is provided, then the length of `edge` is unchanged.
   To ask for a length to be missing, use value `missing` (see example below).

**Note**: If a number outside of [0,1] is used in `lengthratio` or negative values
are used in `lengths`, an error message will appear but the network *will*
still be modified according to the values given.

# examples

```jldoctest
julia> net = readnewick("(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));");

julia> length(net.node)
19

julia> net.edge[4] # edge 4 goes from node -8 to 3
PhyloNetworks.EdgeT{PhyloNetworks.Node}:
 number:4
 length:-1.0
 attached to 2 node(s) (parent first): -8 3


julia> newnode, newedge = PhyloNetworks.breakedge!(net.edge[4], net);

julia> length(net.node) # one more than before
20

julia> newedge # new edge 21 goes from node -8 and 11 (new)
PhyloNetworks.EdgeT{PhyloNetworks.Node}:
 number:21
 length:-1.0
 attached to 2 node(s) (parent first): -8 11


julia> net.edge[4] # original edge 4 now goes from node 11 (new) to 3
PhyloNetworks.EdgeT{PhyloNetworks.Node}:
 number:4
 length:-1.0
 attached to 2 node(s) (parent first): 11 3


julia> writenewick(net) # note extra pair of parentheses around S1
"(((S8,S9),((((S4,(S1)),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));"
```

Another example, this time with edge lengths:
```jldoctest
julia> net = readnewick("((C:1.1,A:2.0)i1,B);"); # edge 1: i1 --1.1--> C

julia> const PN = PhyloNetworks; # just to write less later

julia> PN.breakedge!(net.edge[1],net; lengthratio=.2);

julia> writenewick(net, round=true) # now i1 --0.22--> newnode --0.88--> C
"((A:2.0,(C:0.88):0.22)i1,B);"

julia> PN.breakedge!(net.edge[4],net; lengths=(5,1)); # breaks root-->B

julia> writenewick(net, round=true) # into: root --5--> newnode --1--> B
"((A:2.0,(C:0.88):0.22)i1,(B:1.0):5.0);"

julia> PN.breakedge!(net.edge[4],net; lengths=[missing]); # breaks --1--> B

julia> writenewick(net, round=true) # into: --?--> newnode --1--> B
"((A:2.0,(C:0.88):0.22)i1,((B:1.0)):5.0);"

```

See also: [`fuseedgesat!`](@ref)
"""
function breakedge!(
    edge::Edge,
    net::HybridNetwork;
    lengthratio::Union{Nothing,Real}=nothing,
    lengths=nothing,
)
    useoldlength = isnothing(lengths)
    if !useoldlength
        isnothing(lengthratio) ||
            error("cannot use both options lengthratio=$lengthratio and lengths=$lengths")
        nl = length(lengths)
        nl > 0 || error("'lengths' needs to be of length 1 or more")
        for ell in lengths
            ismissing(ell) && continue
            isa(ell, Real) || error("new length $ell must be a real number")
            ell ≥ 0 || ell == -1 ||
                @error("new length $ell should be ≥ 0 or 'missing' (-1 is interpreted as missing)")
        end
    else
        if isnothing(lengthratio)
            lengthratio = 1/2
        else
            0 ≤ lengthratio ≤ 1 ||
                @error "length ratio $lengthratio should be in [0,1]"
        end
    end
    pn = getparent(edge) # parent node
    # new child edge = old edge, same hybrid attribute
    removeEdge!(pn,edge)
    removeNode!(pn,edge)
    max_edge = maximum(e.number for e in net.edge)
    max_node = maximum(n.number for n in net.node)
    # new parent edge: length=-1,hybrid=false,y=0,z=1,γ=1, ischild1=true,
    newedge = Edge(max_edge+1, -1., false, 0.,1., 1., Node[],true,
                   true,-1, edge.containroot, true,false)
    newnode = Node(max_node+1,false,false,[edge,newedge]) # tree node
    setNode!(edge,newnode) # newnode comes 2nd, and parent node along 'edge'
    edge.ischild1 = true
    setNode!(newedge,newnode) # newnode comes 1st in newedge: ischild1 true
    setEdge!(pn,newedge)
    setNode!(newedge,pn) # pn comes 2nd in newedge
    if useoldlength
        oldlen = edge.length
        if oldlen != -1.0
            len_parent = lengthratio * oldlen
            newedge.length = len_parent
            edge.length = oldlen - len_parent
        end
    else
        ell = lengths[1]
        newedge.length = (ismissing(ell) ? -1 : ell)
        if nl > 1
            ell = lengths[2]
            edge.length = (ismissing(ell) ? -1 : ell)
        end
    end
    pushEdge!(net,newedge)
    pushNode!(net,newnode)
    return newnode, newedge
end

"""
    fuseedgesat!(i::Integer,net::HybridNetwork, multgammas::Bool=false)

Removes `i`th node in net.node, if it is of degree 2.
The parent and child edges of this node are fused.
If either of the edges is hybrid, the hybrid edge is retained. Otherwise, the
edge with the lower edge number is retained.

Reverts the action of [`breakedge!`](@ref).

returns the fused edge.
"""
function fuseedgesat!(i::Integer, net::HybridNetwork, multgammas::Bool=false)
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
        (ce,pe) = (pe,ce) # keep the hybrid edge: keep its ismajor
    end
    isnodeiparent = (nodei ≡ getparent(ce))
    (!ce.hybrid || isnodeiparent) ||
      error("node $(nodei.number) has 1 tree edge ($(pe.number)) and 1 hybrid edge ($(ce.number)), but is child of the hybrid edge.")
    pn = getOtherNode(pe,nodei)
    removeEdge!(nodei,ce) # perhaps useless. in case gc() on ith node affects its edges.
    removeNode!(nodei,ce)
    removeEdge!(pn,pe)
    removeNode!(pn,pe)    # perhaps useless. in case gc() on pe affects its nodes.
    setEdge!(pn,ce)
    setNode!(ce,pn)       # pn comes 2nd in ce now: ce.node is: [original, pn]
    ce.ischild1 = isnodeiparent # to retain same direction as before.
    ce.length = addBL(ce.length, pe.length)
    if multgammas
        ce.gamma = multiplygammas(ce.gamma, pe.gamma)
    end
    if net.rooti==i # isnodeiparent should be true, unless the root and ce's direction were not in sync
        newroot = pn
        if newroot.leaf && !ce.hybrid # then reverse ce's direction. pn.leaf and ce.hybrid should never both occur!
            newroot = ce.node[1] # getOtherNode(ce, pn)
            ce.ischild1 = false
        end
        net.rooti = findfirst(isequal(newroot), net.node)
    end
    deleteNode!(net,nodei)
    deleteEdge!(net,pe,part=false) # do not update partitions. irrelevant for networks of level>1.
    return ce
end

"""
    removedegree2nodes!(net::HybridNetwork, keeproot::Bool=false)

Delete *all* nodes of degree two in `net`, fusing the two adjacent edges
together each time, and return the network.
If the network has a degree-2 root and `keeproot` is false,
then the root is eliminated as well, leaving the network unrooted.
The only exception to this rule is if the root is incident to 2 (outgoing)
hybrid edges. Removing the root should leave a loop-edge (equal end point),
which we don't want to do, to preserve the paths in the original network.
In this case, the root is maintained even if `keeproot` is false.
If `keeproot` is true, then the root is kept even if it's of degree 2.

See [`fuseedgesat!`](@ref).

```jldoctest
julia> net = readnewick("(((((S1,(S2)#H1),(#H1,S3)))#H2),(#H2,S4));");

julia> PhyloNetworks.breakedge!(net.edge[3], net); # create a degree-2 node along hybrid edge

julia> PhyloNetworks.breakedge!(net.edge[3], net); # another one: 2 in a row

julia> PhyloNetworks.breakedge!(net.edge[10], net); # another one, elsewhere

julia> writenewick(net) # extra pairs of parentheses
"((#H2,S4),(((((S1,(((S2)#H1))),(#H1,S3)))#H2)));"

julia> removedegree2nodes!(net);

julia> writenewick(net) # even the root is gone
"(#H2,S4,(((S1,(S2)#H1),(#H1,S3)))#H2);"

julia> net = readnewick("((((C:0.9)I1:0.1)I3:0.1,((A:1.0)I2:0.4)I3:0.6):1.4,(((B:0.2)H1:0.6)I2:0.5)I3:2.1);");

julia> removedegree2nodes!(net, true);

julia> writenewick(net, round=true) # the root was kept
"((C:1.1,A:2.0):1.4,B:3.4);"

```
"""
function removedegree2nodes!(net::HybridNetwork, keeproot::Bool=false)
    rootnode = getroot(net)
    # caution: the root and its incident edges may change when degree-2 nodes
    #          are removed. Indices of nodes to be removed would change too.
    rootin2cycle(nn) = isrootof(nn,net) && all(e.hybrid for e in nn.edge)
    toberemoved(nn) = (keeproot ? length(nn.edge) == 2 && nn !== rootnode :
                                  length(nn.edge) == 2 && !rootin2cycle(nn))
    ndegree2nodes = sum(toberemoved.(net.node))
    for _ in 1:ndegree2nodes # empty if 0 degree-2 nodes
        i = findfirst(toberemoved, net.node)
        # i may be 'nothing' if the root was initially thought to be removed
        # but later its edges turned to be hybrids, so should not be removed
        isnothing(i) || fuseedgesat!(i, net)
    end
    return net
end

"""
    suppressroot!(net::HybridNetwork)

Suppress the root of `net`, in an attempt to convey the interpretation that
its root is in fact unknown, as in unrooted trees and semidirected networks.
Namely:
1. delete the root if has only 1 edge, unless its child is a leaf, and do so
   recursively until the new root has 2 or more edges (or a single child leaf)
2. stop if the root is of degree 3 or more or has 2 outgoing hybrid edges,
3. otherwise fuse its 2 edges, unless its 2 children are leaves
   (in which case a warning is printed).

An error is thrown if the root is a leaf.

See also [`rootonedge!`](@ref) and [`rootatnode!`](@ref), which do the opposite.

# examples

```jldoctest
julia> net = readnewick("(((((a,b))#H2,(#H2,c))));"); # extra root edge

julia> suppressroot!(net); writenewick(net) # root with 3 edges
"(#H2,c,((a,b))#H2);"

julia> deleteleaf!(net, "c"); writenewick(net) # root with 2 hybrid edges
"(#H2,((a,b))#H2);"

julia> suppressroot!(net); writenewick(net) # root cannot be suppressed more
"(#H2,((a,b))#H2);"
```

Below, the resulting network still has a root with only 2 children:

```jldoctest ; filter = r"└ @ PhyloNetworks.*" => s""
julia> net = readnewick("(((a,(b)i1)));");

julia> # plot(net, shownodelabel=true); # extra: root above LSA, and degree-2 node i1 above b

julia> suppressroot!(net); writenewick(net) # root now at i1: still of degree 2
"(b,a)i1;"
```
"""
function suppressroot!(net::HybridNetwork)
    getroot(net).leaf && error("the root is a leaf: it won't be suppressed")
    if length(getroot(net).edge) == 1
        deleteleaf!(net, net.rooti; index=true, unroot=false)
        # recursive. unroot=false to stop when the new root has degree > 1:
        # otherwise deleteleaf! would remove a root if it starts a 2-cycle
    end
    rn = getroot(net) # must be degree > 1, unless net has a single leaf
    if length(rn.edge) == 2 && !all(e.hybrid for e in rn.edge)
        if all(getchild(e).leaf for e in rn.edge)
            @warn "the root's 2 children are leaves: it will not be suppressed"
        else
            fuseedgesat!(net.rooti, net)
        end
    end
    net.isrooted = length(getroot(net).edge) < 3 # .isrooted not used, but hey
    return net
end

"""
    addleaf!(net::HybridNetwork, node::Node, leafname::String, edgelength::Float64=-1.0)
    addleaf!(net::HybridNetwork, edge::Edge, leafname::String, edgelength::Float64=-1.0)

Add a new external edge between `node` or between the "middle" of `edge`
and a newly-created leaf, of name `leafname`.
By default, the new edge length is missing (-1).

output: newly created leaf node.

# examples

```jldoctest
julia> net = readnewick("((S1,(((S2,(S3)#H1),(#H1,S4)))#H2),(#H2,S5));");

julia> net.node[6].name # leaf S4
"S4"

julia> PhyloNetworks.addleaf!(net, net.node[6], "4a"); # adding leaf to a node

julia> writenewick(net, internallabel=true)
"((S1,(((S2,(S3)#H1),(#H1,(4a)S4)))#H2),(#H2,S5));"

julia> PhyloNetworks.addleaf!(net, net.node[6], "4b");

julia> writenewick(net, internallabel=true)
"((S1,(((S2,(S3)#H1),(#H1,(4a,4b)S4)))#H2),(#H2,S5));"
```

```jldoctest
julia> net = readnewick("((S1,(((S2,(S3)#H1),(#H1,S4)))#H2),(#H2,S5));");

julia> [n.name for n in net.edge[7].node] # external edge to S4
2-element Vector{String}:
 "S4"
 ""  

julia> PhyloNetworks.addleaf!(net, net.edge[7], "4a"); # adding leaf to an edge

julia> writenewick(net, internallabel=true)
"((S1,(((S2,(S3)#H1),(#H1,(S4,4a))))#H2),(#H2,S5));"
```
"""
function addleaf!(net::HybridNetwork, speciesnode::Node, leafname::String, edgelength::Float64=-1.0)
    exterioredge = Edge(maximum(e.number for e in net.edge) + 1, edgelength) # ischild1 = true by default in edge creation
    pushEdge!(net, exterioredge)
    setEdge!(speciesnode, exterioredge)
    if speciesnode.hybrid || (!isrootof(speciesnode, net) && !getparentedge(speciesnode).containroot)
        exterioredge.containroot = false
    end
    newleaf = Node(maximum(n.number for n in net.node) + 1, true, false, [exterioredge]) # Node(number, leaf, hybrid, edge array)
    newleaf.name = leafname
    setNode!(exterioredge, [newleaf, speciesnode]) # [child, parent] to match ischild1 = true by default
    if speciesnode.leaf
        deleteat!(net.leaf,findfirst(isequal(speciesnode), net.leaf))
        speciesnode.leaf = false
        net.numtaxa -= 1
    end
    pushNode!(net, newleaf) # push node into network (see auxillary.jl)
    return newleaf
end

function addleaf!(net::HybridNetwork, startingedge::Edge, leafname::String, edgelength::Float64=-1.0)
    newnode, _ = breakedge!(startingedge, net)
    return addleaf!(net, newnode, leafname, edgelength)
end


# Claudia SL & Paul Bastide: November 2015, Cecile: Feb 2016

#################################################
# Direct Edges
#################################################

"""
    directedges!(net::HybridNetwork; checkMajor::Bool=true)

Updates the edges' attribute `ischild1`, according to the root placement.
Also updates edges' attribute `containroot`, for other possible root placements
compatible with the direction of existing hybrid edges.
Relies on hybrid nodes having exactly 1 major hybrid parent edge,
but checks for that if `checkMajor` is true.

Warnings:
1. Assumes that `ischild1` is correct on hybrid edges
   (to avoid changing the identity of which nodes are hybrids and which are not).
2. Does not check for cycles (to maintain a network's DAG status)

Returns the network. Throws a 'RootMismatch' Exception if the root was found to
conflict with the direction of any hybrid edge.
"""
function directedges!(net::HybridNetwork; checkMajor::Bool=true)
    if checkMajor # check each node has 2+ hybrid parent edges (if any), and exactly one major.
        for n in net.node
            nparents = 0 # 0 or 2 normally, but could be >2 if polytomy.
            nmajor = 0   # there should be exactly 1 major parent if nparents>0
            for e in n.edge
                if e.hybrid && n == getchild(e)
                    nparents += 1
                    if (e.ismajor) nmajor +=1; end
                end
            end
            (nparents!=1) || error("node $(n.number) has exactly 1 hybrid parent edge")
            (nparents==0 || nmajor == 1) ||
              error("hybrid node $(n.number) has 0 or 2+ major hybrid parents")
            (nparents!=2 || n.hybrid) ||
              @warn "node $(n.number) has 2 parents but its hybrid attribute is false.
It is not used in directedges!, but might cause an error elsewhere."
            # to fix this: change n.hybrid, net.hybrid, net.numhybrids etc.
            # none of those attributes are used here.
        end
    end
    net.boolg2 = false # attributed used by snaq! Will change ischild1 and containroot
    for e in net.node[net.rooti].edge
        traverseDirectEdges!(net.node[net.rooti],e,true)
    end
    net.isrooted = true
    return net
end

# containroot = true until the path goes through a hybrid node, below which
# containroot is turned to false.
function traverseDirectEdges!(node::Node, edge::Edge, containroot::Bool)
    if edge.hybrid && node==getchild(edge)
        throw(RootMismatch(
"direction (ischild1) of hybrid edge $(edge.number) conflicts with the root.
ischild1 and containroot were updated for a subset of edges in the network only."))
    end
    if node == edge.node[1]
        edge.ischild1 = false
        cn = edge.node[2] # cn = child node
    else
        edge.ischild1 = true
        cn = edge.node[1]
    end
    edge.containroot = containroot
    if !cn.leaf && (!edge.hybrid || edge.ismajor) # continue down recursion
        if edge.hybrid containroot=false; end # changes containroot locally, intentional.
        nchildren=0
        for e in cn.edge
            if e==edge continue; end
            if (e.hybrid && cn == getchild(e)) continue; end
            traverseDirectEdges!(cn,e,containroot)
            nchildren += 1
        end
        if nchildren==0
            throw(RootMismatch("non-leaf node $(cn.number) had 0 children.
Could be a hybrid whose parents' direction conflicts with the root.
ischild1 and containroot were updated for a subset of edges in the network only."))
        end
    end
    return nothing
end

#################################################
## Topological sorting
#################################################

"""
    preorder!(net::HybridNetwork)

Update attribute `net.vec_node` in which the nodes are pre-ordered
(also called topological sorting), such that each node is visited after its parent(s).
The edges' direction needs to be correct before calling `preorder!`, using `directedges!`
"""
function preorder!(net::HybridNetwork)
    net.isrooted || error("net needs to be rooted for preorder!, run root functions or directedges!")
    net.vec_node = Node[] # path of nodes in preorder.
    queue = Node[] # problem with PriorityQueue(): dequeue() takes a
                   # random member if all have the same priority 1.
    # using net.vec_bool track which nodes have already been *visited*
    net.vec_bool = [false for i = 1:size(net.node,1)];
    push!(queue,net.node[net.rooti]) # push root into queue
    while !isempty(queue)
        #println("at this moment, queue is $([n.number for n in queue])")
        curr = pop!(queue); # deliberate choice over shift! for cladewise order
        currind = findfirst(x -> x===curr, net.node)
        # the "curr"ent node may have been already visited: because simple loop (2-cycle)
        !net.vec_bool[currind] || continue
        net.vec_bool[currind] = true # visit curr node
        push!(net.vec_node,curr) #push curr into path
        for e in curr.edge
            if curr == getparent(e)
                other = getchild(e)
                if !e.hybrid
                    push!(queue,other)
                    # print("queuing: "); @show other.number
                else
                    e2 = getpartneredge(e, other)
                    parent = getparent(e2)
                    if net.vec_bool[findfirst(x -> x===parent, net.node)]
                      push!(queue,other)
                      # warning: if simple loop, the same node will be pushed twice: child of "curr" via 2 edges
                    end
                end
            end
        end
    end
    # println("path of nodes is $([n.number for n in net.vec_node])")
end


"""
    cladewiseorder!(net::HybridNetwork)

Update the internal attribute `net.vec_int1`. Used for plotting the network.
In the major tree, all nodes in a given clade are consecutive. On a tree, this function
also provides a pre-ordering of the nodes.
The edges' direction needs to be correct before calling
[`cladewiseorder!`](@ref), using [`directedges!`](@ref)
"""
function cladewiseorder!(net::HybridNetwork)
    net.isrooted || error("net needs to be rooted for cladewiseorder!\n run root functions or directedges!")
    net.vec_int1 = Int[]
    queue = [net.rooti] # index (in net) of nodes in the queue
    # print("queued the root's children's indices: "); @show queue
    while !isempty(queue)
        ni = pop!(queue); # deliberate choice over shift! for cladewise order
        # @show net.node[ni].number
        push!(net.vec_int1, ni)
        for e in net.node[ni].edge
            if net.node[ni] ≡ getparent(e) # net.node[ni] is parent node of e
                if e.ismajor
                    push!(queue, findfirst(isequal(getchild(e)), net.node))
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

Use `plot(net, shownodenumber=true, showedgenumber=false)` to map node and edge numbers
on the network, as shown in the examples below.
(see package [PhyloPlots](https://github.com/juliaphylo/PhyloPlots.jl))

Warning: assumes that edges are correctly directed (ischild1 updated). This is done
by `plot(net)`. Otherwise run `directedges!(net)`.

# Example

```julia
julia> net = readnewick("(A:1.0,((B:1.1,#H1:0.2::0.2):1.2,(((C:0.52,(E:0.5)#H2:0.02::0.7):0.6,(#H2:0.01::0.3,F:0.7):0.8):0.9,(D:0.8)#H1:0.3::0.8):1.3):0.7):0.1;");
julia> using PhyloPlots
julia> plot(net, shownodenumber=true)
julia> rotate!(net, -4)
julia> plot(net)
julia> net=readnewick("(4,((1,(2)#H7:::0.864):2.069,(6,5):3.423):0.265,(3,#H7:::0.136):10.0);");
julia> plot(net, shownodenumber=true, showedgenumber=true)
julia> rotate!(net, -1, orderedEdgeNum=[1,12,9])
julia> plot(net, shownodenumber=true, showedgenumber=true)
julia> rotate!(net, -3)
julia> plot(net)
```

Note that `LinearAlgebra` also exports a function named `rotate!` in Julia v1.5.
If both packages need to be used in Julia v1.5 or higher,
usage of `rotate!` needs to be qualified, such as with `PhyloNetworks.rotate!`.
"""
function rotate!(net::HybridNetwork, nnum::Integer; orderedEdgeNum::Array{Int,1}=Int[])
    nind = 0
    nind = findfirst(n -> n.number == nnum, net.node)
    nind !== nothing || error("cannot find any node with number $nnum in network.")
    n = net.node[nind]
    ci = Int[] # children edge indices
    for i = 1:length(n.edge)
        if n == getparent(n.edge[i])
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
#   inte1, intn1 (cycle number), partition, branch lengths, diamond/triangle types etc.
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

`keeporiginalroot`: if true, keep the root even if it is of degree one
(forcing `unroot` to be false).

Warning: does **not** update edges' `containroot` nor internal attributes
(e.g. those used by SNaQ for level-1 networks).
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
function deleteleaf!(
    net::HybridNetwork,
    nodeNumber::Integer;
    index::Bool=false,
    nofuse::Bool=false,
    simplify::Bool=true,
    unroot::Bool=false,
    multgammas::Bool=false,
    keeporiginalroot::Bool=false
)
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
        @warn "Only 1 node. Removing it: the network will be empty"
        deleteNode!(net,nodei) # empties the network
    elseif nodeidegree == 1
        pe = nodei.edge[1]
        pn = getOtherNode(pe, nodei) # parent node of leaf
        if net.rooti == i && keeporiginalroot
            return nothing
        end
        # keep nodei if pn is a leaf: keep 1 edge for the single remaining leaf
        if pn.leaf
            net.rooti = i # it should have been i before anyway
            length(net.edge)==1 || error("neighbor of degree-1 node $(nodei.name) is a leaf, but network had $(length(net.edge)) edges (instead of 1).")
            length(pn.edge)==1 || error("neighbor of $(nodei.name) is a leaf, incident to $(length(pn.edge)) edges (instead of 1)")
            return nothing
        end
        # remove nodei and pe.
        removeNode!(pn,pe)  # perhaps useless. in case gc() on pe affects pn
        removeEdge!(pn,pe)
        deleteEdge!(net,pe,part=false)
        if net.rooti==i # if node was the root, new root = pn
            net.rooti = findfirst(x -> x===pn, net.node)
        end
        deleteNode!(net,nodei) # this updates the index net.rooti
        deleteleaf!(net, pn.number; nofuse = nofuse, simplify=simplify, unroot=unroot, multgammas=multgammas,
                    keeporiginalroot=keeporiginalroot)
        return nothing
    elseif nodeidegree > 2
        # do nothing: nodei has degree 3+ (through recursive calling)
        return nothing
    end
    # if we get to here, nodei has degree 2 exactly: --e1-- nodei --e2--
    if i==net.rooti && (keeporiginalroot || !unroot)
        return nothing # node = root of degree 2 and we want to keep it
    end
    e1 = nodei.edge[1]
    e2 = nodei.edge[2]
    if e1.hybrid && e2.hybrid
        cn  = getchild(e1)
        cn2 = getchild(e2)
        if !(nodei ≡ cn && nodei ≡ cn2) # nodei *not* the child of both e1 and e2
            # possible at the root, in which case e1,e2 should have same child
            (i==net.rooti && cn ≡ cn2) ||
                error("after removing descendants, node $(nodei.number) has 2 hybrid edges but is not the child of both.")
            # delete e1,e2,nodei and move the root to their child cn
            cn.hybrid || error("child node $(cn.number) of hybrid edges $(e1.number) and $(e2.number) should be a hybrid.")
            # check that cn doesn't have any other parent than e1 and e2
            any(getchild(e) ≡ cn && e !== e1 && e !==e2 for e in cn.edge) &&
                error("root has 2 hybrid edges, but their common child has an extra parent")
            removeEdge!(cn,e1); removeEdge!(cn,e2)
            removeHybrid!(net,cn) # removes cn from net.hybrid, updates net.numhybrids
            cn.hybrid = false # !! allowrootbelow! not called: would require correct ischild1
            empty!(e1.node); empty!(e2.node)
            deleteEdge!(net,e1,part=false); deleteEdge!(net,e2,part=false)
            empty!(nodei.edge)
            deleteNode!(net,nodei)
            net.rooti = findfirst(x -> x ≡ cn, net.node)
            deleteleaf!(net, net.rooti; index=true, nofuse=nofuse, simplify=simplify,
                unroot=unroot, multgammas=multgammas, keeporiginalroot=keeporiginalroot)
            return nothing
        end
        # by now, nodei is the child of both e1 and e2
        p1 = getparent(e1) # find both parents of hybrid leaf
        p2 = getparent(e2)
        # remove node1 and both e1, e2
        sameparent = (p1≡p2) # 2-cycle
        removeNode!(p1,e1);  removeNode!(p2,e2) # perhaps useless
        removeEdge!(p1,e1);  removeEdge!(p2,e2)
        deleteEdge!(net,e1,part=false); deleteEdge!(net,e2,part=false)
        if net.rooti==i net.rooti=getIndex(p1,net); end # should never occur though.
        deleteNode!(net,nodei)
        # recursive call on both p1 and p2.
        deleteleaf!(net, p1.number; nofuse = nofuse, simplify=simplify, unroot=unroot,
                    multgammas=multgammas, keeporiginalroot=keeporiginalroot)
        # p2 may have already been deleted: e.g. if sameparent, or other scenarios
        if !sameparent
          p2idx = findfirst(n -> n.number == p2.number, net.node)
          isnothing(p2idx) ||
            deleteleaf!(net, p2idx; index=true, nofuse=nofuse, simplify=simplify,
                        unroot=unroot, multgammas=multgammas, keeporiginalroot=keeporiginalroot)
        end
    elseif !nofuse
        e1 = fuseedgesat!(i,net, multgammas) # fused edge
        if simplify && e1.hybrid # check for 2-cycle at new hybrid edge
            cn = getchild(e1)
            e2 = getpartneredge(e1, cn) # companion hybrid edge
            pn  = getparent(e1)
            if pn ≡ getparent(e2)
                # e1 and e2 have same child and same parent. Remove e1.
                e2.hybrid = false # assumes bicombining at cn: no third hybrid parent
                e2.ismajor = true
                e2.gamma = addBL(e1.gamma, e2.gamma)
                removeEdge!(pn,e1); removeEdge!(cn,e1)
                deleteEdge!(net,e1,part=false)
                removeHybrid!(net,cn) # removes cn from net.hybrid, updates net.numhybrids
                cn.hybrid = false # !! allowrootbelow! not called: would require correct ischild1
                # call recursion again because pn and/or cn might be of degree 2 (or even 1).
                deleteleaf!(net, cn.number; nofuse = nofuse, simplify=simplify, unroot=unroot,
                            multgammas=multgammas, keeporiginalroot=keeporiginalroot)
                pnidx = findfirst(n -> n.number == pn.number, net.node)
                isnothing(pnidx) ||
                deleteleaf!(net, pnidx; index=true, nofuse=nofuse, simplify=simplify,
                    unroot=unroot, multgammas=multgammas, keeporiginalroot=keeporiginalroot)
            end
        end
    end
    return nothing
end



"""
    deletehybridedge!(net::HybridNetwork, edge::Edge,
                      nofuse=false, unroot=false,
                      multgammas=false, simplify=true, keeporiginalroot=false)

Delete a hybrid `edge` from `net` and return the network.
The network does not have to be of level 1 and may contain polytomies,
although each hybrid node must have exactly 2 parents.
Branch lengths are updated, allowing for missing values.

If `nofuse` is false, when `edge` is removed, its child (hybrid) node is removed
and its partner hybrid edge is removed.
Its child edge is retained (below the hybrid node), fused with the former partner,
with new length: old length + length of `edge`'s old partner.
Any 2-cycle is simplified into a single edge, unless `simplify` is false.

If `nofuse` is true, edges with descendant leaves are kept as is,
and are not fused. Nodes are retained during edge removal,
provided that they have at least one descendant leaf.
The hybrid edge that is partner to `edge` becomes a tree edge,
but has its γ value unchanged (it is not set to 1), since it is not merged
with its child edge after removal of the reticulation.  
Also, 2-cycles are not simplified if `nofuse` is true.
That is, if we get 2 hybrid edges both from the same parent to the same child,
these hybrid edges are retained without being fused into a single tree edge.

If `unroot` is false and if the root is up for deletion during the process,
it will be kept if it's of degree 2 or more.
A root node of degree 1 will be deleted unless `keeporiginalroot` is true.

If `multgammas` is true: inheritance weights are kept by multiplying together
the inheritance γ's of edges that are merged. For example,
if there is a hybrid ladder, the partner hybrid edge remains a hybrid edge
(with a new partner), and its γ is the product of the two hybrid edges
that have been fused. So it won't add up to 1 with its new partner's γ.

If `keeporiginalroot` is true, a root of degree one will not be deleted.

Warnings:

- `containroot` is updated, but this requires correct `ischild1` fields
- if the parent of `edge` is the root and if `nofuse` is false, the root
  is moved to keep the network unrooted with a root of degree two.
- does *not* update attributes needed for snaq! (like inte1, edge.z, edge.y etc.)
"""
function deletehybridedge!(
    net::HybridNetwork,
    edge::Edge,
    nofuse::Bool=false,
    unroot::Bool=false,
    multgammas::Bool=false,
    simplify::Bool=true,
    keeporiginalroot::Bool=false
)
    edge.hybrid || error("edge $(edge.number) has to be hybrid for deletehybridedge!")
    n1 = getchild(edge)  # child of edge, to be deleted unless nofuse
    n1.hybrid || error("child node $(n1.number) of hybrid edge $(edge.number) should be a hybrid.")
    n1degree = length(n1.edge)
    n2 = getparent(edge)  # parent of edge, to be deleted too
    n2degree = length(n2.edge)
    # next: keep hybrid node n1 if it has 4+ edges or if keepNode.
    #       otherwise: detach n1, then delete recursively
    delete_n1_recursively = false
    if n1degree < 3
        error("node $(n1.number) has $(length(n1.edge)) edges instead of 3+");
    # alternatively: error if degree < 2 or leaf,
    #   warning if degree=2 and internal node, then
    #   delete_n1_recursively = true # n1 doesn't need to be detached first
    elseif n1degree == 3 && !nofuse # then fuse 2 of the edges and detach n1
        delete_n1_recursively = true
        pe = nothing # will be other parent (hybrid) edge of n1
        ce = nothing # will be child edge of n1, to be merged with pe
        for e in n1.edge
            if e.hybrid && e!==edge && n1===getchild(e) pe = e; end
            if !e.hybrid || n1===getparent(e)  ce = e; end # does *not* assume correct ischild1 for tree edges :)
        end
        pn = getparent(pe); # parent node of n1, other than n2
        atRoot = (net.node[net.rooti] ≡ n1) # n1 should not be root, but if so, pn will be new root
        # if pe may contain the root, then allow the root on ce and below
        if pe.containroot
            allowrootbelow!(ce) # warning: assumes correct `ischild1` for ce and below
        end
        # next: replace ce by pe+ce, detach n1 from pe & ce, remove pe from network.
        ce.length = addBL(ce.length, pe.length)
        if multgammas
            ce.gamma = multiplygammas(ce.gamma, pe.gamma)
        end
        removeNode!(n1,ce) # ce now has 1 single node cn
        setNode!(ce,pn)    # ce now has 2 nodes in this order: cn, pn
        ce.ischild1 = true
        setEdge!(pn,ce)
        removeEdge!(pn,pe)
        # if (pe.number<ce.number) ce.number = pe.number; end # bad to match edges between networks
        removeEdge!(n1,pe); removeEdge!(n1,ce) # now n1 attached to edge only
        deleteEdge!(net,pe,part=false) # decreases net.numedges   by 1
        removeHybrid!(net,n1) # decreases net.numhybrids by 1
        n1.hybrid = false
        edge.hybrid = false; edge.ismajor = true
        # n1 not leaf, and not added to net.leaf, net.numtaxa unchanged
        if atRoot
            i = findfirst(x -> x===pn, net.node)
            i !== nothing || error("node $(pn.number) not in net!")
            net.rooti = i
        end
        # below: we will need to delete n1 recursively (hence edge)
    else # n1 has 4+ edges (polytomy) or 3 edges but we want to keep it anyway:
        # keep n1 but detach it from 'edge', set its remaining parent to major tree edge
        pe = getpartneredge(edge, n1) # partner edge: keep it this time
        if !pe.ismajor pe.ismajor=true; end
        pe.hybrid = false
        # note: pe.gamma *not* set to 1.0 here
        removeEdge!(n1,edge) # does not update n1.hybrid at this time
        removeHybrid!(net,n1) # removes n1 from net.hybrid, updates net.numhybrids
        n1.hybrid = false
        if pe.containroot
            allowrootbelow!(pe) # warning: assumes correct `ischild1` for pe and below
        end
        # below: won't delete n1, delete edge instead
    end

    formernumhyb = net.numhybrids
    # next: delete n1 recursively, or delete edge and delete n2 recursively.
    # keep n2 if it has 4+ edges (or if nofuse). 1 edge should never occur.
    #       If root, would have no parent: treat network as unrooted and change the root.
    if delete_n1_recursively
        deleteleaf!(net, n1.number; index=false, nofuse=nofuse,
                    simplify=simplify, unroot=unroot, multgammas=multgammas,
                    keeporiginalroot=keeporiginalroot)
    # else: delete "edge" then n2 as appropriate
    elseif n2degree == 1
        error("node $(n2.number) (parent of hybrid edge $(edge.number) to be deleted) has 1 edge only!")
    else
        # fixit: if n2degree == 2 && n2 === net.node[net.rooti] and
        #        if we want to keep original root: then delete edge but keep n2
        # detach n2 from edge, remove hybrid 'edge' from network
        removeEdge!(n2,edge)
        deleteEdge!(net,edge,part=false)
        # remove n2 as appropriate later (recursively)
        deleteleaf!(net, n2.number; index=false, nofuse=nofuse,
                    simplify=simplify, unroot=unroot, multgammas=multgammas,
                    keeporiginalroot=keeporiginalroot)
    end
    if net.numhybrids != formernumhyb # deleteleaf! does not update containroot
        allowrootbelow!(net)
    end
    return net
end


"""
    deleteaboveLSA!(net, preorder=true)

Delete edges and nodes above (ancestral to) the least stable ancestor (LSA)
of the leaves in `net`. See [`leaststableancestor`](@ref) for the definition
of the LSA.
Output: modified network `net`.
"""
function deleteaboveLSA!(net::HybridNetwork, preorder::Bool=true)
    lsa, lsaindex = leaststableancestor(net, preorder)
    for _ in 1:(lsaindex-1)
        # the network may temporarily have multiple "roots"
        nodei = popfirst!(net.vec_node)
        for e in nodei.edge
            # delete all of nodei's edges (which much be outgoing)
            cn = getchild(e)
            removeEdge!(cn, e) # also updates cn.booln1
            empty!(e.node)
            deleteEdge!(net, e; part=false)
        end
        empty!(nodei.edge)
        deleteNode!(net, nodei) # resets net.rooti
        if nodei.name != ""
            j = findfirst(isequal(nodei.name), net.names)
            isnothing(j) || deleteat!(net.names, j)
        end
    end
    net.rooti = findfirst( n -> n===lsa, net.node)
    if lsa.hybrid # edge case: LSA may be hybrid if 1 single leaf in network
        removeHybrid!(net, lsa)
        lsa.hybrid = false
    end
    return net
end

"""
    resetnodenumbers!(net::HybridNetwork; checkpreorder=true, type=:ape)

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
- `checkpreorder`: if false, the `ischild1` edge field and the `net.vec_node`
  network field are supposed to be correct (to get nodes in preorder).
  This is not needed when `type=:internalonly`.

# Examples

```jldoctest
julia> net = readnewick("(A,(B,(C,D)));");

julia> PhyloNetworks.resetnodenumbers!(net)

julia> printnodes(net) # first column "node": root is 5
node leaf  hybrid name i_cycle edges'numbers
1    true  false  A    -1      1   
2    true  false  B    -1      2   
3    true  false  C    -1      3   
4    true  false  D    -1      4   
7    false false       -1      3    4    5   
6    false false       -1      2    5    6   
5    false false       -1      1    6   

julia> net = readnewick("(A,(B,(C,D)));");

julia> PhyloNetworks.resetnodenumbers!(net; type=:postorder)

julia> printnodes(net) # first column "node": root is 7
node leaf  hybrid name i_cycle edges'numbers
1    true  false  A    -1      1   
2    true  false  B    -1      2   
3    true  false  C    -1      3   
4    true  false  D    -1      4   
5    false false       -1      3    4    5   
6    false false       -1      2    5    6   
7    false false       -1      1    6   
```
"""
function resetnodenumbers!(
    net::HybridNetwork;
    checkpreorder::Bool=true,
    type::Symbol=:ape
)
    if checkpreorder
      directedges!(net)
      preorder!(net) # to create/update net.vec_node
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
        nodelist = net.vec_node # pre-order: root first
    elseif type == :postorder
        nodelist = reverse(net.vec_node) # post-order
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
    resetedgenumbers!(net::HybridNetwork, verbose=true)

Check that edge numbers of `net` are consecutive numbers from 1 to the total
number of edges. If not, reset the edge numbers to be so.
"""
function resetedgenumbers!(net::HybridNetwork, verbose::Bool=true)
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

Set `containroot` to `false` for edge `e` and all edges below, recursively.
The traversal stops if `e.containroot` is already `false`, assuming
that `containroot` is already false all the way below down that edge.
"""
function norootbelow!(e::Edge)
    e.containroot || return nothing # if already false: stop
    # if true: turn to false then move down to e's children
    e.containroot = false
    cn = getchild(e) # cn = child node
    for ce in cn.edge
        ce !== e || continue # skip e
        getparent(ce) === cn || continue # skip edges that aren't children of cn
        norootbelow!(ce)
    end
    return nothing
end

"""
    allowrootbelow!(e::Edge)
    allowrootbelow!(n::Node, parent_edge_of_n::Edge)

Set `containroot` to `true` for edge `e` and all edges below, recursively.
The traversal stops whenever a hybrid node is encountered:
if the child of `e` is a hybrid node (that is, if `e` is a hybrid edge)
or if `n` is a hybrid node, then the edges below `e` or `n` are *not* traversed.
"""
function allowrootbelow!(e::Edge)
    e.containroot = true
    e.hybrid && return nothing # e hybrid edge <=> its child hybrid node: stop
    allowrootbelow!(getchild(e), e)
end
function allowrootbelow!(n::Node, pe::Edge)
    # pe assumed to be the parent of n
    for ce in n.edge
        ce !== pe || continue # skip parent edge of n
        getparent(ce) === n || continue # skip edges that aren't children of n
        allowrootbelow!(ce)
    end
    return nothing
end
"""
    allowrootbelow!(net::HybridNetwork)

Set `containroot` to `true` for each edge below the root node, then
traverses `net` in preorder to update `containroot` of all edges (stopping
at hybrid nodes): see the other methods.
Assumes correct `ischild1` edge field.
"""
function allowrootbelow!(net::HybridNetwork)
    rn = net.node[net.rooti]
    for e in rn.edge
        if e.containroot
            allowrootbelow!(getchild(e), e)
        end
    end
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
        ce = getchildedge(h)
        hybchild[h.number] = (ce, ce.length) # original lengths before unzipping
        hybladder[h.number] = (ce.hybrid ? getchild(ce) : nothing)
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
    unzipat_canonical!(hyb::Node, childedge::Edge)

Unzip the reticulation a node `hyb`. See [`unzip_canonical!`](@ref PhyloNetworks.unzip_canonical!).
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
