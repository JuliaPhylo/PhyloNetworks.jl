"""
    biconnectedComponents(network, ignoreTrivial=false)

Calculate biconnected components (aka "blobs") using Tarjan's algorithm:
the output is an array of arrays of edges.
These blobs are returned in post-order, but within a blob,
edges are *not* necessarily sorted in topological order.
If `ignoreTrivial` is true, trivial components (of a single edge)
are not returned.
The network is assumed to be connected.

**Warnings**: for nodes, fields `k`, `inCycle`, and `prev`
are modified during the algorithm. They are used to store the
node's "index" (time of visitation), "lowpoint", and the node's
"parent", as defined by the order in which nodes are visited.

References:
- p. 153 of [Tarjan (1972)](http://epubs.siam.org/doi/pdf/10.1137/0201010).
  Depth first search and linear graph algorithms,
  SIAM Journal on Computing, 1(2):146-160
- on [geeksforgeeks](https://www.geeksforgeeks.org/biconnected-components/),
  there is an error (as of 2018-01-30):
  `elif v != parent[u] and low[u] > disc[v]:` (python version)
  should be replaced by
  `elif v != parent[u] and disc[u] > disc[v]:`
- nice explanation at this
  [url](https://www.cs.cmu.edu/~avrim/451f12/lectures/biconnected.pdf)
"""
function biconnectedComponents(net, ignoreTrivial=false::Bool)
    for n in net.node
        n.inCycle = -1 # -1 for missing, node not visited yet
        n.k = -1       # inCycle = lowpoint
        n.prev = nothing # parent: will be a node later
    end
    S = Edge[] # temporary stack
    blobs = Vector{Edge}[] # will contain the blobs
    biconnectedComponents(net.node[net.root], [0], S, blobs, ignoreTrivial)
    # if stack not empty: create last connected component
    if length(S)>0
        #println("stack not empty at the end.")
        bb = Edge[]
        while length(S)>0
            e2 = pop!(S)
            push!(bb, e2)
        end
        if !ignoreTrivial || length(bb)>1
            push!(blobs, bb)
        end
    end
    return blobs
end

"""
    biconnectedComponents(node, index, S, blobs, ignoreTrivial)

Helper recursive function starting at a node (not a network).
`index` is an array containing a single integer, thus mutable:
order in which nodes are visited.
"""
function biconnectedComponents(node, index, S, blobs, ignoreTrivial)
    #println("\nentering biconnect, index=$(index[1])")
    children = 0
    # set depth index for v to the smallest unused index
    node.k = index[1]       # index / disc
    node.inCycle = index[1] # lowpoint / lowlink
    #println(" for node $(node.number): k=$(node.k), low=$(node.inCycle)")
    index[1] += 1
    #print(" stack: "); @show [e.number for e in S]

    for e in node.edge # each (v, w) in E do
        w = (e.node[1] == node) ? e.node[2] : e.node[1]
        #println(" w: $(w.number) along edge $(e.number)")
        if w.k == -1 # w not yet visited
            w.prev = node # set parent of w to current node
            #println(" w's parent = $(w.prev.number)")
            children += 1
            push!(S, e)
            #println(" edge $(e.number) on stack, intermediate step for node=$(node.number)")
            biconnectedComponents(w, index, S, blobs, ignoreTrivial)
            # check if subtree rooted at w has connection to
            # one of the ancestors of node
            # Case 1 -- per Strongly Connected Components Article
            node.inCycle = min(node.inCycle, w.inCycle) # lowpoint
            #print(" case 1, low=$(node.inCycle) for node $(node.number); ")
            #println("low=$(w.inCycle) for w $(w.number)")
            # if node is an articulation point: pop until e
            if (node.prev === nothing && children > 1) ||
               (node.prev !== nothing && w.inCycle >= node.k)
                # node is either root or an articulation point
                # start a new strongly connected component
                bb = Edge[]
                while length(S)>0
                    e2 = pop!(S)
                    #println(" popping: edge $(e2.number)")
                    push!(bb, e2)
                    e2 != e || break
                end
                if !ignoreTrivial || length(bb)>1
                    push!(blobs, bb)
                end
            end
        elseif w != node.prev && node.k > w.k
            # e is a back edge, not cross edge. Case 2 in article.
            node.inCycle = min(node.inCycle, w.k)
            #println(" case 2, node $(node.number): low=$(node.inCycle)")
            push!(S, e)
            #println(" edge $(e.number) now on stack, final step for node $(node.number)")
        else
            #println(" nothing to do from node $(node.number) to w $(w.number) along edge $(e.number)")
        end
    end
    return nothing
end

"""
    blobInfo(network, ignoreTrivial=true)

Calculate the biconnected components (blobs) using function
[`biconnectedComponents`](@ref) then:
- set node field `isExtBadTriangle` to true at the root of each
  non-trivial blob (and at the network root), false otherwise.
  (a better name for the field would be something like "isBlobRoot".)
- output:
  1. array of nodes that are the roots of each non-trivial blob,
     and the network root. If the root of the full network is
     not part of a non-trivial blob, a corresponding blob is
     added to the list.
  2. array of arrays: for each non-trivial blob,
     array of major hybrid edges in that blob.
  3. array of arrays: same as #2 but for minor hybrid edges,
     with hybrids listed in the same order, for each blob.

Blobs are ordered in reverse topological ordering
(aka post order).
If `ignoreTrivial` is true, trivial components are ignored.

keyword argument: `checkPreorder`, true by default. If false,
the `isChild1` edge field and the `net.nodes_changed` network field
are supposed to be correct.
"""
function blobInfo(net, ignoreTrivial=true::Bool;
    checkPreorder=true::Bool)
    if checkPreorder
      directEdges!(net) # update isChild1, needed for preorder
      preorder!(net) # creates / updates net.nodes_changed
    end
    bcc = biconnectedComponents(net, ignoreTrivial)
    bccRoots = Node[]         # one root node for each blob
    bccMajor = Vector{Edge}[] # one array for each blob
    bccMinor = Vector{Edge}[]
    for bicomp in bcc
        bccMa = Edge[]
        jmin = length(net.node)
        for edge in bicomp
            n = getParent(edge)
            j = something(findfirst(x -> x===n, net.nodes_changed), 0)
            jmin = min(j, jmin)
            if edge.hybrid && edge.isMajor
                push!(bccMa, edge)
            end
        end
        push!(bccRoots, net.nodes_changed[jmin])
        push!(bccMajor, bccMa)
        bccmi = Edge[] # find minor hybrid edges, in same order
        for edge in bccMa
            e = getPartner(edge)
            !e.isMajor || @warn "major edge $(edge.number) has a major partner: edge $(e.number)"
            push!(bccmi, e)
        end
        push!(bccMinor, bccmi)
    end
    # add the network root, if it was in a trivial bi-component (no hybrids)
    # blobs in post-order, so if there the root is there, it's the last blob
    if length(bcc)==0 || bccRoots[end]!=net.node[net.root]
        push!(bccRoots,net.node[net.root])
        push!(bccMajor, Edge[])
        push!(bccMinor, Edge[])
    end
    # update `isExtBadTriangle` to mark nodes that are blob roots:
    # these nodes will serve as dummy leaves when reached from other blobs
    for n in net.node n.isExtBadTriangle = false; end
    for r in bccRoots r.isExtBadTriangle = true;  end
    return(bccRoots, bccMajor, bccMinor)
end

"""
    blobDecomposition!(network)
    blobDecomposition(network)

Find blobs using [`biconnectedComponents`](@ref); find their roots
using [`blobInfo`](@ref); create a forest in the form of a
disconnected network (for efficiency), by deconnecting the
root of each non-trivial blob from its parent.
The root of each blob corresponds to a new leaf
(in another tree of the forest):
the number of the blob's root is given to the newly created leaf.

The first (bang) version modifies the network and returns
the array of blob roots. The second version copies the network
then returns a tuple: the forest and the blob roots.
"""
function blobDecomposition(net)
    net2 = deepcopy(net)
    blobR = blobDecomposition!(net2)
    return net2, blobR
end
function blobDecomposition!(net)
    nextnumber = maximum([n.number for n in net.node])+1
    blobR, tmp, tmp = blobInfo(net, true) # true: ignore trivial single-edge blobs
    for r in blobR
        for e in r.edge
            r == e.node[e.isChild1 ? 1 : 2] || continue
            removeEdge!(r,e) # detach edge e from root r
            removeNode!(r,e) # detach root r from edge e
            dummyleaf = Node(nextnumber, true) # true: leaf
            nextnumber += 1
            dummyleaf.name = string("dummy ", r.number)
            setEdge!(dummyleaf, e) # attach e to new dummy leaf
            setNode!(e,dummyleaf)  # attach new leaf to edge e
            e.isChild1 = false
            pushNode!(net, dummyleaf)
            break
        end
    end
    return blobR
end

"""
    treeedgecomponents(net::HybridNetwork)

Return the tree-edge components of the semidirected network as a `membership`
dictionary `Node => Int`. Nodes with the same membership integer
value are in the same tree-edge component.
The tree-edge components of a network are the connected components of
the network when all hybrid edges are removed.

A `RootMismatch` error is thrown if there exists a cycle in any of the tree-edge
components, or if a tree-edge component has more than one "entry" hybrid node.

Warnings:
- since `Node`s are mutable, the network should not be modified
  until usage of the output `membership` dictionary is over.
- the component IDs are not predicable, but will be consecutive integers
  from 1 to the number of components.
"""
function treeedgecomponents(net::HybridNetwork)
    # partition nodes into tree-edge components (TECs)
    nodes = net.node
    n = length(nodes)
    unvisited = Set(nodes)
    dfs_stack = Vector{Node}()  # stack for iterative depth-first search over tree edges
    dfs_parent = Dict{Node, Node}() # dfs_parent[node] = node's parent in DFS tree
    # membership[node] = id of undirected component that the node belongs to
    membership = Dict{Node, Int}()
    cur_id = 0                  # undirected component id

    while !isempty(unvisited)
        # loop over undirected components
        # start with an unvisited node, DFS with undirected edges
        cur_id += 1
        node = pop!(unvisited) # unpredictable ordering: because unvisited is a set
        push!(dfs_stack, node)
        curnode = node
        dfs_parent[curnode] = curnode
        entrynode = nothing
        while !isempty(dfs_stack)
            # DFS loop over one tree-edge component
            curnode = pop!(dfs_stack)
            delete!(unvisited, curnode)
            membership[curnode] = cur_id
            for e in curnode.edge
                if !e.hybrid
                    # for tree edge, do DFS, check for undirected cycles
                    nextnode = getOtherNode(e, curnode)
                    # if run into visited node (other than parent), then component has cycle
                    if nextnode !== dfs_parent[curnode]
                        if !in(nextnode, unvisited)
                            throw(RootMismatch(
                                "Undirected cycle exists, starting at node number $(nextnode.number)"))
                        else
                            dfs_parent[nextnode] = curnode
                            push!(dfs_stack, nextnode)
                        end
                    end
                else # for hybrid edge, check there is at most one entry node into the TEC
                    if curnode === getChild(e)
                        if isnothing(entrynode)
                            entrynode = curnode
                        elseif entrynode !== curnode
                            throw(RootMismatch(
                                """Multiple entry nodes for tree-edge component, numbered:
                                $(entrynode.number) and $(curnode.number)"""))
                        end
                    end
                end
            end
        end
    end

    return membership
end


"""
    checkroot!(net)
    checkroot!(net::HybridNetwork, membership::Dict{Node, Int})

Set the root of `net` to an appropriate node and update the edges `containRoot`
field appropriately, using the `membership` output by [`treeedgecomponents`](@ref).
A node is appropriate to serve as root if it belongs in the
root tree-edge component, that is, the root of the tree-edge component graph.

- If the current root is appropriate, it is left as is. The direction of
  edges (via `isChild1`) is also left as is, assuming it was in synch with
  the existing root.
- Otherwise, the root is set to the first appropriate node in `net.node`,
  that is not a leaf. Then edges are directed away from this root.

A `RootMismatch` error is thrown if `net` is not a valid semidirected
phylogenetic network (i.e. it is not possible to root the network in a way
compatible with the given hybrid edges).

Output: the `membership` ID of the root component.
The full set of nodes in the root component can be obtained as shown below.
Warning: only use the output component ID after calling the second version
`checkroot!(net, membership)`.

```jldoctest
julia> net = readTopology("(#H1:::0.1,#H2:::0.2,(((b)#H1)#H2,a));");

julia> membership = treeedgecomponents(net);

julia> rootcompID = checkroot!(net, membership);

julia> rootcomp = keys(filter(p -> p.second == rootcompID, membership));

julia> sort([n.number for n in rootcomp]) # number of nodes in the root component
3-element Array{Int64,1}:
 -3
 -2
  4
```
"""
function checkroot!(net::HybridNetwork)
    membership = treeedgecomponents(net) # dict Node => Int
    return checkroot!(net, membership)
end
function checkroot!(net::HybridNetwork, membership::Dict{Node,Int})
    # do not modify nodes or edges until we are done with membership
    nodes = net.node
    ncomp = maximum(values(membership)) # TECs are numbered 1,2,...,ncomp
    #= 1. construct TEC graph: directed graph in which
    vertices = TECs of original network, and
    edge Ci --> Cj for each original hybrid edge: node in Ci --> node in Cj.
    =#
    tecG = [Set{Int}() for comp in 1:ncomp] # tecG[i] = set of children of ith TEC
    noparent = trues(ncomp) # will stay true for TECs with no parent
    for e in net.edge
        e.hybrid || continue # skip tree edges
        up, down = e.isChild1 ? (e.node[2], e.node[1]) : (e.node[1], e.node[2])
        uc_up, uc_down = membership[up], membership[down]
        noparent[uc_down] = false
        push!(tecG[uc_up], uc_down)
    end
    # 2. check that the TEC graph has a single root
    tec_root = findfirst(noparent)
    if isnothing(tec_root) # all TECs have a parent: noparent are all false
        throw(RootMismatch(
            "Semidirected cycle exists: no component has in-degree 0"))
    elseif sum(noparent) > 1
        r1 = findfirst(noparent)
        r2 = findlast(noparent)
        nr1 = findfirst(n -> membership[n] == r1, nodes)
        nr2 = findfirst(n -> membership[n] == r2, nodes)
        throw(RootMismatch("nodes number $(nodes[nr1].number) and $(nodes[nr2].number) have no common ancestor"))
    end

    # 3. topological sort: check the TEC graph has no cycles, that is, is a DAG
    indeg = zeros(Int, ncomp)  # in-degree of vertex in TEC graph
    for uc in tecG
        for i in uc
            indeg[i] += 1
        end
    end
    headstack = [tec_root] # stack of TEC vertices of degree 0 in topological sort
    while !isempty(headstack)
        uc = pop!(headstack)
        indeg[uc] -= 1
        for i in tecG[uc]
            indeg[i] -= 1
            if indeg[i] == 0
                push!(headstack, i)
            end
        end
    end
    cyclehead = findfirst(indeg .!= -1)
    isnothing(cyclehead) || throw(RootMismatch(
        """Semidirected cycle exists, starting at TEC containing node number $(nodes[findfirst(n -> membership[n] == cyclehead, nodes)].number)"""))

    # 4. original network: reset the network root if needed,
    #    and update the edges' containRoot accordingly.
    # NOT done: build the root component to return it. Instead: return tec_root
    #   Set(node for node = nodes if membership[node] == tec_root)
    #   Set(node for (node,comp) in membership if comp == tec_root)
    #rootcomp = keys(filter(p -> p.second == tec_root, membership))
    curroot = nodes[net.root]
    if membership[curroot] == tec_root
        # update containRoot only: true for edges in or out of the root TEC
        for e in net.edge
            e.containRoot = (membership[getParent(e)] == tec_root)
        end
    else
        net.root = findfirst(n -> (!n.leaf && membership[n] == tec_root), nodes)
        directEdges!(net) # also updates containRoot of all edges
    end
    return tec_root # return rootcomp
end

