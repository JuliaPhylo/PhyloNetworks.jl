"""
    biconnectedcomponents(network, ignoreTrivial=false)

Calculate biconnected components (aka "blobs") using Tarjan's algorithm.

Output: array of arrays of edges.
- the length of the array is the number of blobs
- each element is an array of all the edges inside a given blob.

These blobs are returned in post-order, but within a blob,
edges are *not* necessarily sorted in topological order.
If `ignoreTrivial` is true, trivial components (of a single edge)
are not returned.
The network is assumed to be connected.

*Warnings*: for nodes, fields `intn2` and `intn1`
are modified during the algorithm. They are used to store the
node's "index" (time of visitation), "lowpoint", and the node's
"parent", as defined by the order in which nodes are visited.
For edges, field `boole2` is modified, to store whether
the edge has been already been visited or not.

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
function biconnectedcomponents(net, ignoreTrivial::Bool=false)
    for n in net.node
        n.intn1 = -1 # intn1 = lowpoint. -1 for missing, until the node is visited
        n.intn2 = -1 # intn2 = index, order of visit during depth-first search
    end
    for e in net.edge
        e.boole2 = false # true after edge is visited, during depth-first search
    end
    S = Edge[] # temporary stack
    blobs = Vector{Edge}[] # will contain the blobs
    biconnectedcomponents(getroot(net), [0], S, blobs, ignoreTrivial)
    # if stack not empty: create last connected component
    length(S) == 0 || @error("stack of edges not empty at the end: $S")
    return blobs
end

"""
    biconnectedcomponents(node, index, S, blobs, ignoreTrivial)

Helper recursive function starting at a node (not a network).
`index` is an array containing a single integer, thus mutable:
order in which nodes are visited.
"""
function biconnectedcomponents(node, index, S, blobs, ignoreTrivial)
    #println("\nentering biconnect, index=$(index[1])")
    children = 0
    # set depth index for v to the smallest unused index
    node.intn2 = index[1] # index / disc
    node.intn1 = index[1] # lowpoint / lowlink
    #println(" for node $(node.number): index=$(node.intn2), low=$(node.intn1)")
    index[1] += 1
    #print(" stack: "); @show [e.number for e in S]

    for e in node.edge # each (v, w) in E do
        w = (e.node[1] == node) ? e.node[2] : e.node[1]
        #println(" w: $(w.number) along edge $(e.number)")
        if w.intn2 == -1 # w not yet visited, therefore e not yet visited either
            #println(" w's parent = $(node.number)")
            e.boole2 = true
            children += 1
            push!(S, e)
            #println(" edge $(e.number) on stack, intermediate step for node=$(node.number)")
            biconnectedcomponents(w, index, S, blobs, ignoreTrivial)
            # check if subtree rooted at w has connection to
            # one of the ancestors of node
            # Case 1 -- per Strongly Connected Components Article
            node.intn1 = min(node.intn1, w.intn1) # lowpoint
            #print(" case 1, low=$(node.intn1) for node $(node.number); ")
            #println("low=$(w.intn1) for w $(w.number)")
            # if node is an articulation point: pop until e
            if w.intn1 >= node.intn2
                # @info "found root or articulation: node number $(node.number), entry to new blob."
                # start a new strongly connected component
                bb = Edge[]
                while length(S)>0
                    e2 = pop!(S)
                    #println(" popping: edge $(e2.number)")
                    push!(bb, e2)
                    e2 !== e || break
                end
                if !ignoreTrivial || length(bb)>1
                    push!(blobs, bb)
                end
            end
        elseif !e.boole2 && node.intn2 > w.intn2
            # e is a back edge, not cross edge. Case 2 in article.
            e.boole2 = true
            node.intn1 = min(node.intn1, w.intn2)
            #println(" case 2, node $(node.number): low=$(node.intn1)")
            push!(S, e)
            #println(" edge $(e.number) now on stack, final step for node $(node.number)")
        else
            #println(" nothing to do from node $(node.number) to w $(w.number) along edge $(e.number)")
        end
    end
    return nothing
end

"""
    biconnectedcomponent_entrynodes(net, bcc, preorder=true)

Array containing the entry node of the each biconnected component in `bcc`.
`bcc` is supposed to contain the biconnected components as output by
[`biconnectedcomponents`](@ref), that is, an array of array of edges.

These entry nodes depend on the rooting (whereas the BCC only depend on the
unrooted graph). They are either the root of the network or cut node
(articulation points).
"""
function biconnectedcomponent_entrynodes(net, bcc, preorder::Bool=true)
    if preorder
        directedges!(net)
        preorder!(net)
    end
    entrynode = Node[] # one entry node for each blob: cut node or root
    for bicomp in bcc
        jmin = length(net.node)
        for edge in bicomp
            n = getparent(edge)
            j = findfirst(x -> x===n, net.vec_node)
            isnothing(j) && error("node not found in net's pre-ordering 'vec_node'")
            jmin = min(j, jmin)
        end
        push!(entrynode, net.vec_node[jmin])
    end
    return entrynode
end

"""
    biconnectedcomponent_exitnodes(net, bcc, preorder=true)

Array containing an array of the exit node(s) of the each biconnected component
in `bcc`. `bcc` is supposed to contain the biconnected components as output by
[`biconnectedcomponents`](@ref), that is, an array of array of edges.

These exit nodes depend on the rooting (whereas the BCC only depend on the
unrooted graph). The degree of a blob is the number of exit nodes + 1 if
the blob doesn't contain the root (its entry node is a cut node), or + 0 if
the blob contains the root (which enters into the blob but isn't a cut node).

*Warning* (or positive side effect?): the edge `.inte1` attribute is modified.
It stores the index (in `bcc`) of the biconnected component that an edge belongs to.
If an edge doesn't belong in any (e.g. if trivial blobs are ignored),
then its `.inte1` is set to -1.
"""
function biconnectedcomponent_exitnodes(net, bcc, preorder::Bool=true)
    if preorder
        directedges!(net)
        preorder!(net)
    end
    exitnode = Vector{Node}[]  # one array of exit cut nodes for each blob
    for edge in net.edge edge.inte1 = -1; end # in case trivial blobs are ignored
    for (i,bicomp) in enumerate(bcc)
        for edge in bicomp edge.inte1 = i; end
    end
    for (i,bicomp) in enumerate(bcc)
        exitnode_blobi = Node[]
        for edge in bicomp
            edge.ismajor || continue # skip minor edges to avoid duplicating exit node
            n = getchild(edge)
            for e in n.edge
                e !== edge || continue
                if e.inte1 != i # then n is a cut point, incident to another blob
                    push!(exitnode_blobi, n)
                    break
                end
            end
        end
        push!(exitnode, exitnode_blobi)
    end
    return exitnode
end

"""
    blobinfo(network, ignoreTrivial=true)

Calculate the biconnected components (blobs) using function
[`biconnectedcomponents`](@ref) then:
- set node field `booln4` to true at the root of each
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

keyword argument: `checkpreorder`, true by default. If false,
the `ischild1` edge field and the `net.vec_node` network field
are supposed to be correct.

**warning**: see [`biconnectedcomponents`](@ref) for node
attributes modified during the algorithm.
"""
function blobinfo(
    net,
    ignoreTrivial::Bool=true;
    checkpreorder::Bool=true
)
    if checkpreorder
      directedges!(net) # update ischild1, needed for preorder
      preorder!(net) # creates / updates net.vec_node
    end
    bcc = biconnectedcomponents(net, ignoreTrivial)
    bccRoots = biconnectedcomponent_entrynodes(net, bcc, false) # 1 entry node for each blob
    bccMajor = Vector{Edge}[] # one array for each blob
    bccMinor = Vector{Edge}[]
    for bicomp in bcc
        bccMa = Edge[]
        bccmi = Edge[] # find minor hybrid edges, in same order
        for edge in bicomp
            if edge.hybrid && edge.ismajor
                push!(bccMa, edge)
                e = getpartneredge(edge)
                !e.ismajor || @warn "major edge $(edge.number) has a major partner: edge $(e.number)"
                push!(bccmi, e)
            end
        end
        push!(bccMajor, bccMa)
        push!(bccMinor, bccmi)
    end
    # add the network root, if it was in a trivial bi-component (no hybrids)
    rootnode = getroot(net)
    if !any(n -> n === rootnode, bccRoots)
        push!(bccRoots, rootnode)
        push!(bccMajor, Edge[])
        push!(bccMinor, Edge[])
    end
    # update `booln4` to mark nodes that are blob roots:
    # these nodes will serve as dummy leaves when reached from other blobs
    for n in net.node n.booln4 = false; end
    for r in bccRoots r.booln4 = true;  end
    return(bccRoots, bccMajor, bccMinor)
end

"""
    blobdecomposition!(network)
    blobdecomposition(network)

Find blobs using [`biconnectedcomponents`](@ref); find their roots
using [`blobinfo`](@ref); create a forest in the form of a
disconnected network (for efficiency), by deconnecting the
root of each non-trivial blob from its parent.
The root of each blob corresponds to a new leaf
(in another tree of the forest):
the number of the blob's root is given to the newly created leaf.

The first (bang) version modifies the network and returns
the array of blob roots. The second version copies the network
then returns a tuple: the forest and the array of blob roots.

Warnings:
- the forest is represented by a single HybridNetwork object,
  on which most functions don't work (like `writenewick`, plotting etc.)
  because the network is disconnected (to make the forest).
  Revert back to low-level functions, e.g. `printedges` and `printnodes`.
- see [`biconnectedcomponents`](@ref) for node
  attributes modified during the algorithm.
"""
function blobdecomposition(net)
    net2 = deepcopy(net)
    blobR = blobdecomposition!(net2)
    return net2, blobR
end
function blobdecomposition!(net)
    nextnumber = maximum([n.number for n in net.node])+1
    blobR, tmp, tmp = blobinfo(net, true) # true: ignore trivial single-edge blobs
    for r in blobR
        for e in r.edge
            r == e.node[e.ischild1 ? 1 : 2] || continue
            removeEdge!(r,e) # detach edge e from root r
            removeNode!(r,e) # detach root r from edge e
            dummyleaf = Node(nextnumber, true) # true: leaf
            nextnumber += 1
            dummyleaf.name = string("dummy ", r.number)
            setEdge!(dummyleaf, e) # attach e to new dummy leaf
            setNode!(e,dummyleaf)  # attach new leaf to edge e
            e.ischild1 = false
            pushNode!(net, dummyleaf)
            break
        end
    end
    return blobR
end

"""
    process_biconnectedcomponents!(net, preorder=true)

Calculate biconnected components (aka "blobs" in binary networks,
"blocks" more generally) and save them in `net` for re-use later.
The field `net.partition` is reset to store the biconnected components,
by storing them in preorder (a block comes before any of its descendents).
Each one is stored as a `Partition` object, which stores a list of edges in the
biconnected component (which partition edges).

Trivial blocks are formed by one cut-edge and its 2 adjacent nodes.
They are stored.

The edge field `.inte1` is used to store the biconnected component that
the edge belongs to, that is, an edge belongs in `net.partition[edge.inte1]`.
This is an internal coding that could break in the future, and that could
be modified by other functions.

!!! warning
    This preprocessing and `net.partition` become obsolete and incorrect if the
    network's topology is later modified, such as after a semidirected NNI,
    the deletion/addition of a hybrid edge / a leaf. Functions that modify
    a network are **not** required to update the blob decomposition (to save
    time for analyses that don't need the blob decomposition).

    After a re-rooting with `rootatnode!`, the blobs remain the same with the
    same set of edges, but their entry & exit nodes may become incorrect.
    After `rootonedge!`, the number of blobs remains correct but one edge was
    subdivided into 2 edges (and the old root node might have been suppressed),
    so the edge list in some blobs becomes obsolete.

# examples

The network below is binary, so each blocks and blobs are the same.
It has 2 non-trivial blobs: one with 1 hybrid and the other with 2 hybrids.
We first get the blobs, then look at the 2-hybrid blob in more detail.

```jldoctest
julia> net = readnewick("(((D,#H2),((((E)#H2,C),#H1),((B)#H1,A))),((H,#H3),((F)#H3,G)));");

julia> # using PhyloPlots; plot(net, showedgenumber=true, shownodenumber=true);

julia> PhyloNetworks.process_biconnectedcomponents!(net);

julia> length(net.partition) # number of biconnected components
12

julia> findall(!PhyloNetworks.istrivial(blob) for blob in net.partition)
2-element Vector{Int64}:
 3
 7

julia> blob = net.partition[7]; getlevel(blob)
2

julia> [e.number for e in blob.edges]
9-element Vector{Int64}:
 14
  9
 13
 11
  8
  7
  5
  2
  3

julia> i = PhyloNetworks.entrynode_preindex(blob); net.vec_node[i] # entry to the blob
PhyloNetworks.Node:
 number:-3
 attached to 3 edges, numbered: 3 14 15

julia> PhyloNetworks.number_exitnodes(blob) # 5 exit nodes connecting to other blobs below
5

julia> [net.vec_node[i].number for i in PhyloNetworks.exitnodes_preindex(blob)]
5-element Vector{Int64}:
 -9
  5
 -7
  2
 -4
```
"""
function process_biconnectedcomponents!(net::HybridNetwork, preorder=true)
    if preorder
        directedges!(net)
        preorder!(net)
    end
    # isempty(net.partition) || @warn("will erase and re-calculate biconnected components")
    empty!.(net.partition)
    empty!(net.partition)
    bcc = biconnectedcomponents(net, false) # false: don't ignore trivial blobs
    entry = biconnectedcomponent_entrynodes(net, bcc, false)
    entryindex = indexin(entry, net.vec_node)
    exitnodes = biconnectedcomponent_exitnodes(net, bcc, false)
    bloborder = sortperm(entryindex) # pre-ordering for blobs in their own blob tree
    for (ii, ib) in enumerate(bloborder)
        bicomp = bcc[ib]
        for e in bicomp e.inte1 = ii; end
        node_indices = Int[entryindex[ib]]
        append!(node_indices, indexin(exitnodes[ib], net.vec_node))
        push!(net.partition, Partition(node_indices, bicomp))
    end
    return net
end

"""
    getlevel(net, preorder=true, preprocess=false)

Level of `net`: maximum number of edges to delete within a biconnected component
to make it a tree. If the network is bicombining (each hybrid has 2 parents),
then the level is the maximum number of hybrid nodes within a biconnected
component.
Arguments:
- `preorder` to direct edges according to the root and store a pre-ordering of
  nodes in `net.vec_node`, if not done earlier on the network (e.g. after
  a topology modification)
- `preprocess` to recalculate & store the biconnected components with
  [`process_biconnectedcomponents!`](@ref), assuming the preordering was done
"""
function getlevel(
    net::HybridNetwork,
    preorder::Bool=true,
    preprocess::Bool=true
)
    if preorder
        directedges!(net)
        preorder!(net)
    end
    if preprocess
        process_biconnectedcomponents!(net, false)
    elseif isempty(net.partition)
        error("no biconnected components stored in the network: re-run with preprocessing")
    end
    return maximum(getlevel.(net.partition))
end
getlevel(p::Partition) = sum(!e.ismajor for e in p.edges)

"""
    leaststableancestor(net, preorder=true, preprocess=true)

Tuple `(lsa, lsa_index)` where `lsa` is the least stable ancestor (LSA) node
in `net`, and `lsa_index` is the index of `lsa` in `net.vec_node`.
The LSA the lowest node `n` with the following property: *any* path
between *any* leaf and the root must go through `n`. All such nodes with this
property are ancestral to the LSA (and therefore must have an index that is
lower or equal to `lsa_index`).

Exception: if the network has a single leaf, the output `lsa` is the
leaf's parent node, to maintain one external edge between the root and the leaf.

Arguments:
- `preorder` to direct edges according to the root and store a pre-ordering of
  nodes in `net.vec_node`, if not done earlier on the network (e.g. after
  a topology modification)
- `preprocess` to recalculate & store the biconnected components with
  [`process_biconnectedcomponents!`](@ref), assuming the preordering was done

*Warning*:
uses [`biconnectedcomponents`](@ref) and [`biconnectedcomponent_exitnodes`](@ref),
therefore share the same caveats regarding the use of field
`.intn1` and `.intn2` for nodes; `boole2` and `.inte1` for edges.
As a positivie side effect, the biconnected components can be recovered
via the edges' `.inte1` field --including the trivial blobs (cut edges).

See also: [`deleteaboveLSA!`](@ref)
"""
function leaststableancestor(
    net::HybridNetwork,
    preorder::Bool=true,
    preprocess::Bool=true
)
    getroot(net).leaf && error("The root can't be a leaf to find the LSA.")
    if preorder
        directedges!(net)
        preorder!(net)
    end
    if preprocess
        process_biconnectedcomponents!(net, false)
    elseif isempty(net.partition)
        error("no biconnected components stored in the network: re-run with preprocessing")
    end
    bcc = net.partition
    entryindices = entrynode_preindex.(bcc)
    function atlsa(blob) # is blob below the LSA? given that previous blobs are not
        # above LSA if 1 exit and 1 entry that's not an entry to another blob
        # (0 exits: trivial blob (cut-edge) to a leaf)
        return (number_exitnodes(blob) != 1 ||
                sum(isequal(entrynode_preindex(blob)), entryindices) > 1)
    end
    lsablob_i = findfirst(atlsa, bcc) # pre-order important here
    isnothing(lsablob_i) && error("strange: couldn't find the LSA...")
    lsa_index = entrynode_preindex(bcc[lsablob_i])
    return net.vec_node[lsa_index], lsa_index
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
                    if curnode === getchild(e)
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

Set the root of `net` to an appropriate node and update the edges `containroot`
field appropriately, using the `membership` output by [`treeedgecomponents`](@ref).
A node is appropriate to serve as root if it belongs in the
root tree-edge component, that is, the root of the tree-edge component graph.

- If the current root is appropriate, it is left as is. The direction of
  edges (via `ischild1`) is also left as is, assuming it was in synch with
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
julia> net = readnewick("(#H1:::0.1,#H2:::0.2,(((b)#H1)#H2,a));");

julia> membership = treeedgecomponents(net);

julia> rootcompID = checkroot!(net, membership);

julia> rootcomp = keys(filter(p -> p.second == rootcompID, membership));

julia> sort([n.number for n in rootcomp]) # number of nodes in the root component
3-element Vector{Int64}:
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
        up, down = e.ischild1 ? (e.node[2], e.node[1]) : (e.node[1], e.node[2])
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
    #    and update the edges' containroot accordingly.
    # NOT done: build the root component to return it. Instead: return tec_root
    #   Set(node for node = nodes if membership[node] == tec_root)
    #   Set(node for (node,comp) in membership if comp == tec_root)
    #rootcomp = keys(filter(p -> p.second == tec_root, membership))
    curroot = nodes[net.rooti]
    if membership[curroot] == tec_root
        # update containroot only: true for edges in or out of the root TEC
        for e in net.edge
            e.containroot = (membership[getparent(e)] == tec_root)
        end
    else
        net.rooti = findfirst(n -> (!n.leaf && membership[n] == tec_root), nodes)
        directedges!(net) # also updates containroot of all edges
    end
    return tec_root # return rootcomp
end

"""
    treeofblobs(net::HybridNetwork; preorder::Bool=true, preprocess::Bool=true)

Tree of blobs for `net`. This tree is obtained from the network by shrinking any
blob (2-edge connected component) that is not trivial (more than 1 edge) into
a single node. The remaining edges are cut-edges (or bridges): they disconnect
the network if removed. After shrinking blobs, only these cut-edges remain (with
their edge length if any) and the graph forms a tree: the tree-of-blobs.
It typically has polytomies, and may have degree-2 nodes.
See [Allman et al. (2022)](https://doi.org/10.1007/s00285-022-01838-9)
(on binary networks) and
[Xu & Ané (2023)](https://doi.org/10.1007/s00285-022-01847-8).

If the network is binary, each blob is a block (biconnected component).
If the network is not binary, then multiple blocks forming a single blob are
shrunk into the same node.

The network is not modified. `PhyloNetworks.treeofblobs!(net)` modifies `net`,
but is unsafe to use, in the sense that it assumes that the network is already
preordered and preprocessed to get its biconnected component.

In the output tree, the biconnected component information is stored.

Arguments:
- `preorder` to direct edges according to the root and store a pre-ordering of
  nodes in `net.vec_node`, if not done earlier on the network (e.g. after
  a topology modification)
- `preprocess` to recalculate & store the biconnected components with
  [`process_biconnectedcomponents!`](@ref), assuming the preordering was done.

# examples

In this first example, the network is binary (no polytomies of any kind)
and has 2 blobs.

```jldoctest
julia> net = readnewick("((S1:0.1,(((S2,(S3)#H1),(#H1,S4:0.4)):0.2)#H2:0.2),(#H2,S5));");

julia> tob = treeofblobs(net);

julia> writenewick(tob) # 2 polytomies: at root, and to S2,S3,S4
"(S5,S1:0.1,(S4:0.4,S2,S3):0.2);"
```

Now we use the same network but without the edge connecting the 2 blobs:
there are now 2 blocks sharing one node, so 1 blob only.

```jldoctest
julia> net = readnewick("((S1:0.1,((S2,(S3)#H1),(#H1,S4:0.4))#H2:0.2),(#H2,S5));");

julia> # using PhyloPlots; plot(net, style=:majortree, arrowlen=0.1);

julia> tob = treeofblobs(net);

julia> writenewick(tob) # star tree: 1 big polytomy
"(S5,S1:0.1,S4:0.4,S2,S3);"
```
"""
function treeofblobs(
    net::HybridNetwork;
    preorder::Bool=true,
    preprocess::Bool=true
)
    if preorder
        directedges!(net)
        preorder!(net)
    end
    if preprocess
        process_biconnectedcomponents!(net, false)
    elseif isempty(net.partition)
        error("no biconnected components stored in the network: re-run with preprocessing")
    end
    return treeofblobs!(deepcopy(net))
end
function treeofblobs!(net::HybridNetwork)
    bcc = net.partition
    # shrink non-trivial blobs & delete them from blob list
    for ib in reverse(1:length(bcc))
        blob = bcc[ib]
        istrivial(blob) && continue
        # step 1: delete minor hybrid edges
        for i in reverse(1:length(blob.edges))
            e = blob.edges[i]
            e.ismajor && continue
            # check that e is still in net: not already deleted as part of a stack
            ii = findfirst(isequal(e), net.edge)
            isnothing(ii) && continue
            deletehybridedge!(net, e, true, false, false, false, true)
            # nofuse, ...., simplify=false, keeporiginalroot=true
            deleteat!(blob.edges, i)
        end
        # step 2: shrink all other edges, not already removed (non-tree-child)
        for e in blob.edges
            isnothing(findfirst(isequal(e), net.edge)) ||
                shrinkedge!(net, e)
        end
        empty!(blob)
        deleteat!(bcc, ib)
    end
    net.numhybrids == 0 ||
        error("$(net.numhybrids) hybrids after deleting all non-trivial blobs' minor edges")
    # update .inte1 to correct trivial blob index
    for (ib, blob) in enumerate(bcc)
        blob.edges[1].inte1 = ib
        empty!(blob.cycle) # nodes preorder indices are now obsolete
    end
    preorder!(net)
    # update each blob's nodes' preorder indices
    for (i,n) in enumerate(net.vec_node)
        for e in n.edge
            # entry node first in .cycle: achieved by preorder traversal
            push!(bcc[e.inte1].cycle, i)
        end
    end
    return net
end
