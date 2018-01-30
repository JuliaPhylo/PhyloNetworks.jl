"""
    biconnectedComponents(network, ignoreTrivial=false)

Calculate biconnected components (aka "blobs) using Tarjan's algorithm:
the output is an array of arrays of edges.
If `ignoreTrivial` is true, trivial components (of a single edge)
are not returned.

**Warnings**: for nodes, fields `k`, `inCycle`, and `prev``
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
- nice explanation
  [here](https://www.cs.cmu.edu/~avrim/451f12/lectures/biconnected.pdf)
"""
function biconnectedComponents(net, ignoreTrivial=false::Bool)
    for n in net.node
        n.inCycle = -1 # -1 for missing, node not visited yet
        n.k = -1       # inCycle = lowpoint
        n.prev = nothing # parent: will be a node later
    end
    S = Vector{Edge}(0) # temporary stack
    blobs = Vector{Vector{Edge}}(0) # will contain the blobs
    biconnectedComponents(net.node[net.root], [0], S, blobs, ignoreTrivial)
    # if stack not empty: create last connected component
    if length(S)>0
        #println("stack not empty at the end.")
        bb = Vector{Edge}(0)
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
when applied to a node (not a network), index is an array
containing a single integer --thus mutable: order in which
nodes are visited.
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
            if (node.prev == nothing && children > 1) ||
               (node.prev != nothing && w.inCycle >= node.k)
                # node is either root or an articulation point
                # start a new strongly connected component
                bb = Vector{Edge}(0)
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
end
