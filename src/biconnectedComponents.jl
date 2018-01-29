"""
    biconnectedComponents(network, ignoreTrivial=false)

Calculate biconnected components -aka "blobs, using Tarjan's algorithm:
output is a array of arrays of edges. options:

- ignoreTrivial: if true, trivial components are not returned

**warnings**: for nodes, fields `k`, `inCycle` and `isExtBadTriangle`
are modified during the algorithm.
They are used to store "depth", "lowpoint" and "onStack".
"""
function biconnectedComponents(net, ignoreTrivial=false::Bool)
    for n in net.node
        n.inCycle = -1 # -1 for missing, node not visited yet
        n.k = -1
        n.isExtBadTriangle = false
    end
    S = Vector{Node}(0) # temporary stack
    blobs = Vector{Vector{Node}}(0) # will contain the blobs
    biconnectedComponents(net.node[net.root], [0], S, blobs, ignoreTrivial)
    return blobs
end

"""
when applied to a node (not a network), index is an array
containing a single integer --thus mutable. order in which
nodes are found.
"""
function biconnectedComponents(node, index, S, blobs, ignoreTrivial)
    print("\nentering biconnect on: ")
    @show node
    @show index
    # set depth index for v to the smallest unused index
    node.k = index[1]       # index
    node.inCycle = index[1] # lowpoint / lowlink
    println("just updated: index $(node.k), lowlink$(node.inCycle)")
    index[1] += 1
    push!(S, node)
    @show [n.number for n in S]
    node.isExtBadTriangle = true # is on stack

    for e in node.edge # each (v, w) in E do
        w = (e.node[1] == node) ? e.node[2] : e.node[1]
        if w.k == -1 # w not yet visited
            biconnectedComponents(w, index, S, blobs, ignoreTrivial)
            node.inCycle = min(node.inCycle, w.inCycle)
        elseif w.isExtBadTriangle # it is on stack: hence in current biCC
            # if w not on stack: e is cross-edge in DFS tree, ignored
            # the next line may look odd - but is correct.
            # says w.index not w.lowlink; that is deliberate, from original paper
            node.inCycle = min(node.inCycle, w.k)
        end
    end

    # if root node: pop the stack and generate a biCC
    if (node.inCycle == node.k)
      # start a new strongly connected component
      bb = Vector{Node}(0)
      while length(S)>0
        w = pop!(S)
        w.isExtBadTriangle = false
        push!(bb, w)
        w != node || break # break out of while loop if w == node
      end
      push!(blobs, bb)
    end
end
