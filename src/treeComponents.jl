"""
    findrootcomponent!(net::HybridNetwork)

Check that the semidirected graph is a semidirected network (i.e. it
is possible to root the network such that the rooted network, once
unrooted, gives the original semidirected graph), throw a
`RootMismatch` exception if this is not the case. 

If the graph is indeed a semidirected network, return the set of nodes
in the root tree component.  The function also modifies the edge
attribute `containRoot` appropriately along the way.
"""
function findrootcomponent!(net::HybridNetwork)
    # partition nodes into undirected components (UCs)
    nodes = net.node
    n = length(nodes)
    unvisited = Set(nodes)
    dfs_stack = Vector{Node}()
    dfs_parent = Dict{Node, Node}()
    membership = Dict{Node, Int}()
    cur_id = 0

    while !isempty(unvisited)
        # loop over undirected components
        # start with an unvisited node, DFS with undirected edges
        cur_id += 1
        node = pop!(unvisited)
        push!(dfs_stack, node)
        curnode = node
        dfs_parent[curnode] = curnode
        entrynode = nothing
        while !isempty(dfs_stack)
            # DFS loop over one UC
            curnode = pop!(dfs_stack)
            delete!(unvisited, curnode)
            membership[curnode] = cur_id
            for e in curnode.edge
                if !e.hybrid
                    # for undirected edge, do DFS, check for undirected cycles
                    nextnode = curnode == e.node[1] ? e.node[2] : e.node[1]
                    # if run into visited node (that is not the
                    # parent), then the UC has cycle
                    if nextnode != dfs_parent[curnode]
                        if !in(nextnode, unvisited)
                            throw(RootMismatch("Undirected cycle exists"))
                        else
                            dfs_parent[nextnode] = curnode
                            push!(dfs_stack, nextnode)
                        end
                    end
                else
                    # for directed edge, check there is at most one
                    # entry node for the UC
                    child = e.node[e.isChild1 ? 1 : 2]
                    if curnode == child 
                        if isnothing(entrynode)
                            entrynode = curnode
                        elseif entrynode != curnode
                            throw(RootMismatch("Multiple entry nodes for one component"))
                        end
                    end
                end
            end
        end
    end
    
    # construct UC graph and check it's rooted DAG
    ucg = [Set{Int}() for comp in 1:cur_id] # ucg[i] = set of children of ith UC
    noparent = trues(cur_id)
    for e in net.edge
        if e.hybrid
            edge = e.isChild1 ? (e.node[2], e.node[1]) : (e.node[1], e.node[2])
            ucedge = (membership[edge[1]], membership[edge[2]])
            noparent[ucedge[2]] = false
            push!(ucg[ucedge[1]], ucedge[2])
        end
    end
    if sum(noparent) == 0
        throw(RootMismatch("Semidirected cycle exists"))
    elseif sum(noparent) > 1
        throw(RootMismatch("No possible common root"))
    end
    root = findfirst(noparent)
    # DFS from root
    stack = [root]
    compvisited = falses(cur_id)
    while !isempty(stack)
        uc = pop!(stack)
        !compvisited[uc] || throw(RootMismatch("Semidirected cycle exists"))
        compvisited[uc] = true
        for comp in ucg[uc]
            if !compvisited[comp]
                push!(stack, comp)
            end
        end
    end

    rootcomp = Set(node for node = nodes if membership[node] == root)
    # mark all edges in first component as contain root
    for e in net.edge
        up, down = e.isChild1 ? (e.node[2], e.node[1]) : (e.node[1], e.node[2])
        e.containRoot = in(up, root) && (e.hybrid || in(down, root))
    end
    return rootcomp
end

