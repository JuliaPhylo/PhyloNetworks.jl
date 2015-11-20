# functions for trait evolution on network
# Claudia November 2015

function traverseDirectEdges!(node::Node, edge::Edge)
    if(isEqual(node,edge.node[1]))
        DEBUG && println("changing edge $(edge.number)")
        edge.isChild1 = false
    else
        DEBUG && println("changing edge $(edge.number)")
        edge.isChild1 = true
    end
    node2 = getOtherNode(edge,node)
    if(!node2.leaf)
        for(e in node2.edge)
            if(!isEqual(e,edge) && !e.hybrid)
                traverseDirectEdges!(node2,e)
            end
        end
    end
end

function directEdges!(net::HybridNetwork)
    for(e in net.node[net.root].edge)
        traverseDirectEdges!(net.node[net.root],e)
    end
    net.isRooted = true
end

# fixit: add this to root! function

# function to order nodes in topological sorting
# saves vector of nodes in the right order in net.nodes_changed
function topSorting!(net::HybridNetwork)
    net.isRooted || error("net needs to be rooted for topological sorting, run root! or directEdges")
    println("starting traversal of network from root: $(net.node[net.root].number)")
    net.nodes_changed = Node[] #path of nodes
    queue = PriorityQueue();
    net.visited = [false for i = 1:size(net.node,1)];
    push!(net.nodes_changed,net.node[net.root]) #push root into path
    net.visited[getIndex(net.node[net.root],net)] = true # visit root
    for(e in net.node[net.root].edge)
        enqueue!(queue,getOtherNode(e,net.node[net.root]),1) #enqueue child of root
    end
    while(!isempty(queue))
        #println("at this moment, queue is $([n.number for n in queue])")
        curr = dequeue!(queue);
        #println("$(typeof(curr))")
        println("curr $(curr.number)")
        net.visited[getIndex(curr,net)] = true # visit curr node
        push!(net.nodes_changed,curr) #push curr into path
        for(e in curr.edge)
            if(e.isMajor) #only traverse major edges
                if(isEqual(curr,e.node[e.isChild1 ? 2 : 1])) # e is child edge of curr
                    other = getOtherNode(e,curr)
                    if(!e.hybrid)
                        enqueue!(queue,other,1)
                    else
                        for(e2 in other.edge)
                            if(e2.hybrid && !isEqual(e,e2)) #find the other hybrid parent edge for other
                                parent = getOtherNode(e2,other)
                                if(net.visited[getIndex(parent,net)])
                                    enqueue!(queue,other,1)
                                end
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    println("path of nodes is $([n.number for n in net.nodes_changed])")
end



