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

# function to get the parent node of a given node
# it assumes the isChild1 attributes are correct
function getParent(node::Node)
    parents = Node[]
    if(node.hybrid)
        for(e in node.edge)
            if(e.hybrid)
                push!(parents,getOtherNode(e,node))
            end
        end
    else
        for(e in node.edge)
            if(isEqual(node,e.isChild1 ? e.node[1] : e.node[2])) #node is child of e
                push!(parents,getOtherNode(e,node))
            end
        end
    end
    return parents
end

# function to order nodes in topological sorting
# saves vector of nodes in the right order in net.nodes_changed
function preorder!(net::HybridNetwork)
    net.isRooted || error("net needs to be rooted for preorder, run root! or directEdges")
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
                if(isEqual(curr,e.node[e.isChild1 ? 2 : 1])) # curr is the parent node if e
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

# function to compute the variance-covariance between Node and its parents
function updateSharedPathMatrix!(i::Int,nodes::Vector{Node},V::Matrix)
    parent = getParent(nodes[i]) #array of nodes (empty, size 1 or 2)
    if(isempty(parent)) #nodes[i] is root
        return
    elseif(length(parent) == 1) #nodes[i] is tree
        parentIndex = getIndex(parent[1],nodes)
        for(j in 1:(i-1))
            V[i,j] = V[j,parentIndex]
            V[j,i] = V[j,parentIndex]
        end
        V[i,i] = V[parentIndex,parentIndex] + getConnectingEdge(nodes[i],parent[1]).length
    elseif(length(parent) == 2) #nodes[i] is hybrid
        parentIndex1 = getIndex(parent[1],nodes)
        parentIndex2 = getIndex(parent[2],nodes)
        edge1 = getConnectingEdge(nodes[i],parent[1])
        edge2 = getConnectingEdge(nodes[i],parent[2])
        edge1.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[1].number) should be a hybrid egde")
        edge2.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[2].number) should be a hybrid egde")
        for(j in 1:(i-1))
            V[i,j] = V[j,parentIndex1]*edge1.gamma + V[j,parentIndex2]*edge2.gamma
            V[j,i] = V[i,j]
        end
        V[i,i] = edge1.gamma*edge1.gamma*(V[parentIndex1,parentIndex1] + edge1.length) + edge2.gamma*edge2.gamma*(V[parentIndex2,parentIndex2] + edge2.length) + 2*edge1.gamma*edge2.gamma*V[parentIndex1,parentIndex2]
    end
end


function sharedPathMatrix(net::HybridNetwork; checkPreorder=true::Bool) #maybe we only need to input
    net.isRooted || error("net needs to be rooted to get matrix of shared path lengths")
    if(checkPreorder)
        preorder!(net)
    end
    sharedPathMatrix(net.nodes_changed)
end

function sharedPathMatrix(nodes::Vector{Node})
    n = length(net.nodes_changed)
    V = zeros(Float64,n,n)
    for(i in 1:n) #sorted list of nodes
        updateSharedPathMatrix!(i,net.nodes_changed,V)
    end
    return V
end






