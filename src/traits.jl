# functions for trait evolution on network
# Claudia November 2015

#################################################
# Direct Edges
#################################################

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

#################################################
## Topological sorting
#################################################

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

#################################################
## Function to traverse the network in the pre-order, updating a matrix
#################################################

# This function takes an init and update funtions as arguments
function recursionPreOrder(net::HybridNetwork, checkPreorder=true::Bool, init=identity::Function, updateRoot=identity::Function, updateTree=identity::Function, updateHybrid=identity::Function, params...)
	net.isRooted || error("net needs to be rooted to get matrix of shared path lengths")
	if(checkPreorder)
		preorder!(net)
	end
	recursionPreOrder(net.nodes_changed, init, updateRoot, updateTree, updateHybrid, params)
end

function recursionPreOrder(nodes::Vector{Node}, init::Function, updateRoot::Function, updateTree::Function, updateHybrid::Function, params)
    n = length(nodes)
    M = init(nodes, params)
    for(i in 1:n) #sorted list of nodes
        updatePreOrder!(i, nodes, M, updateRoot, updateTree, updateHybrid, params)
    end
    return M
end

# Update on the network
# Takes three function as arguments : updateRoot, updateTree, updateHybrid
function updatePreOrder!(i::Int,nodes::Vector{Node},V::Matrix, updateRoot::Function, updateTree::Function, updateHybrid::Function, params)
    parent = getParent(nodes[i]) #array of nodes (empty, size 1 or 2)
    if(isempty(parent)) #nodes[i] is root
        updateRoot(V, i, params)
    elseif(length(parent) == 1) #nodes[i] is tree
        parentIndex = getIndex(parent[1],nodes)
	edge = getConnectingEdge(nodes[i],parent[1])
	updateTree(V, i, parentIndex, edge, params)
    elseif(length(parent) == 2) #nodes[i] is hybrid
        parentIndex1 = getIndex(parent[1],nodes)
        parentIndex2 = getIndex(parent[2],nodes)
        edge1 = getConnectingEdge(nodes[i],parent[1])
        edge2 = getConnectingEdge(nodes[i],parent[2])
        edge1.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[1].number) should be a hybrid egde")
        edge2.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[2].number) should be a hybrid egde")
	updateHybrid(V, i, parentIndex1, parentIndex2, edge1, edge2, params)
    end
end

# Function to get the indexes of the tips. Returns a mask.
function getTipsIndexes(net::HybridNetwork)
	tipsNumbers = [n.number for n in net.leaf]
	nodesOrder = [n.number for n in net.nodes_changed]
	mask = BitArray(length(nodesOrder)) ## Function Match ??
	for tip in tipsNumbers
		mask = mask | (tip .== nodesOrder)
	end
	return(mask)
end


#################################################
## Functions to compute the variance-covariance between Node and its parents
#################################################

function sharedPathMatrix(net::HybridNetwork; checkPreorder=true::Bool)
	recursionPreOrder(net, checkPreorder, initsharedPathMatrix, updateRootSharedPathMatrix!, updateTreeSharedPathMatrix!, updateHybridSharedPathMatrix!)
end

function updateRootSharedPathMatrix!(V::Matrix, i::Int, params)
	return	
end


function updateTreeSharedPathMatrix!(V::Matrix, i::Int, parentIndex::Int, edge::Edge, params)
	for(j in 1:(i-1))
            V[i,j] = V[j,parentIndex]
            V[j,i] = V[j,parentIndex]
        end
        V[i,i] = V[parentIndex,parentIndex] + edge.length
end

function updateHybridSharedPathMatrix!(V::Matrix, i::Int, parentIndex1::Int, parentIndex2::Int, edge1::Edge, edge2::Edge, params)
        for(j in 1:(i-1))
            V[i,j] = V[j,parentIndex1]*edge1.gamma + V[j,parentIndex2]*edge2.gamma
            V[j,i] = V[i,j]
        end
        V[i,i] = edge1.gamma*edge1.gamma*(V[parentIndex1,parentIndex1] + edge1.length) + edge2.gamma*edge2.gamma*(V[parentIndex2,parentIndex2] + edge2.length) + 2*edge1.gamma*edge2.gamma*V[parentIndex1,parentIndex2]
end


#function updateSharedPathMatrix!(i::Int,nodes::Vector{Node},V::Matrix, params)
#    parent = getParent(nodes[i]) #array of nodes (empty, size 1 or 2)
#    if(isempty(parent)) #nodes[i] is root
#        return
#    elseif(length(parent) == 1) #nodes[i] is tree
#        parentIndex = getIndex(parent[1],nodes)
#        for(j in 1:(i-1))
#            V[i,j] = V[j,parentIndex]
#            V[j,i] = V[j,parentIndex]
#        end
#        V[i,i] = V[parentIndex,parentIndex] + getConnectingEdge(nodes[i],parent[1]).length
#    elseif(length(parent) == 2) #nodes[i] is hybrid
#        parentIndex1 = getIndex(parent[1],nodes)
#        parentIndex2 = getIndex(parent[2],nodes)
#        edge1 = getConnectingEdge(nodes[i],parent[1])
#        edge2 = getConnectingEdge(nodes[i],parent[2])
#        edge1.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[1].number) should be a hybrid egde")
#        edge2.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[2].number) should be a hybrid egde")
#        for(j in 1:(i-1))
#            V[i,j] = V[j,parentIndex1]*edge1.gamma + V[j,parentIndex2]*edge2.gamma
#            V[j,i] = V[i,j]
#        end
#        V[i,i] = edge1.gamma*edge1.gamma*(V[parentIndex1,parentIndex1] + edge1.length) + edge2.gamma*edge2.gamma*(V[parentIndex2,parentIndex2] + edge2.length) + 2*edge1.gamma*edge2.gamma*V[parentIndex1,parentIndex2]
#    end
#end

function initsharedPathMatrix(nodes::Vector{Node}, params)
	n = length(nodes)
	return(zeros(Float64,n,n))
end

# Extract the variance at the tips
function extractVarianceTips(V::Matrix, net::HybridNetwork)
	mask = getTipsIndexes(net)
	return(V[mask, mask])
end

#function sharedPathMatrix(net::HybridNetwork; checkPreorder=true::Bool) #maybe we only need to input
#    net.isRooted || error("net needs to be rooted to get matrix of shared path lengths")
#    if(checkPreorder)
#        preorder!(net)
#    end
#    sharedPathMatrix(net.nodes_changed)
#end

#function sharedPathMatrix(nodes::Vector{Node})
#    n = length(net.nodes_changed)
#    V = zeros(Float64,n,n)
#    for(i in 1:n) #sorted list of nodes
#        updateSharedPathMatrix!(i,net.nodes_changed,V)
#    end
#    return V
#end


#################################################
## Functions for Phylgenetic Network regression
#################################################

# New type for phyloNetwork regression
type phyloNetworkRegression
    coefficients::Vector
    sigma2::Real
    loglik::Real
    V::Matrix
    Vy::Matrix
    fittedValues::Vector
    residuals::Vector
end

# Function for lm with net residuals
function phyloNetorklm(Y::Vector, X::Matrix, net::HybridNetwork, model="BM"::AbstractString)
	# Geting variance covariance
	V = sharedPathMatrix(net)
	Vy = extractVarianceTips(V, net)
	# Needed quantities (naive)
	ntaxa = length(Y)
	Vyinv = inv(Vy)
	XtVyinv = X' * Vyinv
	logdetVy = logdet(Vy)
	# beta hat
	betahat = inv(XtVyinv * X) * XtVyinv * Y
	# sigma2 hat
	fittedValues =  X * betahat
	residuals = Y - fittedValues
	sigma2hat = 1/ntaxa * (residuals' * Vyinv * residuals)
	# log likelihood
	loglik = - 1 / 2 * (ntaxa + ntaxa * log(2 * pi) + ntaxa * log(sigma2hat) + logdetVy)
	# Result
	res = phyloNetworkRegression(betahat, sigma2hat[1], loglik[1], V, Vy, fittedValues, residuals)
	return(res)
end

# Add methods on type phyloNetworkRegression



#################################################
## Types for params process
#################################################

# Abstract type of all the (future) types (BM, OU, ...)
abstract paramsProcess

# BM type
type paramsBM <: paramsProcess
    mu::Real # Ancestral value or mean
    sigma2::Real # variance
    randomRoot::Bool # Root is random ? default false
    varRoot::Real # root variance. Default NaN
end
# Constructor
paramsBM(mu, sigma2) = paramsBM(mu, sigma2, false, NaN) # default values



#################################################
## Simulation Function
#################################################

# Uses recursion on the network.
# Takes params of type paramsProcess as an entry
# Returns a matrix with two lines:
# - line one = expectations at all the nodes
# - line two = simulated values at all the nodes
# The nodes are ordered as given by topological sorting
function simulate(net::HybridNetwork, params::paramsProcess, model="BM"::AbstractString, checkPreorder=true::Bool)
	recursionPreOrder(net, checkPreorder, initSimulateBM, updateRootSimulateBM!, updateTreeSimulateBM!, updateHybridSimulateBM!, params)
end

function initSimulateBM(nodes::Vector{Node}, params::Tuple{paramsBM})
	return(zeros(2, length(nodes)))
end

function updateRootSimulateBM!(M::Matrix, i::Int, params::Tuple{paramsBM})
	params = params[1]
	if (params.randomRoot)
		M[1, i] = params.mu # expectation
		M[2, i] = params.mu + sqrt(params.varRoot) * randn() # random value
	else
		M[1, i] = params.mu # expectation
		M[2, i] = params.mu # random value (root fixed)
	end
end


function updateTreeSimulateBM!(M::Matrix, i::Int, parentIndex::Int, edge::Edge, params::Tuple{paramsBM})
	params = params[1]
	M[1, i] = params.mu  # expectation
	M[2, i] = M[2, parentIndex] + sqrt(params.sigma2 * edge.length) * randn() # random value
end

function updateHybridSimulateBM!(M::Matrix, i::Int, parentIndex1::Int, parentIndex2::Int, edge1::Edge, edge2::Edge, params::Tuple{paramsBM})
	params = params[1]
       	M[1, i] = params.mu  # expectation
	M[2, i] =  edge1.gamma * (M[2, parentIndex1] + sqrt(params.sigma2 * edge1.length) * randn()) + edge2.gamma * (M[2, parentIndex2] + sqrt(params.sigma2 * edge2.length) * randn()) # random value
end


# function updateSimulateBM!(i::Int, nodes::Vector{Node}, M::Matrix, params::Tuple{paramsBM})
#     params = params[1]
#     parent = getParent(nodes[i]) #array of nodes (empty, size 1 or 2)
#     if(isempty(parent)) #nodes[i] is root
#         if (params.randomRoot)
# 		M[1, i] = params.mu # expectation
# 		M[2, i] = params.mu + sqrt(params.varRoot) * randn() # random value
# 	else
# 		M[1, i] = params.mu # expectation
# 		M[2, i] = params.mu # random value (root fixed)
# 	end
# 
#     elseif(length(parent) == 1) #nodes[i] is tree
#         parentIndex = getIndex(parent[1],nodes)
# 	l = getConnectingEdge(nodes[i],parent[1]).length
# 	M[1, i] = params.mu  # expectation
# 	M[2, i] = M[2, parentIndex] + sqrt(params.sigma2 * l) * randn() # random value
# 
#     elseif(length(parent) == 2) #nodes[i] is hybrid
#         parentIndex1 = getIndex(parent[1],nodes)
#         parentIndex2 = getIndex(parent[2],nodes)
#         edge1 = getConnectingEdge(nodes[i],parent[1])
#         edge2 = getConnectingEdge(nodes[i],parent[2])
#         edge1.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[1].number) should be a hybrid egde")
#         edge2.hybrid || error("connecting edge between node $(nodes[i].number) and $(parent[2].number) should be a hybrid egde")
# 	M[1, i] = params.mu  # expectation
# 	M[2, i] =  edge1.gamma * (M[2, parentIndex1] + sqrt(params.sigma2 * edge1.length) * randn()) + edge2.gamma * (M[2, parentIndex2] + sqrt(params.sigma2 * edge2.length) * randn()) # random value
#     end
# end

# Extract the vector of simulated values at the tips
function extractSimulateTips(sim::Matrix, net::HybridNetwork)
	mask = getTipsIndexes(net)
	return(squeeze(sim[2, mask], 1))
end


