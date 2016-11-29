# functions for trait evolution on network
# Claudia & Paul Bastide: November 2015

###############################################################################
###############################################################################
## Function to traverse the network in the pre-order, updating a matrix
###############################################################################
###############################################################################

# Matrix with rows and/or columns in topological order of the net.
"""
`matrixTopologicalOrder`

Matrix associated to an `HybridNetwork` sorted in topological order.

The following functions and extractors can be applied to it: `tipLabels`, `obj[:Tips]`, `obj[:InternalNodes]`, `obj[:TipsNodes]` (see documentation for function `getindex(obj, d,[ indTips, msng])`).

Functions `sharedPathMatrix` and `simulate` return objects of this type.

Has fields: `V`, `nodesNumbersTopOrder`, `internalNodesNumbers`, `tipsNumbers`, `tipsNames`, `indexation`.
Type in "?matrixTopologicalOrder.field" to get documentation on a specific field.

"""
type matrixTopologicalOrder
  "V: the matrix per se"
	V::Matrix # Matrix in itself
  "nodesNumbersTopOrder: vector of nodes numbers in the topological order, used for the matrix"
	nodesNumbersTopOrder::Vector{Int} # Vector of nodes numbers for ordering of the matrix
  "internalNodesNumbers: vector of internal nodes number, in the original net order"
	internalNodesNumbers::Vector{Int} # Internal nodes numbers (original net order)
  "tipsNumbers: vector of tips numbers, in the origial net order"
	tipsNumbers::Vector{Int} # Tips numbers (original net order)
  "tipsNames: vector of tips names, in the original net order"
	tipsNames::Vector # Tips Names (original net order)
  """
	indexation: a string giving the type of matrix `V`:
		-"r": rows only are indexed by the nodes of the network
		-"c": columns only are indexed by the nodes of the network
	  -"b": both rows and columns are indexed by the nodes of the network
	"""
	indexation::AbstractString # Are rows ("r"), columns ("c") or both ("b") indexed by nodes numbers in the matrix ?
end

function Base.show(io::IO, obj::matrixTopologicalOrder)
	println(io, "$(typeof(obj)):\n$(obj.V)")
end

function tipLabels(obj::matrixTopologicalOrder)
	return obj.tipsNames
end

# This function takes an init and update funtions as arguments
# It does the recursion using these functions on a preordered network.
function recursionPreOrder(
	net::HybridNetwork,
	checkPreorder=true::Bool,
	init=identity::Function,
	updateRoot=identity::Function,
	updateTree=identity::Function,
	updateHybrid=identity::Function,
	indexation="b"::AbstractString,
	params...
	)
	net.isRooted || error("net needs to be rooted to get matrix of shared path lengths")
	if(checkPreorder)
		preorder!(net)
	end
	M = recursionPreOrder(net.nodes_changed, init, updateRoot, updateTree, updateHybrid, params)
	# Find numbers of internal nodes
	nNodes = [n.number for n in net.node]
	nleaf = [n.number for n in net.leaf]
	deleteat!(nNodes, indexin(nleaf, nNodes))
	matrixTopologicalOrder(M, [n.number for n in net.nodes_changed], nNodes, nleaf, [n.name for n in net.leaf], indexation)
end

function recursionPreOrder(
	nodes::Vector{Node},
	init::Function,
	updateRoot::Function,
	updateTree::Function,
	updateHybrid::Function,
	params
	)
	n = length(nodes)
	M = init(nodes, params)
	for i in 1:n #sorted list of nodes
		updatePreOrder!(i, nodes, M, updateRoot, updateTree, updateHybrid, params)
	end
	return M
end

# Update on the network
# Takes three function as arguments : updateRoot, updateTree, updateHybrid
function updatePreOrder!(
	i::Int,
	nodes::Vector{Node},
	V::Matrix, updateRoot::Function,
	updateTree::Function,
	updateHybrid::Function,
	params
	)
	parent = getParents(nodes[i]) #array of nodes (empty, size 1 or 2)
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
# function getTipsIndexes(net::HybridNetwork)
# 	tipsNumbers = [n.number for n in net.leaf]
# 	nodesOrder = [n.number for n in net.nodes_changed]
#     getTipsIndexes(nodesOrder, tipsNumbers)
# end

# function getTipsIndexes(nodesOrder::Vector{Int64}, tipsNumbers::Vector{Int64})
# 	mask = BitArray(length(nodesOrder)) ## Function Match ??
# 	for tip in tipsNumbers
# 		mask = mask | (tip .== nodesOrder)
# 	end
# 	return(mask)
# end

# Extract the right part of a matrix in topological order
# Tips : submatrix corresponding to tips
# InternalNodes : submatrix corresponding to internal nodes
# TipsNodes : submatrix nTips x nNodes of interactions
# !! Extract sub-matrices in the original net nodes numbers !!
# function Base.getindex(obj::matrixTopologicalOrder, d::Symbol)
# 	if d == :Tips # Extract rows and/or columns corresponding to the tips
# 		mask = indexin(obj.tipsNumbers, obj.nodesNumbersTopOrder)
# 		obj.indexation == "b" && return obj.V[mask, mask] # both columns and rows are indexed by nodes
# 		obj.indexation == "c" && return obj.V[:, mask] # Only the columns
# 		obj.indexation == "r" && return obj.V[mask, :] # Only the rows
# 	end
# 	if d == :InternalNodes # Idem, for internal nodes
# 		mask = indexin(obj.internalNodesNumbers, obj.nodesNumbersTopOrder)
# 		obj.indexation == "b" && return obj.V[mask, mask]
# 		obj.indexation == "c" && return obj.V[:, mask] 
# 		obj.indexation == "r" && return obj.V[mask, :] 
# 	end
# 	if d == :TipsNodes
# 		maskNodes = indexin(obj.internalNodesNumbers, obj.nodesNumbersTopOrder)
# 		maskTips = indexin(obj.tipsNumbers, obj.nodesNumbersTopOrder,)
# 		obj.indexation == "b" && return obj.V[maskTips, maskNodes]
# 		obj.indexation == "c" && error("Both rows and columns must be net
# 		ordered to take the submatrix tips vs internal nodes.")
# 		obj.indexation == "r" && error("Both rows and columns must be net
# 		ordered to take the submatrix tips vs internal nodes.")
# 	end
# 	d == :All && return obj.V
# end

# If some tips are missing, treat them as "internal nodes"
"""
`getindex(obj, d,[ indTips, msng])`

Getting submatrices of an object of type `matrixTopologicalOrder`.

# Arguments
* `obj::matrixTopologicalOrder`: the matrix from which to extract.
* `d::Symbol`: a symbol precising which sub-matrix to extract. Can be:
  * `:Tips` columns and/or rows corresponding to the tips
  * `:InternalNodes` columns and/or rows corresponding to the internal nodes
  * `:TipsNodes` columns corresponding to internal nodes, and row to tips (works only is indexation="b")
* `indTips::Vector{Int}`: optional argument precising a specific order for the tips (internal use).
* `msng::BitArray{1}`: optional argument precising the missing tips (internal use).

"""
function Base.getindex(
	obj::matrixTopologicalOrder,
	d::Symbol,
	indTips=collect(1:length(obj.tipsNumbers))::Vector{Int},
	msng=trues(length(obj.tipsNumbers))::BitArray{1}
	)
	if d == :Tips # Extract rows and/or columns corresponding to the tips with data
		maskTips = indexin(obj.tipsNumbers, obj.nodesNumbersTopOrder)
		maskTips = maskTips[indTips]
		maskTips = maskTips[msng]
		obj.indexation == "b" && return obj.V[maskTips, maskTips] # both columns and rows are indexed by nodes
		obj.indexation == "c" && return obj.V[:, maskTips] # Only the columns
		obj.indexation == "r" && return obj.V[maskTips, :] # Only the rows
	end
	if d == :InternalNodes # Idem, for internal nodes
		maskNodes = indexin(obj.internalNodesNumbers, obj.nodesNumbersTopOrder)
		maskTips = indexin(obj.tipsNumbers, obj.nodesNumbersTopOrder)
		maskTips = maskTips[indTips]
		maskNodes = [maskNodes; maskTips[!msng]]
		obj.indexation == "b" && return obj.V[maskNodes, maskNodes]
		obj.indexation == "c" && return obj.V[:, maskNodes] 
		obj.indexation == "r" && return obj.V[maskNodes, :] 
	end
	if d == :TipsNodes
		maskNodes = indexin(obj.internalNodesNumbers, obj.nodesNumbersTopOrder)
		maskTips = indexin(obj.tipsNumbers, obj.nodesNumbersTopOrder)
		maskTips = maskTips[indTips]
		maskNodes = [maskNodes; maskTips[!msng]]
		maskTips = maskTips[msng]
		obj.indexation == "b" && return obj.V[maskTips, maskNodes]
		obj.indexation == "c" && error("Both rows and columns must be net
		ordered to take the submatrix tips vs internal nodes.")
		obj.indexation == "r" && error("Both rows and columns must be net
		ordered to take the submatrix tips vs internal nodes.")
	end
	d == :All && return obj.V
end

###############################################################################
###############################################################################
## Functions to compute the variance-covariance between Node and its parents
###############################################################################
###############################################################################
"""
`sharedPathMatrix(net::HybridNetwork; checkPreorder=true::Bool)`

This function computes the shared path matrix between all the nodes of a
network. It assumes that the network is in the pre-order. If checkPreorder is
true (default), then it runs function 'preoder' on the network beforehand.

Returns an object of type 'matrixTopologicalOrder'.

"""
function sharedPathMatrix(
	net::HybridNetwork;
	checkPreorder=true::Bool
	)
	recursionPreOrder(
		net,
		checkPreorder,
		initsharedPathMatrix,
		updateRootSharedPathMatrix!,
		updateTreeSharedPathMatrix!,
		updateHybridSharedPathMatrix!,
		"b")
end

function updateRootSharedPathMatrix!(V::Matrix, i::Int, params)
	return
end


function updateTreeSharedPathMatrix!(
	V::Matrix,
	i::Int,
	parentIndex::Int,
	edge::Edge,
	params
	)
	for j in 1:(i-1)
		V[i,j] = V[j,parentIndex]
		V[j,i] = V[j,parentIndex]
	end
	V[i,i] = V[parentIndex,parentIndex] + edge.length
end

function updateHybridSharedPathMatrix!(
	V::Matrix,
	i::Int,
	parentIndex1::Int,
	parentIndex2::Int,
	edge1::Edge,
	edge2::Edge,
	params
	)
	for j in 1:(i-1)
		V[i,j] = V[j,parentIndex1]*edge1.gamma + V[j,parentIndex2]*edge2.gamma
		V[j,i] = V[i,j]
	end
	V[i,i] = edge1.gamma*edge1.gamma*(V[parentIndex1,parentIndex1] + edge1.length) + edge2.gamma*edge2.gamma*(V[parentIndex2,parentIndex2] + edge2.length) + 2*edge1.gamma*edge2.gamma*V[parentIndex1,parentIndex2]
end


#function updateSharedPathMatrix!(i::Int,nodes::Vector{Node},V::Matrix, params)
#    parent = getParents(nodes[i]) #array of nodes (empty, size 1 or 2)
#    if(isempty(parent)) #nodes[i] is root
#        return
#    elseif(length(parent) == 1) #nodes[i] is tree
#        parentIndex = getIndex(parent[1],nodes)
#        for j in 1:(i-1)
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
#        for j in 1:(i-1)
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
# function extractVarianceTips(V::Matrix, net::HybridNetwork)
# 	mask = getTipsIndexes(net)
# 	return(V[mask, mask])
# end

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
#    for i in 1:n #sorted list of nodes
#        updateSharedPathMatrix!(i,net.nodes_changed,V)
#    end
#    return V
#end


###############################################################################
###############################################################################
## Types for params process
###############################################################################
###############################################################################

# Abstract type of all the (future) types (BM, OU, ...)
abstract paramsProcess

"""
`paramsBM <: paramsProcess`

Type for a BM process on a network. Fields are `mu` (expectation),
`sigma2` (variance), `randomRoot` (whether the root is random, default to `false`),
and `varRoot` (if the root is random, the variance of the root, defalut to `NaN`).

"""
# BM type
type paramsBM <: paramsProcess
	mu::Real # Ancestral value or mean
	sigma2::Real # variance
	randomRoot::Bool # Root is random ? default false
	varRoot::Real # root variance. Default NaN
end
# Constructor
paramsBM(mu, sigma2) = paramsBM(mu, sigma2, false, NaN) # default values

function Base.show(io::IO, obj::paramsBM)
	disp =  "$(typeof(obj)):\n"
	pt = paramstable(obj)
	if obj.randomRoot
		disp = disp * "Parameters of a BM with random root:\n" * pt
	else
		disp = disp * "Parameters of a BM with fixed root:\n" * pt
	end
	println(io, disp)
end

function paramstable(obj::paramsBM)
	disp = "mu: $(obj.mu)\nSigma2: $(obj.sigma2)"
	if obj.randomRoot
		disp = disp * "\nvarRoot: $(obj.varRoot)"
	end
	return(disp)
end


###############################################################################
###############################################################################
## Simulation Function
###############################################################################
###############################################################################

"""
`traitSimulation`

Result of a trait simulation on an `HybridNetwork` with function `simulate`.

The following functions and extractors can be applied to it: `tipLabels`, `obj[:Tips]`, `obj[:InternalNodes]` (see documentation for function `getindex(obj, d,[ indTips, msng])`).

Has fields: `M`, `params`, `model`.
"""
type traitSimulation
	M::matrixTopologicalOrder
	params::paramsProcess
	model::AbstractString
end

function Base.show(io::IO, obj::traitSimulation)
	disp = "$(typeof(obj)):\n"
	disp = disp * "Trait simulation results on a network with $(length(obj.M.tipsNames)) tips, using a using a $(obj.model) model, with parameters:\n"
	disp = disp * paramstable(obj.params)
	println(io, disp)
end

function tipLabels(obj::traitSimulation)
	return tipLabels(obj.M)
end


"""
`simulate(net::HybridNetwork, params::paramsProcess, checkPreorder=true::Bool)`

Simualte some traits on `net` using the parameters `params`. For now, only
parameters of type `paramsBM` (Brownian Motion) are accepted.

Assumes that the network is in the pre-order. If `checkPreorder=true` (default),
then it runs function `preoder` on the network beforehand.

Returns an object of type `traitSimulation`.

# Examples
```julia
julia> phy = readTopology("examples/carnivores_tree.txt");
julia> par = paramsBM(1, 0.1); # BM with expectation 1 and variance 0.1.
julia> par
julia> sim = simulate(phy, par); # Simulate on the tree.
julia> sim
julia> traits = sim[:Tips] # Extract simulated values at the tips.
```
"""
# Uses recursion on the network.
# Takes params of type paramsProcess as an entry
# Returns a matrix with two lines:
# - line one = expectations at all the nodes
# - line two = simulated values at all the nodes
# The nodes are ordered as given by topological sorting
function simulate(
	net::HybridNetwork,
	params::paramsProcess,
	checkPreorder=true::Bool)
	if isa(params, paramsBM)
		model = "BM"
	else
		error("The 'simulate' function only works for a BM process (for now).")
	end
	M = recursionPreOrder(
				net,
				checkPreorder,
				initSimulateBM,
				updateRootSimulateBM!,
				updateTreeSimulateBM!,
				updateHybridSimulateBM!,
				"c",
				params
				)
	traitSimulation(M, params, model)
end

# Initialization of the structure
function initSimulateBM(nodes::Vector{Node}, params::Tuple{paramsBM})
	return(zeros(2, length(nodes)))
end

# Initialization of the root
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

# Going down to a tree node
function updateTreeSimulateBM!(
	M::Matrix,
	i::Int,
	parentIndex::Int,
	edge::Edge,
	params::Tuple{paramsBM}
	)
	params = params[1]
	M[1, i] = params.mu  # expectation
	M[2, i] = M[2, parentIndex] + sqrt(params.sigma2 * edge.length) * randn() # random value
end

# Going down to an hybrid node
function updateHybridSimulateBM!(
	M::Matrix,
	i::Int,
	parentIndex1::Int, 
	parentIndex2::Int,
	edge1::Edge,
	edge2::Edge,
	params::Tuple{paramsBM}
	)
	params = params[1]
  M[1, i] = params.mu  # expectation
	M[2, i] =  edge1.gamma * (M[2, parentIndex1] + sqrt(params.sigma2 * edge1.length) * randn()) + edge2.gamma * (M[2, parentIndex2] + sqrt(params.sigma2 * edge2.length) * randn()) # random value
end


# function updateSimulateBM!(i::Int, nodes::Vector{Node}, M::Matrix, params::Tuple{paramsBM})
#     params = params[1]
#     parent = getParents(nodes[i]) #array of nodes (empty, size 1 or 2)
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
function Base.getindex(obj::traitSimulation, d::Symbol)
#    if d == :Tips
#       res = obj.M[:Tips]
#       squeeze(res[2, :], 1)
#    end
#	squeeze(getindex(obj.M, d)[2, :], 1)
 getindex(obj.M, d)[2, :]
end

# function extractSimulateTips(sim::Matrix, net::HybridNetwork)
# 	mask = getTipsIndexes(net)
# 	return(squeeze(sim[2, mask], 1))
# end

###############################################################################
###############################################################################
## Functions for Phylgenetic Network regression
###############################################################################
###############################################################################

# New type for phyloNetwork regression
"""
`phyloNetworkLinearModel<:LinPredModel`

Regression object for a phylogenetic regression. Result of fitting function `phyloNetworklm`.
Dominated by the `LinPredModel` class, from package `GLM`.

The following StatsBase functions can be applied to it:
`coef`, `nobs`, `vcov`, `stderr`, `confint`, `coeftable`, `dof_residual`, `dof`, `deviance`,
`residuals`, `model_response`, `predict`, `loglikelihood`, `nulldeviance`, `nullloglikelihood`,
`r2`, `adjr2`, `aic`, `aicc`, `bic`.

The following DataFrame functions can also be applied to it:
`ModelFrame`, `ModelMatrix`, `Formula`.

Estimated variance and mean of the BM process used can be retrieved with 
functions `sigma2_estim` and `mu_estim`.

If a Pagel's lambda model is fitted, the parameter can be retrieved with function 
	`lambda_estim`.

An ancestral state reconstruction can be performed from this fitted object using function:
	`ancestralStateReconstruction`.

Has fields: `lm`, `V`, `Vy`, `RL`, `Y`, `X`, `logdetVy`, `ind`, `msng`, `model`, `lambda`.
Type in "?phyloNetworkLinearModel.field" to get help on a specific field.
"""
type phyloNetworkLinearModel <: LinPredModel 
  "lm: a GLM.LinearModel object, fitted on the cholesky-tranformend problem"
	lm::GLM.LinearModel # result of a lm on a matrix
  "V: a matrixTopologicalOrder object of the network-induced correlations"
	V::matrixTopologicalOrder
  "Vy: the sub matrix corresponding to the tips and actually used for the correction"
	Vy::Matrix
  "RL: a LowerTriangular matrix, Cholesky transform of Vy=RL*RL'"
	RL::LowerTriangular
  "Y: the vector of data"
	Y::Vector
  "X: the matrix of regressors"
	X::Matrix
  "logdetVy: the log-determinent of Vy"
	logdetVy::Real
  "ind: vector matching the tips of the network against the names of the dataframe provided. 0 if the match could not be performed."
	ind::Vector{Int} 
  "msng: vector indicating which of the tips are missing"
	msng::BitArray{1} # Which tips are not missing
	"model: the model used for the fit"
	model::AbstractString
	"If applicable, value of lambda (default to 1)."
	lambda::Real
end

phyloNetworkLinearModel(lm_fit, V, Vy, RL, Y, X, logdetVy, ind, msng, model) = phyloNetworkLinearModel(lm_fit, V, Vy, RL, Y, X, logdetVy, ind, msng, model, 1.0) 

# Function for lm with net residuals
function phyloNetworklm(
	X::Matrix,
	Y::Vector,
	net::HybridNetwork,
	msng=trues(length(Y))::BitArray{1}, # Which tips are not missing ?
	model="BM"::AbstractString,
	ind=[0]::Vector{Int}
	)
	# Geting variance covariance
	V = sharedPathMatrix(net)
	# Fit
	phyloNetworklm(X, Y, V, msng, model, ind)
end

# Same function, but when the matrix V is already known.
function phyloNetworklm(
	X::Matrix,
	Y::Vector,
	V::matrixTopologicalOrder,
	msng=trues(length(Y))::BitArray{1}, # Which tips are not missing ?
	model="BM"::AbstractString,
	ind=[0]::Vector{Int}
	)
	## Choose Model
	if (model == "BM")
		return phyloNetworklm_BM(X, Y, V, msng, ind) 
	end
	if (model == "lambda")
		return phyloNetworklm_lambda(X, Y, V, msng, ind) 
	end
end

###############################################################################
## Fit BM

function phyloNetworklm_BM(
	X::Matrix,
	Y::Vector,
	V::matrixTopologicalOrder,
	msng=trues(length(Y))::BitArray{1}, # Which tips are not missing ?
	ind=[0]::Vector{Int}
	)
   	# Extract tips matrix
	Vy = V[:Tips]
	# Re-order if necessary
	if (ind != [0]) Vy = Vy[ind, ind] end
	# Keep only not missing values
	Vy = Vy[msng, msng]
	# Cholesky decomposition
	R = cholfact(Vy)
	RL = R[:L]
	# Fit
	phyloNetworkLinearModel(lm(RL\X, RL\Y), V, Vy, RL, Y, X, logdet(Vy), ind, msng, "BM")
end

###############################################################################
## Fit Pagel's Lambda

function transform_matrix_lambda!{T <: AbstractFloat}(V::matrixTopologicalOrder, lam::T)
	# WARNING : This transformation is not the expected one if branch length are modified.
	# Need a function for node heigh computation.
	for i in 1:size(V.V, 1)
		for j in 1:size(V.V, 2)
			if i != j
				V.V[i,j] *= lam
			end
		end
	end
#	V_diag = diagm(diag(V.V))
#	V.V = lam * V.V + (1 - lam) * V_diag
end

function logLik_lam{T <: AbstractFloat}(
	lam::T,
	X::Matrix,
	Y::Vector,
	V::matrixTopologicalOrder,
	msng=trues(length(Y))::BitArray{1}, # Which tips are not missing ?
	ind=[0]::Vector{Int}
	)
	# Transform V according to lambda
	transform_matrix_lambda!(V, lam)
	# Fit and take likelihood
	fit_lam = phyloNetworklm_BM(X, Y, V, msng, ind)
	res = - loglikelihood(fit_lam)
	# Go back to original V
	transform_matrix_lambda!(V, 1/lam)
	return res
end

# Code for optim taken from PhyloNetworks.jl/src/optimization.jl, lines 276 - 331
const fAbsTr = 1e-10 
const fRelTr = 1e-10
const xAbsTr = 1e-10 
const xRelTr = 1e-10

function phyloNetworklm_lambda(
	X::Matrix,
	Y::Vector,
	V::matrixTopologicalOrder,
	msng=trues(length(Y))::BitArray{1}, # Which tips are not missing ?
	ind=[0]::Vector{Int},
	ftolRel=fRelTr::AbstractFloat,
	xtolRel=xRelTr::AbstractFloat,
	ftolAbs=fAbsTr::AbstractFloat,
	xtolAbs=xAbsTr::AbstractFloat,
	startingValue=0.5::Real
	)
	# Find Best lambda using optimize from package NLopt
	opt = NLopt.Opt(:LN_BOBYQA, 1)
	NLopt.ftol_rel!(opt, ftolRel) # relative criterion
  NLopt.ftol_abs!(opt, ftolAbs) # absolute critetion 
	NLopt.xtol_rel!(opt, xtolRel) # criterion on parameter value changes
	NLopt.xtol_abs!(opt, xtolAbs) # criterion on parameter value changes
	NLopt.maxeval!(opt, 1000) # max number of iterations
	NLopt.lower_bounds!(opt, 1e-100) # Lower bound  
	NLopt.upper_bounds!(opt, 1.0)
	count = 0
  function fun(x::Vector{Float64}, g::Vector{Float64})
		x = convert(AbstractFloat, x[1])
		res = logLik_lam(x, X, Y, V, msng, ind)
		count =+ 1
		return res
	end
	NLopt.min_objective!(opt, fun)
	fmin, xmin, ret = NLopt.optimize(opt, [startingValue])
	# Best value dans result
	transform_matrix_lambda!(V, xmin[1])
	res = phyloNetworklm_BM(X, Y, V, msng, ind)
	res.lambda = xmin[1]
	res.model = "lambda"
	return res
end

"""
    phyloNetworklm(f, fr, net, model="BM",
		fTolRel=1e^-10, fTolAbs=1e^-10, xTolRel=1e^-10, xTolAbs=1e^-10,
		startingValue=0.5)`

Phylogenetic regression, using the correlation structure induced by the network.

Returns an object of class `phyloNetworkLinearModel`. See documentation for this type and
example to see all the functions that can be applied to it.

# Arguments
* `f::Formula`: formula to use for the regression (see the `DataFrame` package)
* `fr::AbstractDataFrame`: DataFrame containing the data and regressors at the tips. It should have an extra column labelled "tipsNames", that gives the names of the taxa for each observation.
* `net::HybridNetwork`: phylogenetic network to use. Should have labelled tips.
* `model::AbstractString="BM"`: the model to use, "BM" being the default and only available model for now. If the entry is a TREE, then "lambda" can fit a Pagel's lambda model.
* `no_names::Bool=false`: if `true`, force the function to ignore the tips names. The data is then assumed to be in the same order as the tips of the network. Default to false, setting it to true is dangerous, and strongly discouraged.
If `model="lambda"`, there are a few more parameters to control the optimization in the parameter:
* `fTolRel::AbstractFloat=1e-10`: relative tolerance on the likelihood value for the optimization in lambda.
* `fTolAbs::AbstractFloat=1e-10`: absolute tolerance on the likelihood value for the optimization in lambda.
* `xTolRel::AbstractFloat=1e-10`: relative tolerance on the parameter value for the optimization in lambda.
* `xTolAbs::AbstractFloat=1e-10`: absolute tolerance on the parameter value for the optimization in lambda.
* `startingValue::Real=0.5`: the starting value for the parameter in the optimization in lambda.

# Example 
```julia
julia> phy = readTopology("examples/caudata_tree.txt");
julia> dat = readtable("examples/caudata_trait.txt");
julia> fitBM = phyloNetworklm(trait ~ 1, dat, phy);
julia> fitBM # Shows a summary
julia> sigma2_estim(fitBM)
julia> mu_estim(fitBM)
julia> loglikelihood(fitBM) 
julia> aic(fitBM)
julia> aicc(fitBM)
julia> bic(fitBM)
julia> coef(fitBM)
julia> confint(fitBM)
julia> r2(fitBM)
julia> adjr2(fitBM)
julia> vcov(fitBM)
julia> residuals(fitBM)
julia> model_response(fitBM)
julia> predict(fitBM)
```

# See also 
Type `phyloNetworkLinearModel`, Function `ancestralStatesReconstruction`
""" #"
# Deal with formulas
function phyloNetworklm(
	f::Formula,
	fr::AbstractDataFrame,
	net::HybridNetwork;
	model="BM"::AbstractString,
	no_names=false::Bool,
	ftolRel=fRelTr::AbstractFloat,
	xtolRel=xRelTr::AbstractFloat,
	ftolAbs=fAbsTr::AbstractFloat,
	xtolAbs=xAbsTr::AbstractFloat,
	startingValue=0.5::Real
	)
	# Match the tips names: make sure that the data provided by the user will
	# be in the same order as the ordered tips in matrix V.
	V = sharedPathMatrix(net)
	if no_names # The names should not be taken into account.
		ind = [0]
		info("As requested (no_names=true), I am ignoring the tips names on the
		network and in the dataframe.")
	elseif (any(V.tipsNames == "") || !any(DataFrames.names(fr) .== :tipsNames))
		if (any(V.tipsNames == "") && !any(DataFrames.names(fr) .== :tipsNames))
			error("The network provided has no tip names, and the input dataframe has
			no column labelled tipsNames, so I can't match the data on the network
			unambiguously. If you are sure that the tips of the network are in the
			same order as the values of the dataframe provided, then please re-run
			this function with argument no_name=true.")
		end
		if any(V.tipsNames == "")
			error("The network provided has no tip names, so I can't match the data
			on the network unambiguously. If you are sure that the tips of the
			network are in the same order as the values of the dataframe provided,
			then please re-run this function with argument no_name=true.")
		end
		if !any(DataFrames.names(fr) .== :tipsNames)
			error("The input dataframe has no column labelled tipsNames, so I can't
			match the data on the network unambiguously. If you are sure that the
			tips of the network are in the same order as the values of the dataframe
			provided, then please re-run this function with argument no_name=true.")
		end
	else
#        ind = indexin(V.tipsNames, fr[:tipsNames])
		ind = indexin(fr[:tipsNames], V.tipsNames)
		if any(ind == 0) || length(unique(ind)) != length(ind)
				error("Tips names of the network and names provided in column tipsNames
				of the dataframe do not match.")
		end
#	fr = fr[ind, :]
	end
	# Find the regression matrix and answer vector
	mf = ModelFrame(f,fr)
	mm = ModelMatrix(mf)
	Y = convert(Vector{Float64},DataFrames.model_response(mf))
	# Fit the model (Method copied from DataFrame/src/statsmodels/statsmodels.jl, lines 47-58)
	DataFrames.DataFrameRegressionModel(phyloNetworklm(mm.m, Y, V, mf.msng, model, ind), mf, mm)
#    # Create the object
#    phyloNetworkLinPredModel(DataFrames.DataFrameRegressionModel(fit, mf, mm),
#    fit.V, fit.Vy, fit.RL, fit.Y, fit.X, fit.logdetVy, ind, mf.msng)
end

### Methods on type phyloNetworkRegression

## Un-changed Quantities
# Coefficients of the regression
StatsBase.coef(m::phyloNetworkLinearModel) = coef(m.lm) 
# Number of observations
StatsBase.nobs(m::phyloNetworkLinearModel) = nobs(m.lm) 
# vcov matrix
StatsBase.vcov(m::phyloNetworkLinearModel) = vcov(m.lm)
# Standart error
StatsBase.stderr(m::phyloNetworkLinearModel) = stderr(m.lm)
# Confidence Intervals
StatsBase.confint(m::phyloNetworkLinearModel; level=0.95::Real) = confint(m.lm, level)
# coef table (coef, stderr, confint)
StatsBase.coeftable(m::phyloNetworkLinearModel) = coeftable(m.lm)
# Degrees of freedom for residuals
StatsBase.dof_residual(m::phyloNetworkLinearModel) =  nobs(m) - length(coef(m))
# Degrees of freedom consumed in the model
function StatsBase.dof(m::phyloNetworkLinearModel)
	res = length(coef(m)) + 1 # (+1: dispersion parameter)
	if (m.model == "lambda")
		res += 1 # lambda is one parameter
	end
	return res
end
# Deviance (sum of squared residuals with metric V)
StatsBase.deviance(m::phyloNetworkLinearModel) = deviance(m.lm)

## Changed Quantities
# Compute the residuals
# (Rescaled by cholesky of variance between tips)
StatsBase.residuals(m::phyloNetworkLinearModel) = m.RL * residuals(m.lm)
# Tip data
StatsBase.model_response(m::phyloNetworkLinearModel) = m.Y
# Predicted values at the tips
# (rescaled by cholesky of tips variances)
StatsBase.predict(m::phyloNetworkLinearModel) = m.RL * predict(m.lm)
#Log likelihood of the fitted linear model
StatsBase.loglikelihood(m::phyloNetworkLinearModel) =  loglikelihood(m.lm) - 1/2 * m.logdetVy
# Null  Deviance (sum of squared residuals with metric V)
# REMARK Not just the null deviance of the cholesky regression
# Might be something better to do than this, though.
function StatsBase.nulldeviance(m::phyloNetworkLinearModel)
	vo = ones(length(m.Y), 1)
	vo = m.RL \ vo
	bo = inv(vo'*vo)*vo'*model_response(m.lm)
	ro = model_response(m.lm) - vo*bo
	return sum(ro.^2)
end
# Null Log likelihood (null model with only the intercept)
# Same remark
function StatsBase.nullloglikelihood(m::phyloNetworkLinearModel)
	n = length(m.Y)
	return -n/2 * (log(2*pi * nulldeviance(m)/n) + 1) - 1/2 * m.logdetVy
end
# coefficient of determination (1 - SS_res/SS_null)
# Copied from GLM.jl/src/lm.jl, line 139
StatsBase.r2(m::phyloNetworkLinearModel) = 1 - deviance(m)/nulldeviance(m)
# adjusted coefficient of determination 
# Copied from GLM.jl/src/lm.jl, lines 141-146
function StatsBase.adjr2(obj::phyloNetworkLinearModel)
	n = nobs(obj)
	# dof() includes the dispersion parameter
	p = dof(obj) - 1
	1 - (1 - r2(obj))*(n-1)/(n-p)
end

## REMARK
# As phyloNetworkLinearModel <: LinPredModel, the following functions are automatically defined:
# aic, aicc, bic

## New quantities
# ML estimate for variance of the BM
sigma2_estim(m::phyloNetworkLinearModel) = deviance(m.lm) / nobs(m)
# Need to be adapted manually to DataFrameRegressionModel beacouse it's a new function
sigma2_estim(m::DataFrames.DataFrameRegressionModel) = sigma2_estim(m.model)
# ML estimate for ancestral state of the BM
function mu_estim(m::phyloNetworkLinearModel)
	warn("You fitted the data against a custom matrix, so I have no way of knowing which column is your intercept (column of ones). I am using the first coefficient for ancestral mean mu by convention, but that might not be what you are looking for.")
	return coef(m)[1]
end
# Need to be adapted manually to DataFrameRegressionModel beacouse it's a new function
function mu_estim(m::DataFrames.DataFrameRegressionModel)#{PhyloNetworks.phyloNetworkLinearModel,Float64})
	if (!m.mf.terms.intercept)
		error("The fit was done without intercept, so I cannot estimate mu")
	end
	return coef(m)[1]
end
# Lambda estim
lambda_estim(m::phyloNetworkLinearModel) = m.lambda
lambda_estim(m::DataFrames.DataFrameRegressionModel) = lambda_estim(m.model)

### Functions specific to DataFrameRegressionModel
DataFrames.ModelFrame(m::DataFrames.DataFrameRegressionModel) = m.mf
DataFrames.ModelMatrix(m::DataFrames.DataFrameRegressionModel) = m.mm
DataFrames.Formula(m::DataFrames.DataFrameRegressionModel) = Formula(m.mf.terms)

### Print the results
# Variance
function paramstable(m::phyloNetworkLinearModel)
	Sig = sigma2_estim(m)
	res = "Sigma2: $(Sig)"
	if (m.model == "lambda")
	  Lamb = lambda_estim(m)
		res = res*"\nLambda: $(Lamb)"
	end
	return(res)
end
function Base.show(io::IO, obj::phyloNetworkLinearModel)
	println(io, "$(typeof(obj)):\n\nParameter(s) Estimates:\n", paramstable(obj), "\n\nCoefficients:\n", coeftable(obj))
end
# For DataFrameModel. Copied from DataFrames/jl/src/statsmodels/statsmodels.jl, lines 101-118
function Base.show(io::IO, model::DataFrames.DataFrameRegressionModel)#{PhyloNetworks.phyloNetworkLinearModel,Float64})
	ct = coeftable(model)
	println(io, "$(typeof(model))")
	println(io)
	println(io, Formula(model.mf.terms))
	println(io)
	println(io,"Parameter(s) Estimates:")
	println(io, paramstable(model.model))
	println(io)
	println(io,"Coefficients:")
	show(io, ct)
	println(io)
	println(io, "Log Likelihood: "*"$(loglikelihood(model))")
	println(io, "AIC: "*"$(aic(model))")
end


## Deprecated
# function StatsBase.vcov(obj::phyloNetworkLinearModel)
#    sigma2_estim(obj) * inv(obj.X' * obj.X) 
# end
#function StatsBase.vcov(obj::phyloNetworkLinPredModel)
#   sigma2_estim(obj) * inv(obj.X' * obj.X) 
#end
#StatsBase.stderr(m::phyloNetworkLinPredModel) = sqrt(diag(vcov(m)))
# Confidence intervals on coeficients
# function StatsBase.confint(obj::phyloNetworkLinearModel, level=0.95::Real)
#     hcat(coef(obj),coef(obj)) + stderr(obj) *
#     quantile(TDist(dof_residual(obj)), (1. - level)/2.) * [1. -1.]
# end
# Log likelihood of the fitted BM
# StatsBase.loglikelihood(m::phyloNetworkLinearModel) = - 1 / 2 * (nobs(m) + nobs(m) * log(2 * pi) + nobs(m) * log(sigma2_estim(m)) + m.logdetVy)
#StatsBase.loglikelihood(m::phyloNetworkLinPredModel) = - 1 / 2 * (nobs(m) + nobs(m) * log(2 * pi) + nobs(m) * log(sigma2_estim(m)) + m.logdetVy)
# Coefficients
# function StatsBase.coeftable(mm::phyloNetworkLinearModel)
#     cc = coef(mm)
#     se = stderr(mm)
#     tt = cc ./ se
#     CoefTable(hcat(cc,se,tt,ccdf(FDist(1, dof_residual(mm)), abs2(tt))),
#               ["Estimate","Std.Error","t value", "Pr(>|t|)"],
#               ["x$i" for i = 1:size(mm.lm.pp.X, 2)], 4)
# end

###############################################################################
###############################################################################
## Ancestral State Reconstruction
###############################################################################
###############################################################################
# Class for reconstructed states on a network
"""
`reconstructedStates`

Type containing the inferred information about the law of the ancestral states
given the observed tips values. The missing tips are considered as ancestral states.

The following functions can be applied to it:
`expectations` (vector of expectations at all nodes), `stderr` (the standard error),
`predint` (the prediction interval), `plot`.

Has fields: `traits_nodes`, `variances_nodes`, `NodesNumbers`, `traits_tips`, `TipsNumbers`, `model`.
Type in "?reconstructedStates.field" to get help on a specific field.
"""
type reconstructedStates
  "traits_nodes: the infered expectation of 'missing' values (ancestral nodes and missing tips)"
	traits_nodes::Vector # Nodes are actually "missing" data (including tips)
  "variances_nodes: the variance covariance matrix between all the 'missing' nodes"
	variances_nodes::Matrix
  "NodesNumbers: vector of the nodes numbers, in the same order as `traits_nodes`"
	NodesNumbers::Vector{Int}
  "traits_tips: the observed traits values at the tips"
	traits_tips::Vector # Observed values at tips
  "TipsNumbers: vector of tips numbers, in the same order as `traits_tips`"
	TipsNumbers::Vector # Observed tips only
  "model (Nullable): if not null, the `phyloNetworkLinearModel` used for the computations."
	model::Nullable{phyloNetworkLinearModel} # If empirical, the corresponding fitted object.
end

function expectations(obj::reconstructedStates)
	return DataFrame(nodeNumber = [obj.NodesNumbers; obj.TipsNumbers], condExpectation = [obj.traits_nodes; obj.traits_tips])
end

StatsBase.stderr(obj::reconstructedStates) = sqrt(diag(obj.variances_nodes))

function predint(obj::reconstructedStates, level=0.95::Real)
	if isnull(obj.model)
		qq = quantile(Normal(), (1. - level)/2.)
	else
		qq = quantile(TDist(dof_residual(get(obj.model))), (1. - level)/2.)
		warn("As the variance is estimated, the predictions intervals are not exact, and should probably be larger.")
	end
	tmpnode = hcat(obj.traits_nodes, obj.traits_nodes) + stderr(obj) * qq * [1. -1.]
	return vcat(tmpnode, hcat(obj.traits_tips, obj.traits_tips))
end

function Base.show(io::IO, obj::reconstructedStates)
	println(io, "$(typeof(obj)):\n",
	CoefTable(hcat(vcat(obj.NodesNumbers, obj.TipsNumbers), vcat(obj.traits_nodes, obj.traits_tips), predint(obj)),
  					["Node index", "Pred.", "Min.", "Max. (95%)"],
						fill("", length(obj.NodesNumbers)+length(obj.TipsNumbers))))
end

function predintPlot(obj::reconstructedStates, level=0.95::Real)
	pri = predint(obj, level)
	pritxt = Array{AbstractString}(size(pri, 1))
	for i=1:length(obj.NodesNumbers)
		pritxt[i] = "[" * string(round(pri[i, 1], 2)) * ", " * string(round(pri[i, 2], 2)) * "]"
	end
	for i=(length(obj.NodesNumbers)+1):size(pri, 1)
		pritxt[i] = string(round(pri[i, 1], 2))
	end
	return DataFrame(nodeNumber = [obj.NodesNumbers; obj.TipsNumbers], PredInt = pritxt)
end

"""
'plot(net::HybridNetwork, obj::reconstructedStates; kwargs...)

Plot the reconstructed states computed by function `ancestralStateReconstruction`
on a network.

# Arguments
* `net::HybridNetwork`: a phylogenetic network.
* `obj::reconstructedStates`: the reconstructed states on the network. 
* `kwargs...`: further arguments to be passed to the netwotk `plot` function.

See documentation for function `ancestralStateReconstruction(obj::phyloNetworkLinearModel[, X_n::Matrix])` for examples.

"""
function Gadfly.plot(net::HybridNetwork, obj::reconstructedStates; kwargs...)
	plot(net, nodeLabel = predintPlot(obj); kwargs...)
end

"""
`ancestralStateReconstruction(net::HybridNetwork, Y::Vector, params::paramsBM)`

Compute the conditional expectations and variances of the ancestral (un-observed)
traits values at the internal nodes of the phylogenetic network (`net`), 
given the values of the traits at the tips of the network (`Y`) and some
known parameters of the process used for trait evolution (`params`, only BM with fixed root
works for now).

This function assumes that the parameters of the process are known. For a more general
function, see `ancestralStateReconstruction(obj::phyloNetworkLinearModel[, X_n::Matrix])`.

"""
# Reconstruction from known BM parameters
function ancestralStateReconstruction(
	net::HybridNetwork,
  Y::Vector,
  params::paramsBM
	)
	V = sharedPathMatrix(net)
	ancestralStateReconstruction(V, Y, params)
end

function ancestralStateReconstruction(
	V::matrixTopologicalOrder,
	Y::Vector,
	params::paramsBM
	)
	# Variances matrices
	Vy = V[:Tips]
	Vz = V[:InternalNodes]
	Vyz = V[:TipsNodes]
	R = cholfact(Vy)
	RL = R[:L]
	temp = RL \ Vyz
	# Vectors of means
	m_y = ones(size(Vy)[1]) .* params.mu # !! correct only if no predictor. 
	m_z = ones(size(Vz)[1]) .* params.mu # !! works if BM no shift.
	# Actual computation
	ancestralStateReconstruction(
	Vz, temp, RL,
	Y, m_y, m_z,
	V.internalNodesNumbers,
	V.tipsNumbers,
	params.sigma2
	)
end

# Reconstruction from all the needed quantities
function ancestralStateReconstruction(
	Vz::Matrix,
  VyzVyinvchol::Matrix,
  RL::LowerTriangular,
  Y::Vector, m_y::Vector, m_z::Vector,
	NodesNumbers::Vector,
	TipsNumbers::Vector,
	sigma2::Real,
	add_var=zeros(size(Vz))::Matrix, # Additional variance for BLUP
	model=Nullable{phyloNetworkLinearModel}()::Nullable{phyloNetworkLinearModel}
	) 
	m_z_cond_y = m_z + VyzVyinvchol' * (RL \ (Y - m_y))
	V_z_cond_y = sigma2 .* (Vz - VyzVyinvchol' * VyzVyinvchol)
	reconstructedStates(m_z_cond_y, V_z_cond_y + add_var, NodesNumbers, Y, TipsNumbers, model)
end

# """
# `ancestralStateReconstruction(obj::phyloNetworkLinearModel, X_n::Matrix)`
# Function to find the ancestral traits reconstruction on a network, given an
# object fitted by function phyloNetworklm, and some predictors expressed at all the nodes of the network.
# 
# - obj: a phyloNetworkLinearModel object, or a
# DataFrameRegressionModel{phyloNetworkLinearModel}, if data frames were used.
# - X_n a matrix with as many columns as the number of predictors used, and as
# many lines as the number of unknown nodes or tips.
# 
# Returns an object of type ancestralStateReconstruction.
# """

# Empirical reconstruciton from a fitted object
# TO DO: Handle the order of internal nodes for matrix X_n
function ancestralStateReconstruction(obj::phyloNetworkLinearModel, X_n::Matrix)
	if (size(X_n)[2] != length(coef(obj)))
		error("The number of predictors for the ancestral states (number of columns of X_n) do not match the number of predictors at the tips.")
	end
	if (size(X_n)[1] != length(obj.V.internalNodesNumbers) + sum(!obj.msng))
		error("The number of lines of the predictors do not match the number of nodes plus the number of missing tips.")
	end
	m_y = predict(obj)
	m_z = X_n * coef(obj) 
	# If the tips were re-organized, do the same for Vyz
	if (obj.ind != [0])
# 		iii = indexin(1:length(obj.msng), obj.ind[obj.msng])
# 		iii = iii[iii .> 0]
# 		jjj = [1:length(obj.V.internalNodesNumbers); indexin(1:length(obj.msng), obj.ind[!obj.msng])]
# 		jjj = jjj[jjj .> 0]
# 		Vyz = Vyz[iii, jjj]
		Vyz = obj.V[:TipsNodes, obj.ind, obj.msng]
		missingTipsNumbers = obj.V.tipsNumbers[obj.ind][!obj.msng]
		nmTipsNumbers = obj.V.tipsNumbers[obj.ind][obj.msng]
	else
		warn("There were no indication for the position of the tips on the network. Assuming that they are given in the same order. Please check that this is what you intended.")
		Vyz = obj.V[:TipsNodes, collect(1:length(obj.V.tipsNumbers)), obj.msng]
		missingTipsNumbers = obj.V.tipsNumbers[!obj.msng]
		nmTipsNumbers = obj.V.tipsNumbers[obj.msng]
	end
	temp = obj.RL \ Vyz
	U = X_n - temp' * (obj.RL \ obj.X)
	add_var = U * vcov(obj) * U'
	ancestralStateReconstruction(
		obj.V[:InternalNodes, obj.ind, obj.msng],
    temp,
    obj.RL,
    obj.Y,
		m_y,
		m_z,
		[obj.V.internalNodesNumbers; missingTipsNumbers],
		nmTipsNumbers, 
		sigma2_estim(obj),
		add_var,
		obj
		)
end

"""
`ancestralStateReconstruction(obj::phyloNetworkLinearModel[, X_n::Matrix])`

Function to find the ancestral traits reconstruction on a network, given an
object fitted by function `phyloNetworklm`. By default, the function assumes
that the regressor is just an intercept. If the value of the regressor for 
all the ancestral states is known, it can be entered in X_n, a matrix with as
many columns as the number of predictors used, and as many lines as the number
of unknown nodes or tips.

Returns an object of type `reconstructedStates`.
See documentation for this type and examples for functions that can be applied to it.

# Example 
```julia
julia> phy = readTopology("examples/carnivores_tree.txt");
julia> dat = readtable("examples/carnivores_trait.txt");
julia> fitBM = phyloNetworklm(trait ~ 1, dat, phy);
julia> ancStates = ancestralStateReconstruction(fitBM); 
julia> ancStates # Should produce a warning, as variance is unknown.
julia> plot(phy, ancStates)
julia> expectations(ancStates)
julia> predint(ancStates)
julia> ## Some tips may also be missing
julia> dat[[2, 5], :trait] = NA
julia> fitBM = phyloNetworklm(trait ~ 1, dat, phy);
julia> ancStates = ancestralStateReconstruction(fitBM); 
julia> plot(phy, ancStates)
julia> expectations(ancStates)
julia> predint(ancStates)
```
"""
# Default reconstruction for a simple BM (known predictors)
function ancestralStateReconstruction(obj::phyloNetworkLinearModel)
	if ((size(obj.X)[2] != 1) || !any(obj.X .== 1)) # Test if the regressor is just an intercept.
		error("As the predictor is not reduced to the intercept, I can't guess the ancestral predictors values. If you know all the ancestral predictors, please provide them as a matrix argument to the function. Otherwise, you might consider doing a multivariate linear regression (not implemented yet).")
	end
  X_n = ones((length(obj.V.internalNodesNumbers) + sum(!obj.msng), 1))
	ancestralStateReconstruction(obj, X_n)
end
# For a DataFrameRegressionModel
function ancestralStateReconstruction{T<:Union{Float32,Float64}}(obj::DataFrames.DataFrameRegressionModel{PhyloNetworks.phyloNetworkLinearModel, Array{T,2}})
	ancestralStateReconstruction(obj.model)
end
function ancestralStateReconstruction{T<:Union{Float32,Float64}}(obj::DataFrames.DataFrameRegressionModel{PhyloNetworks.phyloNetworkLinearModel, Array{T,2}}, X_n::Matrix)
	ancestralStateReconstruction(obj.model, X_n::Matrix)
end

"""
`ancestralStateReconstruction(fr::AbstractDataFrame, net::HybridNetwork; kwargs...)`

Function to find the ancestral traits reconstruction on a network, given some data at the tips.
Uses function `phyloNetworklm` to perform a phylogenetic regression of the data against an
intercept (amounts to fitting an evolutionary model on the network, BM being the only option 
available for now).

See documentation on `phyloNetworklm` and `ancestralStateReconstruction(obj::phyloNetworkLinearModel[, X_n::Matrix])`
for further details.

Returns an object of type `reconstructedStates`.
"""
# Deal with formulas
function ancestralStateReconstruction(
	fr::AbstractDataFrame,
	net::HybridNetwork;
	kwargs...
	)
	nn = names(fr)
	datpos = nn .!= :tipsNames
  if sum(datpos) > 1
		error("Beside one column labelled 'tipsNames', the dataframe fr should have only one column, corresponding to the data at the tips of the network.")
	end
	f = Formula(nn[datpos][1], 1)
	reg = phyloNetworklm(f, fr, net; kwargs...)
	return ancestralStateReconstruction(reg)
end

# # Default reconstruction for a simple BM
# function ancestralStateReconstruction(obj::phyloNetworkLinearModel, mu::Real)
# 	m_y = predict(obj)
#   m_z = ones(size(obj.V.internalNodesNumbers)[1]) .* mu # !! works if BM no shift.
# 	ancestralStateReconstruction(obj, m_y, m_z)
# end
# # Handling the ancestral mean
# function ancestralStateReconstruction(obj::phyloNetworkLinearModel)
# 	mu = mu_estim(obj) # Produce a warning if intercept is unknown.
# 	ancestralStateReconstruction(obj, mu)
# end
# # For a DataFrameRegressionModel
# function ancestralStateReconstruction(obj::DataFrames.DataFrameRegressionModel{PhyloNetworks.phyloNetworkLinearModel,Float64})
# 	mu = mu_estim(obj) # Produces an error if intercept is not here.
# 	ancestralStateReconstruction(obj.model, mu)
# end

#################################################
## Old version of phyloNetworklm (naive) 
#################################################

# function phyloNetworklmNaive(X::Matrix, Y::Vector, net::HybridNetwork, model="BM"::AbstractString)
# 	# Geting variance covariance
# 	V = sharedPathMatrix(net)
# 	Vy = extractVarianceTips(V, net)
# 	# Needed quantities (naive)
# 	ntaxa = length(Y)
# 	Vyinv = inv(Vy)
# 	XtVyinv = X' * Vyinv
# 	logdetVy = logdet(Vy)
#        # beta hat
# 	betahat = inv(XtVyinv * X) * XtVyinv * Y
#        # sigma2 hat
# 	fittedValues =  X * betahat
# 	residuals = Y - fittedValues
# 	sigma2hat = 1/ntaxa * (residuals' * Vyinv * residuals)
#        # log likelihood
# 	loglik = - 1 / 2 * (ntaxa + ntaxa * log(2 * pi) + ntaxa * log(sigma2hat) + logdetVy)
# 	# Result
# #	res = phyloNetworkRegression(betahat, sigma2hat[1], loglik[1], V, Vy, fittedValues, residuals)
# 	return((betahat, sigma2hat[1], loglik[1], V, Vy, logdetVy, fittedValues, residuals))
# end
