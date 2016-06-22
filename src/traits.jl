# functions for trait evolution on network
# Claudia & Paul Bastide: November 2015

#################################################
## Function to traverse the network in the pre-order, updating a matrix
#################################################

# Matrix with rows and/or columns in topological order of the net.
type matrixTopologicalOrder
    V::Matrix # Matrix in itself
    nodesNumbersTopOrder::Vector{Int64} # Vector of nodes numbers for ordering of the matrix
    internalNodesNumbers::Vector{Int64} # Internal nodes numbers (original net order)
    tipsNumbers::Vector{Int64} # Tips numbers (original net order)
    tipsNames::Vector # Tips Names (original net order)
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
    for(i in 1:n) #sorted list of nodes
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
function Base.getindex(obj::matrixTopologicalOrder, d::Symbol)
    if d == :Tips # Extract rows and/or columns corresponding to the tips
        mask = indexin(obj.tipsNumbers, obj.nodesNumbersTopOrder)
        obj.indexation == "b" && return obj.V[mask, mask] # both columns and rows are indexed by nodes
        obj.indexation == "c" && return obj.V[:, mask] # Only the columns
        obj.indexation == "r" && return obj.V[mask, :] # Only the rows
    end
    if d == :InternalNodes # Idem, for internal nodes
        mask = indexin(obj.internalNodesNumbers, obj.nodesNumbersTopOrder)
        obj.indexation == "b" && return obj.V[mask, mask]
        obj.indexation == "c" && return obj.V[:, mask] 
        obj.indexation == "r" && return obj.V[mask, :] 
    end
    if d == :TipsNodes
        maskNodes = indexin(obj.internalNodesNumbers, obj.nodesNumbersTopOrder)
        maskTips = indexin(obj.tipsNumbers, obj.nodesNumbersTopOrder,)
        obj.indexation == "b" && return obj.V[maskTips, maskNodes]
        obj.indexation == "c" && error("Both rows and columns must be net
        ordered to take the submatrix tips vs internal nodes.")
        obj.indexation == "r" && error("Both rows and columns must be net
        ordered to take the submatrix tips vs internal nodes.")
    end
    d == :All && return obj.V
end

#################################################
## Functions to compute the variance-covariance between Node and its parents
#################################################
"""
`sharedPathMatrix(net::HybridNetwork; checkPreorder=true::Bool)`

This function computes the shared path matrix between all the nodes of a
network. It assumes that the network is in the pre-order. If checkPreorder is
true, then it runs function 'preoder' on the network beforehand. Returns an
object of type 'matrixTopologicalOrder'.
"""
function sharedPathMatrix(
	net::HybridNetwork;
	checkPreorder=true::Bool
	)
	recursionPreOrder(net, checkPreorder, initsharedPathMatrix, updateRootSharedPathMatrix!, updateTreeSharedPathMatrix!, updateHybridSharedPathMatrix!, "b")
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
	for(j in 1:(i-1))
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
        for(j in 1:(i-1))
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
#    for(i in 1:n) #sorted list of nodes
#        updateSharedPathMatrix!(i,net.nodes_changed,V)
#    end
#    return V
#end


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


#################################################
## Functions for Phylgenetic Network regression
#################################################

# New type for phyloNetwork regression
type phyloNetworkLinearModel <: LinPredModel 
    lm::GLM.LinearModel # result of a lm on a matrix
    V::matrixTopologicalOrder
    Vy::Matrix
    RL::LowerTriangular
    Y::Vector
    X::Matrix
    logdetVy::Real
    ind::Vector{Int} # vector matching the tips of the network against the names of the data frame provided. 0 if the match could not be preformed.
    msng::BitArray{1} # Which tips are not missing
end

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
   	phyloNetworkLinearModel(lm(RL\X, RL\Y), V, Vy, RL, Y, X, logdet(Vy), ind, msng)
end


"""
`phyloNetworklm(f::Formula, fr::AbstractDataFrame, net::HybridNetwork)`

Performs a regression according to the formula provided by the user, using
the correlation structure induced by the network.
The data frame fr should have an extra column labelled "tipsNames" that gives
the names of the taxa for each observation.
"""
# Deal with formulas
function phyloNetworklm(
	f::Formula,
	fr::AbstractDataFrame,
	net::HybridNetwork,
	model="BM"::AbstractString
	)
    # Match the tips names: make sure that the data provided by the user will
    # be in the same order as the ordered tips in matrix V.
    V = sharedPathMatrix(net)
    if any(V.tipsNames == "")
        warn("The network provided has no tip names. The tips are assumed te be is the same order than the data. You'd better know what you're doing.")
        ind = [0]
    elseif !any(DataFrames.names(fr) .== :tipsNames)
        warn("The entry data frame has no column labelled tipsNames. Please add such a column to match the tips against the network. Otherwise the tips are assumed to be is the same order than the data and you'd better know what you're doing.")
        ind = [0]
    else
#        ind = indexin(V.tipsNames, fr[:tipsNames])
        ind = indexin(fr[:tipsNames], V.tipsNames)
        if any(ind == 0) || length(unique(ind)) != length(ind)
            error("Tips names of the network and names provided in column tipsNames of the dataframe do not match.")
        end
#	fr = fr[ind, :]
    end
    # Find the regression matrix and answer vector
    mf = ModelFrame(f,fr)
    mm = ModelMatrix(mf)
    Y = convert(Vector{Float64},DataFrames.model_response(mf))
    # Fit the model
    DataFrames.DataFrameRegressionModel(phyloNetworklm(mm.m, Y, V, mf.msng, model, ind), mf, mm)
		# (Method copied from DataFrame/src/statsmodels/statsmodels.jl, lines 47-58)
    # Create the object
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
StatsBase.confint(m::phyloNetworkLinearModel) = confint(m.lm)
# coef table (coef, stderr, confint)
StatsBase.coeftable(m::phyloNetworkLinearModel) = coeftable(m.lm)
# Degrees of freedom for residuals
StatsBase.df_residual(m::phyloNetworkLinearModel) =  nobs(m) - length(coef(m))
# Degrees of freedom consumed in the model (+1: dispersion parameter)
StatsBase.df(m::phyloNetworkLinearModel) = length(coef(m)) + 1 
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
StatsBase.R2(m::phyloNetworkLinearModel) = 1 - deviance(m)/nulldeviance(m)
# adjusted coefficient of determination 
# Copied from GLM.jl/src/lm.jl, lines 141-146
function StatsBase.adjR2(obj::phyloNetworkLinearModel)
    n = nobs(obj)
    # df() includes the dispersion parameter
    p = df(obj) - 1
    1 - (1 - R2(obj))*(n-1)/(n-p)
end

## REMARK
# As phyloNetworkLinearModel <: LinPredModel, the following functions are automatically defined:
# AIC, AICc, BIC

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
function mu_estim(m::DataFrames.DataFrameRegressionModel{PhyloNetworks.phyloNetworkLinearModel,Float64})
	if (!m.mf.terms.intercept) error("The fit was not without intercept, so I cannot estimate mu") end
	return coef(m)[1]
end

### Functions specific to DataFrameRegressionModel
DataFrames.ModelFrame(m::DataFrames.DataFrameRegressionModel) = m.mf
DataFrames.ModelMatrix(m::DataFrames.DataFrameRegressionModel) = m.mm
DataFrames.Formula(m::DataFrames.DataFrameRegressionModel) = Formula(m.mf.terms)

### Print the results
# Variance
function paramstable(m::phyloNetworkLinearModel)
    Sig = sigma2_estim(m)
    "Sigma2: $(Sig)"
end
function Base.show(io::IO, obj::phyloNetworkLinearModel)
    println(io, "$(typeof(obj)):\n\nParameter(s) Estimates:\n", paramstable(obj), "\n\nCoefficients:\n", coeftable(obj))
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
#     quantile(TDist(df_residual(obj)), (1. - level)/2.) * [1. -1.]
# end
# Log likelihood of the fitted BM
# StatsBase.loglikelihood(m::phyloNetworkLinearModel) = - 1 / 2 * (nobs(m) + nobs(m) * log(2 * pi) + nobs(m) * log(sigma2_estim(m)) + m.logdetVy)
#StatsBase.loglikelihood(m::phyloNetworkLinPredModel) = - 1 / 2 * (nobs(m) + nobs(m) * log(2 * pi) + nobs(m) * log(sigma2_estim(m)) + m.logdetVy)
# Coefficients
# function StatsBase.coeftable(mm::phyloNetworkLinearModel)
#     cc = coef(mm)
#     se = stderr(mm)
#     tt = cc ./ se
#     CoefTable(hcat(cc,se,tt,ccdf(FDist(1, df_residual(mm)), abs2(tt))),
#               ["Estimate","Std.Error","t value", "Pr(>|t|)"],
#               ["x$i" for i = 1:size(mm.lm.pp.X, 2)], 4)
# end

#################################################
## Ancestral State Reconstruction
#################################################
# Class for reconstructed states on a network
type reconstructedStates
    traits_nodes::Vector
    variances_nodes::Matrix
    internalNodesNumbers::Vector{Int}
end

function Base.show(io::IO, obj::reconstructedStates)
    println(io, "$(typeof(obj)):\n",
#    		"Conditional Law of Ancestral traits:\n",
	   	 CoefTable(hcat(obj.internalNodesNumbers, obj.traits_nodes, diag(obj.variances_nodes)),
	    			["Node index", "Cond. Exp.", "Cond. Var."],
				fill("", size(obj.internalNodesNumbers))))
end


"""
`ancestralStateReconstruction(net::HybridNetwork, Y::Vector, params::paramsBM)`

Computes the conditional expectations and variances of the ancestral (un-observed)
traits values at the internal nodes of the phylogenetic network (net), 
given the values of the traits at the tips of the network (Y) and some
known parameters of the process used for trait evolution (params, only BM with fixed root
works for now).
"""
# Reconstruction from known BM parameters
function ancestralStateReconstruction(net::HybridNetwork,
                                      Y::Vector,
                                      params::paramsBM)
    V = sharedPathMatrix(net)
    ancestralStateReconstruction(V, Y, params)
end

function ancestralStateReconstruction(V::matrixTopologicalOrder,
				      Y::Vector,
				      params::paramsBM)
    Vy = V[:Tips]
    Vz = V[:InternalNodes]
    Vyz = V[:TipsNodes]
    R = cholfact(Vy)
    RL = R[:L]
    ancestralStateReconstruction(V, Vy, Vz, Vyz, RL, Y, params.mu,
    params.sigma2)
end

# Reconstruction from all the needed quantities
function ancestralStateReconstruction(V::matrixTopologicalOrder,
                                      Vy::Matrix, Vz::Matrix,
                                      Vyz::Matrix,
                                      RL::LowerTriangular{Float64,Array{Float64,2}},
                                      Y::Vector, mu::Real, sigma2::Real)
    # Vectors of means
    m_y = ones(size(Vy)[1]) .* mu # !! correct only of no predictor. 
    m_z = ones(size(Vz)[1]) .* mu # !! works because BM no shift.
    temp = RL \ Vyz
    m_z_cond_y = m_z + temp' * (RL \ (Y - m_y))
    V_z_cond_y = sigma2 .* (Vz - temp' * temp)
    reconstructedStates(m_z_cond_y, V_z_cond_y, V.internalNodesNumbers)
end

# Empirical reconstruciton from a fitted object
function ancestralStateReconstruction(obj::phyloNetworkLinearModel, mu::Real)
  Vyz = obj.V[:TipsNodes]
	# If the tips were re-organized, do the same for Vyz
	if (obj.ind != [0])
		Vyz = Vyz[obj.ind, :]
	else
		warn("There were no indication for the position of the tips on the network. Assuming that they are given in the same order. Please check that this is what you intended.")
	end
	ancestralStateReconstruction(obj.V,
															 obj.Vy,
															 obj.V[:InternalNodes],
                               Vyz,
                               obj.RL,
                               obj.Y,
															 mu,
															 sigma2_estim(obj))
end
# Handling the ancestral mean
function ancestralStateReconstruction(obj::phyloNetworkLinearModel)
	mu = mu_estim(obj) # Produce a warning if intercept is unknown.
	ancestralStateReconstruction(obj, mu)
end
# For a DataFrameRegressionModel
function ancestralStateReconstruction(obj::DataFrames.DataFrameRegressionModel{PhyloNetworks.phyloNetworkLinearModel,Float64})
	mu = mu_estim(obj)
	ancestralStateReconstruction(obj.model, mu)
end

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




#################################################
## Simulation Function
#################################################

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


# Uses recursion on the network.
# Takes params of type paramsProcess as an entry
# Returns a matrix with two lines:
# - line one = expectations at all the nodes
# - line two = simulated values at all the nodes
# The nodes are ordered as given by topological sorting
function simulate(net::HybridNetwork, params::paramsProcess, model="BM"::AbstractString, checkPreorder=true::Bool)
	M = recursionPreOrder(net, checkPreorder, initSimulateBM, updateRootSimulateBM!, updateTreeSimulateBM!, updateHybridSimulateBM!, "c", params)
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
function updateTreeSimulateBM!(M::Matrix, i::Int, parentIndex::Int, edge::Edge, params::Tuple{paramsBM})
	params = params[1]
	M[1, i] = params.mu  # expectation
	M[2, i] = M[2, parentIndex] + sqrt(params.sigma2 * edge.length) * randn() # random value
end

# Going down to an hybrid node
function updateHybridSimulateBM!(M::Matrix, i::Int, parentIndex1::Int, parentIndex2::Int, edge1::Edge, edge2::Edge, params::Tuple{paramsBM})
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
    squeeze(obj.M[:d][2, :], 1)
end

# function extractSimulateTips(sim::Matrix, net::HybridNetwork)
# 	mask = getTipsIndexes(net)
# 	return(squeeze(sim[2, mask], 1))
# end

