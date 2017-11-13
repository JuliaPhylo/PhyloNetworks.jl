"""
    discrete_tree_corelikelihood(tree, tips, logtrans, forwardlik, directlik)

Calculate likelihood of discrete character states on a phylogenetic tree. Designed for 
use inside the `discrete_optimlikelihood` function.

"""
function discrete_tree_corelikelihood(tree::HybridNetwork, tips::Dict{String,Int64},
    logtrans::AbstractArray, forwardlik::AbstractArray, directlik::AbstractArray)
    k=size(logtrans)[1]
    @show k
    @show logtrans
    for ni in reverse(1:length(tree.nodes_changed)) # post-order
        n = tree.nodes_changed[ni]
        @show ni
        @show n
        if n.leaf
            # assumes that forwardlik was initialized at 0.
            # keep it at 0. = log(1.) if no data
            tiplabel = n.name
            if haskey(tips, tiplabel)
                for i in 1:k
                    forwardlik[i,ni] = -Inf64 # log(0) = -Inf if i != observed state
                end
                forwardlik[tips[n.name], ni] = 0.
            end
        else # foward likelihood = product of direct likelihood over all children edges
            for e in n.edge
                n == getParent(e) || continue # to next edge if n is not parent of e
                forwardlik[:,ni] += directlik[:,e.number]
            end
        end
        @show forwardlik[:,ni]
        if ni==1 # root is first index in nodes changed
            logprior = [-log(k) for i in 1:k] # uniform prior; could be changed later based on user input
            @show logprior[1]
            @show forwardlik[1,ni]
            loglik = logprior[1] + forwardlik[1,ni] # log of prob of data AND root in state 1
            @show loglik
            for i in 2:k
                @show logprior[i]
                @show forwardlik[i,ni]
                @show logprior[i] + forwardlik[i,ni]
                loglik = logsumexp(loglik, logprior[i] + forwardlik[i,ni])
                @show loglik
            end
            @show forwardlik
            return loglik
        end
        # if we keep going: n is not the root
        # calculate direct likelihood on the parent edge of n
        for e in n.edge
            if n == getChild(e)
                lt = view(logtrans, :,:,e.number)
                @show lt
                directlik[:,e.number] = lt[:,1] + forwardlik[1,ni]
                for i in 1:k # state at parent node
                    for j in 2:k # j = state at node n
                        tmp = lt[i,j] + forwardlik[j,ni]
                        @show tmp
                        directlik[i,e.number] = logsumexp(directlik[i,e.number],tmp)
                    end
                end
                @show e
                @show directlik[:,e.number]
                break # we visited the parent edge: break out of for loop
            end
        end
    end
end

"""
    function discrete_tree_core_ancestralstate(tree, tips, logtrans, forwardlik, 
        directlik, backwardlik)

Estimates ancestral states for discrete characters for each tree in a 
reticulated network. Designed to be used within the `discrete_optimliklihood function.`

# Examples

"""

function discrete_tree_core_ancestralstate(tree::HybridNetwork, tips::Dict{String,Int64},
    logtrans::AbstractArray, backwardlik::AbstractArray)
    k=size(logtrans)[1]
    logprior = [-log(k) for i in 1:k]
    for n in tree.nodes_changed
        if n.root
            backwardlik[:,n] = logprior
        else
            pn = n.isChild1 ? 1 : 2
            for e in n.edge
                pe = e.isChild1 ? 1 : 2
            end 
            for i in 1:k
                for j in 1:k
                    tmp = backwardlik[j,pn] +logtrans[j,i,pn] + # Fixit: sum child of pn: directlik[j,e]
                    if j==1
                        backwardlik[i,n] = tmp
                    elseif j > 1
                        backwardlik[i,n] = logsumexp(backwardlik[i,n],tmp)
                    end
                end
            end
        end
    end
    return backwardlik
end

"""
    discrete_corelikelihood(tips, mod, trees, ltw,
        logtrans,forwardlik,directlik,backwardlik)

Calculate likelihood for discrete characters on a network,
given the likelihood of each tree in the network. Designed for use inside
the `discrete_optimlikelihood` function.

"""

function discrete_corelikelihood(tips::Dict{String,Int64}, mod::TraitSubstitutionModel,
    trees::Array{HybridNetwork}, ltw::AbstractVector, logtrans::AbstractArray,
    forwardlik::AbstractArray, directlik::AbstractArray, backwardlik::AbstractArray)
    ll = Array{Float64,1}(length(trees))
    for t in 1:length(trees)
        ll[t] = discrete_tree_corelikelihood(trees[t],tips,logtrans,
                   view(forwardlik, :,:,t),view(directlik, :,:,t))
    end
    #f(t) = discrete_tree_corelikelihood(trees[t],tips,logtrans,view(forwardlik,:,:,t),
    #         view(directlik,:,:,t),view(backwardlik,:,:,t))
    #ll = pmap(f, 1:length(trees)) # ll = loglikelihood given each tree
    @show ll
    res = ll[1] + ltw[1] # result: loglikelihood given the network
    @show res
    for t in 2:length(trees)
        res = logsumexp(res, ll[t] + ltw[t])
    end
    return res
end

"""
    discrete_optimlikelihood(tips, mod, net, NLoptMethod=:LD_MMA,
        ftolRel=1e-12, ftolAbs=1e-8, xtolRel=1e-10, xtolAbs=1e-10)

Calculate optimized likelihood of a reticulate network given a
dictionary of discrete character states at the tips, 
a continous time Markov model, and a network.

Assumes discrete character states are consecutive integers starting from 1.

optional arguments (default):
- checkPreorder (true)
- NLoptMethod (:LD_MMA) for the optimization algorithm.
Other options include :LN_COBYLA (derivative-free); see NLopt package.
- tolerance values to control when the optimization is stopped:
  ftolRel (1e-12), ftolAbs (1e-10) on the criterion, and
  xtolRel (1e-10), xtolAbs (1e-10)

# Examples:

julia-repl> net = readTopology("(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);")
julia-repl> tips = Dict("A" => 1, "B" => 1, "C" => 2, "D" => 2)
julia-repl> PhyloNetworks.discrete_optimlikelihood(tips, m1, net)
res = -2.6638637960257583
"""

function discrete_optimlikelihood(tips::Dict{String,Int64}, mod::TraitSubstitutionModel, net::HybridNetwork;
    NLoptMethod=:LD_MMA::Symbol, ftolRel=1e-12::Float64, ftolAbs=1e-8::Float64,
    xtolRel=1e-10::Float64, xtolAbs=1e-10::Float64)
    # fixit new function to detect number of states, make model labels compatible with julia indexing
    trees = displayedTrees(net, 0.0)
    @show trees
    for tree in trees
        preorder!(tree)
        directEdges!(tree)
    end  
    #tips::Dict{Int64,Set{T}}
    ntrees = length(trees)
    @show ntrees
    k = nStates(mod)
    @show k
    #both mlik and logtrans should be 3-d arrays
    #mlik[i,n,t] = log P{data below n in tree t given state i above n}
    # fixit: re-number edges to be consective, positive numbers; check edges are positive numbers
    #keep track of largest edge number
    #initialize 3d array: Array{Float64}((i,j,e))
    # fixit: later change the arrays below into SharedArray
    forwardlik = zeros(Float64, k,length(net.node),ntrees)
    directlik  = zeros(Float64, k,length(net.edge),ntrees)
    backwardlik= zeros(Float64, k,length(net.node),ntrees)
    #logtrans[i,j,e]; i = start_state, j = end_state, e = edge.number
    #Step 1
    ltw = Array{Float64,1}(length(trees))
    t = 0
    for tree in trees
        t+=1
        ltw[t] = 0.0
        for e in tree.edge
            if e.gamma != 1.0
                ltw[t] += log(e.gamma)
            end
        end
    end
    @show ltw
    #Step 2
    logtrans = zeros(Float64, k,k,length(net.edge))
    for edge in net.edge
        @show edge
        @show logtrans[:,:,edge.number]
        logtrans[:,:,edge.number] = log.(P(mod,edge.length)) # element-wise
        @show logtrans[:,:,edge.number]
    end
    #Step 3
    #NLoptMethod=:LD_MMA::Symbol,
    #ftolRel=1e-12::Float64, ftolAbs=1e-8::Float64,
    #xtolRel=1e-10::Float64, xtolAbs=1e-10::Float64)
    
    #k=length(values[tips])
    #opt = NLopt.Opt(NLoptMethod,k) # :LD_MMA to use gradient
    #NLopt.ftol_rel!(opt,ftolRel) # relative criterion
    #NLopt.ftol_abs!(opt,ftolAbs) # absolute criterion
    #NLopt.xtol_rel!(opt,xtolRel)
    #NLopt.xtol_abs!(opt,xtolAbs)
    #NLopt.maxeval!(opt,1000) # max number of iterations
    #NLopt.lower_bounds!(opt, zeros(k))
    #NLopt.upper_bounds!(opt,upper(net))
    #count = 0

    #function obj(x::Vector{Int64},g::Vector{Float64})
    #    count += 1
    #    val = discrete_corelikelihood(tips,mod,trees,ltw,logtrans,forwardlik,directlik,backwardlik)
    #    return val
    #end
    #NLopt.min_objective!(opt,obj)
    #fmin, xmin, ret = NLopt.optimize(opt,tips) # optimization here!
    #return fmin,xmin,ret
    discrete_corelikelihood(tips,mod,trees,ltw,logtrans,forwardlik,directlik,backwardlik)
end

"""
consecutive_tips(tips)

Modify dictionary states into consecutive positive integers starting from 1.

# Examples

julia-repl> tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1)
julia-repl> consecutive_tips(tips)
Dict{Any,Int64} with 4 entries:
  "B" => 1
  "A" => 1
  "C" => 2
  "D" => 2
"""

function consecutive_tips(tips::Dict)
    newLabels = Dict{Any,Int64}()
    orderedTips = sort(collect(zip(values(tips),keys(tips))))
    t = 1 #consecutive integer
    for i in 1:length(keys(tips))
        if i==1
            newLabels[orderedTips[i][2]] = t
        else
            if orderedTips[i][1] == orderedTips[i-1][1]
                newLabels[orderedTips[i][2]] = t
            else
                t+=1
                newLabels[orderedTips[i][2]] = t
            end
        end
    end
    return newLabels
end

function readDataFrametoDict(filename::String)
    tmp = CSV.read(filename)
    DataFrameToDict(tmp)
end

"""
DataFrameToDict(tmp)

Convert DataFrame of character states for leaves of a given phylogenetic 
network into a dictionary of leaf names as keys and leaf states as consecutive integers
beginning from 1.

Assumes data frame for a single trait, in which case the taxon names
are to appear in column 1 or in a column named "taxon" or "species", and
trait values are to appear in column 2 or in a column named "trait".

# Examples
julia-repl> tmp = DataFrame(species = ["A","B","C","D"], trait = [1,1,2,2,])
4×2 DataFrames.DataFrame
│ Row │ species │ trait │
├─────┼─────────┼───────┤
│ 1   │ A       │ 1     │
│ 2   │ B       │ 1     │
│ 3   │ C       │ 2     │
│ 4   │ D       │ 2     │
julia-repl> DataFrameToDict(tmp)
Dict{Any,Int64} with 4 entries:
  "B" => 1
  "A" => 1
  "C" => 2
  "D" => 2

"""


function DataFrameToDict(tmp::DataFrame)
    tipLabels = Array{String}(0)
    tipStates = Array{Any}(0)
    tips = Dict{Any,Any}()
    for i in tmp[1]
        push!(tipLabels, i)
    end
    for i in tmp[2]
        push!(tipStates, i)
    end
    length(tipLabels) == length(tipStates) || 
        error("number of tips does not match number of states")
    for i in 1:length(tipLabels)
        tips[tipLabels[i]] = tipStates[i]
    end
    tips = PhyloNetworks.consecutive_tips(tips)
    return tips
end