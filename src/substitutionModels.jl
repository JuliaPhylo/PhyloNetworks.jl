using StaticArrays #move to PhyloNetworks.jl and add to REQUIRE

abstract type TraitSubstitutionModel end
const SM = TraitSubstitutionModel
const BTSM = BinaryTraitSubstitutionModel
const Bmatrix = SMatrix{2, 2, Float64}

"""
`BinaryTraitSubstitutionModel` is an abstract type that contains all models
describing a substitution process impacting biological characters with binary states
with continous time Markov models.
"""

struct BinaryTraitSubstitutionModel <: TraitSubstitutionModel
    α::Float64
    β::Float64
    π0::Float64
    π1::Float64
    label::SVector{2, String}
    function BinaryTraitSubstitutionModel(α::Float64, β::Float64)
	    α >= 0. || error("parameter α must be non-negative")
	    β >= 0. || error("parameter β must be non-negative")
	    ab = α+β
	    ab > 0. || error("α+β must be positive")
        new(α, β, β/ab, α/ab, SVector("0", "1"))
    end
    # Fixit function with Strings[ label0, label1] in SVector
end

function nStates(mod::BTSM)
    return 2::Int
end

function show(io::IO, object::BinaryTraitSubstitutionModel)
    str = "Binary Trait Substitution Model:\n"
    str *= "rate 0→ 1 α=$(object.α)\n"
    str *= "rate 1→ 0 β =$(object.β)\n"
    print(io, str)
end

struct EqualRatesSubstitutionModel <: TraitSubstitutionModel
    k::Int
    α::Float64
    function EqualRatesSubstitutionModel(k::Int, α::Float64)
        k >= 2 || error("parameter k must be greater than or equal to 2")
        α > 0 || error("parameter α must be positive")
        new(k, α)
    end
end

function show(io::IO, object::EqualRatesSubstitutionModel)
    str = "Equal Rates Substitution Model:\n"
    str *= "all rates α=$(object.α)\n"
    str *= "number of states, k=$(object.k)\n"
    print(io, str)
end

function nStates(mod::EqualRatesSubstitutionModel)
    return mod.k
end

"""
Generate a Q matrix for a `BinaryTraitSubstitutionModel`, of the form:
    ```math
    Q = \begin{bmatrix}
        Q_{0, 0} & Q_{0, 1} \\
        Q_{1, 0} & Q_{1, 1} \end{bmatrix}
    ```
"""

@inline function Q(mod::BTSM)
    return Bmatrix(-mod.α, mod.α,
                   mod.β, -mod.β)
end

function Q(mod::EqualRatesSubstitutionModel)
    M = fill(mod.α,(mod.k,mod.k))
    d = -(mod.k-1)*mod.α
    for i in 1:mod.k
        M[i,i] = d
    end
    return M
end

"""
Generate a P matrix for a `BinaryTraitSubstitutionModel`, of the form:
    ```math
    P = \begin{bmatrix}
        P_{0, 0} & P_{0, 1} \\
        P_{1, 0} & P_{1, 1} \end{bmatrix}.
    ```
    for specified time
"""

@inline function P(mod::SM, t::Float64)
    if t < 0.0
        error("t must be positive")
    end
    return expm(Q(mod) * t)
end

function P(mod::SM, t::Array{Float64})
    if any(t .< 0.0)
        error("t must be positive")
    end
    eig_vals, eig_vecs = eig(Q(mod))
    return [eig_vecs * expm(diagm(eig_vals)*i) * eig_vecs' for i in t]
end

@inline function P(mod::BTSM, t::Float64)
    if t < 0.0
        error("t must be positive")
    end
    e1 = exp(-(mod.α+mod.β)t)
    p0 = mod.π0
    p1 = mod.π1
    a0= p0 *e1
    a1= p1*e1
    return BMatrix (p0+a1, p1-a1, #fixit check order of matrix elements
                    p0-a0, p1+a0)
end

"""
    randomTrait(model, t, start) #fixit document both options
simulate traits along one edge of length t.
`start` must be a vector of integers, each representing the starting value of one trait.
"""

function randomTrait!(endTrait::Vector{Int}, mod::SM, t::Float64, start::Vector{Int})
    Pt = P(mod, t)
    w = [aweights(Pt[1,:]), aweights(Pt[2,:])]
    k = size(Pt, 1) # number of states
    for i in 1:length(start)
        endTrait[i] = sample(1:k, w[start[i]]) #fixit use sample! instead
    end
    return endTrait
end

function randomTrait(mod::SM, t::Float64, start::Vector{Int})
    res = Vector{Int}(length(start))
    randomTrait!(res, mod, t, start)
end

"""
#Fixit adject name and options
    discreteTraitDistanceMatrix(net::HybridNetwork;
        keepInternal=false::Bool, checkPreorder=true::Bool)

Compute the shared path matrix between all the nodes of a
network. Assumes that the network is in pre-order. If checkPreorder is
true (default), then it runs function `preorder` on the network beforehand.

# Examples
"""

function randomTrait(mod::SM, net::HybridNetwork;
    ntraits=1::Int, keepInternal=true::Bool, checkPreorder=true::Bool)
    net.isRooted || error("net needs to be rooted for preorder recursion")
    if(checkPreorder)
        preorder!(net)
    end
    nnodes = net.numNodes
    M = Matrix{Int}(ntraits, nnodes) # M[i,j]= trait i for node j
    randomTrait!(M,mod,net)
    if !keepInternal
        M = getTipSubmatrix(M, net)
    end
    return M
end

function getTipSubmatrix(M::Matrix, net::HybridNetwork)
    nodenames = [n.name for n in net.nodes_changed]
    tipind = Int[]
    for l in net.leaf
        push!(tipind, findfirst(nodenames, l.name))
    end
    return M[tipind, tipind]
end

function randomTrait!(M::Matrix{Float64}, mod::SM, net::HybridNetwork)
    recursionPreOrder!(net.nodes_changed, M, # updates M in place
            updateRootRandomTrait!,
            updateTreeRandomTrait!,
            updateHybridRandomTrait!,
            mod)
end

function updateRootRandomTrait!(V::AbstractArray, i::Int, mod)
    sample!(1:nStates(mod), view(V, :, i))
    return
end

function updateTreeRandomTrait!(V::Matrix,
    i::Int,parentIndex::Int,edge::PhyloNetworks.Edge,
    mod)
    for j in 1:(i-1)
        V[i,j] = mod*(V[i,parentIndex] + edge.length)
    end
end

function updateHybridRandomTrait!(V::Matrix,
        i::Int, parentIndex1::Int, parentIndex2::Int,
        edge1::PhyloNetworks.Edge, edge2::PhyloNetworks.Edge)
    γ1 = edge1.gamma
    γ2 = edge2.gamma
    r = :random in [0,1]
    if r <= γ1
        V[i,j] = V[i,parentIndex1]
    else
        V[i,j] = V[i,parentIndex2]
    end
    return
end
