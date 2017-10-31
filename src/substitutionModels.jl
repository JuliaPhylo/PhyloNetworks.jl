abstract type TraitSubstitutionModel end
const SM = TraitSubstitutionModel
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
    label0::String
    label1::String
    function BinaryTraitSubstitutionModel(α::Float64, β::Float64, label0::String, label1::String)
	    α >= 0. || error("parameter α must be non-negative")
	    β >= 0. || error("parameter β must be non-negative")
	    ab = α+β
        ab > 0. || error("α+β must be positive")
        new(α, β, β/ab, α/ab, label0, label1)
    end
end

BinaryTraitSubstitutionModel(α, β) = BinaryTraitSubstitutionModel(α, β, "0", "1")

#struct BinaryTraitSubstitutionModel <: TraitSubstitutionModel
#    α::Float64
#    β::Float64
#    π0::Float64
#    π1::Float64
#    label::SVector{2, String}
#    function BinaryTraitSubstitutionModel(α::Float64, β::Float64)
#	    α >= 0. || error("parameter α must be non-negative")
#	    β >= 0. || error("parameter β must be non-negative")
#	    ab = α+β
#	    ab > 0. || error("α+β must be positive")
#        new(α, β, β/ab, α/ab, SVector("0", "1"))
#    end
#    # Fixit function with Strings[ label0, label1] in SVector
#end

const BTSM = BinaryTraitSubstitutionModel

function nStates(mod::BTSM)
    return 2::Int
end

function show(io::IO, object::BinaryTraitSubstitutionModel)
    str = "Binary Trait Substitution Model:\n"
    str *= "rate $(object.label0)→$(object.label1) α=$(object.α)\n"
    str *= "rate $(object.label1)→$(object.label0) β=$(object.β)\n"
    print(io, str)
end

#function show(io::IO, object::BinaryTraitSubstitutionModel)
#    str = "Binary Trait Substitution Model:\n"
#    str *= "rate 0→1 α=$(object.α)\n"
#    str *= "rate 1→0 β=$(object.β)\n"
#    print(io, str)
#end

"""
`EqualRatesSubstitutionModel` is an abstract type that contains all models
describing a substitution process impacting biological characters with equal rates
of transition between all character states with continous time Markov models.
"""

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
    Q(mod)

Generate a Q matrix for a `TraitSubstitutionModel`, of the form:
    ```math
    Q = \begin{bmatrix}
        Q_{0, 0} & Q_{0, 1} \\
        Q_{1, 0} & Q_{1, 1} \end{bmatrix}
    ```
"""

@inline function Q(mod::BTSM)
    return Bmatrix(-mod.α, mod.β, mod.α, -mod.β)
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
    function P(mod, t)

Generate a P matrix for a `TraitSubstitutionModel`, of the form:
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
    try
        eig_vals, eig_vecs = eig(Q(mod)) # Only hermitian matrices are diagonalizable by 
        # *StaticArrays*. Non-Hermitian matrices should be converted to `Array`first.
        return [eig_vecs * expm(diagm(eig_vals)*i) * eig_vecs' for i in t]
    catch
        eig_vals, eig_vecs = eig(Array(Q(mod)))
        k = nStates(mod)
        return [SMatrix{k,k}(eig_vecs * expm(diagm(eig_vals)*i) * inv(eig_vecs)) for i in t]
    end
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
    return Bmatrix(p0+a1, p0-a0, p1-a1, p1+a0) # By columns
end

function randomTrait!(endTrait::AbstractVector{Int}, mod::SM, t::Float64, start::AbstractVector{Int})
    Pt = P(mod, t)
    k = size(Pt, 1) # number of states
    w = [aweights(Pt[i,:]) for i in 1:k]
    for i in 1:length(start)
        endTrait[i] =sample(1:k, w[start[i]])
    end
    return endTrait
end

"""
    randomTrait(model, t, start) 

Simulate traits along one edge of length t.
`start` must be a vector of integers, each representing the starting value of one trait.

# Examples
```julia-repl
julia> m1 = PhyloNetworks.BinaryTraitSubstitutionModel(1.0, 2.0) 
julia> srand(12345);
julia> randomTrait(m1, 0.2, [1,2,1,2,2])
 5-element Array{Int64,1}:
 1
 2
 1
 2
 2
```
"""

function randomTrait(mod::SM, t::Float64, start::AbstractVector{Int})
    res = Vector{Int}(length(start))
    randomTrait!(res, mod, t, start)
end

"""
    randomTrait(mod, net; ntraits=1, keepInternal=true, checkPreorder=true)

Simulates evolution of discrete traits on a rooted evolutionary network based on 
the supplied evolutionary model and returns an array of character states for each character
at each node. Trait sampling is uniform at the root.

# Arguments
- ntraits: the number of traits to be evaluated. Default value of 1.
- keepInternal: if true, export character states at internal nodes

# Examples
```julia-repl
julia> m1 = PhyloNetworks.BinaryTraitSubstitutionModel(1.0, 2.0) 
julia> net = readTopology("(A:1.0,(B:1.0,(C:1.0,D:1.0):1.0):1.0);")
julia> srand(12345);
julia> a,b = randomTrait(m1, net)
 ([1 2 … 1 2], String["-2", "-3", "-4", "D", "C", "B", "A"])
julia> a
 1×7 Array{Int64,2}:
 1  2  1  1  1  1  2
julia> b
 7-element Array{String,1}:
 "-2"
 "-3"
 "-4"
 "D"
 "C"
 "B"
 "A"
```
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
        nodeLabels = [n.name for n in net.nodes_changed if n.leaf]
    else
        nodeLabels = [n.name == "" ? string(n.number) : n.name for n in net.nodes_changed]    
    end
    return M, nodeLabels
end

function getTipSubmatrix(M::Matrix, net::HybridNetwork)
    nodenames = [n.name for n in net.nodes_changed]
    tipind = Int[]
    for l in net.leaf
        push!(tipind, findfirst(nodenames, l.name))
    end
    return M[:, tipind]
end

function randomTrait!(M::Matrix{Int}, mod::SM, net::HybridNetwork)
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
    randomTrait!(view(V, :, i), mod, edge.length, view(V, :, parentIndex))
end

function updateHybridRandomTrait!(V::Matrix,
        i::Int, parentIndex1::Int, parentIndex2::Int,
        edge1::PhyloNetworks.Edge, edge2::PhyloNetworks.Edge, mod)
    randomTrait!(view(V, :, i), mod, edge1.length, view(V, :, parentIndex1))
    tmp = randomTrait(mod, edge2.length, view(V, :, parentIndex2))
    for j in 1:size(V,1)
        if V[j,i] == tmp[j]
            continue
        end
        if rand() > edge1.gamma
            V[j,i] = tmp[j]
        end
    end
end