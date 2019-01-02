"""
    TraitSubstitutionModel

Abstract type for discrete trait substitution models,
using a continous time Markov model on a phylogeny.
Adapted from the substitutionModels module in BioJulia.
The same [`Q`](@ref) and [`P`](@ref) function names are used for the
transition rates and probabilities.

see [`BinaryTraitSubstitutionModel`](@ref),
[`EqualRatesSubstitutionModel`](@ref),
[`TwoBinaryTraitSubstitutionModel`](@ref)
"""
abstract type TraitSubstitutionModel{T} end
const SM = TraitSubstitutionModel{T} where T
const Bmatrix = SMatrix{2, 2, Float64}

"""
    nStates(model)

Number of character states for a given trait evolution model.
"""
nStates(mod::SM) = error("nStates not defined for $(typeof(mod)).")

"""
    nparams(model)

Number of parameters for a given trait evolution model
(length of field `model.rate`).
"""
nparams(mod::SM) = error("nparams not defined for $(typeof(mod)).")

"""
    Q(model)

Substitution rate matrix for a given substitution model:
Q[i,j] is the rate of transitioning from state i to state j.
"""
Q(mod::SM) = error("rate matrix Q not defined for $(typeof(mod)).")

"""
    showQ(IO, model)

Print the Q matrix to the screen, with trait states as labels on rows and columns.
adapted from prettyprint function by mcreel, found 2017/10 at
https://discourse.julialang.org/t/display-of-arrays-with-row-and-column-names/1961/6
"""
function showQ(io::IO, object::SM)
    M = Q(object)
    pad = max(8,maximum(length.(object.label))+1)
    for i = 1:size(M,2) # print the header
        print(io, lpad(object.label[i],(i==1 ? 2*pad : pad), " "))
    end
    print(io, "\n")
    for i = 1:size(M,1) # print one row per state
        if object.label != ""
            print(io, lpad(object.label[i],pad," "))
        end
        for j = 1:size(M,2)
            if j == i
                print(io, lpad("*",pad," "))
            else
                fmt = "%$(pad).4f"
                @eval(@printf($io,$fmt,$(M[i,j])))
            end
        end
        print(io, "\n")
    end
end

"""
    P(mod, t)

Probability transition matrix for a [`TraitSubstitutionModel`](@ref), of the form

    P[1,1] ... P[1,k]
       .          .
       .          .
    P[k,1] ... P[k,k]

where P[i,j] is the probability of ending in state j after time t,
given that the process started in state i.
"""
@inline function P(mod::SM, t::Float64)
    t >= 0.0 || error("substitution model: >=0 branch lengths are needed")
    return exp(Q(mod) * t)
end

"""
    P(mod, t::Array{Float64})

When applied to a general substitution model, matrix exponentiation is used.
The time argument `t` can be an array.
"""
function P(mod::SM, t::Array{Float64})
    all(t .>= 0.0) || error("t's must all be positive")
    try
        eig_vals, eig_vecs = eigen(Q(mod)) # Only hermitian matrices are diagonalizable by
        # *StaticArrays*. Non-Hermitian matrices should be converted to `Array`first.
        return [eig_vecs * Matrix(Diagonal(exp.(eig_vals*i))) * eig_vecs' for i in t]
    catch
        eig_vals, eig_vecs = eigen(Array(Q(mod)))
        k = nStates(mod)
        return [SMatrix{k,k}(eig_vecs * Matrix(Diagonal(exp.(eig_vals*i))) * inv(eig_vecs)) for i in t]
    end
end

"""
    BinaryTraitSubstitutionModel(α, β [, label])

[`TraitSubstitutionModel`](@ref) for binary traits (with 2 states).
Default labels are "0" and "1".
α is the rate of transition from "0" to "1", and β from "1" to "0".
"""
mutable struct BinaryTraitSubstitutionModel{T} <: TraitSubstitutionModel{T}
    rate::Vector{Float64}
    label::Vector{T} # most often: T = String, but could be BioSymbols.DNA
    function BinaryTraitSubstitutionModel{T}(rate, label::Vector{T}) where T
        @assert length(rate) == 2 "binary state: need 2 rates"
        rate[1] >= 0. || error("parameter α must be non-negative")
        rate[2] >= 0. || error("parameter β must be non-negative")
        ab = rate[1] + rate[1]
        ab > 0. || error("α+β must be positive")
        @assert length(label) == 2 "need 2 labels exactly"
        new(rate, label)
    end
end
const BTSM = BinaryTraitSubstitutionModel{T} where T
BinaryTraitSubstitutionModel(r::AbstractVector, label::AbstractVector) = BinaryTraitSubstitutionModel{eltype(label)}(r, label)
BinaryTraitSubstitutionModel(α::Float64, β::Float64, label) = BinaryTraitSubstitutionModel([α,β], label)
BinaryTraitSubstitutionModel(α::Float64, β::Float64) = BinaryTraitSubstitutionModel(α, β, ["0", "1"])

"""
# Examples

```jldoctest
julia> m1 = BinaryTraitSubstitutionModel([1.0,2.0], ["low","high"])
Binary Trait Substitution Model:
rate low→high α=1.0
rate high→low β=2.0


julia> nStates(m1)
2
```
"""
function nStates(::BTSM)
    return 2::Int
end

nparams(::BTSM) = 2::Int

"""
For a BinaryTraitSubstitutionModel, the rate matrix Q is of the form:

    -α  α
     β -β
"""
@inline function Q(mod::BTSM)
    return Bmatrix(-mod.rate[1], mod.rate[2], mod.rate[1], -mod.rate[2])
end

function Base.show(io::IO, object::BTSM)
    str = "Binary Trait Substitution Model:\n"
    str *= "rate $(object.label[1])→$(object.label[2]) α=$(round(object.rate[1], digits=5))\n"
    str *= "rate $(object.label[2])→$(object.label[1]) β=$(round(object.rate[2], digits=5))\n"
    print(io, str)
end

@inline function P(mod::BTSM, t::Float64)
    t >= 0.0 || error("substitution model: >=0 branch lengths are needed")
    ab = mod.rate[1] + mod.rate[2]
    e1 = exp(-ab*t)
    p0 = mod.rate[2]/ab # asymptotic frequency of state "0"
    p1 = mod.rate[1]/ab # asymptotic frequency of state "1"
    a0= p0 *e1
    a1= p1*e1
    return Bmatrix(p0+a1, p0-a0, p1-a1, p1+a0) # by columns
end

"""
    TwoBinaryTraitSubstitutionModel(rate [, label])

[`TraitSubstitutionModel`](@ref) for two binary traits, possibly correlated.
Default labels are "x0", "x1" for trait 1, and "y0", "y1" for trait 2.
If provided, `label` should be a vector of size 4, listing labels for
trait 1 first then labels for trait 2.
`rate` should be a vector of substitution rates of size 8.
rate[1],...,rate[4] describe rates of changes in trait 1.
rate[5],...,rate[8] describe rates of changes in trait 2.

In the transition matrix, trait combinations are listed in the following order:
x0-y0, x0-y1, x1-y0, x1-y1.

# example

```julia
model = TwoBinaryTraitSubstitutionModel([2.0,1.2,1.1,2.2,1.0,3.1,2.0,1.1],
        ["carnivory", "noncarnivory", "wet", "dry"]);
model
using PhyloPlots
plot(model) # to visualize states and rates
```
"""
mutable struct TwoBinaryTraitSubstitutionModel <: TraitSubstitutionModel{String}
    rate::Vector{Float64}
    label::Vector{String}
    function TwoBinaryTraitSubstitutionModel(α, label)
        all( x -> x >= 0., α) || error("rates must be non-negative")
        @assert length(α)==8 "need 8 rates"
        @assert length(label)==4 "need 4 labels for all combinations of 2 binary traits"
        new(α, [string(label[1], "-", label[3]), # warning: original type of 'label' lost here
                string(label[1], "-", label[4]),
                string(label[2], "-", label[3]),
                string(label[2], "-", label[4])])
    end
end
const TBTSM = TwoBinaryTraitSubstitutionModel
TwoBinaryTraitSubstitutionModel(α::Vector{Float64}) = TwoBinaryTraitSubstitutionModel(α, ["x0", "x1", "y0", "y1"])

nStates(::TBTSM) = 4::Int
nparams(::TBTSM) = 8::Int

function Q(mod::TBTSM)
    M = fill(0.0,(4,4))
    a = mod.rate
    M[1,3] = a[1]
    M[3,1] = a[2]
    M[2,4] = a[3]
    M[4,2] = a[4]
    M[1,2] = a[5]
    M[2,1] = a[6]
    M[3,4] = a[7]
    M[4,3] = a[8]
    M[1,1] = -M[1,2] - M[1,3]
    M[2,2] = -M[2,1] - M[2,4]
    M[3,3] = -M[3,4] - M[3,1]
    M[4,4] = -M[4,3] - M[4,2]
    return M
end

function Base.show(io::IO, object::TBTSM)
    print(io, "Substitution model for 2 binary traits, with rate matrix:\n")
    showQ(io, object)
end

"""
    EqualRatesSubstitutionModel(numberStates, α, labels)

[`TraitSubstitutionModel`](@ref) for traits with any number of states
and equal substitution rates α between all states.
Default labels are "1","2",...
"""
mutable struct EqualRatesSubstitutionModel{T} <: TraitSubstitutionModel{T}
    k::Int
    rate::Vector{Float64}
    label::Vector{T}
    function EqualRatesSubstitutionModel{T}(k, rate, label::Vector{T}) where T
        k >= 2 || error("parameter k must be greater than or equal to 2")
        @assert length(rate)==1 "rate must be a vector of length 1"
        rate[1] > 0 || error("parameter α (rate) must be positive")
        @assert length(label)==k "label vector of incorrect length"
        new(k, rate, label)
    end
end
const ERSM = EqualRatesSubstitutionModel{T} where T
EqualRatesSubstitutionModel(k::Int, α::Float64, label::AbstractVector) = EqualRatesSubstitutionModel{eltype(label)}(k,[α],label)
EqualRatesSubstitutionModel(k::Int, α::Float64) = EqualRatesSubstitutionModel{String}(k, [α], string.(1:k))

function nStates(mod::ERSM)
    return mod.k
end
nparams(::ERSM) = 1::Int

function Base.show(io::IO, object::ERSM)
    str = "Equal Rates Substitution Model with k=$(object.k),\n"
    str *= "all rates equal to α=$(object.rate[1]).\n"
    str *= "rate matrix Q:\n"
    print(io, str)
    showQ(io, object)
end

function Q(mod::ERSM)
    α = mod.rate[1]
    M = fill(α, (mod.k,mod.k))
    d = -(mod.k-1) * α
    for i in 1:mod.k
        M[i,i] = d
    end
    return M
end

"""
    randomTrait(model, t, start)
    randomTrait!(end, model, t, start)

Simulate traits along one edge of length t.
`start` must be a vector of integers, each representing the starting value of one trait.
The bang version (ending with !) uses the vector `end` to store the simulated values.

# Examples
```jldoctest
julia> m1 = BinaryTraitSubstitutionModel(1.0, 2.0)
Binary Trait Substitution Model:
rate 0→1 α=1.0
rate 1→0 β=2.0


julia> using Random; Random.seed!(12345);

julia> randomTrait(m1, 0.2, [1,2,1,2,2])
5-element Array{Int64,1}:
 1
 2
 1
 1
 2
```
"""
function randomTrait(mod::SM, t::Float64, start::AbstractVector{Int})
    res = Vector{Int}(undef, length(start))
    randomTrait!(res, mod, t, start)
end

@doc (@doc randomTrait) randomTrait!
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
    randomTrait(model, net; ntraits=1, keepInternal=true, checkPreorder=true)

Simulate evolution of discrete traits on a rooted evolutionary network based on
the supplied evolutionary model. Trait sampling is uniform at the root.

optional arguments:

- `ntraits`: number of traits to be simulated (default: 1 trait).
- `keepInternal`: if true, export character states at all nodes, including
  internal nodes. if false, export character states at tips only.

output:

- matrix of character states with one row per trait, one column per node;
  these states are *indices* in `model.label`, not the trait labels themselves.
- vector of node labels (for tips) or node numbers (for internal nodes)
  in the same order as columns in the character state matrix

# examples

```jldoctest
julia> m1 = BinaryTraitSubstitutionModel(1.0, 2.0, ["low","high"]);

julia> net = readTopology("(((A:4.0,(B:1.0)#H1:1.1::0.9):0.5,(C:0.6,#H1:1.0::0.1):1.0):3.0,D:5.0);");

julia> using Random; Random.seed!(1234);

julia> trait, lab = randomTrait(m1, net)
([1 2 … 1 1], ["-2", "D", "-3", "-6", "C", "-4", "#H1", "B", "A"])

julia> trait
1×9 Array{Int64,2}:
 1  2  1  1  2  2  1  1  1

julia> lab
9-element Array{String,1}:
 "-2" 
 "D"  
 "-3" 
 "-6" 
 "C"  
 "-4" 
 "#H1"
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
    M = Matrix{Int}(undef, ntraits, nnodes) # M[i,j]= trait i for node j
    randomTrait!(M,mod,net)
    if !keepInternal
        M = getTipSubmatrix(M, net, indexation=:cols) # subset columns only. rows=traits
        nodeLabels = [n.name for n in net.nodes_changed if n.leaf]
    else
        nodeLabels = [n.name == "" ? string(n.number) : n.name for n in net.nodes_changed]    
    end
    return M, nodeLabels
end

function randomTrait!(M::Matrix, mod::SM, net::HybridNetwork)
    recursionPreOrder!(net.nodes_changed, M, # updates M in place
            updateRootRandomTrait!,
            updateTreeRandomTrait!,
            updateHybridRandomTrait!,
            mod)
end

function updateRootRandomTrait!(V::AbstractArray, i::Int, mod)
    sample!(1:nStates(mod), view(V, :, i)) # uniform at the root
    return
end

function updateTreeRandomTrait!(V::Matrix,
    i::Int,parentIndex::Int,edge::Edge,
    mod)
    randomTrait!(view(V, :, i), mod, edge.length, view(V, :, parentIndex))
end

function updateHybridRandomTrait!(V::Matrix,
        i::Int, parentIndex1::Int, parentIndex2::Int,
        edge1::Edge, edge2::Edge, mod)
    randomTrait!(view(V, :, i), mod, edge1.length, view(V, :, parentIndex1))
    tmp = randomTrait(mod, edge2.length, view(V, :, parentIndex2))
    for j in 1:size(V,1) # loop over traits
        if V[j,i] == tmp[j] # both parents of the hybrid node have the same trait
            continue # skip the rest: go to next trait
        end
        if rand() > edge1.gamma
            V[j,i] = tmp[j] # switch to inherit trait of parent 2
        end
    end
end
