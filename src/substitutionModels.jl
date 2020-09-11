"""
    SubstitutionModel

Abstract type for substitution models,
using a continous time Markov model on a phylogeny.
Adapted from [SubstitutionModels.jl](https://github.com/BioJulia/SubstitutionModels.jl/)
in BioJulia.

For variable rates, see [`RateVariationAcrossSites`](@ref)

For sub types, see [`NucleicAcidSubstitutionModel`](@ref), [`TraitSubstitutionModel`](@ref)

All these models are supposed to have fields `rate` and `eigeninfo`.
"""
abstract type SubstitutionModel end #ideally, we'd like this to be SubstitutionModels.SubstitionModel
const SM = SubstitutionModel
const Qmatrix = StaticArrays.SMatrix{4, 4, Float64}
const Pmatrix = StaticArrays.MMatrix{4, 4, Float64}
const Bmatrix = StaticArrays.MMatrix{2, 2, Float64}

"""
    TraitSubstitutionModel

For subtypes, see [`BinaryTraitSubstitutionModel`](@ref),
[`EqualRatesSubstitutionModel`](@ref),
[`TwoBinaryTraitSubstitutionModel`](@ref)
"""
abstract type TraitSubstitutionModel{T} <: SubstitutionModel end #this accepts labels
const TSM = TraitSubstitutionModel{T} where T #T is type of labels

"""
    NucleicAcidSubstitutionModel

Adapted from [SubstitutionModels.jl](https://github.com/BioJulia/SubstitutionModels.jl/)
in BioJulia. The same [`Q`](@ref) and [`P`](@ref) function names are used for the
transition rates and probabilities.

For subtypes, see [`JC69`](@ref), [`HKY85`](@ref)
"""
abstract type NucleicAcidSubstitutionModel <: SubstitutionModel end
const NASM = NucleicAcidSubstitutionModel

"""
    nparams(model)

Number of parameters for a given trait evolution model
(length of field `model.rate`).
"""
nparams(obj::SM) = error("nparams not defined for $(typeof(obj)).")

"""
    setrates!(model, rates)

update rates then call [`seteigeninfo!`](@ref) to update a model's eigeninfo
"""
function setrates!(obj::SM, rates::AbstractVector)
    obj.rate[:] = rates
    seteigeninfo!(obj)
end

"""
    seteigeninfo!(obj)

Calculate eigenvalue & eigenfector information for a substitution model (SM) object
(as needed to calculate transition rate matrices) and store this info within the object.
"""
function seteigeninfo!(obj::SM)
    error("seteigeninfo! not finalized yet for generic substitution models")
    # eig_vals, eig_vecs = eigen(Q(obj)) #this assumes that every SM has an eigeninfo and has
    # obj.eigeninfo[:] = eig_vals, eig_vecs, inv(eig_vecs)
    # then make P! use seteigeninfo! by for a generic SM
end

"""
    getlabels(model)

State labels of a substitution model.
"""
function getlabels(obj::SM)
    error("Model must be of type TraitSubstitutionModel or NucleicAcidSubstitutionModel. Got $(typeof(obj))")
end

"""
    nstates(model)

Number of character states for a given evolution model.
"""
nstates(obj::SM) = error("nstates not defined for $(typeof(obj)).")

"""
For example, this is 4 for a `NucleicAcidSubstitutionModel`.

```jldoctest
julia> nstates(JC69([0.03], false))
4

julia> nstates(HKY85([.5], [0.25, 0.25, 0.25, 0.25]))
4
```
"""
function nstates(obj::NASM)
    return 4::Int
end

function getlabels(obj::TSM)
    return obj.label
end

"""
for a given [`NucleicAcidSubstitutionModel`](@ref), labels are symbols
from [BioSymbols](https://github.com/BioJulia/BioSymbols.jl).
For now, only ACGTs are allowed. (When fitting data, any ambiguity code
in the data would be treated as missing value).

# examples

```jldoctest
julia> getlabels(JC69([0.03], false))
4-element Array{BioSymbols.DNA,1}:
 DNA_A
 DNA_C
 DNA_G
 DNA_T

julia> getlabels(HKY85([.5], repeat([0.25], 4)))
4-element Array{BioSymbols.DNA,1}:
 DNA_A
 DNA_C
 DNA_G
 DNA_T

```
"""
function getlabels(::NASM)
    return [BioSymbols.DNA_A, BioSymbols.DNA_C, BioSymbols.DNA_G, BioSymbols.DNA_T]
end

"""
    Q(model)

Substitution rate matrix for a given substitution model:
Q[i,j] is the rate of transitioning from state i to state j.
"""
Q(obj::SM) = error("rate matrix Q not defined for $(typeof(obj)).")

"""
    showQ(IO, model)

Print the Q matrix to the screen, with trait states as labels on rows and columns.
adapted from prettyprint function by mcreel, found 2017/10 at
https://discourse.julialang.org/t/display-of-arrays-with-row-and-column-names/1961/6
"""
function showQ(io::IO, obj::SM)
    M = Q(obj)
    pad = max(8,maximum(length(getlabels(obj))+1))
    for i = 1:size(M,2) # print the header
        print(io, lpad(getlabels(obj)[i],(i==1 ? 2*pad : pad), " "))
    end
    for i = 1:size(M,1) # print one row per state
        print(io, "\n")
        if getlabels(obj) != ""
            print(io, lpad(getlabels(obj)[i],pad," "))
        end
        for j = 1:size(M,2)
            if j == i
                print(io, lpad("*",pad," "))
            else
                fmt = "%$(pad).4f"
                @eval(@printf($io,$fmt,$(M[i,j])))
            end
        end
    end
end

"""
    P(obj, t)

Probability transition matrix for a [`TraitSubstitutionModel`](@ref), of the form

    P[1,1] ... P[1,k]
       .          .
       .          .
    P[k,1] ... P[k,k]

where P[i,j] is the probability of ending in state j after time t,
given that the process started in state i.
see also: [`P!`](@ref).

HKY example:
```jldoctest
julia> m1 = HKY85([0.5], [0.20, 0.30, 0.30, 0.20])
HKY85 Substitution Model base frequencies: [0.2, 0.3, 0.3, 0.2]
relative rate version with transition/tranversion ratio kappa = 0.5,
 scaled so that there is one substitution per unit time
rate matrix Q:
               A       C       G       T
       A       *  0.4839  0.2419  0.3226
       C  0.3226       *  0.4839  0.1613
       G  0.1613  0.4839       *  0.3226
       T  0.3226  0.2419  0.4839       *

julia> PhyloNetworks.P(m1, 0.2)
4×4 StaticArrays.MArray{Tuple{4,4},Float64,2,16} with indices SOneTo(4)×SOneTo(4):
 0.81592    0.0827167  0.0462192  0.0551445
 0.0551445  0.831326   0.0827167  0.0308128
 0.0308128  0.0827167  0.831326   0.0551445
 0.0551445  0.0462192  0.0827167  0.81592  
```

Juke-Cantor example:
```jldoctest
julia> m1 = JC69([1.]);

julia> PhyloNetworks.P(m1, 0.2)
4×4 StaticArrays.MArray{Tuple{4,4},Float64,2,16} with indices SOneTo(4)×SOneTo(4):
 0.824446   0.0585179  0.0585179  0.0585179
 0.0585179  0.824446   0.0585179  0.0585179
 0.0585179  0.0585179  0.824446   0.0585179
 0.0585179  0.0585179  0.0585179  0.824446 
```
"""
@inline function P(obj::SM, t::Float64)
    t >= 0.0 || error("substitution model: >=0 branch lengths are needed")
    k = nstates(obj)
    Pmat = MMatrix{k,k,Float64}(undef)
    return P!(Pmat, obj, t)
end

"""
    P!(Pmat::AbstractMatrix, obj::SM, t::Float64)

Fill in the input matrix `Pmat` with the transition rates
to go from each state to another in time `t`, according to rates in `Q`.
see also: [`P`](@ref).

```jldoctest
julia> m1 = BinaryTraitSubstitutionModel([1.0,2.0], ["low","high"])
Binary Trait Substitution Model:
rate low→high α=1.0
rate high→low β=2.0

julia> PhyloNetworks.P!(Matrix{Float64}(undef,2,2), m1, 0.3) # fills an uninitialized 2x2 matrix of floats
2×2 Array{Float64,2}:
 0.80219  0.19781
 0.39562  0.60438

julia> m2 = JC69([1.]);

julia> PhyloNetworks.P!(Matrix{Float64}(undef,4,4), m2, 0.2)
4×4 Array{Float64,2}:
 0.824446   0.0585179  0.0585179  0.0585179
 0.0585179  0.824446   0.0585179  0.0585179
 0.0585179  0.0585179  0.824446   0.0585179
 0.0585179  0.0585179  0.0585179  0.824446 
```
"""
@inline function P!(Pmat::AbstractMatrix, obj::SM, t::Float64)
    Pmat[:] = exp(Q(obj) * t)
    return Pmat
end

"""
    BinaryTraitSubstitutionModel(α, β [, label])

Model for binary traits, that is, with 2 states. Default labels are "0" and "1".
α is the rate of transition from "0" to "1", and β from "1" to "0".
"""
struct BinaryTraitSubstitutionModel{T} <: TraitSubstitutionModel{T} # fixit: back to mutable struct?
    rate::Vector{Float64}
    label::Vector{T} # most often: T = String, but could be BioSymbols.DNA
    eigeninfo::Vector{Float64}
    function BinaryTraitSubstitutionModel{T}(rate::Vector{Float64}, label::Vector{T}, eigeninfo::Vector{Float64}) where T
        #Warning: this constructor should not be used directly. Use it with the constructors below.
        @assert length(rate) == 2 "binary state: need 2 rates"
        rate[1] >= 0. || error("parameter α must be non-negative")
        rate[2] >= 0. || error("parameter β must be non-negative")
        @assert length(label) == 2 "need 2 labels exactly"
        new(rate, label, eigeninfo)
    end
end
const BTSM = BinaryTraitSubstitutionModel{T} where T
function BinaryTraitSubstitutionModel(r::AbstractVector, label::AbstractVector)
    obj = BinaryTraitSubstitutionModel{eltype(label)}(r::AbstractVector, label::AbstractVector, zeros(3)) # Vector{Float64}(undef,3) for julia v1.0
    seteigeninfo!(obj)
    return obj
end
BinaryTraitSubstitutionModel(α::Float64, β::Float64, label::AbstractVector) = BinaryTraitSubstitutionModel([α,β], label)
BinaryTraitSubstitutionModel(α::Float64, β::Float64) = BinaryTraitSubstitutionModel(α, β, ["0", "1"])

"""
for a [`BinaryTraitSubstitutionModel`]: store eigenvalue (q_01+q_10) and stationary distribution
"""
function seteigeninfo!(obj::BTSM)
    ab = obj.rate[1] + obj.rate[2] #eigenvalue = -(a+b)
    ab > 0. || error("α+β must be positive")
    p0 = obj.rate[2]/ab # asymptotic frequency of state "0"
    p1 = obj.rate[1]/ab # asymptotic frequency of state "1"
    obj.eigeninfo[1] = ab
    obj.eigeninfo[2] = p0
    obj.eigeninfo[3] = p1
end

"""
for a `BinaryTraitSubstitutionModel`, this is 2:

```jldoctest
julia> m1 = BinaryTraitSubstitutionModel([1.0,2.0], ["low","high"])
Binary Trait Substitution Model:
rate low→high α=1.0
rate high→low β=2.0

julia> nstates(m1)
2
```
"""
function nstates(::BTSM)
    return 2::Int
end

nparams(::BTSM) = 2::Int

"""
For a BinaryTraitSubstitutionModel, the rate matrix Q is of the form:

    -α  α
     β -β
"""
@inline function Q(obj::BTSM)
    return Bmatrix(-obj.rate[1], obj.rate[2], obj.rate[1], -obj.rate[2])
end

function Base.show(io::IO, obj::BTSM)
    str = "Binary Trait Substitution Model:\n"
    str *= "rate $(obj.label[1])→$(obj.label[2]) α=$(round(obj.rate[1], digits=5))\n"
    str *= "rate $(obj.label[2])→$(obj.label[1]) β=$(round(obj.rate[2], digits=5))"
    print(io, str)
end

function P!(Pmat::AbstractMatrix, obj::BTSM, t::Float64)
    e1 = exp(-obj.eigeninfo[1]*t)
    a0= obj.eigeninfo[2]*e1
    a1= obj.eigeninfo[3]*e1
    Pmat[1,1] = obj.eigeninfo[2]+a1
    Pmat[2,1] = obj.eigeninfo[2]-a0
    Pmat[1,2] = obj.eigeninfo[3]-a1
    Pmat[2,2] = obj.eigeninfo[3]+a0
    return Pmat
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
struct TwoBinaryTraitSubstitutionModel <: TraitSubstitutionModel{String}
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

nstates(::TBTSM) = 4::Int
nparams(::TBTSM) = 8::Int

function Q(obj::TBTSM)
    M = fill(0.0,(4,4))
    a = obj.rate
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

function Base.show(io::IO, obj::TBTSM)
    print(io, "Substitution model for 2 binary traits, with rate matrix:\n")
    showQ(io, obj)
end

"""
    EqualRatesSubstitutionModel(numberStates, α, labels)

[`TraitSubstitutionModel`](@ref) for traits with any number of states
and equal substitution rates α between all states.
Default labels are "1","2",...

# example

```jldoctest
julia> m1 = EqualRatesSubstitutionModel(2, [1.0], ["low","high"])
Equal Rates Substitution Model with k=2,
all rates equal to α=1.0.
rate matrix Q:
             low    high
     low       *  1.0000
    high  1.0000       *
```
"""
struct EqualRatesSubstitutionModel{T} <: TraitSubstitutionModel{T}
    k::Int
    rate::Vector{Float64}
    label::Vector{T}
    eigeninfo::Vector{Float64}
    function EqualRatesSubstitutionModel{T}(k::Int, rate::Vector{Float64},
        label::Vector{T}, eigeninfo::Vector{Float64}) where T
        k >= 2 || error("parameter k must be greater than or equal to 2")
        @assert length(rate)==1 "rate must be a vector of length 1"
        rate[1] > 0 || error("parameter α (rate) must be positive")
        @assert length(label)==k "label vector of incorrect length"
        new(k, rate, label, eigeninfo)
    end
end
const ERSM = EqualRatesSubstitutionModel{T} where T
function EqualRatesSubstitutionModel(k::Int, rate::Vector{Float64}, label::AbstractVector)
    obj = EqualRatesSubstitutionModel{eltype(label)}(k::Int, rate::Vector{Float64}, label::AbstractVector, zeros(1)) # Vector{Float64}(undef,3) for julia v1.0
    seteigeninfo!(obj)
    return obj
end
EqualRatesSubstitutionModel(k::Int, α::Float64, label::AbstractVector) = EqualRatesSubstitutionModel(k,[α],label)
EqualRatesSubstitutionModel(k::Int, α::Float64) = EqualRatesSubstitutionModel(k, [α], string.(1:k))


"""
for a [`EqualRatesSubstitutionModel`]: store lambda = k/(k-1), where k is the number of states
"""
function seteigeninfo!(obj::ERSM)
    q = (obj.k/(obj.k-1.0))
    obj.eigeninfo[1] = q
end

function nstates(obj::ERSM)
    return obj.k
end
nparams(::ERSM) = 1::Int

function Base.show(io::IO, obj::ERSM)
    str = "Equal Rates Substitution Model with k=$(obj.k),\n"
    str *= "all rates equal to α=$(round(obj.rate[1], digits=5)).\n"
    str *= "rate matrix Q:\n"
    print(io, str)
    showQ(io, obj)
end
function Q(obj::ERSM)
    #this might be wrong, doesnt match ERSM proof
    α = obj.rate[1]
    M = fill(α, (obj.k,obj.k))
    d = -(obj.k-1) * α
    for i in 1:obj.k
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
function randomTrait(obj::TSM, t::Float64, start::AbstractVector{Int})
    res = Vector{Int}(undef, length(start))
    randomTrait!(res, obj, t, start)
end

@doc (@doc randomTrait) randomTrait!
function randomTrait!(endTrait::AbstractVector{Int}, obj::TSM, t::Float64, start::AbstractVector{Int})
    Pt = P(obj, t)
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

julia> using Random; Random.seed!(47);

julia> trait, lab = randomTrait(m1, net)
([1 2 … 1 1], ["-2", "D", "-3", "-6", "C", "-4", "H1", "B", "A"])

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
 "H1"
 "B"  
 "A"  
```
"""
function randomTrait(obj::TSM, net::HybridNetwork;
    ntraits=1::Int, keepInternal=true::Bool, checkPreorder=true::Bool)
    net.isRooted || error("net needs to be rooted for preorder recursion")
    if(checkPreorder)
        preorder!(net)
    end
    nnodes = net.numNodes
    M = Matrix{Int}(undef, ntraits, nnodes) # M[i,j]= trait i for node j
    randomTrait!(M,obj,net)
    if !keepInternal
        M = getTipSubmatrix(M, net, indexation=:cols) # subset columns only. rows=traits
        nodeLabels = [n.name for n in net.nodes_changed if n.leaf]
    else
        nodeLabels = [n.name == "" ? string(n.number) : n.name for n in net.nodes_changed]
    end
    return M, nodeLabels
end

function randomTrait!(M::Matrix, obj::TSM, net::HybridNetwork)
    recursionPreOrder!(net.nodes_changed, M, # updates M in place
            updateRootRandomTrait!,
            updateTreeRandomTrait!,
            updateHybridRandomTrait!,
            obj)
end

function updateRootRandomTrait!(V::AbstractArray, i::Int, obj)
    sample!(1:nstates(obj), view(V, :, i)) # uniform at the root
    return
end

function updateTreeRandomTrait!(V::Matrix,
    i::Int,parentIndex::Int,edge::Edge,
    obj)
    randomTrait!(view(V, :, i), obj, edge.length, view(V, :, parentIndex))
end

function updateHybridRandomTrait!(V::Matrix,
        i::Int, parentIndex1::Int, parentIndex2::Int,
        edge1::Edge, edge2::Edge, obj)
    randomTrait!(view(V, :, i), obj, edge1.length, view(V, :, parentIndex1))
    tmp = randomTrait(obj, edge2.length, view(V, :, parentIndex2))
    for j in 1:size(V,1) # loop over traits
        if V[j,i] == tmp[j] # both parents of the hybrid node have the same trait
            continue # skip the rest: go to next trait
        end
        if rand() > edge1.gamma
            V[j,i] = tmp[j] # switch to inherit trait of parent 2
        end
    end
end

"""
    JC69(rate, relative)

Jukes Cantor (1969) nucleic acid substitution model, which has a single rate parameter.
`rate` corresponds to the absolute diagonal elements, that is, the rate of change
(to any of the other 2 states). Individual rates are `rate`/3.
If `relative` is true (default), the transition matrix [`Q`](@ref) is normalized
to an average of 1 transition per unit of time: in which case `rate` is set to 1.0.

# examples

```jldoctest
julia> m1 = JC69([0.25], false)
Jukes and Cantor 69 Substitution Model,
absolute rate version
off-diagonal rates equal to 0.25/3.
rate matrix Q:
               A       C       G       T
       A       *  0.0833  0.0833  0.0833
       C  0.0833       *  0.0833  0.0833
       G  0.0833  0.0833       *  0.0833
       T  0.0833  0.0833  0.0833       *

julia> nstates(m1)
4

julia> nparams(m1)
1

julia> m2 = JC69([0.5])
Jukes and Cantor 69 Substitution Model,
relative rate version
off-diagonal rates equal to 1/3
rate matrix Q:
               A       C       G       T
       A       *  0.3333  0.3333  0.3333
       C  0.3333       *  0.3333  0.3333
       G  0.3333  0.3333       *  0.3333
       T  0.3333  0.3333  0.3333       *

julia> nparams(m2)
0
```
"""
struct JC69 <: NucleicAcidSubstitutionModel
    rate::Vector{Float64}
    relative::Bool
    eigeninfo::Vector{Float64} #TODO change to MVector, see if its faster
    function JC69(rate::Vector{Float64}, relative::Bool, eigeninfo::Vector{Float64})
        #Warning: this constructor should not be used directly. Use the constructor below
        #   which will call this.
        0 <= length(rate) <= 1 || error("rate not a valid length for a JC69 model")
        all(rate .> 0.0) || error("All elements of rate must be positive for a JC69 model")
        new(rate, relative, eigeninfo)
    end
end
function JC69(rate::AbstractVector, relative=true::Bool)
    obj = JC69(rate, relative, zeros(1))
    seteigeninfo!(obj)
    return obj
end
JC69(rate::Float64, relative=true::Bool) = JC69([rate], relative)


"""
for [`JC69`](@ref): store lambda = 4/3 (if relative) or rate * 4/3 (absolute).
"""
function seteigeninfo!(obj::JC69)
    if obj.relative
        lambda = (4.0/3.0) #eigen value of the Q matrix
    else
        lambda = (4.0/3.0)*obj.rate[1] #eigen value of the Q matrix
    end
    obj.eigeninfo[1] = lambda
end

function Base.show(io::IO, obj::JC69)
    str = "Jukes and Cantor 69 Substitution Model,\n"
    if obj.relative == true
        str *= "relative rate version\n"
        str *= "off-diagonal rates equal to 1/3\n"
    else
        str *= "absolute rate version\n"
        str *= "off-diagonal rates equal to $(round(obj.rate[1], digits=5))/3.\n"
    end
    str *= "rate matrix Q:\n"
    print(io, str)
    showQ(io, obj)
end

"""
    HKY85(rate, pi, relative)

A nucleic acid substitution model based on Hasegawa et al. 1985 substitution model.
`rate` should be a vector of 1 or 2 rates, and `pi` a vector of 4 probabilities summing to 1.

If `relative` is false, the 2 rates represent the transition rate and the transversion rate,
α and β. If `relative` is true (default), only the first rate is used and represents the
transition/transversion ratio: κ=α/β. The rate transition matrix Q is normalized to have
1 change / unit of time on average, i.e. the absolute version of Q is divided by
`2(piT*piC + piA*piG)α + 2(piY*piR)β`.

`nparams` returns 1 or 2.
In other words: the stationary distribution is not counted in the number of parameters
(and `fitdiscrete` does not optimize the pi values at the moment).

# examples

```jldoctest
julia> m1 = HKY85([.5], [0.20, 0.30, 0.30, 0.20])
HKY85 Substitution Model base frequencies: [0.2, 0.3, 0.3, 0.2]
relative rate version with transition/tranversion ratio kappa = 0.5,
 scaled so that there is one substitution per unit time
rate matrix Q:
               A       C       G       T
       A       *  0.4839  0.2419  0.3226
       C  0.3226       *  0.4839  0.1613
       G  0.1613  0.4839       *  0.3226
       T  0.3226  0.2419  0.4839       *

julia> nstates(m1)
4

julia> m2 = HKY85([0.5, 0.5], [0.20, 0.30, 0.30, 0.20], false)
HKY85 Substitution Model base frequencies: [0.2, 0.3, 0.3, 0.2]
absolute rate version with transition/transversion ratio kappa = a/b = 1.0
 with rates a = 0.5 and b = 0.5
rate matrix Q:
               A       C       G       T
       A       *  0.1500  0.1500  0.1000
       C  0.1000       *  0.1500  0.1000
       G  0.1000  0.1500       *  0.1000
       T  0.1000  0.1500  0.1500       *

```
"""
struct HKY85 <: NucleicAcidSubstitutionModel
    rate::Vector{Float64}
    pi::Vector{Float64}
    relative::Bool
    eigeninfo::Vector{Float64}
    function HKY85(rate::Vector{Float64}, pi::Vector{Float64}, relative::Bool, eigeninfo::Vector{Float64})
        #Warning: this constructor should not be used directly. Use the constructor below,
        #   which will call this.
        all(rate .> 0.) || error("All elements of rate must be positive")
        1 <= length(rate) <= 2 || error("rate has invalid length for HKY85 model")
        if relative length(rate) == 1 || error("the relative version of HKY85 takes a single rate")
        else        length(rate) == 2 || error("the absolute version of HKY85 takes 2 rates")
        end
        length(pi) == 4 || error("pi must be of length 4")
        all(0. .< pi.< 1.) || error("All base proportions must be between 0 and 1")
        isapprox(sum(pi), 1.; atol = 1e-12) || error("Base proportions must sum to 1.")
        new(rate, pi, relative, eigeninfo)
    end
end
function HKY85(rate::AbstractVector, pi::Vector{Float64}, relative=true::Bool)
    obj = HKY85(rate, pi, relative, zeros(5))
    seteigeninfo!(obj)
    return obj
end
HKY85(rate::Float64, pi::Vector{Float64}) = HKY85([rate], pi, true)

"""
for [`HKY85`](@ref): store piR, piY, the 2 non-zero eigenvalues and a scaling factor
"""
function seteigeninfo!(obj::HKY85)
    piA = obj.pi[1]; piC = obj.pi[2]; piG = obj.pi[3]; piT = obj.pi[4]
    piR = piA + piG
    piY = piT + piC
    obj.eigeninfo[4] = piR
    obj.eigeninfo[5] = piY

    a = obj.rate[1] # relative: a = kappa. absolute: a = kappa*b
    if obj.relative
        # b=1/lambda: scaling factor to have branch lengths in substitutions/site
        b = 1.0 / (2*a*(piT*piC + piA*piG) + 2*(piY*piR))
        obj.eigeninfo[2] = - (piR * a + piY) * b # lambda_R
        obj.eigeninfo[3] = - (piY * a + piR) * b # lambda_Y
    else
        b = obj.rate[2]
        obj.eigeninfo[2] = -((piR * a) + (piY * b))
        obj.eigeninfo[3] = -((piY * a) + (piR * b))
    end
    obj.eigeninfo[1] = -b
end

function Base.show(io::IO, obj::HKY85)
    str = "HKY85 Substitution Model base frequencies: $(obj.pi)\n"
    if obj.relative
        str *= "relative rate version with transition/tranversion ratio kappa = $(round(obj.rate[1], digits=5)),"
        str *= "\n scaled so that there is one substitution per unit time\n"
    else
        str *= "absolute rate version with transition/transversion ratio kappa = a/b = "
        str *= "$(round(obj.rate[1]/obj.rate[2], digits=5))"
        str *= "\n with rates a = $(round(obj.rate[1], digits=5)) and b = $(round(obj.rate[2], digits=5))\n"
    end
    str *= "rate matrix Q:\n"
    print(io, str)
    showQ(io, obj)
end

"""
for `JC69` model: 0 if relative, 1 if absolute
"""
function nparams(obj::JC69)
    return (obj.relative ? 0 : 1)
end

"""
for `HKY85` model: 1 if relative, 2 if absolute
"""
function nparams(obj::HKY85)
    return (obj.relative ? 1 : 2)
end

@inline function Q(obj::JC69)
    if obj.relative
        Q0 = 1.0/3.0
        #we multiply by 1/3 to make branch lengths interpretable as avg(#substitutions) (pi-weighted avg of diagonal is -1)
    else
        Q0 = (1.0/3.0)*obj.rate[1]
    end

    return Qmatrix(-3*Q0, Q0, Q0, Q0,
                    Q0, -3*Q0, Q0, Q0,
                    Q0, Q0, -3*Q0, Q0,
                    Q0, Q0, Q0, -3*Q0)
end

@inline function Q(obj::HKY85)
    piA = obj.pi[1]; piC = obj.pi[2]; piG = obj.pi[3]; piT = obj.pi[4]
    piR = piA + piG
    piY = piT + piC
    if obj.relative
        k = obj.rate[1]
        lambda = (2*k*(piT*piC + piA*piG) + 2*(piY*piR)) # for sub/time interpretation. see HKY docstring
        Q₁  = piA/lambda
        Q₂  = k * Q₁
        Q₃  = piC/lambda
        Q₄  = k * Q₃
        Q₆  = piG/lambda
        Q₅  = k * Q₆
        Q₇  = piT/lambda
        Q₈  = k * Q₇
        Q₉  = -(Q₃ + Q₅ + Q₇)
        Q₁₀ = -(Q₁ + Q₆ + Q₈)
        Q₁₁ = -(Q₂ + Q₃ + Q₇)
        Q₁₂ = -(Q₁ + Q₄ + Q₆)
    else #GT, AT, CG, and AC are less likely (transversions)
        a = obj.rate[1]; b = obj.rate[2]
        Q₁  = b * piA
        Q₂  = a * piA
        Q₃  = b * piC
        Q₄  = a * piC
        Q₅  = a * piG
        Q₆  = b * piG
        Q₇  = b * piT
        Q₈  = a * piT
        Q₉  = -(Q₃ + Q₅ + Q₇) #AA
        Q₁₀ = -(Q₁ + Q₆ + Q₈) #CC
        Q₁₁ = -(Q₂ + Q₃ + Q₇) #GG
        Q₁₂ = -(Q₁ + Q₄ + Q₆) #TT
    end
    return Qmatrix(Q₉,  Q₁,  Q₂,  Q₁, #ACGT
                    Q₃,  Q₁₀, Q₃,  Q₄,
                    Q₅,  Q₆,  Q₁₁, Q₆,
                    Q₇,  Q₈,  Q₇,  Q₁₂)
end

function P!(Pmat::AbstractMatrix, obj::JC69, t::Float64)
    P0 = 0.25 + 0.75 * exp(-t * obj.eigeninfo[1]) #lambda
    P1 = 0.25 - 0.25 * exp(-t * obj.eigeninfo[1])
    # P1 off-diagonal and P0 on diagonal
    Pmat[1,1] = P0; Pmat[2,1] = P1; Pmat[3,1] = P1; Pmat[4,1] = P1
    Pmat[1,2] = P1; Pmat[2,2] = P0; Pmat[3,2] = P1; Pmat[4,2] = P1
    Pmat[1,3] = P1; Pmat[2,3] = P1; Pmat[3,3] = P0; Pmat[4,3] = P1
    Pmat[1,4] = P1; Pmat[2,4] = P1; Pmat[3,4] = P1; Pmat[4,4] = P0
    return Pmat
end

function P!(Pmat::AbstractMatrix, obj::HKY85, t::Float64)
    #returns ACGT
    piA = obj.pi[1]; piC = obj.pi[2]; piG = obj.pi[3]; piT = obj.pi[4]
    piR = obj.eigeninfo[4]
    piY = obj.eigeninfo[5]
    # expm1(x) = e^x-1 accurately. important when t is small
    ebm1 = expm1(obj.eigeninfo[1] * t) # -b eigevalue
    eRm1 = expm1(obj.eigeninfo[2] * t) # lambda_R
    eYm1 = expm1(obj.eigeninfo[3] * t) # lambda_Y
    # transversions
    P_TA_CA = - piA * ebm1 # C or T -> A: P{A|T} = P{A|C}
    P_AC_GC = - piC * ebm1 # P{C|A} = P{C|G}
    P_TG_CG = - piG * ebm1 # P{G|T} = P{G|C}
    P_AT_GT = - piT * ebm1 # P{T|A} = P{T|G}
    # transitions: A<->G (R) or C<->T (Y)
    tmp = (ebm1 * piY - eRm1) / piR
    P_GA = piA * tmp # P{A|G}
    P_AG = piG * tmp # P{G|A}
    tmp = (ebm1 * piR - eYm1) / piY
    P_TC = piC * tmp # P{C|T}
    P_CT = piT * tmp # P{T|C}
    # no change
    transv = P_AC_GC + P_AT_GT # transversion to Y: - piY * ebm1
    P_AA = 1.0 - P_AG - transv
    P_GG = 1.0 - P_GA - transv
    transv = P_TA_CA + P_TG_CG # transversion to R
    P_CC = 1.0 - P_CT - transv
    P_TT = 1.0 - P_TC - transv
    # fill in the P matrix
    Pmat[1,1] = P_AA
    Pmat[2,2] = P_CC
    Pmat[3,3] = P_GG
    Pmat[4,4] = P_TT
    # to A: j=1
    Pmat[2,1] = P_TA_CA
    Pmat[3,1] = P_GA
    Pmat[4,1] = P_TA_CA
    # to C: j=2
    Pmat[1,2] = P_AC_GC
    Pmat[3,2] = P_AC_GC
    Pmat[4,2] = P_TC
    # to G: j=3
    Pmat[1,3] = P_AG
    Pmat[2,3] = P_TG_CG
    Pmat[4,3] = P_TG_CG
    # to T: j=4
    Pmat[1,4] = P_AT_GT
    Pmat[2,4] = P_CT
    Pmat[3,4] = P_AT_GT
    return Pmat
end

abstract type RateVariationAcrossSites end

"""
    RateVariationAcrossSites(; pinv=0.0, alpha=Inf, ncat=4)

Model for variable substitution rates across sites (or across traits) using
the discrete Gamma model (+G, Yang 1994, Journal of Molecular Evolution) or
the invariable-sites model (+I, Hasegawa, Kishino & Yano 1985 J Mol Evol).
Both types of rate variation can be combined (+G+I, Gu, Fu & Li 1995, Mol Biol Evol)
but this is discouraged (Jia, Lo & Ho 2014 PLOS One).
Using rate variation increases the number of parameters by one (+G or +I)
or by two (+G+I).

Because the mean of the desired distribution or rates is 1, we use a Gamma
distribution with shape α and scale θ=1/α (rate β=α) if no invariable sites,
or scale θ=1/(α(1-pinv)), that is rate β=α(1-pinv) with a proportion pinv
of invariable sites.
The shape parameter is referred to as alpha here.
The Gamma distribution is discretized into `ncat` categories.
In each category, the category's rate multiplier is a normalized quantile of the gamma distribution.

```jldoctest
julia> rv = RateVariationAcrossSites()
Rate variation across sites: discretized Gamma
categories for Gamma discretization: 1
rates: [1.0]

julia> nparams(rv)
0

julia> typeof(rv)
PhyloNetworks.RVASGamma{1}

julia> rv = RateVariationAcrossSites(alpha=1.0, ncat=4)
Rate variation across sites: discretized Gamma
alpha: 1.0
categories for Gamma discretization: 4
rates: [0.146, 0.513, 1.071, 2.27]

julia> typeof(rv)
PhyloNetworks.RVASGamma{4}

julia> PhyloNetworks.setalpha!(rv, 2.0)
Rate variation across sites: discretized Gamma
alpha: 2.0
categories for Gamma discretization: 4
rates: [0.319, 0.683, 1.109, 1.889]

julia> nparams(rv)
1

julia> rv = RateVariationAcrossSites(pinv=0.3)
Rate variation across sites: +I (invariable sites)
pinv: 0.3
rates: [0.0, 1.429]

julia> nparams(rv)
1

julia> typeof(rv)
PhyloNetworks.RVASInv

julia> PhyloNetworks.setpinv!(rv, 0.05)
Rate variation across sites: +I (invariable sites)
pinv: 0.05
rates: [0.0, 1.053]

julia> rv = RateVariationAcrossSites(pinv=0.3, alpha=2.0, ncat=4)
Rate variation across sites: discretized Gamma+I
pinv: 0.3
alpha: 2.0
categories for Gamma discretization: 4
rates: [0.0, 0.456, 0.976, 1.584, 2.698]
probabilities: [0.3, 0.175, 0.175, 0.175, 0.175]

julia> nparams(rv)
2

julia> typeof(rv)
PhyloNetworks.RVASGammaInv{5}

julia> PhyloNetworks.setalpha!(rv, 3.0)
Rate variation across sites: discretized Gamma+I
pinv: 0.3
alpha: 3.0
categories for Gamma discretization: 4
rates: [0.0, 0.6, 1.077, 1.584, 2.454]
probabilities: [0.3, 0.175, 0.175, 0.175, 0.175]

julia> PhyloNetworks.setpinv!(rv, 0.05)
Rate variation across sites: discretized Gamma+I
pinv: 0.05
alpha: 3.0
categories for Gamma discretization: 4
rates: [0.0, 0.442, 0.793, 1.167, 1.808]
probabilities: [0.05, 0.238, 0.238, 0.238, 0.238]

julia> PhyloNetworks.setpinvalpha!(rv, 0.1, 5.0)
Rate variation across sites: discretized Gamma+I
pinv: 0.1
alpha: 5.0
categories for Gamma discretization: 4
rates: [0.0, 0.593, 0.91, 1.221, 1.721]
probabilities: [0.1, 0.225, 0.225, 0.225, 0.225]
```
"""
function RateVariationAcrossSites(; pinv=0.0::Float64, alpha=Inf64::Float64, ncat=1::Int)
    ncat>1 && alpha == Inf && error("please specify ncat=1 or alpha<Inf")
    # ncat is 1 or α is finite here
    if pinv==0.0 # if α is infinite: no rate variation, RVASGamma with ncat=1
        return RVASGamma(alpha, ncat)
    end
    # pinv>0 here
    if ncat == 1
        return RVASInv(pinv)
    end
    # pinv>0, ncat>1 and α is finite here
    return RVASGammaInv(pinv, alpha, ncat)
end
"""
    RateVariationAcrossSites(rvsymbol::Symbol, ncategories=4::Int)

Default model for rate variation across site, specified by a symbol:
- `:noRV` for no rate variation
- `:G` or `:Gamma` for gamma-distributed rates
- `:I` or `:Inv` for two categories: invariable and variable
- `:GI` or `:GI` for both.
"""
function RateVariationAcrossSites(rvsymbol::Symbol, ncategories=4::Int)
    if rvsymbol == :noRV
        rvas = RateVariationAcrossSites()
    elseif rvsymbol == :Gamma || rvsymbol == :G
        rvas = RateVariationAcrossSites(alpha=1.0, ncat=ncategories)
    elseif rvsymbol == :GammaInv || rvsymbol == :GI
        rvas = RateVariationAcrossSites(pinv=0.05, alpha=1.0, ncat=ncategories)
    elseif rvsymbol == :Inv || rvsymbol == :I
        rvas = RateVariationAcrossSites(pinv=0.05)
    else
        error("model $rvsymbol unknown or not implemented yet:\nrate variation model needs to be :Gamma or :Inv or :GammaInv")
    end
end

struct RVASGamma{S} <: RateVariationAcrossSites
    # S = ncat, and size of vectors
    alpha::StaticArrays.MVector{1,Float64} # mutable
    ncat::Int
    ratemultiplier::StaticArrays.MVector{S,Float64}
    lograteweight::StaticArrays.SVector{S,Float64} # will be uniform: log(1/ncat)
end
function RVASGamma(alpha=1.0::Float64, ncat=4::Int)
    @assert ncat > 0 "ncat must be 1 or greater"
    uniflw = -log(ncat) # = log(1/ncat)
    obj = RVASGamma{ncat}(
            StaticArrays.MVector{1,Float64}(alpha), ncat,
            StaticArrays.MVector{ncat,Float64}(undef), # rates
            StaticArrays.SVector{ncat,Float64}([uniflw for i in 1:ncat]))
    if ncat == 1
        obj.ratemultiplier[1] = 1.0
    else
        setalpha!(obj, alpha) # checks for alpha >= 0
    end
    return obj
end

struct RVASInv <: RateVariationAcrossSites
    pinv::StaticArrays.MVector{1,Float64} # mutable
    ratemultiplier::StaticArrays.MVector{2,Float64}
    lograteweight::StaticArrays.MVector{2,Float64}
end
function RVASInv(pinv=0.05::Float64)
    r = StaticArrays.MVector{2,Float64}(undef) # rates
    r[1] = 0.0 # invariable category
    obj = RVASInv(StaticArrays.MVector{1,Float64}(pinv),
            r,
            StaticArrays.MVector{2,Float64}(undef)) # log weights
    setpinv!(obj, pinv) # checks for 0 <= pinv < 1
    return obj
end

struct RVASGammaInv{S} <: RateVariationAcrossSites
    # S = ncat+1, and size of vectors
    pinv::StaticArrays.MVector{1,Float64}  # mutable
    alpha::StaticArrays.MVector{1,Float64} # mutable
    ncat::Int
    ratemultiplier::StaticArrays.MVector{S,Float64}
    lograteweight::StaticArrays.MVector{S,Float64}
end
function RVASGammaInv(pinv::Float64, alpha::Float64, ncat::Int)
    @assert ncat > 1 "ncat must be 2 or more for the Gamma+I model"
    s = 1+ncat
    r = StaticArrays.MVector{s,Float64}(undef) # rates
    r[1] = 0.0 # invariable category
    obj = RVASGammaInv{s}(
            StaticArrays.MVector{1,Float64}(pinv),
            StaticArrays.MVector{1,Float64}(alpha), ncat,
            r,
            StaticArrays.MVector{s,Float64}(undef)) # log weights
    setpinvalpha!(obj, pinv, alpha) # checks for α >= 0 and 0 <= pinv < 1
    return obj
end

"""
    setalpha!(obj, alpha)

Set the shape parameter `alpha` in a RateVariationAcrossSites model `obj`,
and update the rate multipliers accordingly.
Return the modified object.
"""
function setalpha!(obj::RVASGamma{S}, alpha::Float64) where S
    @assert alpha >= 0 "alpha must be >= 0"
    obj.alpha[1] = alpha
    gammadist = Distributions.Gamma(alpha, 1/alpha)
    cumprob = 1/(2obj.ncat) .+ (0:(obj.ncat-1))/obj.ncat # cumulative prob to discretize Gamma
    obj.ncat > 1 || return obj
    rv = obj.ratemultiplier
    for (i,cp) in enumerate(cumprob)
        @inbounds rv[i] = quantile(gammadist, cp)
    end
    rv ./= mean(rv)
    return obj
end
function setalpha!(obj::RVASGammaInv{S}, alpha::Float64) where S
    @assert alpha >= 0 "alpha must be >= 0"
    obj.alpha[1] = alpha
    ncat = obj.ncat
    pvar = 1.0 - obj.pinv[1]
    gammadist = Distributions.Gamma(alpha, 1/alpha)
    r0 = quantile.(gammadist, 1/(2ncat) .+ (0:(ncat-1))/ncat)
    r0 ./= mean(r0)
    for i in 2:(ncat+1)
        @inbounds obj.ratemultiplier[i] = r0[i-1]/pvar
    end
    return obj
end

"""
    setpinv!(obj, pinv)

Set the proportion of invariable sites `pinv` in a RateVariationAcrossSites
model `obj`, and update the rate multipliers & weights accordingly.
For `RVASInvGamma` objects, the original rate multipliers are assumed correct,
according to the original `pinv` value.
Return the modified object.
"""
function setpinv!(obj::RVASInv, pinv::Float64)
    @assert 0.0 <= pinv < 1.0 "pinv must be in [0,1)"
    obj.pinv[1] = pinv
    pvar = 1.0-pinv # 0 not okay here: ratemultiplier would be infinite
    obj.lograteweight[1] = log(pinv) # -Inf is okay
    obj.lograteweight[2] = log(pvar)
    obj.ratemultiplier[2] = 1.0/pvar # to get an average rate = 1
    return obj
end
function setpinv!(obj::RVASGammaInv{S}, pinv::Float64) where S
    @assert 0.0 <= pinv < 1.0 "pinv must be in [0,1)"
    ncat = obj.ncat
    pvar = 1.0-pinv # 0 not okay here: ratemultiplier would be infinite
    pvarratio = (1.0-obj.pinv[1])/pvar # old p_variable / new p_variable
    obj.pinv[1] = pinv
    obj.lograteweight[1] = log(pinv) # -Inf is okay
    uniflw = -log(ncat)+log(pvar)
    for i in 2:(ncat+1)
        @inbounds obj.lograteweight[i] = uniflw
        @inbounds obj.ratemultiplier[i] *= pvarratio # gamma rate / p_variable
    end
    return obj
end

"""
    setpinvalpha!(obj, pinv, alpha)

Set the proportion of invariable sites `pinv` and the `alpha` parameter for
the discretized gamma distribution in a model `obj` of type `RVASGammaInv{S}`.
Update the rate multipliers & weights accordingly.
The mean of the distribution is constrained to 1.

Return the modified object.
"""
function setpinvalpha!(obj::RVASGammaInv{S}, pinv::Float64, alpha::Float64) where S
    @assert 0.0 <= pinv < 1.0 "pinv must be in [0,1)"
    @assert alpha >= 0 "alpha must be >= 0"
    obj.alpha[1] = alpha
    obj.pinv[1] = pinv
    obj.lograteweight[1] = log(pinv) # -Inf is okay
    ncat = obj.ncat
    gammadist = Distributions.Gamma(alpha, 1/alpha)
    r0 = quantile.(gammadist, 1/(2ncat) .+ (0:(ncat-1))/ncat)
    r0 ./= mean(r0)
    pvar = 1.0-pinv # 0 not okay here: ratemultiplier would be infinite
    uniflw = -log(ncat)+log(pvar)
    for i in 2:(ncat+1)
        @inbounds obj.lograteweight[i] = uniflw
        @inbounds obj.ratemultiplier[i] = r0[i-1]/pvar
    end
    return obj
end

function Base.show(io::IO, obj::RVASGamma{S})  where S
    str = "Rate variation across sites: discretized Gamma\n"
    if length(obj.ratemultiplier)>1
        str *= "alpha: $(round(obj.alpha[1], digits=5))\n"
    end
    str *= "categories for Gamma discretization: $(obj.ncat)\n"
    str *= "rates: $(round.(obj.ratemultiplier, digits=3))"
    print(io, str)
end
function Base.show(io::IO, obj::RVASInv)
    str = "Rate variation across sites: +I (invariable sites)\n"
    str *= "pinv: $(round(obj.pinv[1], digits=5))\n"
    str *= "rates: $(round.(obj.ratemultiplier, digits=3))"
    print(io, str)
end
function Base.show(io::IO, obj::RVASGammaInv{S})  where S
    str = "Rate variation across sites: discretized Gamma+I\n"
    str *= "pinv: $(round(obj.pinv[1], digits=5))\n"
    str *= "alpha: $(round(obj.alpha[1], digits=5))\n"
    str *= "categories for Gamma discretization: $(obj.ncat)\n"
    str *= "rates: $(round.(obj.ratemultiplier, digits=3))\n"
    str *= "probabilities: $(round.(exp.(obj.lograteweight), digits=3))"
    print(io, str)
end

function nparams(obj::RVASGamma{S}) where S
    return (obj.ncat == 1 ? 0 : 1)
end
nparams(::RVASInv) = 1::Int
nparams(::RVASGammaInv{S}) where S = 2::Int # ncat must be >1

"""
    getparameters(obj::RateVariationAcrossSites)

Return a copy of the alpha and/or pinv parameters of model `obj`,
in a single vector.
"""
getparameters(obj::RVASInv) = copy(obj.pinv)
getparameters(obj::RVASGamma{S}) where S = copy(obj.alpha)
getparameters(obj::RVASGammaInv{S}) where S = [obj.pinv[1], obj.alpha[1]]

"""
    setparameters!(obj::RateVariationAcrossSites, par::AbstractVector)

Set the values of the alpha and/or pinv parameters of model `obj`.
See also [`setalpha!`](@ref), [`setpinv!`](@ref) and [`setpinvalpha!`](@ref)
"""
setparameters!(obj::RVASInv, par::AbstractVector) = setpinv!(obj, par[1])
setparameters!(obj::RVASGamma{S}, par::AbstractVector) where S =
    setalpha!(obj, par[1])
setparameters!(obj::RVASGammaInv{S}, par::AbstractVector) where S =
    setpinvalpha!(obj, par[1], par[2])

"""
    getparamindex(obj::RateVariationAcrossSites)

Indices of parameters in (p_invariable, alpha).
"""
getparamindex(::RVASInv) = [1]
getparamindex(::RVASGamma{S}) where S = [2]
getparamindex(::RVASGammaInv{S}) where S = [1,2]

"""
    empiricalDNAfrequencies(DNAdata::AbstractDataFrame, DNAweights,
                            correction=true, useambiguous=true)

Estimate base frequencies in DNA data `DNAdata`, ordered ACGT.

- `DNAdata`: data frame. All columns are used. If the first column
  gives species names, find a way to ignore it before calculating empirical
  frequencies, e.g. `empiricalDNAfrequencies(view(DNAdata, :, 2:size(DNAdata, 2)))`.
  Data type must be `BioSymbols.DNA` or `Char` or `String`.
  WARNING: this is checked on the first column only.
- `DNAweights`: vector of weights, to weigh each column in `DNAdata`.
- `correction`: if `true`, add 1 to each count and 4 to the denominator
  for a more stable estimator, similar to Bayes prior of 1/4 and
  the Agresti-Coull interval in binomial estimation.
- `useambiguous`: if `true`, ambiguous bases are used (except gaps and Ns).
  For example, `Y` adds 0.5 weight to `C` and 0.5 weight to `T`.
"""
function empiricalDNAfrequencies(dnaDat::AbstractDataFrame, dnaWeights::Vector,
    correctedestimate=true::Bool, useambiguous=true::Bool)

    # warning: checking first column and first row only
    dnadat1type = eltype(dnaDat[!,1])
    dnadat1type == BioSymbols.DNA || dnadat1type == Char ||
      dnaDat[1,1] ∈ string.(BioSymbols.alphabet(DNA)) ||
        error("empiricalDNAfrequencies requires data of type String, Char, or BioSymbols.DNA")

    # initialize counts: keys same as BioSymbols.ACGT (which are ordered)
    prior = correctedestimate ? 1.0 : 0.0
    dnacounts = Dict(DNA_A=>prior, DNA_C=>prior, DNA_G=>prior, DNA_T=>prior)

    convert2dna = dnadat1type != BioSymbols.DNA
    for j in 1:size(dnaDat, 2) # for each column
        col = dnaDat[!,j]
        wt = dnaWeights[j]
        for nuc in col      # for each row
            if convert2dna
                nuc = convert(DNA, nuc[1]) # if nuc is a string, nuc[1] = 1st character
            end
            if nuc ∈ BioSymbols.ACGT
                dnacounts[nuc] += wt
            elseif nuc == DNA_Gap || nuc == DNA_N || !useambiguous
                continue # to next row
            # else: ambiguity, uses BioSequences definitions
            elseif nuc == DNA_M # A or C
                dnacounts[DNA_A] += 0.5*wt
                dnacounts[DNA_C] += 0.5*wt
            elseif nuc == DNA_R # A or G
                dnacounts[DNA_A] += 0.5*wt
                dnacounts[DNA_G] += 0.5*wt
            elseif nuc == DNA_W # A or T/U
                dnacounts[DNA_A] += 0.5*wt
                dnacounts[DNA_T] += 0.5*wt
            elseif nuc == DNA_S # C or G
                dnacounts[DNA_C] += 0.5*wt
                dnacounts[DNA_G] += 0.5*wt
            elseif nuc == DNA_Y # C or T/U
                dnacounts[DNA_C] += 0.5*wt
                dnacounts[DNA_T] += 0.5*wt
            elseif nuc == DNA_K # G or T/U
                dnacounts[DNA_G] += 0.5*wt
                dnacounts[DNA_T] += 0.5*wt
            elseif nuc == DNA_V # A or C or G
                dnacounts[DNA_A] += wt/3
                dnacounts[DNA_C] += wt/3
                dnacounts[DNA_G] += wt/3
            elseif nuc == DNA_H # A or C or T
                dnacounts[DNA_A] += wt/3
                dnacounts[DNA_C] += wt/3
                dnacounts[DNA_T] += wt/3
            elseif nuc == DNA_D # A or G or T/U
                dnacounts[DNA_A] += wt/3
                dnacounts[DNA_G] += wt/3
                dnacounts[DNA_T] += wt/3
            elseif nuc == DNA_B # C or G or T/U
                dnacounts[DNA_C] += wt/3
                dnacounts[DNA_G] += wt/3
                dnacounts[DNA_T] += wt/3
            end
        end
    end
    totalweight = sum(values(dnacounts))
    res = [dnacounts[key]/totalweight for key in BioSymbols.ACGT] # to control the order
    all(0. .<= res .<= 1.) || error("""weird: empirical base frequency < 0 or > 1""")
    return res
end

"""
    stationary(substitutionmodel)

Stationary distribution of a Markov model
"""
stationary(mod::SM) = error("stationary not defined for $(typeof(mod)).")

function stationary(mod::JC69) #for JC, stationary = uniform
    return [0.25,0.25,0.25,0.25]
end

function stationary(mod::HKY85)
    return mod.pi
end

function stationary(mod::ERSM) #for ERSM, uniform = stationary
    return [1/mod.k for i in 1:mod.k]
end

function stationary(mod::BTSM)
    return [mod.eigeninfo[2], mod.eigeninfo[3]]
end
