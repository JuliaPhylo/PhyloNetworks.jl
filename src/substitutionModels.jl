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
function Base.show(io::IO, obj::SM)
    error("show not defined for $(typeof(obj)).")
end

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
    print(io, "\n")
    for i = 1:size(M,1) # print one row per state
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
        print(io, "\n")
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



julia> PhyloNetworks.P(m1, 3.0)
4×4 StaticArrays.MArray{Tuple{4,4},Float64,2,16} with indices SOneTo(4)×SOneTo(4):
 0.217509  0.198417  0.190967  0.198417
 0.297625  0.312992  0.297625  0.28645
 0.28645   0.297625  0.312992  0.297625
 0.198417  0.190967  0.198417  0.217509
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
    Pmat = MMatrix{k,k}(repeat([0.0], k*k))
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

julia> PhyloNetworks.P!(Matrix{Float64}(undef,2,2), m1, 3.0) # fills an uninitialized 2x2 matrix of floats
2×2 Array{Float64,2}:
 0.666708  0.333292
 0.666584  0.333416

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
    str *= "rate $(obj.label[2])→$(obj.label[1]) β=$(round(obj.rate[2], digits=5))\n"
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

julia> using Random; Random.seed!(1234);

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
        sum(pi) == 1. || error("Base proportions must sum to 1")
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

    a = obj.rate[1]
    if obj.relative
        lambda = 2*a*(piT*piC + piA*piG) + 2*(piY*piR) #for sub/year interpretation
        b = 1.0/lambda #scaling factor so that branch lengths can be interpreted as substitutions
        obj.eigeninfo[2] = -((piR * a*b) + (piY * b))
        obj.eigeninfo[3] = -((piY * a*b) + (piR * b))
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
    Pmat[:,:] .= P1 # puts P1 on the off-diagonals (in particular)
    for i in 1:4 Pmat[i,i]=P0; end # diagonal
    return Pmat
end

function P!(Pmat::AbstractMatrix, obj::HKY85, t::Float64)
    #returns ACGT
    piA = obj.pi[1]; piC = obj.pi[2]; piG = obj.pi[3]; piT = obj.pi[4]
    piR = obj.eigeninfo[4]
    piY = obj.eigeninfo[5]

    e2 = exp(obj.eigeninfo[1] * t)
    e3 = exp(obj.eigeninfo[2] * t)
    e4 = exp(obj.eigeninfo[3] * t)

    P_AA = piA + ((piA * piY)/piR) * e2 + (piG/piR) * e3
    P_CC = piC + ((piC * piR)/piY) * e2 + (piT/piY) * e4
    P_GG = piG + ((piG * piY)/piR) * e2 + (piA/piR) * e3
    P_TT = piT + ((piT * piR)/piY) * e2 + (piC/piY) * e4

    P_TA_CA = piA * (1 - e2) #double name because TA and CA are equal
    P_GA = piA + (piA*piY/piR)*e2 - (piA/piR)*e3 #P6

    P_AC_GC = piC * (1 - e2) #double name because AC and GC are equal #p7
    P_TC = piC + ((piC * piR)/piY) * e2 - (piC/piY) * e4

    P_AG = piG + ((piG * piY)/piR) * e2 - (piG/piR) * e3 #p9
    P_TG_CG = piG * (1 - e2) #double name because TG and CG are equal #p10

    P_AT_GT = piT * (1 - e2) #double name because AT and GT are equal #P11
    P_CT = piT + ((piT * piR)/piY) * e2 - (piT/piY) * e4

    Pmat[1,1] = P_AA
    Pmat[2,2] = P_CC
    Pmat[3,3] = P_GG
    Pmat[4,4] = P_TT

    Pmat[1,2] = P_TA_CA
    Pmat[1,3] = P_GA
    Pmat[1,4] = P_TA_CA

    Pmat[2,1] = P_AC_GC
    Pmat[2,3] = P_AC_GC
    Pmat[2,4] = P_TC

    Pmat[3,1] = P_AG
    Pmat[3,2] = P_TG_CG
    Pmat[3,4] = P_TG_CG

    Pmat[4,1] = P_AT_GT
    Pmat[4,2] = P_CT
    Pmat[4,3] = P_AT_GT
    return Pmat
end

"""
    RateVariationAcrossSites(α=1.0, ncat=4)

Model for variable substitution rates across sites (or across traits) using
the discrete Gamma model (Yang 1994, Journal of Molecular Evolution).
Turn any SM model to SM + Gamma.
Using this model increases the number of parameters by one.

Because the mean of the desired continuous Gamma distribution is 1, we use
shape α and scale θ=1/α (e.g. rate β=α).
The shape parameter is referred to as alpha here.
The Gamma distribution is discretized into `ncat` categories.
In each category, the category's rate multiplier is a normalized quantile of the gamma distribution.

```jldoctest
julia> rv = RateVariationAcrossSites()
Rate Variation Across Sites using Discretized Gamma Model
alpha: 1.0
categories for Gamma discretization: 4
ratemultiplier: [0.14578, 0.51313, 1.07083, 2.27025]

julia> PhyloNetworks.setalpha!(rv, 2.0)

julia> rv
Rate Variation Across Sites using Discretized Gamma Model
alpha: 2.0
categories for Gamma discretization: 4
ratemultiplier: [0.31907, 0.68336, 1.10898, 1.8886]

julia> RateVariationAcrossSites(2.0, 4)
Rate Variation Across Sites using Discretized Gamma Model
alpha: 2.0
categories for Gamma discretization: 4
ratemultiplier: [0.31907, 0.68336, 1.10898, 1.8886]
```
"""
mutable struct RateVariationAcrossSites
    alpha::Float64
    ncat::Int
    ratemultiplier::Array{Float64}
    function RateVariationAcrossSites(alpha = 1.0::Float64, ncat = 4::Int)
        @assert alpha >= 0.0 "alpha must be >= 0"
        @assert ncat > 0 "ncat must be 1 or greater"
        if ncat == 1
            ratemultiplier = [1.0]
        else
            cuts = Vector((0:(ncat-1))/ncat) + repeat([1/(2ncat)], ncat)
            ratemultiplier = quantile.(Distributions.Gamma(alpha, 1/alpha), cuts)
            ratemultiplier ./= mean(ratemultiplier)
        end
        new(alpha, ncat, ratemultiplier)
    end
end
const RVAS = RateVariationAcrossSites

"""
    setalpha!(obj, alpha)

Set the shape parameter `alpha` in a RateVariationAcrossSites model `obj`,
and update the rate multipliers accordingly.
"""
function setalpha!(obj::RateVariationAcrossSites, alpha::Float64)
    @assert alpha >= 0 "alpha must be a float >= 0"
    obj.alpha = alpha
    if obj.ncat > 1
        rv = obj.ratemultiplier
        cuts = Vector((0:(obj.ncat-1))/obj.ncat) + repeat([1/(2*obj.ncat)], obj.ncat)
        rv[:] = quantile.(Distributions.Gamma(alpha, 1/alpha), cuts)
        rv ./= mean(rv)
    end
    return nothing
end

function Base.show(io::IO, obj::RateVariationAcrossSites)
    str = "Rate Variation Across Sites using Discretized Gamma Model\n"
    str *= "alpha: $(round(obj.alpha, digits=5))\n"
    str *= "categories for Gamma discretization: $(obj.ncat)\n"
    str *= "ratemultiplier: $(round.(obj.ratemultiplier, digits=5))\n"
    print(io, str)
end

function nparams(obj::RateVariationAcrossSites)
    return (obj.ncat == 1 ? 0 : 1)
end

"""
    empiricalDNAfrequencies(DNAdata::AbstractDataFrame, DNAweights,
                            correction=true, useambiguous=true)

Estimate base frequencies in DNA data `DNAdata`, ordered ACGT.

- `DNAdata`: data frame. All columns are used. If the first column
  gives species names, find a way to ignore it before calculating empirical
  frequencies, e.g. `empiricalDNAfrequencies(view(DNAdata, :, 2:ncol(DNAdata)))`.
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
    eltypes(dnaDat)[1] == BioSymbols.DNA || eltypes(dnaDat)[1] == Char ||
      dnaDat[1,1] ∈ string.(BioSymbols.alphabet(DNA)) ||
        error("empiricalDNAfrequencies requires data of type String, Char, or BioSymbols.DNA")

    # initialize counts: keys same as BioSymbols.ACGT (which are ordered)
    prior = correctedestimate ? 1.0 : 0.0
    dnacounts = Dict(DNA_A=>prior, DNA_C=>prior, DNA_G=>prior, DNA_T=>prior)

    convert2dna = eltypes(dnaDat)[1] != BioSymbols.DNA
    for j in 1:ncol(dnaDat) # for each column
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

function stationary(mod::JC69)
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
