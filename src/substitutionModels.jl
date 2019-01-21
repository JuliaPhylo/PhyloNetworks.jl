"""
    SubstitutionModel

Abstract type for substitution models, 
using a continous time Markov model on a phylogeny.
Adapted from the SubstitutionModels module in BioJulia.
The same [`Q`](@ref) and [`P`](@ref) function names are used for the
transition rates and probabilities.

For variable rates, see [`RateVariationAcrossSites`](@ref)

For sub types, see [`NucleicAcidSubstitutionModel`](@ref), [`TraitSubstitutionModel`](@ref)
"""
abstract type SubstitutionModel end #ideally, we'd like this to be SubstitutionModels.SubstitionModel
const SM = SubstitutionModel
const Qmatrix = StaticArrays.SMatrix{4, 4, Float64}
const Pmatrix = StaticArrays.MMatrix{4, 4, Float64}
const Bmatrix = StaticArrays.MMatrix{2, 2, Float64}

"""
    TraitSubstitutionModel

Adapted from the SubstitutionModels module in BioJulia.
The same [`Q`](@ref) and [`P`](@ref) function names are used for the
transition rates and probabilities.

For subtypes, see [`BinaryTraitSubstitutionModel`](@ref),
[`EqualRatesSubstitutionModel`](@ref),
[`TwoBinaryTraitSubstitutionModel`](@ref)
"""
abstract type TraitSubstitutionModel{T} <: SubstitutionModel end #this accepts labels
const TSM = TraitSubstitutionModel{T} where T #T is type of labels

"""
    NucleicAcidSubstitutionModel

Adapted from the SubstitutionModels module in BioJulia and 
depricated PhyloModels.jl by Justin Angevaare
The same [`Q`](@ref) and [`P`](@ref) function names are used for the
transition rates and probabilities.

For subtypes, see [`JC69`](@ref), [`HKY85`](@ref)
"""
abstract type NucleicAcidSubstitutionModel <: SubstitutionModel end
const NASM = NucleicAcidSubstitutionModel

function Base.show(io::IO, obj::SM)
    error("show not defined for $(typeof(obj)).") #TODO 1.0 change to @error
end

"""
    nparams(model)

Number of parameters for a given trait evolution model
(length of field `model.rate`).
"""
nparams(obj::SM) = error("nparams not defined for $(typeof(obj)).")

"""
    setrates!(model, rates)

update rates then call seteigeninfo!() to update a model's eigeninfo
"""
#? if seteigeninfo isnt set for SM, setrates shouldnt work either, right?
function setrates!(obj::SM, rates::AbstractVector) #TODO see pull request for SM packages
    obj.rate = rates
    # seteigeninfo!(obj)
end

"""
    seteigeninfo!(model)

calculate eigeninfo for a model and stores it within the object. Constructor allocate memory then call this function.
"""
seteigeninfo!(obj::SM)  = error("seteigeninfo! not defined for $(typeof(obj)).")
# seteigeninfo!(obj::TSM)  = error("seteigeninfo! not defined for $(typeof(obj)).") #? needed?

#TODO change constructors

"""
    getlabels(model)

return labels for a model
"""
function getlabels(obj::SM)
    error("Model must be of type TraitSubstitutionModel or NucleicAcidSubstitutionModel. Got $(typeof(obj))")
end

"""
    nstates(model)

return number of character states for a given trait evolution model.
"""
nstates(obj::TSM) = error("nstates not defined for $(typeof(obj)).")

"""
    nstates(model::NASM)

return number of character states for a NucleicAcidSubstitutionModel

# Examples

```julia-repl
julia> nstates(JC69())
4
julia> nstates(HKY85([.5], [0.25, 0.25, 0.25, 0.25]))
4
```
"""
function nstates(obj::NASM)
    return 4::Int
end

"""
    getlabels(model::TSM)

Return labels for trait substitution model
"""
function getlabels(obj::TSM)
    return obj.label
end

"""
    getlabels(mod::NASM)

return labels for a given NuceleicAcidSubstitutionModel using BioSymbols symbols.

# examples

```julia-repl 
julia> getlabels(JC69())
[DNA_A, DNA_C, DNA_G, DNA_T]
julia> getlabels(HKY85([.5], [0.25, 0.25, 0.25, 0.25]))
[DNA_A, DNA_C, DNA_G, DNA_T]
````
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
    pad = max(8,maximum(length(getlabels(obj))+1)) #? Error here with length. Removed . Is this okay?
    for i = 1:size(M,2) # print the header
        print(io, lpad(getlabels(obj)[i],(i==1? 2*pad : pad), " "))
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
faster when t is scalar (not an array)
"""
@inline function P(obj::SM, t::Float64)
    t >= 0.0 || error("substitution model: >=0 branch lengths are needed")
    return expm(Q(obj) * t)
end

"""
    P(obj, t::Array{Float64})

When applied to a general substitution model, matrix exponentiation is used.
The time argument `t` can be an array.
#TODO 
"""
function P(obj::SM, t::Array{Float64})
    all(t .>= 0.0) || error("t's must all be positive")
    try
        eig_vals, eig_vecs = eig(Q(obj)) # Only hermitian matrices are diagonalizable by #TODO change for 1.0 to eigen
        # *StaticArrays*. Non-Hermitian matrices should be converted to `Array`first.
        return [eig_vecs * expm(diagm(eig_vals)*i) * eig_vecs' for i in t]
    catch
        eig_vals, eig_vecs = eig(Array(Q(obj)))
        k = nstates(obj)
        return [(eig_vecs * expm(diagm(eig_vals)*i) * inv(eig_vecs)) for i in t] 
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
    seteigeninfo!(obj::BinaryTraitSubstitutionModel)

calculate eigeninfo for a model and store it within the object.
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
# Examples

```julia-repl
julia> m1 = BinaryTraitSubstitutionModel([1.0,2.0], ["low","high"])
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
    str *= "rate $(obj.label[1])→$(obj.label[2]) α=$(obj.rate[1])\n"
    str *= "rate $(obj.label[2])→$(obj.label[1]) β=$(obj.rate[2])\n"
    print(io, str)
end

"""
    P(obj::BTSM, t::Float64)

```julia
julia> m1 = BinaryTraitSubstitutionModel([1.0,2.0], ["low","high"])
julia> P(m1, 1.0)
2×2 StaticArrays.SArray{Tuple{2,2},Float64,2,4}:
 0.683262  0.316738
 0.633475  0.366525
```
"""
@inline function P(obj::BTSM, t::Float64)
    t >= 0.0 || error("substitution model: >=0 branch lengths are needed")
    e1 = exp(-obj.eigeninfo[1]*t)
    a0= obj.eigeninfo[2]*e1
    a1= obj.eigeninfo[3]*e1
    return Bmatrix(obj.eigeninfo[2]+a1, obj.eigeninfo[2]-a0, obj.eigeninfo[3]-a1, obj.eigeninfo[3]+a0) # by columns
end

"""
    P!(Pmat::AbstractMatrix, obj::BTSM, t::Float64)
modifies P rate matrix (see traitsLikDiscrete.jl) during optimization. 

```julia-repl 
julia> m1 = BinaryTraitSubstitutionModel([1.0,2.0], ["low","high"])
julia> P!(P(m1, 1.0), m1, 3.0)
2×2 StaticArrays.MArray{Tuple{2,2},Float64,2,4}:
 0.666708  0.333292
 0.666584  0.333416
````
"""
function P!(Pmat::AbstractMatrix, obj::BTSM, t::Float64)
    e1 = exp(-obj.eigeninfo[1]*t)
    a0= obj.eigeninfo[2]*e1
    a1= obj.eigeninfo[3]*e1
    Pmat[1,1] = obj.eigeninfo[2]+a1
    Pmat[2,2] = obj.eigeninfo[3]+a0
    Pmat[1,2] = obj.eigeninfo[3]-a1
    Pmat[2,1] = obj.eigeninfo[2]-a0
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

[`TraitSubstitutionModel`](@ref) for traits with any k number of states
and equal substitution rates α between all states.
Default labels are "1","2",...

# example
```julia
julia> m1 = EqualRatesSubstitutionModel(2, [1.0], ["low","high"])
Equal Rates Substitution Model with k=2,
all rates equal to α=1.0.
rate matrix Q:
             low    high
     low       *  1.0000
    high  1.0000       *
````
#? Is rate alpha equal to 1/k? Doesnt seem so so do we need rate input here? could cause conflicts
"""
mutable struct EqualRatesSubstitutionModel{T} <: TraitSubstitutionModel{T}
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
    #obj = BinaryTraitSubstitutionModel{eltype(label)}(r::AbstractVector, label::AbstractVector, zeros(3)) # Vector{Float64}(undef,3) for julia v1.0
    obj = EqualRatesSubstitutionModel{eltype(label)}(k::Int, rate::Vector{Float64}, label::AbstractVector, zeros(1)) # Vector{Float64}(undef,3) for julia v1.0
    seteigeninfo!(obj)
    return obj
end
#EqualRatesSubstitutionModel(k::Int, α::Float64, label::AbstractVector) = EqualRatesSubstitutionModel{eltype(label)}(k,[α],label)
#EqualRatesSubstitutionModel(k::Int, α::Float64) = EqualRatesSubstitutionModel{String}(k, [α], string.(1:k))
#this isnt working ^
#TODO check this
"""
    seteigeninfo!(obj::EqualRatesSubstitutionModel)

calculate eigeninfo for a model and store it within the object.
"""
#? obj.rate[1] = alpha = 1/k, right? we store this multiple times.
function seteigeninfo!(obj::ERSM) 
    q = (obj.k/(obj.k-1)) 
    obj.eigeninfo[1] = 1
end

function nstates(obj::ERSM)
    return obj.k
end
nparams(::ERSM) = 1::Int

function Base.show(io::IO, obj::ERSM)
    str = "Equal Rates Substitution Model with k=$(obj.k),\n"
    str *= "all rates equal to α=$(obj.rate[1]).\n"
    str *= "rate matrix Q:\n"
    print(io, str)
    showQ(io, obj)
end

#? alpha = 1/k, right? Why do we store obj.rate = alpha and k if this is the case? could cause conflicts
function Q(obj::ERSM)
    α = obj.rate[1]
    M = fill(α, (obj.k,obj.k))
    d = -(obj.k-1) * α
    for i in 1:obj.k
        M[i,i] = d
    end
    return M
end
"""
    P(obj::ERSM)
```julia
julia> m1 = EqualRatesSubstitutionModel(2, [1.0], ["low","high"]);
julia> P(m1, 2.0)
2×2 StaticArrays.MArray{Tuple{2,2},Float64,2,4}:
 1.0       0.864665
 0.864665  1.0
```
"""
function P(obj::ERSM, t::Float64)
    p1 = exp(-obj.eigeninfo[1]*t) + obj.rate[1]*(1-exp(-obj.eigeninfo[1]*t))
    p0 = obj.rate[1]*(1-exp(-obj.eigeninfo[1]*t))
    return Bmatrix(p1, p0, 
                    p0, p1)
end

"""
    P!(Pmat::AbstractMatrix, obj::ERSM, t::Float64)
modifies P rate matrix (see traitsLikDiscrete.jl) during optimization. 

```julia-repl 
julia> m1 = EqualRatesSubstitutionModel(2, [1.0], ["low","high"]);
julia> P!(P(m1, 1.0), m1, 3.0)
2×2 StaticArrays.MArray{Tuple{2,2},Float64,2,4}:
 1.0       0.950213
 0.950213  1.0
````
"""
function P!(Pmat::AbstractMatrix, obj::ERSM, t::Float64)
    e1 = exp(-obj.eigeninfo[1]*t)
    p1 = e1 + obj.rate[1]*(1-e1)
    p0 = obj.rate[1]*(1-e1)
    Pmat[1,1] = p1
    Pmat[2,2] = p1
    Pmat[1,2] = p0
    Pmat[2,1] = p0
    return Pmat
end

"""
    randomTrait(model, t, start)
    randomTrait!(end, model, t, start)

Simulate traits along one edge of length t.
`start` must be a vector of integers, each representing the starting value of one trait.
The bang version (ending with !) uses the vector `end` to store the simulated values.

# Examples
```julia-repl
julia> m1 = BinaryTraitSubstitutionModel(1.0, 2.0)

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
function randomTrait(obj::TSM, t::Float64, start::AbstractVector{Int})
    res = Vector{Int}(length(start)) #this should be an abstract array?
    randomTrait!(res, obj, t, start)
end

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

```julia-repl
julia> m1 = BinaryTraitSubstitutionModel(1.0, 2.0, ["low","high"]);
julia> net = readTopology("(((A:4.0,(B:1.0)#H1:1.1::0.9):0.5,(C:0.6,#H1:1.0::0.1):1.0):3.0,D:5.0);");
julia> srand(1234);
julia> trait, lab = randomTrait(m1, net)
([1 2 … 1 1], String["-2", "D", "-3", "-6", "C", "-4", "#H1", "B", "A"])
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

function randomTrait(obj::TSM, net::HybridNetwork;
    ntraits=1::Int, keepInternal=true::Bool, checkPreorder=true::Bool)
    net.isRooted || error("net needs to be rooted for preorder recursion")
    if(checkPreorder)
        preorder!(net)
    end
    nnodes = net.numNodes
    M = Matrix{Int}(ntraits, nnodes) # M[i,j]= trait i for node j
    randomTrait!(M,obj,net)
    if !keepInternal
        M = getTipSubmatrix(M, net, indexation=:cols) # subset columns only. rows=traits
        nodeLabels = [n.name for n in net.nodes_changed if n.leaf]
    else
        nodeLabels = [n.name == "" ? string(n.number) : n.name for n in net.nodes_changed]    
    end
    return M, nodeLabels
end

@doc (@doc randomTrait) randomTrait!
function randomTrait!(M::Matrix{Int}, obj::TSM, net::HybridNetwork)
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
    JC69()

[`JC69`](@ref) nucleic acid substitution model based on Jukes and Cantor's 1969 substitution model. 
Default is relative rate model. Based on depricated PhyloModels.jl by Justin Angevaare

# examples

```julia-repl
julia> m1 = JC69([0.25], false)
julia> nstates(m1)
4
julia> m2 = JC69([0.5])
julia> nparams(m2)
1
```
"""
struct JC69 <: NucleicAcidSubstitutionModel
    rate::Vector{Float64}
    relative::Bool
    eigeninfo::Vector{Float64} #TODO change to MVector, see if its faster
    
    function JC69(rate::Vector{Float64}, relative::Bool, eigeninfo::Vector{Float64})
        #Warning: this constructor should not be used directly. Use the constructor below
        #   which will call this.        
        if !(0 <= length(rate) <= 1)
            error("rate not a valid length for a JC69 model")
        elseif any(rate .<= 0.)
            error("All elements of rate must be positive for a JC69 model")
        end
        new(rate, relative, eigeninfo)
    end
end
function JC69(rate::AbstractVector, relative=true::Bool)
    obj = JC69(rate, relative, zeros(1)) # Vector{Float64}(undef,1) for julia v1.0 mutable static array
    seteigeninfo!(obj)
    return obj
end
JC69(rate::Float64, relative=true::Bool) = JC69([rate], relative)


"""
    seteigeninfo!(obj::JC69)

calculate eigeninfo for a model and stores it within the object.
"""
function seteigeninfo!(obj::JC69) 
    if obj.relative
        lambda = (4.0/3.0) #to make branch lengths interpretable. lambda = substitutions/year
    else
        lambda = (4.0/3.0)*obj.rate[1]
    end
    obj.eigeninfo[1] = lambda
end

function Base.show(io::IO, obj::JC69)
    str = "Jukes and Cantor 69 Substitution Model,\n"
    if obj.relative == true
        str *= "relative rate version\n"
        str *= "all rates equal to 0.25.\n"
    else
        str *= "absolute rate version\n"
        str *= "all rates equal to $(obj.rate).\n"
    end
    str *= "rate matrix Q:\n"
    print(io, str)
    showQ(io, obj)
end

"""
    HKY85(rate, pi, relative)

[`HKY85`](@ref) nucleic acid substitution model based on Hasegawa et al. 1984 substitution model. 
Default is relative rate model. Based on depricated PhyloModels.jl by Justin Angevaare

For the relative HKY85 model, there are 
# examples

```julia-repl
julia> m1 = HKY85([.5], [0.25, 0.25, 0.25, 0.25])
julia> nstates(m1)
4
julia> m2 = HKY85([0.5, 0.5], [0.25, 0.25, 0.25, 0.25], false)
julia> nstates(m2)
4
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
        if any(rate .<= 0.)
            error("All elements of rate must be positive")
        elseif !(1 <= length(rate) <= 2)
            error("rate is not a valid length for HKY85 model")
        elseif length(pi) !== 4
            error("pi must be of length 4")
        elseif !all(0. .< pi.< 1.)
            error("All base proportions must be between 0 and 1")
        elseif sum(pi) !== 1.
            error("Base proportions must sum to 1")
        end
      
        if length(rate) == 1
            new(rate, pi, true, eigeninfo) #check this because of defaults
        else
            new(rate, pi, false, eigeninfo)
        end
    end
end
function HKY85(rate::AbstractVector, pi::Vector{Float64}, relative=true::Bool)
    obj = HKY85(rate, pi, relative, zeros(5)) # Vector{Float64}(undef,3) for julia v1.0
    seteigeninfo!(obj)
    return obj
end
HKY85(rate::Float64, pi::Vector{Float64}, relative=true::Bool) = HKY85([rate], pi, relative)

"""
    seteigeninfo!(obj::HKY85)

calculate eigeninfo for a model and stores it within the object.
"""
function seteigeninfo!(obj::HKY85)
    piA = obj.pi[1]; piC = obj.pi[2]; piG = obj.pi[3]; piT = obj.pi[4]
    piR = piA + piG
    piY = piT + piC

    a = obj.rate[1]
    if obj.relative
        b = 1.0
    else
        b = obj.rate[2]
    end
    obj.eigeninfo[1] = -b
    obj.eigeninfo[2] = -((piR * a) + (piY * b))
    obj.eigeninfo[3] = -((piY * a) + (piR * b))
    obj.eigeninfo[4] = piR
    obj.eigeninfo[5] = piY
end

function Base.show(io::IO, obj::HKY85)
    str = "HKY85 Substitution Model with rate=$(obj.rate),\n"
    if obj.relative
        str *= "relative rate version\n"
    else
        str *= "absolute rate version\n"
    end
    str *= "rate matrix Q:\n"
    print(io, str)
    showQ(io, obj)
end

"""
    nparams(model::JC69)

return number of parameters for JC69

"""
@inline function nparams(obj::JC69)
    (obj.relative ? 0 : 1)
end

"""
    nparams(model::HKY85)

Return number of parameters for HKY85
"""
@inline function nparams(obj::HKY85)
    (obj.relative ? 4 : 5)
end


"""
    Q(obj::JC69)
return Q rate matrix for the given model. Mutable version modeled after BioJulia/SubstitionModel.jl
and PyloModels.jl

Substitution rate matrix for a given substitution model:
Q[i,j] is the rate of transitioning from state i to state j.
"""

@inline function Q(obj::JC69)
    if obj.relative
        lambda = 1.0/3.0 #this allows branch lengths to be interpreted as substitutions/site (pi-weighted avg of diagonal is -1)
    else
        lambda = (1.0/3.0)*obj.rate[1]
    end

    return Qmatrix(-3*lambda, lambda, lambda, lambda,
            lambda, -3*lambda, lambda, lambda,
            lambda, lambda, -3*lambda, lambda,
            lambda, lambda, lambda, -3*lambda)
end

"""
    Q(obj::HKY85)
    return Q rate matrix. Mutable version modeled after BioJulia/SubstitionModel.jl
    and PhyloModels.jl
    
    Substitution rate matrix for a given substitution model:
    Q[i,j] is the rate of transitioning from state i to state j.
    For absolute, 
    average number of substitutions per site for time 1. lambda (or d) = (2*(piT*piC + piA*piG)*alpha + 2*(piY+piR)*beta)
    For relative,
    average number of substitutions per site is 1 (lambda = 1). kappa = alpha/beta. 
    ```julia-repl 
    julia> m1 = HKY85([.5], [0.25, 0.25, 0.25, 0.25])
    ````
"""
@inline function Q(obj::HKY85) #we end up not using this much -- we use this mostly to eamine our model
    piA = obj.pi[1]; piC = obj.pi[2]; piG = obj.pi[3]; piT = obj.pi[4]
    piR = piA + piG
    piY = piT + piC
    if obj.relative
        k = obj.rate[1]
        lambda = (2*(piT*piC + piA*piG)k + 2*(piY+piR)) #for sub/year interpretation

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
    else 
        a = obj.rate[1]; b = obj.rate[2]      
        Q₁  = a * piA
        Q₂  = a * piA
        Q₃  = b * piC
        Q₄  = a * piC
        Q₅  = a * piG
        Q₆  = b * piG
        Q₇  = b * piT
        Q₈  = a * piT
        Q₉  = -(Q₃ + Q₅ + Q₇)
        Q₁₀ = -(Q₁ + Q₆ + Q₈)
        Q₁₁ = -(Q₂ + Q₃ + Q₇)
        Q₁₂ = -(Q₁ + Q₄ + Q₆)
    end
    return Qmatrix(Q₉,  Q₁,  Q₂,  Q₁,
                    Q₃,  Q₁₀, Q₃,  Q₄,
                    Q₅,  Q₆,  Q₁₁, Q₆,
                    Q₇,  Q₈,  Q₇,  Q₁₂)

end

"""
    P(JC69, t)

Probability transition matrix for a [`NucleicAcidSubstitutionModel`](@ref), of the form

P[1,1] ... P[1,k]
   .          .
   .          .
P[k,1] ... P[k,k]

where P[i,j] is the probability of ending in state j after time t,
given that the process started in state i.
```julia-repl 
julia> m1 = JC69([0.25])
julia> P(m1, 3)
TODO
````
"""

@inline function P(obj::JC69, t::Float64)
    if t < 0
        error("Time must be positive")
    end
    P0 = 0.25 + 0.75 * exp(-t * obj.eigeninfo[1]) #lambda
    P1 = 0.25 - 0.25 * exp(-t * obj.eigeninfo[1])
    return Pmatrix(P0, P1, P1, P1, 
            P1, P0, P1, P1,
            P1, P1, P0, P1,
            P1, P1, P1, P0)
end

"""
    P(HKY85, t)

Probability transition matrix for a [`NucleicAcidSubstitutionModel`](@ref), of the form

P[1,1] ... P[1,k]
   .          .
   .          .
P[k,1] ... P[k,k]

where P[i,j] is the probability of ending in state j after time t,
given that the process started in state i.

rate = [k] (relative rate model)
or
rate = [a, b] (absolute rate model)
Both models according to Molecular Evolution (Yang 2014)

```julia-repl 
julia> m1 = HKY85([.5], [0.25, 0.25, 0.25, 0.25])
julia> P(m1, 3.0)
TODO
````
"""
@inline function P(obj::HKY85, t::Float64)
    if t < 0.0
        error("t must be positive")
    end
    piA = obj.pi[1]; piC = obj.pi[2]; piG = obj.pi[3]; piT = obj.pi[4]
    piR = obj.eigeninfo[4]
    piY = obj.eigeninfo[5]

    e2 = exp(obj.eigeninfo[1] * t)
    e3 = exp(obj.eigeninfo[2] * t)
    e4 = exp(obj.eigeninfo[3] * t)

    P_TT = piT + ((piT * piR)/piY) * e2 + (piC/piY) * e4
    P_TC = piC + ((piT * piR)/piY) * e2 - (piC/piY) * e4
    P_TA_CA = piA * (1 - e2) #double name because TA and CA are equal
    P_TG_CG = piG * (1 - e2) #double name because TG and CG are equal

    P_CT = piT + ((piT * piR)/piY) * e2 - (piT/piY) * e4
    P_CC = piC + ((piT * piR)/piY) * e2 + (piT/piY) * e4

    P_AT_GT = piT * (1 - e2) #double name because AT and GT are equal
    P_AC_GC = piC * (1 - e2) #double name because AC and GC are equal
    P_AA = piA + ((piA * piY)/piR) * e2 + (piG/piR) * e3
    P_AG = piG + ((piG * piY)/piR) * e2 - (piG/piR) * e3

    P_GA = piA + ((piA * piY)/piR) * e2 - (piA/piR) * e3
    P_GG = piG + ((piG * piY)/piR) * e2 + (piA/piR) * e3

    return Pmatrix(P_TT, P_TC, P_TA_CA, P_TG_CG,
                    P_CT, P_CC, P_TA_CA, P_TG_CG,
                    P_AT_GT, P_AC_GC, P_AA, P_AG,
                    P_AT_GT, P_AC_GC, P_GA, P_GG)
end

"""
    P!(Pmat::AbstractMatrix, obj::JC69, t::Float64)
modifies P rate matrix (see traitsLikDiscrete.jl) during optimization. 

```julia-repl 
julia> m1 = JC69([0.25])
julia> P!(P(m1, 1.0), m1, 3.0)
````
"""
function P!(Pmat::AbstractMatrix, obj::JC69, t::Float64)
    # if t < 0 
    #This requirement should be handled by the optimization function 
    #   and when we fit it to a tree we'll check that all branch lengths are positive
    #     error("Time must be positive")
    # end
    P0 = 0.25 + 0.75 * exp(-t * obj.eigeninfo[1]) #lambda
    P1 = 0.25 - 0.25 * exp(-t * obj.eigeninfo[1])
    Pmat[:,:] = P1 #puts P1 on the off-diagonals
    for i in 1:4 Pmat[i,i]=P0;end #diagonal
    return Pmat
end

"""
    P!(Pmat::AbstractMatrix, obj::HKY85, t::Float64)
modifies P rate matrix (see traitsLikeDiscrete)
The model will have optimized the rate using NLopt in loglikfun 
in traitsLikDiscrete.jl.
```julia-repl 
julia> m1 = HKY85([.5], [0.25, 0.25, 0.25, 0.25])
julia> P!(P(m1, 1.0), m1, 3.0)
````
"""
function P!(Pmat::AbstractMatrix, obj::HKY85, t::Float64)
    piA = obj.pi[1]; piC = obj.pi[2]; piG = obj.pi[3]; piT = obj.pi[4]
    piR = obj.eigeninfo[4]
    piY = obj.eigeninfo[5]

    e2 = exp(obj.eigeninfo[1] * t)
    e3 = exp(obj.eigeninfo[2] * t)
    e4 = exp(obj.eigeninfo[3] * t)

    P_TT = piT + ((piT * piR)/piY) * e2 + (piC/piY) * e4
    P_TC = piC + ((piT * piR)/piY) * e2 - (piC/piY) * e4
    P_TA_CA = piA * (1 - e2) #double name because TA and CA are equal
    P_TG_CG = piG * (1 - e2) #double name because TG and CG are equal

    P_CT = piT + ((piT * piR)/piY) * e2 - (piT/piY) * e4
    P_CC = piC + ((piT * piR)/piY) * e2 + (piT/piY) * e4

    P_AT_GT = piT * (1 - e2) #double name because AT and GT are equal
    P_AC_GC = piC * (1 - e2) #double name because AC and GC are equal
    P_AA = piA + ((piA * piY)/piR) * e2 + (piG/piR) * e3
    P_AG = piG + ((piG * piY)/piR) * e2 - (piG/piR) * e3

    P_GA = piA + ((piA * piY)/piR) * e2 - (piA/piR) * e3
    P_GG = piG + ((piG * piY)/piR) * e2 + (piA/piR) * e3

    Pmat[1,1] = P_TT
    Pmat[2,2] = P_CC
    Pmat[3,3] = P_AA
    Pmat[4,4] = P_GG
    Pmat[1,2] = P_TC
    Pmat[[1,2],3] = P_TA_CA
    Pmat[[1,2],4] = P_TG_CG
    Pmat[2,1] = P_CT
    Pmat[[3,4],1] = P_AT_GT
    Pmat[[3,4],2] = P_AC_GC
    P[3,4] = P_AG
    P[4,3] = P_GA
    return Pmat
end

"""
    RateVariationAcrossSites ()

[`RateVariationAcrossSites`](@ref) Allow variable substitution rates across sites using the discrete gamma model
(Yang 1994, Journal of Molecular Evolution). Turn any NASM to NASM + gamma. 

Using this model increases nparams by one.

Because mean(gamma) should equal 1, alpha = beta. Refer to this parameter as alpha here.

```julia-repl
julia> rv = RateVariationAcrossSites()
julia> show(rv)
#TODO
julia> setalpha(rv)
#TODO
```
"""
mutable struct RateVariationAcrossSites
    alpha::Float64
    ncat::Int
    ratemultiplier::Array{Float64}
    function RateVariationAcrossSites(alpha = 1.0::Float64, ncat = 4::Int)
        @assert alpha >= 0 "alpha must be >= 0"
        if ncat == 1
            ratemultiplier = [1.0]
        else
            cuts = (0:(ncat-1))/ncat + 1/2ncat
            ratemultiplier = quantile.(Distributions.Gamma(alpha, alpha), cuts)
        end
        new(alpha, ncat, ratemultiplier)
    end
end
const RVAS = RateVariationAcrossSites

"""
    setalpha!(obj, alpha)
Set alpha in RateVariationAcrossSites model
"""
function setalpha!(obj::RateVariationAcrossSites, alpha::Int)
    @assert alpha >= 0 "alpha must be >= 0"
    obj.alpha = alpha
    cuts = (0:(obj.ncat-1))/obj.ncat + 1/2obj.ncat
    obj.ratemultiplier[:] = quantile.(Distributions.Gamma(alpha, alpha), cuts)
end

function Base.show(io::IO, obj::RateVariationAcrossSites)
    str = "$(obj) Rate Variation Across Sites using Discretized Gamma Model\n"
    str *= "alpha: $(obj.alpha)\n"
    str *= "categories for Gamma discretization: $(obj.ncat)\n"
    str *= "ratemultiplier: $(obj.ratemultiplier)\n"
    print(io, str)
    showQ(io, obj)
end