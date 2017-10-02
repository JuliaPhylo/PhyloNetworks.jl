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
  
function show(io::IO, object::BinaryTraitSubstitutionModel)
	print(io, "Binary Trait Substitution Model: \n rate 0→ 1 α=$(object.α) rate 1→ 0 β =$(object.β) ") #fixit separate into multiple strings; print all at end; * = concatenate
end

# fixit: function Q, function P

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


#fixit generic Q function

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
    randomTrait(model, t, start)
simulate traits along one edge of length t.  
`start` must be a vector of integers, each representing the starting value of one trait.
"""
function randomTrait(mod::SM, t::Float64, start::Vector{Int})
    Pt = P(mod, t)
    w = [aweights(Pt[1,:]), aweights(Pt[2,:])]
    res = Int[]
    k = size(Pt, 1) # number of states
    for i in 1:length(start)
        push!(res, sample(1:k, w[start[i]]))
    end
    return res
end


