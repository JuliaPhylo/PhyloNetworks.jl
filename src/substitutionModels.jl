abstract type TraitSubstitutionModel end
const SM = TraitSubstitutionModel
const Bmatrix = SMatrix{2, 2, Float64}

"""
`BinaryTraitSubstitutionModel` is an abstract type that contains all models
describing a substitution process impacting biological characters with binary states
with continous time Markov models.
"""

mutable struct BinaryTraitSubstitutionModel <: TraitSubstitutionModel
    α::Float64
    β::Float64
    π0::Float64
    π1::Float64
    label::SVector{2, String}
    function BinaryTraitSubstitutionModel(α::Float64, β::Float64, label::SVector{2, String})
	    α >= 0. || error("parameter α must be non-negative")
	    β >= 0. || error("parameter β must be non-negative")
	    ab = α+β
        ab > 0. || error("α+β must be positive")
        new(α, β, β/ab, α/ab, label)
    end
end

BinaryTraitSubstitutionModel(α, β) = BinaryTraitSubstitutionModel(α, β, SVector("0", "1"))

"""
    TwoBinaryTraitSubstitutionModel(α [, label])

    α1 = rate 0->1
    α2 = rate 
"""

mutable struct TwoBinaryTraitSubstitutionModel <: TraitSubstitutionModel
    α::Vector{Float64}
    label::Vector{String}
    function TwoBinaryTraitSubstitutionModel(α::AbstractVector{Float64}, label::AbstractVector{String})
	    all( x -> x >= 0., α) || error("rates in α must be non-negative")
        new(α, label)
    end
end

TwoBinaryTraitSubstitutionModel(α) = TwoBinaryTraitSubstitutionModel(α, ["0", "1", "0", "1"])

function show(io::IO, object::TwoBinaryTraitSubstitutionModel)
    R"""
    signif<-3
    par(mfrow=c(2,1))
    plot.new()
    par(mar=c(1.1,2.1,3.1,2.1))
    plot.window(xlim=c(0,2),ylim=c(0,1),asp=1)
    """
    R"""
    mtext("Two Binary Trait Substitution Model",side=3,adj=0,line=1.2,cex=1.2)
    arrows(x0=0.15,y0=0.15,y1=0.85,lwd=2,length=0.1)
    arrows(x0=0.2,y0=0.85,y1=0.15,lwd=2,length=0.1)
    arrows(x0=1.6,y0=0.05,x1=0.4,lwd=2,length=0.1)
    arrows(x0=0.4,y0=0.1,x1=1.6,lwd=2,length=0.1)
    arrows(x0=1.8,y0=0.15,y1=0.85,lwd=2,length=0.1)
    arrows(x0=1.85,y0=0.85,y1=0.15,lwd=2,length=0.1)
    arrows(x0=1.6,y0=0.9,x1=0.4,lwd=2,length=0.1)
    arrows(x0=0.4,y0=0.95,x1=1.6,lwd=2,length=0.1)
    text(x=0.175,y=0.95,paste($(object.label[1]), ",", $(object.label[1])))
    text(x=1.825,y=0.95,paste($(object.label[1]), ",", $(object.label[2])))
    text(x=1.825,y=0.05,paste($(object.label[2]), ",", $(object.label[2])))
    text(x=0.175,y=0.05,paste($(object.label[2]), ",", $(object.label[1])))
    """
    R"""
    text(x=1,y=1,round($(object.α[1]),signif),cex=0.8)
    """
    R"""
    text(x=1,y=0.85,round($(object.α[2]),signif),cex=0.8)
    text(x=1.9,y=0.5,round($(object.α[3]),signif),cex=0.8,srt=90)
    text(x=1.75,y=0.5,round($(object.α[4]),signif),cex=0.8,srt=90)
    """
    R"""
    text(x=1,y=0,round($(object.α[5]),signif),cex=0.8)
    text(x=1,y=0.15,round($(object.α[6]),signif),cex=0.8)
    text(x=0.1,y=0.5,round($(object.α[7]),signif),cex=0.8,srt=90)
    text(x=0.25,y=0.5,round($(object.α[8]),signif),cex=0.8,srt=90)
    """
end

function Q(mod::TwoBinaryTraitSubstitutionModel)
    M = fill(0.0,(4,4))
    a = mod.α
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
        
const BTSM = BinaryTraitSubstitutionModel

const TBTSM = TwoBinaryTraitSubstitutionModel

"""
    nStates(mod)

Show number of character states in a given model.

# Examples

```julia-repl
julia> m1 = BinaryTraitSubstitutionModel(1.0, 2.0)
julia> nStates(m1)
 2
```
"""

function nStates(mod::BTSM)
    return 2::Int
end

function show(io::IO, object::BinaryTraitSubstitutionModel)
    str = "Binary Trait Substitution Model:\n"
    str *= "rate $(object.label[1])→$(object.label[2]) α=$(object.α)\n"
    str *= "rate $(object.label[2])→$(object.label[1]) β=$(object.β)\n"
    print(io, str)
end

"""
`EqualRatesSubstitutionModel` is an abstract type that contains all models
describing a substitution process impacting biological characters with equal rates
of transition between all character states with continous time Markov models.
"""

mutable struct EqualRatesSubstitutionModel <: TraitSubstitutionModel
    k::Int
    α::Float64
    label::Vector{String}
    function EqualRatesSubstitutionModel(k::Int, α::Float64, label::Vector{String})
        k >= 2 || error("parameter k must be greater than or equal to 2")
        α > 0 || error("parameter α must be positive")
        new(k, α, label)
    end
end

function show(io::IO, object::EqualRatesSubstitutionModel)
    str = "Equal Rates Substitution Model:\n"
    str *= "all rates α=$(object.α)\n"
    str *= "number of states, k=$(object.k)\n"  
    print(io, str)
    M = fill(object.α, object.k, object.k)
    for i = 1:size(M,2)
        pad = 8
        if object.label != "" && i==1
            pad = 2*8
        end    
        @printf("%s", lpad(object.label[i],pad," "))
    end
    @printf("\n")
    # print the rows
    for i = 1:size(M,1)
        if object.label != ""
            @printf("%s", lpad(object.label[i],8," "))
        end
        for j = 1:size(M,2)
            # TBD: use fmt defined above to print array contents
            @printf("%8.4f",(M[i,j]))
        end
        @printf("\n")
    end  
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