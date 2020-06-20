# Tests to simulate multivariate traits

## Get an ultrametric network ('net' already in global scope)
net = readTopology("(((Ag:5,(#H1:1::0.056,((Ak:2,(E:1,#H2:1::0.004):1):1,(M:2)#H2:1::0.996):1):1):1,(((((Az:1,Ag2:1):1,As:2):1)#H1:1::0.944,Ap:4):1,Ar:5):1):1,(P:4,20:4):3,165:7);");

@testset "Simulate data and check means and dimensions" begin

# one single trait
pars = ParamsMultiBM([5.0], ones(1,1))
sim = simulate(net, pars)
@test sim[:Tips, :Exp] == [5.0 5.0 5.0 5.0 5.0 5.0 5.0 5.0 5.0 5.0 5.0 5.0]

trait_dim = 3

## Simulate a MBD
Random.seed!(17920921); # fix the seed

μ = randn(trait_dim)
Σ = randn(trait_dim, trait_dim)
Σ = Σ * Σ' # needs to be positive definite

@test_throws ErrorException ParamsMultiBM(μ[1:2], Σ)

pars = ParamsMultiBM(μ, Σ); # params of a MBD
@test_logs show(devnull, pars)

sim = simulate(net, pars); # simulate according to a BM
@test_logs show(devnull, sim)
@test_throws ErrorException sim[:Tips, :Broken]

# Extract simulated values
traitsTips = sim[:Tips];
traitsNodes = sim[:InternalNodes];

# Check dimensions
@test size(traitsTips) == (trait_dim, net.numTaxa)
@test size(traitsNodes) == (trait_dim, net.numNodes - net.numTaxa)

# Check means (no shifts)
@test sim[:All, :Exp] ≈ μ * ones(net.numNodes)'

end

###############################################################################
## Test of distibution (fixed root)
###############################################################################
@testset "Simulate test distribution (fixed root)" begin
trait_dim = 3

## Generate some values
Random.seed!(18480224); # fix the seed
μ = randn(trait_dim)
Σ = randn(trait_dim, trait_dim)
Σ = Σ * Σ' # needs to be positive definite
pars = ParamsMultiBM(μ, Σ); # params of a MBD

N = 10000
S = length(tipLabels(net));
μ_sim = zeros(trait_dim, S)
Σ_sim = zeros(trait_dim * S, trait_dim * S)
for i = 1:N
    tips = simulate(net, pars)[:Tips]
    μ_sim .+= tips
    v_sim = vec(tips)
    Σ_sim .+= v_sim * v_sim'
end
μ_sim ./= N
Σ_sim ./= N
Σ_sim = Σ_sim - vec(μ_sim) * vec(μ_sim)'

## Check means
μ_true = μ * ones(S)'
@test isapprox(μ_sim, μ_true, atol=0.2)

## Check covariance
Ψ = Matrix(vcv(net))
Σ_true = kron(Ψ, Σ)
@test isapprox(Σ_sim, Σ_true, atol=3.0) # norm L2 of 36x36 matrix

end

###############################################################################
## Test of distibution (random root)
###############################################################################
@testset "Simulate test distribution (random root)" begin
trait_dim = 3

## Generate some values
Random.seed!(24452384); # fix the seed
μ = randn(trait_dim)
Σ = randn(trait_dim, trait_dim)
Σ = Σ * Σ' # needs to be positive definite
pars = ParamsMultiBM(μ, Σ); # params of a MBD

Σ_root = randn(trait_dim, trait_dim)
Σ_root = Σ_root * Σ_root' / 10.0 # needs to be positive definite (and not too big to reduce variance in test)
pars.varRoot = Σ_root
pars.randomRoot = true

@test_logs show(devnull, pars)
show(devnull, pars)

N = 10000
S = length(tipLabels(net));
μ_sim = zeros(trait_dim, S)
Σ_sim = zeros(trait_dim * S, trait_dim * S)
for i = 1:N
    tips = simulate(net, pars)[:Tips]
    μ_sim .+= tips
    v_sim = vec(tips)
    Σ_sim .+= v_sim * v_sim'
end
μ_sim ./= N
Σ_sim ./= N
Σ_sim = Σ_sim - vec(μ_sim) * vec(μ_sim)'

## Check means
@test isapprox(μ_sim, μ * ones(S)', atol=0.2)

## Check covariance
Ψ = Matrix(vcv(net))
Σ_true = kron(Ψ, Σ) + kron(ones(S, S), pars.varRoot)
@test isapprox(Σ_sim, Σ_true, atol=3.6) # norm L2 of 36x36 matrix

end


################################################################################
## With shifts
################################################################################

net = readTopology("(A:2.5,((B:1,#H1:0.5::0.4):1,(C:1,(D:0.5)#H1:0.5::0.6):1):0.5);")

@testset "Simulate with Shifts" begin
trait_dim = 3

Random.seed!(275698234545); # fix the seed

## Test construction function
@test_throws ErrorException ShiftNet(net.edge[7], [1.0, 2.0],  net) # can't put a shift on hybrid branch
@test_throws ErrorException ShiftNet(net.node[6], [1.0, 2.0],  net) # can't put a shift on hybrid branch
@test_throws ErrorException ShiftNet([net.node[7], net.node[6]], [1.0 2.0], net) # dimensions don't match
@test ShiftNet(net.edge[8], [1.0, 2.0],  net).shift ≈ ShiftNet([net.edge[8]], [1.0 2.0],  net).shift
@test ShiftNet(net.edge[8], [1.0, 2.0],  net).shift ≈ ShiftNet(net.node[7], [1.0, 2.0],  net).shift
@test ShiftNet(net.node[7], [1.0, 2.0],  net).shift ≈ ShiftNet([net.node[7]], [1.0 2.0],  net).shift

## Concatenate function
sh1 = ShiftNet(net.node[7], [1.0, 2.0],  net)*ShiftNet(net.node[9], [3.0, -1.5],  net)
@test_logs show(devnull, sh1)
@test sh1.shift ≈ ShiftNet([net.node[7], net.node[9]], [1.0 2.0; 3.0 -1.5],  net).shift
@test_throws ErrorException sh1*ShiftNet(net.edge[7], [4.0, 3.5],  net) # can't concatenate if the two affect the same edges
@test sh1.shift ≈ (sh1*ShiftNet([net.node[7]], [1.0 2.0],  net)).shift
@test_throws ErrorException sh1*ShiftNet(net.edge[8], [4.0, 3.5, 5.0],  net) # can't concatenate if the two affect the same edges

## Values and edge numbers functions
@test getShiftEdgeNumber(sh1) == [-1, 8]
@test getShiftValue(sh1) == [3.0 -1.5; 1.0 2.0]

## Hybrid shifts
@test shiftHybrid([4.5 2.0], net).shift ≈ ShiftNet(net.edge[6], [4.5, 2.0], net).shift
@test_throws ErrorException shiftHybrid([4.5 2.0; 3.0 5.0], net) # dimension mismatch

## Distributions
μ = randn(trait_dim)
Σ = randn(trait_dim, trait_dim)
Σ = Σ * Σ' # needs to be positive definite

@test ParamsMultiBM(μ, Σ, net).shift.shift ≈ ParamsMultiBM(μ, Σ, ShiftNet(net, trait_dim)).shift.shift
@test_throws ErrorException ParamsMultiBM(μ, Σ, ShiftNet(net, 1))
# below: the shift at the root, net.node[9], is not used by the simulation
# It doesn't even affect the mean μ
sh = ShiftNet(net.node[7], [1.0, 2.0, -1.0],  net)*ShiftNet(net.node[9], [3.0, -1.5, 4.2],  net)

pars = ParamsMultiBM(μ, Σ, sh)

@test_logs show(devnull, pars)
show(devnull, pars)

N = 10000
S = length(tipLabels(net));
μ_sim = zeros(trait_dim, S)
Σ_sim = zeros(trait_dim * S, trait_dim * S)
for i = 1:N
    tips = simulate(net, pars)[:Tips]
    μ_sim .+= tips
    v_sim = vec(tips)
    Σ_sim .+= v_sim * v_sim'
end
μ_sim ./= N
Σ_sim ./= N
Σ_sim = Σ_sim - vec(μ_sim) * vec(μ_sim)'

## Check means
sim = simulate(net, pars)
μ_true = sim[:Tips, :Exp]
@test μ_true ≈ μ .+ [0 0 1 0.6; 0 0 2 1.2; 0 0 -1 -.6] # from shift on edge 8 only
@test isapprox(μ_sim, μ_true, atol=0.2)

## Check covariance
Ψ = Matrix(vcv(net))
Σ_true = kron(Ψ, Σ)
@test isapprox(Σ_sim, Σ_true, norm = x -> maximum(abs.(x)), atol=0.4)
end
