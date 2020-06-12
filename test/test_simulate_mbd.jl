# Tests to simulate multivariate traits

using PhyloNetworks, LinearAlgebra, Test
import Random

## Get a network and make it ultrametric
net = readTopology("(((Ag:5,(#H1:1::0.056,((Ak:2,(E:1,#H2:1::0.004):1):1,(M:2)#H2:1::0.996):1):1):1,(((((Az:1,Ag2:1):1,As:2):1)#H1:1::0.944,Ap:4):1,Ar:5):1):1,(P:4,20:4):3,165:7);");
trait_dim = 3


@testset "Simulate data and check means and dimensions" begin

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

# Extract simulated values
traitsTips = sim[:Tips];
traitsNodes = sim[:InternalNodes];

# Check dimensions
@test size(traitsTips) == (trait_dim, net.numTaxa)
@test size(traitsNodes) == (trait_dim, net.numNodes - net.numTaxa)

expectations = sim[:All, :Exp]

# Check means (no shifts)
max_diff = maximum(expectations - μ * ones(net.numNodes)')
@test max_diff ≈ 0.0 atol=1e-10

end

###############################################################################
## Test of distibution (fixed root)
###############################################################################
@testset "Simulate test distribution (fixed root)" begin

## Generate some values
Random.seed!(18480224); # fix the seed
μ = randn(trait_dim)
Σ = randn(trait_dim, trait_dim)
Σ = Σ * Σ' # needs to be positive definite
pars = ParamsMultiBM(μ, Σ); # params of a MBD

N = 50000
S = length(tipLabels(net));
μ_sim = zeros(trait_dim, S)
Σ_sim = zeros(trait_dim * S, trait_dim * S)
for i = 1:N
    tips = simulate(net, pars)[:Tips]
    μ_sim .+= tips
    v_sim = vec(tips)
    Σ_sim += v_sim * v_sim'
end

μ_sim ./= N
Σ_sim ./= N
Σ_sim = Σ_sim - vec(μ_sim) * vec(μ_sim)'

## Check means

μ_true = μ * ones(S)'

μ_max = maximum(abs.(μ_true - μ_sim))
@test μ_max < 1e-1

## Check covariance

Ψ = Matrix(vcv(net))
Σ_true = kron(Ψ, Σ)
Σ_max = maximum(abs.(Σ_true - Σ_sim))
@test Σ_max < 2e-1


end

###############################################################################
## Test of distibution (random root)
###############################################################################
@testset "Simulate test distribution (random root)" begin

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

N = 50000
S = length(tipLabels(net));
μ_sim = zeros(trait_dim, S)
Σ_sim = zeros(trait_dim * S, trait_dim * S)
for i = 1:N
    tips = simulate(net, pars)[:Tips]
    μ_sim .+= tips
    v_sim = vec(tips)
    Σ_sim += v_sim * v_sim'
end

μ_sim ./= N
Σ_sim ./= N
Σ_sim = Σ_sim - vec(μ_sim) * vec(μ_sim)'

## Check means

μ_true = μ * ones(S)'

μ_max = maximum(abs.(μ_true - μ_sim))
@test μ_max < 1e-1

## Check covariance

Ψ = Matrix(vcv(net))
Σ_true = kron(Ψ, Σ) + kron(ones(S, S), pars.varRoot)
Σ_max = maximum(abs.(Σ_true - Σ_sim))
@show Σ_max
@test Σ_max < 2e-1


end


################################################################################
## With shifts
################################################################################

@testset "Simulate with Shifts" begin

Random.seed!(275698234545); # fix the seed

global net
net = readTopology("(A:2.5,((B:1,#H1:0.5::0.4):1,(C:1,(D:0.5)#H1:0.5::0.6):1):0.5);")

## Test construction function
@test_throws ErrorException ShiftNet(net.edge[7], [1.0, 2.0],  net) # can't put a shift on hybrid branch
@test_throws ErrorException ShiftNet(net.node[6], [1.0, 2.0],  net) # can't put a shift on hybrid branch
@test ShiftNet(net.edge[8], [1.0, 2.0],  net).shift ≈ ShiftNet([net.edge[8]], [1.0 2.0],  net).shift
@test ShiftNet(net.edge[8], [1.0, 2.0],  net).shift ≈ ShiftNet(net.node[7], [1.0, 2.0],  net).shift
@test ShiftNet(net.node[7], [1.0, 2.0],  net).shift ≈ ShiftNet([net.node[7]], [1.0 2.0],  net).shift

## Concatenate function
sh1 = ShiftNet(net.node[7], [1.0, 2.0],  net)*ShiftNet(net.node[9], [3.0, -1.5],  net)
@test sh1.shift ≈ ShiftNet([net.node[7], net.node[9]], [1.0 2.0; 3.0 -1.5],  net).shift
@test_throws ErrorException sh1*ShiftNet(net.edge[7], [4.0, 3.5],  net) # can't concatenate if the two affect the same edges
@test sh1.shift ≈ (sh1*ShiftNet([net.node[7]], [1.0 2.0],  net)).shift
@test_throws ErrorException sh1*ShiftNet(net.edge[8], [4.0, 3.5, 5.0],  net) # can't concatenate if the two affect the same edges


## Values and edge numbers functions
@test getShiftEdgeNumber(sh1) == [-1, 8]
@test getShiftValue(sh1) == [3.0 -1.5; 1.0 2.0]

## Hybrid shifts
@test shiftHybrid([4.5 2.0], net).shift ≈ ShiftNet(net.edge[6], [4.5, 2.0], net).shift

## Distributions

μ = randn(trait_dim)
Σ = randn(trait_dim, trait_dim)
Σ = Σ * Σ' # needs to be positive definite

@test ParamsMultiBM(μ, Σ, net).shift.shift ≈ ParamsMultiBM(μ, Σ, ShiftNet(net, trait_dim)).shift.shift
@test_throws ErrorException ParamsMultiBM(μ, Σ, ShiftNet(net, 1))

pars = ParamsMultiBM(μ, Σ, net)

N = 50000
S = length(tipLabels(net));
μ_sim = zeros(trait_dim, S)
Σ_sim = zeros(trait_dim * S, trait_dim * S)
for i = 1:N
    tips = simulate(net, pars)[:Tips]
    μ_sim .+= tips
    v_sim = vec(tips)
    Σ_sim += v_sim * v_sim'
end

μ_sim ./= N
Σ_sim ./= N
Σ_sim = Σ_sim - vec(μ_sim) * vec(μ_sim)'

## Check means
sim = simulate(net, pars)
μ_true = sim[:Tips, :Exp]

μ_max = maximum(abs.(μ_true - μ_sim))
@show μ_max
@test μ_max < 1e-1

## Check covariance

Ψ = Matrix(vcv(net))
Σ_true = kron(Ψ, Σ)
Σ_max = maximum(abs.(Σ_true - Σ_sim))
@test Σ_max < 2e-1


end
