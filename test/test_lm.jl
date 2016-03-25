## Test of PhyloNetworklm

using PhyloNetworks
using GLM
include("../src/types.jl")
include("../src/functions.jl")
include("../src/traits.jl")

tree= "(A,((B,#H1),(C,(D)#H1)));"
net=readTopologyLevel1(tree)
printEdges(net)

# Re-root the tree so that it matches my example
root!(net, "A")
printEdges(net)
# directEdges!(net) ## I am forced to do thi step here, because root sends a net with net.isRooted = false. Expected behavior ? 
preorder!(net)
plot(net, useEdgeLength = true,  showEdgeNumber=true)

# Make the network ultrametric
net.edge[1].length = 2.5
net.edge[6].length = 0.5
net.edge[7].length = 0.5
net.edge[3].length = 0.5
plot(net, useEdgeLength = true)
# Rk: Is there a way to check that the branch length are coherents with 
# one another (Especialy for hybrids) ?

# Ancestral state reconstruction with ready-made matrices
params = paramsBM(10, 1)
sim = simulate(net, params)
Y = sim[:Tips]
X = ones(4, 1)
fit = phyloNetworklm(Y, X, net)

fitNaive = phyloNetworklmNaive(Y, X, net) # old naive version

# Simulate correlated data in data frames
b0 = 1
b1 = 10
sim = simulate(net, paramsBM(1, 1))
A = sim[:Tips]
B = b0 + b1 * A + simulate(net,  paramsBM(0, 0.1))[:Tips]
# perfect user using right format and formula
df = DataFrame(trait = B, pred = A, tipsNames = sim.M.tipsNames)
fit = phyloNetworklm(trait ~ pred, df, net)

# With Matrices
X = hcat(ones(4), A)
fit_mat = phyloNetworklm(B, X, net)

# unordered data
df = df[[2, 1, 4, 3], :]
df
fitbis = phyloNetworklm(trait ~ pred, df, net)

# unnamed ordered data
df = DataFrame(trait = B, pred = A)
fitter = phyloNetworklm(trait ~ pred, df, net)

# unnamed un-ordered data
df = df[[2, 1, 4, 3], :]
fitter = phyloNetworklm(trait ~ pred, df, net) # Wrong pred


fit
loglikelihood(fit)

### Add NAs
b0 = 1
b1 = 10
sim = simulate(net, paramsBM(1, 1))
A = sim[:Tips]
B = b0 + b1 * A + simulate(net,  paramsBM(0, 0.1))[:Tips]
df = DataFrame(trait = B, pred = A, tipsNames = sim.M.tipsNames)
df[2, :pred] = NA

fit = phyloNetworklm(trait ~ pred, df, net)
