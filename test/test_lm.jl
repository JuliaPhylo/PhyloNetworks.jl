## Test of PhyloNetworklm
## still not an automatic test function, needs work

using PhyloNetworks
using GLM
using DataFrames
using Base.Test
#include("../src/types.jl")
#include("../src/functions.jl")
#include("../src/traits.jl")

tree= "(A,((B,#H1),(C,(D)#H1)));"
net=readTopologyLevel1(tree)
printEdges(net)

# Re-root the tree so that it matches my example
rootatnode!(net, "A")
printEdges(net)
preorder!(net)
# plot(net, useEdgeLength = true,  showEdgeNumber=true)

# Make the network ultrametric
net.edge[1].length = 2.5
net.edge[6].length = 0.5
net.edge[7].length = 0.5
net.edge[3].length = 0.5
# plot(net, useEdgeLength = true)
# Rk: Is there a way to check that the branch length are coherents with 
# one another (Especialy for hybrids) ?

# Ancestral state reconstruction with ready-made matrices
params = paramsBM(10, 1)
sim = simulate(net, params)
Y = sim[:Tips]
X = ones(4, 1)
fit = phyloNetworklm(X, Y, net)
# fitNaive = phyloNetworklmNaive(X, Y, net) # old naive version
# @test_approx_eq coef(fit) coef(fitNaive)
# @test_approx_eq vcov(fit) vcov(fitNaive)
# @test_approx_eq nobs(fit) nobs(fitNaive)
# @test_approx_eq residuals(fit) residuals(fitNaive)
# @test_approx_eq model_response(fit) model_response(fitNaive)
# @test_approx_eq predict(fit) predict(fitNaive)
# @test_approx_eq df_residual(fit) df_residual(fitNaive)
# @test_approx_eq sigma2_estim(fit) sigma2_estim(fitNaive)
# @test_approx_eq stderr(fit) stderr(fitNaive)
# @test_approx_eq confint(fit) confint(fitNaive)
# @test_approx_eq loglikelihood(fit) loglikelihood(fitNaive)

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
fit_mat = phyloNetworklm(X, B, net)

@test_approx_eq coef(fit) coef(fit_mat)
@test_approx_eq vcov(fit) vcov(fit_mat)
@test_approx_eq nobs(fit) nobs(fit_mat)
@test_approx_eq residuals(fit) residuals(fit_mat)
@test_approx_eq model_response(fit) model_response(fit_mat)
@test_approx_eq predict(fit) predict(fit_mat)
@test_approx_eq df_residual(fit) df_residual(fit_mat)
@test_approx_eq sigma2_estim(fit) sigma2_estim(fit_mat)
@test_approx_eq stderr(fit) stderr(fit_mat)
@test_approx_eq confint(fit) confint(fit_mat)
@test_approx_eq loglikelihood(fit) loglikelihood(fit_mat)


# unordered data
df = df[[2, 1, 4, 3], :]
fitbis = phyloNetworklm(trait ~ pred, df, net)

@test_approx_eq coef(fit) coef(fitbis)
@test_approx_eq vcov(fit) vcov(fitbis)
@test_approx_eq nobs(fit) nobs(fitbis)
@test_approx_eq residuals(fit)[fitbis.ind] residuals(fitbis)
@test_approx_eq model_response(fit)[fitbis.ind] model_response(fitbis)
@test_approx_eq predict(fit)[fitbis.ind] predict(fitbis)
@test_approx_eq df_residual(fit) df_residual(fitbis)
@test_approx_eq sigma2_estim(fit) sigma2_estim(fitbis)
@test_approx_eq stderr(fit) stderr(fitbis)
@test_approx_eq confint(fit) confint(fitbis)
@test_approx_eq loglikelihood(fit) loglikelihood(fitbis)

# unnamed ordered data
df = DataFrame(trait = B, pred = A)
fitter = phyloNetworklm(trait ~ pred, df, net)

@test_approx_eq coef(fit) coef(fitter)
@test_approx_eq vcov(fit) vcov(fitter)
@test_approx_eq nobs(fit) nobs(fitter)
@test_approx_eq residuals(fit) residuals(fitter)
@test_approx_eq model_response(fit) model_response(fitter)
@test_approx_eq predict(fit) predict(fitter)
@test_approx_eq df_residual(fit) df_residual(fitter)
@test_approx_eq sigma2_estim(fit) sigma2_estim(fitter)
@test_approx_eq stderr(fit) stderr(fitter)
@test_approx_eq confint(fit) confint(fitter)
@test_approx_eq loglikelihood(fit) loglikelihood(fitter)

# unnamed un-ordered data
df = df[[2, 1, 4, 3], :]
fitter = phyloNetworklm(trait ~ pred, df, net) # Wrong pred


### Add NAs
df = DataFrame(trait = B, pred = A, tipsNames = tipLabels(sim))
df[2, :pred] = NA
fitna = phyloNetworklm(trait ~ pred, df, net)

df = df[[2, 1, 4, 3], :]
fitnabis = phyloNetworklm(trait ~ pred, df, net)

@test_approx_eq coef(fitna) coef(fitnabis)
@test_approx_eq vcov(fitna) vcov(fitnabis)
@test_approx_eq nobs(fitna) nobs(fitnabis)
@test_approx_eq sort(residuals(fitna)) sort(residuals(fitnabis))
@test_approx_eq sort(model_response(fitna)) sort(model_response(fitnabis))
@test_approx_eq sort(predict(fitna)) sort(predict(fitnabis))
@test_approx_eq df_residual(fitna) df_residual(fitnabis)
@test_approx_eq sigma2_estim(fitna) sigma2_estim(fitnabis)
@test_approx_eq stderr(fitna) stderr(fitnabis)
@test_approx_eq confint(fitna) confint(fitnabis)
@test_approx_eq loglikelihood(fitna) loglikelihood(fitnabis)


### BLUP
params = paramsBM(3, 1)
sim = simulate(net, params)
Y = sim[:Tips]
ancestral_traits = ancestralStateReconstruction(net, Y, params)

