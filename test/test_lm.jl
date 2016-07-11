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
#printEdges(net)

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
phynetlm = phyloNetworklm(X, Y, net)
@show phynetlm
# Naive version (GLS)
ntaxa = length(Y)
Vy = phynetlm.Vy
Vyinv = inv(Vy)
XtVyinv = X' * Vyinv
logdetVy = logdet(Vy)
betahat = inv(XtVyinv * X) * XtVyinv * Y
fittedValues =  X * betahat
resids = Y - fittedValues
sigma2hat = 1/ntaxa * (resids' * Vyinv * resids)
 # log likelihood
loglik = - 1 / 2 * (ntaxa + ntaxa * log(2 * pi) + ntaxa * log(sigma2hat) + logdetVy)
# null version
nullX = ones(ntaxa, 1)
nullXtVyinv = nullX' * Vyinv
nullresids = Y - nullX * inv(nullXtVyinv * nullX) * nullXtVyinv * Y 
nullsigma2hat = 1/ntaxa * (nullresids' * Vyinv * nullresids)
nullloglik = - 1 / 2 * (ntaxa + ntaxa * log(2 * pi) + ntaxa * log(nullsigma2hat) + logdetVy)


@test_approx_eq coef(phynetlm) betahat
@test_approx_eq nobs(phynetlm) ntaxa
@test_approx_eq residuals(phynetlm) resids
@test_approx_eq model_response(phynetlm) Y
@test_approx_eq predict(phynetlm) fittedValues
@test_approx_eq df_residual(phynetlm) ntaxa-length(betahat)
@test_approx_eq sigma2_estim(phynetlm) sigma2hat 
@test_approx_eq loglikelihood(phynetlm) loglik
@test_approx_eq vcov(phynetlm) sigma2hat*ntaxa/(ntaxa-length(betahat))*inv(XtVyinv * X)
@test_approx_eq stderr(phynetlm) sqrt(diag(sigma2hat*ntaxa/(ntaxa-length(betahat))*inv(XtVyinv * X)))
@test_approx_eq df(phynetlm)  length(betahat)+1
@test_approx_eq deviance(phynetlm) sigma2hat * ntaxa
@test_approx_eq nulldeviance(phynetlm) nullsigma2hat * ntaxa 
@test_approx_eq nullloglikelihood(phynetlm) nullloglik
@test_approx_eq loglikelihood(phynetlm) nullloglikelihood(phynetlm)
@test_approx_eq deviance(phynetlm) nulldeviance(phynetlm)
#@test_approx_eq R2(phynetlm) 1-sigma2hat / nullsigma2hat 
#@test_approx_eq adjR2(phynetlm) 1 - (1 - (1-sigma2hat/nullsigma2hat))*(ntaxa-1)/(ntaxa-length(betahat)) 
@test_approx_eq AIC(phynetlm) -2*loglik+2*(length(betahat)+1)
@test_approx_eq AICc(phynetlm) -2*loglik+2*(length(betahat)+1)+2(length(betahat)+1)*((length(betahat)+1)+1)/(ntaxa-(length(betahat)+1)-1)
@test_approx_eq BIC(phynetlm) -2*loglik+(length(betahat)+1)*log(ntaxa)


# Simulate correlated data in data frames
b0 = 1
b1 = 10
sim = simulate(net, paramsBM(1, 1))
A = sim[:Tips]
B = b0 + b1 * A + simulate(net,  paramsBM(0, 0.1))[:Tips]

# With Matrices
X = hcat(ones(4), A)
fit_mat = phyloNetworklm(X, B, net)
@show fit_mat

# Naive version (GLS)
ntaxa = length(B)
Vy = fit_mat.Vy
Vyinv = inv(Vy)
XtVyinv = X' * Vyinv
logdetVy = logdet(Vy)
betahat = inv(XtVyinv * X) * XtVyinv * B
fittedValues =  X * betahat
resids = B - fittedValues
sigma2hat = 1/ntaxa * (resids' * Vyinv * resids)
# log likelihood
loglik = - 1 / 2 * (ntaxa + ntaxa * log(2 * pi) + ntaxa * log(sigma2hat) + logdetVy)
# null version
nullX = ones(ntaxa, 1)
nullXtVyinv = nullX' * Vyinv
nullresids = B - nullX * inv(nullXtVyinv * nullX) * nullXtVyinv * B
nullsigma2hat = 1/ntaxa * (nullresids' * Vyinv * nullresids)
nullloglik = - 1 / 2 * (ntaxa + ntaxa * log(2 * pi) + ntaxa * log(nullsigma2hat) + logdetVy)
@test_approx_eq coef(fit_mat) betahat
@test_approx_eq nobs(fit_mat) ntaxa
@test_approx_eq residuals(fit_mat) resids
@test_approx_eq model_response(fit_mat) B
@test_approx_eq predict(fit_mat) fittedValues
@test_approx_eq df_residual(fit_mat) ntaxa-length(betahat)
@test_approx_eq sigma2_estim(fit_mat) sigma2hat
@test_approx_eq loglikelihood(fit_mat) loglik
@test_approx_eq vcov(fit_mat) sigma2hat*ntaxa/(ntaxa-length(betahat)).*inv(XtVyinv * X)
@test_approx_eq stderr(fit_mat) sqrt(diag(sigma2hat*ntaxa/(ntaxa-length(betahat)).*inv(XtVyinv * X)))
@test_approx_eq df(fit_mat)  length(betahat)+1
@test_approx_eq deviance(fit_mat) sigma2hat * ntaxa
@test_approx_eq nulldeviance(fit_mat) nullsigma2hat * ntaxa
@test_approx_eq nullloglikelihood(fit_mat) nullloglik
@test_approx_eq R2(fit_mat) 1-sigma2hat / nullsigma2hat
@test_approx_eq adjR2(fit_mat) 1 - (1 - (1-sigma2hat/nullsigma2hat))*(ntaxa-1)/(ntaxa-length(betahat))
@test_approx_eq AIC(fit_mat) -2*loglik+2*(length(betahat)+1)
@test_approx_eq AICc(fit_mat) -2*loglik+2*(length(betahat)+1)+2(length(betahat)+1)*((length(betahat)+1)+1)/(ntaxa-(length(betahat)+1)-1)
@test_approx_eq BIC(fit_mat) -2*loglik+(length(betahat)+1)*log(ntaxa)

## perfect user using right format and formula
dfr = DataFrame(trait = B, pred = A, tipsNames = sim.M.tipsNames)
phynetlm = phyloNetworklm(trait ~ pred, dfr, net)
@show phynetlm

@test_approx_eq coef(phynetlm) coef(fit_mat)
@test_approx_eq vcov(phynetlm) vcov(fit_mat)
@test_approx_eq nobs(phynetlm) nobs(fit_mat)
@test_approx_eq residuals(phynetlm) residuals(fit_mat)
@test_approx_eq model_response(phynetlm) model_response(fit_mat)
@test_approx_eq predict(phynetlm) predict(fit_mat)
@test_approx_eq df_residual(phynetlm) df_residual(fit_mat)
@test_approx_eq sigma2_estim(phynetlm) sigma2_estim(fit_mat)
@test_approx_eq stderr(phynetlm) stderr(fit_mat)
@test_approx_eq confint(phynetlm) confint(fit_mat)
@test_approx_eq loglikelihood(phynetlm) loglikelihood(fit_mat)
@test_approx_eq df(phynetlm)  df(fit_mat)
@test_approx_eq deviance(phynetlm)  deviance(fit_mat)
@test_approx_eq nulldeviance(phynetlm)  nulldeviance(fit_mat)
@test_approx_eq nullloglikelihood(phynetlm)  nullloglikelihood(fit_mat)
@test_approx_eq R2(phynetlm)  R2(fit_mat)
@test_approx_eq adjR2(phynetlm)  adjR2(fit_mat)
@test_approx_eq AIC(phynetlm)  AIC(fit_mat)
@test_approx_eq AICc(phynetlm)  AICc(fit_mat)
@test_approx_eq BIC(fit_mat)  BIC(fit_mat)



# unordered data
dfr = dfr[[2, 1, 4, 3], :]
fitbis = phyloNetworklm(trait ~ pred, dfr, net)
@show fitbis

@test_approx_eq coef(phynetlm) coef(fitbis)
@test_approx_eq vcov(phynetlm) vcov(fitbis)
@test_approx_eq nobs(phynetlm) nobs(fitbis)
@test_approx_eq residuals(phynetlm)[fitbis.model.ind] residuals(fitbis)
@test_approx_eq model_response(phynetlm)[fitbis.model.ind] model_response(fitbis)
@test_approx_eq predict(phynetlm)[fitbis.model.ind] predict(fitbis)
@test_approx_eq df_residual(phynetlm) df_residual(fitbis)
@test_approx_eq sigma2_estim(phynetlm) sigma2_estim(fitbis)
@test_approx_eq stderr(phynetlm) stderr(fitbis)
@test_approx_eq confint(phynetlm) confint(fitbis)
@test_approx_eq loglikelihood(phynetlm) loglikelihood(fitbis)
@test_approx_eq df(phynetlm)  df(fitbis)
@test_approx_eq deviance(phynetlm)  deviance(fitbis)
@test_approx_eq nulldeviance(phynetlm)  nulldeviance(fitbis)
@test_approx_eq nullloglikelihood(phynetlm)  nullloglikelihood(fitbis)
@test_approx_eq R2(phynetlm)  R2(fitbis)
@test_approx_eq adjR2(phynetlm)  adjR2(fitbis)
@test_approx_eq AIC(phynetlm)  AIC(fitbis)
@test_approx_eq AICc(phynetlm)  AICc(fitbis)
@test_approx_eq BIC(fitbis)  BIC(fitbis)
@test_approx_eq mu_estim(fitbis)  mu_estim(fitbis)


# unnamed ordered data
dfr = DataFrame(trait = B, pred = A)
fitter = phyloNetworklm(trait ~ pred, dfr, net, no_names=true)
@show fitter

@test_approx_eq coef(phynetlm) coef(fitter)
@test_approx_eq vcov(phynetlm) vcov(fitter)
@test_approx_eq nobs(phynetlm) nobs(fitter)
@test_approx_eq residuals(phynetlm) residuals(fitter)
@test_approx_eq model_response(phynetlm) model_response(fitter)
@test_approx_eq predict(phynetlm) predict(fitter)
@test_approx_eq df_residual(phynetlm) df_residual(fitter)
@test_approx_eq sigma2_estim(phynetlm) sigma2_estim(fitter)
@test_approx_eq stderr(phynetlm) stderr(fitter)
@test_approx_eq confint(phynetlm) confint(fitter)
@test_approx_eq loglikelihood(phynetlm) loglikelihood(fitter)
@test_approx_eq df(phynetlm)  df(fitter)
@test_approx_eq deviance(phynetlm)  deviance(fitter)
@test_approx_eq nulldeviance(phynetlm)  nulldeviance(fitter)
@test_approx_eq nullloglikelihood(phynetlm)  nullloglikelihood(fitter)
@test_approx_eq R2(phynetlm)  R2(fitter)
@test_approx_eq adjR2(phynetlm)  adjR2(fitter)
@test_approx_eq AIC(phynetlm)  AIC(fitter)
@test_approx_eq AICc(phynetlm)  AICc(fitter)
@test_approx_eq BIC(fitter)  BIC(fitter)


# unnamed un-ordered data
dfr = dfr[[2, 1, 4, 3], :]
@test_throws ErrorException fitter = phyloNetworklm(trait ~ pred, dfr, net) # Wrong pred


### Add NAs
dfr = DataFrame(trait = B, pred = A, tipsNames = tipLabels(sim))
dfr[2, :pred] = NA
fitna = phyloNetworklm(trait ~ pred, dfr, net)
@show fitna

dfr = dfr[[2, 1, 4, 3], :]
fitnabis = phyloNetworklm(trait ~ pred, dfr, net)
@show fitnabis

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
@test_approx_eq df(fitna)  df(fitnabis)
@test_approx_eq deviance(fitna)  deviance(fitnabis)
@test_approx_eq nulldeviance(fitna)  nulldeviance(fitnabis)
@test_approx_eq nullloglikelihood(fitna)  nullloglikelihood(fitnabis)
@test_approx_eq R2(fitna)  R2(fitnabis)
@test_approx_eq adjR2(fitna)  adjR2(fitnabis)
@test_approx_eq AIC(fitna)  AIC(fitnabis)
@test_approx_eq AICc(fitna)  AICc(fitnabis)
@test_approx_eq BIC(fitnabis)  BIC(fitnabis)



### Ancestral State Reconstruction
params = paramsBM(3, 1)
sim = simulate(net, params)
Y = sim[:Tips]
# From known parameters
ancestral_traits = ancestralStateReconstruction(net, Y, params)
# BLUP
dfr = DataFrame(trait = Y, tipsNames = tipLabels(sim))
phynetlm = phyloNetworklm(trait~1, dfr, net)
eblup = ancestralStateReconstruction(phynetlm)



