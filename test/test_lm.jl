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
# Naive version (GLS)
ntaxa = length(Y)
Vy = fit.Vy
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


@test_approx_eq coef(fit) betahat
@test_approx_eq nobs(fit) ntaxa
@test_approx_eq residuals(fit) resids
@test_approx_eq model_response(fit) Y
@test_approx_eq predict(fit) fittedValues
@test_approx_eq df_residual(fit) ntaxa-length(betahat)
@test_approx_eq sigma2_estim(fit) sigma2hat 
@test_approx_eq loglikelihood(fit) loglik
@test_approx_eq vcov(fit) sigma2hat*ntaxa/(ntaxa-length(betahat))*inv(XtVyinv * X)
@test_approx_eq stderr(fit) sqrt(diag(sigma2hat*ntaxa/(ntaxa-length(betahat))*inv(XtVyinv * X)))
@test_approx_eq df(fit)  length(betahat)+1
@test_approx_eq deviance(fit) sigma2hat * ntaxa
@test_approx_eq nulldeviance(fit) nullsigma2hat * ntaxa 
@test_approx_eq nullloglikelihood(fit) nullloglik
@test_approx_eq loglikelihood(fit) nullloglikelihood(fit)
@test_approx_eq deviance(fit) nulldeviance(fit)
#@test_approx_eq R2(fit) 1-sigma2hat / nullsigma2hat 
#@test_approx_eq adjR2(fit) 1 - (1 - (1-sigma2hat/nullsigma2hat))*(ntaxa-1)/(ntaxa-length(betahat)) 
@test_approx_eq AIC(fit) -2*loglik+2*(length(betahat)+1)
@test_approx_eq AICc(fit) -2*loglik+2*(length(betahat)+1)+2(length(betahat)+1)*((length(betahat)+1)+1)/(ntaxa-(length(betahat)+1)-1)
@test_approx_eq BIC(fit) -2*loglik+(length(betahat)+1)*log(ntaxa)


# Simulate correlated data in data frames
b0 = 1
b1 = 10
sim = simulate(net, paramsBM(1, 1))
A = sim[:Tips]
B = b0 + b1 * A + simulate(net,  paramsBM(0, 0.1))[:Tips]

# With Matrices
X = hcat(ones(4), A)
fit_mat = phyloNetworklm(X, B, net)

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
fit = phyloNetworklm(trait ~ pred, dfr, net)

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
@test_approx_eq df(fit)  df(fit_mat)
@test_approx_eq deviance(fit)  deviance(fit_mat)
@test_approx_eq nulldeviance(fit)  nulldeviance(fit_mat)
@test_approx_eq nullloglikelihood(fit)  nullloglikelihood(fit_mat)
@test_approx_eq R2(fit)  R2(fit_mat)
@test_approx_eq adjR2(fit)  adjR2(fit_mat)
@test_approx_eq AIC(fit)  AIC(fit_mat)
@test_approx_eq AICc(fit)  AICc(fit_mat)
@test_approx_eq BIC(fit_mat)  BIC(fit_mat)



# unordered data
dfr = dfr[[2, 1, 4, 3], :]
fitbis = phyloNetworklm(trait ~ pred, dfr, net)

@test_approx_eq coef(fit) coef(fitbis)
@test_approx_eq vcov(fit) vcov(fitbis)
@test_approx_eq nobs(fit) nobs(fitbis)
@test_approx_eq residuals(fit)[fitbis.model.ind] residuals(fitbis)
@test_approx_eq model_response(fit)[fitbis.model.ind] model_response(fitbis)
@test_approx_eq predict(fit)[fitbis.model.ind] predict(fitbis)
@test_approx_eq df_residual(fit) df_residual(fitbis)
@test_approx_eq sigma2_estim(fit) sigma2_estim(fitbis)
@test_approx_eq stderr(fit) stderr(fitbis)
@test_approx_eq confint(fit) confint(fitbis)
@test_approx_eq loglikelihood(fit) loglikelihood(fitbis)
@test_approx_eq df(fit)  df(fitbis)
@test_approx_eq deviance(fit)  deviance(fitbis)
@test_approx_eq nulldeviance(fit)  nulldeviance(fitbis)
@test_approx_eq nullloglikelihood(fit)  nullloglikelihood(fitbis)
@test_approx_eq R2(fit)  R2(fitbis)
@test_approx_eq adjR2(fit)  adjR2(fitbis)
@test_approx_eq AIC(fit)  AIC(fitbis)
@test_approx_eq AICc(fit)  AICc(fitbis)
@test_approx_eq BIC(fitbis)  BIC(fitbis)
@test_approx_eq mu_estim(fitbis)  mu_estim(fitbis)


# unnamed ordered data
dfr = DataFrame(trait = B, pred = A)
fitter = phyloNetworklm(trait ~ pred, dfr, net)

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
@test_approx_eq df(fit)  df(fitter)
@test_approx_eq deviance(fit)  deviance(fitter)
@test_approx_eq nulldeviance(fit)  nulldeviance(fitter)
@test_approx_eq nullloglikelihood(fit)  nullloglikelihood(fitter)
@test_approx_eq R2(fit)  R2(fitter)
@test_approx_eq adjR2(fit)  adjR2(fitter)
@test_approx_eq AIC(fit)  AIC(fitter)
@test_approx_eq AICc(fit)  AICc(fitter)
@test_approx_eq BIC(fitter)  BIC(fitter)


# unnamed un-ordered data
dfr = dfr[[2, 1, 4, 3], :]
fitter = phyloNetworklm(trait ~ pred, dfr, net) # Wrong pred


### Add NAs
dfr = DataFrame(trait = B, pred = A, tipsNames = tipLabels(sim))
dfr[2, :pred] = NA
fitna = phyloNetworklm(trait ~ pred, dfr, net)

dfr = dfr[[2, 1, 4, 3], :]
fitnabis = phyloNetworklm(trait ~ pred, dfr, net)

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
fit = phyloNetworklm(trait~1, dfr, net)
eblup = ancestralStateReconstruction(fit)



