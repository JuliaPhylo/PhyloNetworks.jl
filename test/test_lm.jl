## tests of PhyloNetworklm

tree_str= "(A:2.5,((B:1,#H1:0.5::0.1):1,(C:1,(D:0.5)#H1:0.5::0.9):1):0.5);"
net = readTopology(tree_str)
preorder!(net)
# printEdges(net)
# plot(net, useEdgeLength = true,  showEdgeLength=true, showGamma=true)

# Rk: Is there a way to check that the branch length are coherent with
# one another (Especialy for hybrids) ?
# Not yet (CA, 2016-12-01).
# Would be great to add functions to calculate distance node - root.
# several such distances depending on path: 2 parent choices at each hybrid

# Ancestral state reconstruction with ready-made matrices
params = ParamsBM(10, 1)
srand(2468) # sets the seed for reproducibility, to debug potential error
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
@test_approx_eq dof_residual(phynetlm) ntaxa-length(betahat)
@test_approx_eq sigma2_estim(phynetlm) sigma2hat
@test_approx_eq loglikelihood(phynetlm) loglik
@test_approx_eq vcov(phynetlm) sigma2hat*ntaxa/(ntaxa-length(betahat))*inv(XtVyinv * X)
@test_approx_eq stderr(phynetlm) sqrt(diag(sigma2hat*ntaxa/(ntaxa-length(betahat))*inv(XtVyinv * X)))
@test_approx_eq dof(phynetlm)  length(betahat)+1
@test_approx_eq deviance(phynetlm) sigma2hat * ntaxa
@test_approx_eq nulldeviance(phynetlm) nullsigma2hat * ntaxa
@test_approx_eq nullloglikelihood(phynetlm) nullloglik
@test_approx_eq loglikelihood(phynetlm) nullloglikelihood(phynetlm)
@test_approx_eq deviance(phynetlm) nulldeviance(phynetlm)
@test_approx_eq_eps r2(phynetlm) 1-sigma2hat / nullsigma2hat 1e-15
@test_approx_eq_eps adjr2(phynetlm) 1 - (1 - (1-sigma2hat/nullsigma2hat))*(ntaxa-1)/(ntaxa-length(betahat)) 1e-15
@test_approx_eq aic(phynetlm) -2*loglik+2*(length(betahat)+1)
@test_approx_eq aicc(phynetlm) -2*loglik+2*(length(betahat)+1)+2(length(betahat)+1)*((length(betahat)+1)+1)/(ntaxa-(length(betahat)+1)-1)
@test_approx_eq bic(phynetlm) -2*loglik+(length(betahat)+1)*log(ntaxa)

# with data frames
dfr = DataFrame(trait = Y, tipNames = sim.M.tipNames)
fitbis = phyloNetworklm(trait ~ 1, dfr, net)
@show fitbis

@test_approx_eq coef(phynetlm) coef(fitbis)
@test_approx_eq vcov(phynetlm) vcov(fitbis)
@test_approx_eq nobs(phynetlm) nobs(fitbis)
@test_approx_eq residuals(phynetlm)[fitbis.model.ind] residuals(fitbis)
@test_approx_eq model_response(phynetlm)[fitbis.model.ind] model_response(fitbis)
@test_approx_eq predict(phynetlm)[fitbis.model.ind] predict(fitbis)
@test_approx_eq dof_residual(phynetlm) dof_residual(fitbis)
@test_approx_eq sigma2_estim(phynetlm) sigma2_estim(fitbis)
@test_approx_eq stderr(phynetlm) stderr(fitbis)
@test_approx_eq confint(phynetlm) confint(fitbis)
@test_approx_eq loglikelihood(phynetlm) loglikelihood(fitbis)
@test_approx_eq dof(phynetlm)  dof(fitbis)
@test_approx_eq deviance(phynetlm)  deviance(fitbis)
@test_approx_eq nulldeviance(phynetlm)  nulldeviance(fitbis)
@test_approx_eq nullloglikelihood(phynetlm)  nullloglikelihood(fitbis)
@test_approx_eq_eps r2(phynetlm)  r2(fitbis) 1e-15
@test_approx_eq_eps adjr2(phynetlm)  adjr2(fitbis) 1e-15
@test_approx_eq aic(phynetlm)  aic(fitbis)
@test_approx_eq aicc(phynetlm)  aicc(fitbis)
@test_approx_eq bic(phynetlm)  bic(fitbis)
@test_approx_eq mu_estim(phynetlm)  mu_estim(fitbis)

## Pagel's Lambda
fitlam = phyloNetworklm(trait ~ 1, dfr, net, model = "lambda", fixedValue=1.0)
@show fitlam

@test_approx_eq lambda_estim(fitlam) 1.0
@test_approx_eq coef(fitlam) coef(fitbis)
@test_approx_eq vcov(fitlam) vcov(fitbis)
@test_approx_eq nobs(fitlam) nobs(fitbis)
@test_approx_eq residuals(fitlam)[fitbis.model.ind] residuals(fitbis)
@test_approx_eq model_response(fitlam)[fitbis.model.ind] model_response(fitbis)
@test_approx_eq predict(fitlam)[fitbis.model.ind] predict(fitbis)
@test_approx_eq dof_residual(fitlam) dof_residual(fitbis)
@test_approx_eq sigma2_estim(fitlam) sigma2_estim(fitbis)
@test_approx_eq stderr(fitlam) stderr(fitbis)
@test_approx_eq confint(fitlam) confint(fitbis)
@test_approx_eq loglikelihood(fitlam) loglikelihood(fitbis)
@test_approx_eq dof(fitlam)  dof(fitbis) + 1
@test_approx_eq deviance(fitlam)  deviance(fitbis)
@test_approx_eq nulldeviance(fitlam)  nulldeviance(fitbis)
@test_approx_eq nullloglikelihood(fitlam)  nullloglikelihood(fitbis)
@test_approx_eq_eps r2(fitlam)  r2(fitbis) 1e-15
@test_approx_eq_eps adjr2(fitlam)  adjr2(fitbis) - 0.5 1e-15
@test_approx_eq aic(fitlam)  aic(fitbis) + 2
#@test_approx_eq aicc(fitlam)  aicc(fitbis)
@test_approx_eq bic(fitlam)  bic(fitbis) + log(nobs(fitbis))
@test_approx_eq mu_estim(fitlam)  mu_estim(fitbis)

## Pagel's Lambda
fitlam = phyloNetworklm(trait ~ 1, dfr, net, model = "lambda")
@show fitlam
@test_approx_eq lambda_estim(fitlam) 1.24875

###############################################################################
#### Other Network
###############################################################################
# originally: "(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);"
# followed by changes in net.edge[?].length values to make the network ultrametric
net = readTopology("(((Ag:5,(#H1:1::0.056,((Ak:2,(E:1,#H2:1::0.004):1):1,(M:2)#H2:1::0.996):1):1):1,(((((Az:1,Ag2:1):1,As:2):1)#H1:1::0.944,Ap:4):1,Ar:5):1):1,(P:4,20:4):3,165:7);");
# plot(net, useEdgeLength = true,  showEdgeNumber=true)

#### Simulate correlated data in data frames ####
b0 = 1
b1 = 10
srand(5678)
sim = simulate(net, ParamsBM(1, 1))
A = sim[:Tips]
B = b0 + b1 * A + simulate(net,  ParamsBM(0, 0.1))[:Tips]

# With Matrices
X = hcat(ones(12), A)
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
@test_approx_eq dof_residual(fit_mat) ntaxa-length(betahat)
@test_approx_eq sigma2_estim(fit_mat) sigma2hat
@test_approx_eq loglikelihood(fit_mat) loglik
@test_approx_eq vcov(fit_mat) sigma2hat*ntaxa/(ntaxa-length(betahat)).*inv(XtVyinv * X)
@test_approx_eq stderr(fit_mat) sqrt(diag(sigma2hat*ntaxa/(ntaxa-length(betahat)).*inv(XtVyinv * X)))
@test_approx_eq dof(fit_mat)  length(betahat)+1
@test_approx_eq deviance(fit_mat) sigma2hat * ntaxa
@test_approx_eq nulldeviance(fit_mat) nullsigma2hat * ntaxa
@test_approx_eq nullloglikelihood(fit_mat) nullloglik
@test_approx_eq_eps r2(fit_mat) 1-sigma2hat / nullsigma2hat 1e-15
@test_approx_eq_eps adjr2(fit_mat) 1 - (1 - (1-sigma2hat/nullsigma2hat))*(ntaxa-1)/(ntaxa-length(betahat)) 1e-15
@test_approx_eq aic(fit_mat) -2*loglik+2*(length(betahat)+1)
@test_approx_eq aicc(fit_mat) -2*loglik+2*(length(betahat)+1)+2(length(betahat)+1)*((length(betahat)+1)+1)/(ntaxa-(length(betahat)+1)-1)
@test_approx_eq bic(fit_mat) -2*loglik+(length(betahat)+1)*log(ntaxa)

## perfect user using right format and formula
dfr = DataFrame(trait = B, pred = A, tipNames = sim.M.tipNames)
phynetlm = phyloNetworklm(trait ~ pred, dfr, net)
@show phynetlm

@test_approx_eq coef(phynetlm) coef(fit_mat)
@test_approx_eq vcov(phynetlm) vcov(fit_mat)
@test_approx_eq nobs(phynetlm) nobs(fit_mat)
@test_approx_eq residuals(phynetlm) residuals(fit_mat)
@test_approx_eq model_response(phynetlm) model_response(fit_mat)
@test_approx_eq predict(phynetlm) predict(fit_mat)
@test_approx_eq dof_residual(phynetlm) dof_residual(fit_mat)
@test_approx_eq sigma2_estim(phynetlm) sigma2_estim(fit_mat)
@test_approx_eq stderr(phynetlm) stderr(fit_mat)
@test_approx_eq confint(phynetlm) confint(fit_mat)
@test_approx_eq loglikelihood(phynetlm) loglikelihood(fit_mat)
@test_approx_eq dof(phynetlm)  dof(fit_mat)
@test_approx_eq deviance(phynetlm)  deviance(fit_mat)
@test_approx_eq nulldeviance(phynetlm)  nulldeviance(fit_mat)
@test_approx_eq nullloglikelihood(phynetlm)  nullloglikelihood(fit_mat)
@test_approx_eq r2(phynetlm)  r2(fit_mat)
@test_approx_eq adjr2(phynetlm)  adjr2(fit_mat)
@test_approx_eq aic(phynetlm)  aic(fit_mat)
@test_approx_eq aicc(phynetlm)  aicc(fit_mat)
@test_approx_eq bic(phynetlm)  bic(fit_mat)

# unordered data
srand(1234)
dfr = dfr[sample(1:12, 12, replace=false), :]
fitbis = phyloNetworklm(trait ~ pred, dfr, net)

@test_approx_eq coef(phynetlm) coef(fitbis)
@test_approx_eq vcov(phynetlm) vcov(fitbis)
@test_approx_eq nobs(phynetlm) nobs(fitbis)
@test_approx_eq residuals(phynetlm)[fitbis.model.ind] residuals(fitbis)
@test_approx_eq model_response(phynetlm)[fitbis.model.ind] model_response(fitbis)
@test_approx_eq predict(phynetlm)[fitbis.model.ind] predict(fitbis)
@test_approx_eq dof_residual(phynetlm) dof_residual(fitbis)
@test_approx_eq sigma2_estim(phynetlm) sigma2_estim(fitbis)
@test_approx_eq stderr(phynetlm) stderr(fitbis)
@test_approx_eq confint(phynetlm) confint(fitbis)
@test_approx_eq loglikelihood(phynetlm) loglikelihood(fitbis)
@test_approx_eq dof(phynetlm)  dof(fitbis)
@test_approx_eq deviance(phynetlm)  deviance(fitbis)
@test_approx_eq nulldeviance(phynetlm)  nulldeviance(fitbis)
@test_approx_eq nullloglikelihood(phynetlm)  nullloglikelihood(fitbis)
@test_approx_eq r2(phynetlm)  r2(fitbis)
@test_approx_eq adjr2(phynetlm)  adjr2(fitbis)
@test_approx_eq aic(phynetlm)  aic(fitbis)
@test_approx_eq aicc(phynetlm)  aicc(fitbis)
@test_approx_eq bic(phynetlm)  bic(fitbis)
@test_approx_eq mu_estim(phynetlm)  mu_estim(fitbis)

# unnamed ordered data
dfr = DataFrame(trait = B, pred = A)
fitter = phyloNetworklm(trait ~ pred, dfr, net, no_names=true)

@test_approx_eq coef(phynetlm) coef(fitter)
@test_approx_eq vcov(phynetlm) vcov(fitter)
@test_approx_eq nobs(phynetlm) nobs(fitter)
@test_approx_eq residuals(phynetlm) residuals(fitter)
@test_approx_eq model_response(phynetlm) model_response(fitter)
@test_approx_eq predict(phynetlm) predict(fitter)
@test_approx_eq dof_residual(phynetlm) dof_residual(fitter)
@test_approx_eq sigma2_estim(phynetlm) sigma2_estim(fitter)
@test_approx_eq stderr(phynetlm) stderr(fitter)
@test_approx_eq confint(phynetlm) confint(fitter)
@test_approx_eq loglikelihood(phynetlm) loglikelihood(fitter)
@test_approx_eq dof(phynetlm)  dof(fitter)
@test_approx_eq deviance(phynetlm)  deviance(fitter)
@test_approx_eq nulldeviance(phynetlm)  nulldeviance(fitter)
@test_approx_eq nullloglikelihood(phynetlm)  nullloglikelihood(fitter)
@test_approx_eq r2(phynetlm)  r2(fitter)
@test_approx_eq adjr2(phynetlm)  adjr2(fitter)
@test_approx_eq aic(phynetlm)  aic(fitter)
@test_approx_eq aicc(phynetlm)  aicc(fitter)
@test_approx_eq bic(phynetlm)  bic(fitter)

# unnamed un-ordered data
dfr = dfr[sample(1:12, 12, replace=false), :]
@test_throws ErrorException fitter = phyloNetworklm(trait ~ pred, dfr, net) # Wrong pred

### Add NAs
dfr = DataFrame(trait = B, pred = A, tipNames = tipLabels(sim))
dfr[[2, 8, 11], :pred] = NA
fitna = phyloNetworklm(trait ~ pred, dfr, net)
@show fitna

dfr = dfr[sample(1:12, 12, replace=false), :]
fitnabis = phyloNetworklm(trait ~ pred, dfr, net)

@test_approx_eq coef(fitna) coef(fitnabis)
@test_approx_eq vcov(fitna) vcov(fitnabis)
@test_approx_eq nobs(fitna) nobs(fitnabis)
@test_approx_eq sort(residuals(fitna)) sort(residuals(fitnabis))
@test_approx_eq sort(model_response(fitna)) sort(model_response(fitnabis))
@test_approx_eq sort(predict(fitna)) sort(predict(fitnabis))
@test_approx_eq dof_residual(fitna) dof_residual(fitnabis)
@test_approx_eq sigma2_estim(fitna) sigma2_estim(fitnabis)
@test_approx_eq stderr(fitna) stderr(fitnabis)
@test_approx_eq confint(fitna) confint(fitnabis)
@test_approx_eq loglikelihood(fitna) loglikelihood(fitnabis)
@test_approx_eq dof(fitna)  dof(fitnabis)
@test_approx_eq deviance(fitna)  deviance(fitnabis)
@test_approx_eq nulldeviance(fitna)  nulldeviance(fitnabis)
@test_approx_eq nullloglikelihood(fitna)  nullloglikelihood(fitnabis)
@test_approx_eq r2(fitna)  r2(fitnabis)
@test_approx_eq adjr2(fitna)  adjr2(fitnabis)
@test_approx_eq aic(fitna)  aic(fitnabis)
@test_approx_eq aicc(fitna)  aicc(fitnabis)
@test_approx_eq bic(fitna)  bic(fitnabis)

## Pagel's Lambda
fitlam = phyloNetworklm(trait ~ pred, dfr, net, model = "lambda", fixedValue = 1.0)
@show fitlam

@test_approx_eq lambda_estim(fitlam) 1.0
@test_approx_eq coef(fitlam) coef(fitnabis)
@test_approx_eq vcov(fitlam) vcov(fitnabis)
@test_approx_eq nobs(fitlam) nobs(fitnabis)
@test_approx_eq residuals(fitlam) residuals(fitnabis)
@test_approx_eq model_response(fitlam) model_response(fitnabis)
@test_approx_eq predict(fitlam) predict(fitnabis)
@test_approx_eq dof_residual(fitlam) dof_residual(fitnabis)
@test_approx_eq sigma2_estim(fitlam) sigma2_estim(fitnabis)
@test_approx_eq stderr(fitlam) stderr(fitnabis)
@test_approx_eq confint(fitlam) confint(fitnabis)
@test_approx_eq loglikelihood(fitlam) loglikelihood(fitnabis)
@test_approx_eq dof(fitlam)  dof(fitnabis) + 1
@test_approx_eq deviance(fitlam)  deviance(fitnabis)
@test_approx_eq nulldeviance(fitlam)  nulldeviance(fitnabis)
@test_approx_eq nullloglikelihood(fitlam)  nullloglikelihood(fitnabis)
@test_approx_eq_eps r2(fitlam)  r2(fitnabis) 1e-15
@test_approx_eq_eps adjr2(fitlam)-1  (adjr2(fitnabis)-1)*(nobs(fitnabis)-dof(fitnabis)+1)/(nobs(fitnabis)-dof(fitlam)+1) 1e-15
@test_approx_eq aic(fitlam)  aic(fitnabis) + 2
#@test_approx_eq aicc(fitlam)  aicc(fitnabis)
@test_approx_eq bic(fitlam)  bic(fitnabis) + log(nobs(fitnabis))
@test_approx_eq mu_estim(fitlam)  mu_estim(fitnabis)

## Pagel's Lambda
fitlam = phyloNetworklm(trait ~ pred, dfr, net, model = "lambda")
@show fitlam
@test_approx_eq_eps lambda_estim(fitlam) 1.1135518305 1e-10


### Ancestral State Reconstruction
params = ParamsBM(3, 1)
sim = simulate(net, params)
Y = sim[:Tips]
# From known parameters
ancestral_traits = ancestralStateReconstruction(net, Y, params)
# BLUP
dfr = DataFrame(trait = Y, tipNames = tipLabels(sim))
phynetlm = phyloNetworklm(trait~1, dfr, net)
blup = ancestralStateReconstruction(phynetlm)
# plot(net, blup)
@show blup

# BLUP same, using the function dirrectly
blup_bis = ancestralStateReconstruction(dfr, net)

@test_approx_eq expectations(blup)[:condExpectation] expectations(blup_bis)[:condExpectation]
@test_approx_eq expectations(blup)[:nodeNumber] expectations(blup_bis)[:nodeNumber]
@test_approx_eq blup.traits_tips blup_bis.traits_tips
@test_approx_eq blup.TipNumbers blup_bis.TipNumbers
@test_approx_eq predint(blup) predint(blup_bis)

dfr = DataFrame(trait = Y, tipNames = tipLabels(sim), reg = Y)
@test_throws ErrorException fitter = ancestralStateReconstruction(dfr, net) # cannot handle a predictor

# Unordered
dfr2 = dfr[sample(1:12, 12, replace=false), :]
phynetlm = phyloNetworklm(trait~1, dfr2, net)
blup2 = ancestralStateReconstruction(phynetlm)

@test_approx_eq expectations(blup)[:condExpectation][1:length(blup.NodeNumbers)] expectations(blup2)[:condExpectation][1:length(blup.NodeNumbers)]
@test_approx_eq blup.traits_tips[phynetlm.model.ind] blup2.traits_tips
@test_approx_eq blup.TipNumbers[phynetlm.model.ind] blup2.TipNumbers
@test_approx_eq predint(blup)[1:length(blup.NodeNumbers), :] predint(blup2)[1:length(blup.NodeNumbers), :]

# With unknown tips
dfr[[2, 4], :trait] = NA
phynetlm = phyloNetworklm(trait~1, dfr, net)
blup = ancestralStateReconstruction(phynetlm)
# plot(net, blup)

# Unordered
dfr2 = dfr[[1, 2, 5, 3, 4, 6, 7, 8, 9, 10, 11, 12], :]
phynetlm = phyloNetworklm(trait~1, dfr, net)
blup2 = ancestralStateReconstruction(phynetlm)

@test_approx_eq expectations(blup)[:condExpectation][1:length(blup.NodeNumbers)] expectations(blup2)[:condExpectation][1:length(blup.NodeNumbers)]
@test_approx_eq predint(blup)[1:length(blup.NodeNumbers), :] predint(blup2)[1:length(blup.NodeNumbers), :]

#################
## Data with no phylogenetic signal
#################

net = readTopology("(((Ag:5,(#H1:1::0.056,((Ak:2,(E:1,#H2:1::0.004):1):1,(M:2)#H2:1::0.996):1):1):1,(((((Az:1,Ag2:1):1,As:2):1)#H1:1::0.944,Ap:4):1,Ar:5):1):1,(P:4,20:4):3,165:7);");
# plot(net, useEdgeLength = true,  showEdgeNumber=true)

#### Simulate correlated data in data frames ####
b0 = 1
b1 = 10
srand(5678)
A = randn(size(tipLabels(net), 1))
B = b0 + b1 * A + randn(size(tipLabels(net), 1))
dfr = DataFrame(trait = B, pred = A, tipNames = tipLabels(net))

## Network
phynetlm = phyloNetworklm(trait ~ pred, dfr, net, model = "lambda")

@test_approx_eq_eps lambda_estim(phynetlm) 0.5894200143 1e-8

## Major Tree
tree = majorTree(net)
phynetlm = phyloNetworklm(trait ~ pred, dfr, tree, model = "lambda")

@test_approx_eq_eps lambda_estim(phynetlm) 0.5903394415 1e-6

############################
## Against no regressor
###########################
params = ParamsBM(10, 1)
srand(2468) # sets the seed for reproducibility, to debug potential error
sim = simulate(net, params)
Y = sim[:Tips]
phynetlm = phyloNetworklm(zeros(length(Y),0), Y, net)
@show phynetlm
# Naive version (GLS)
ntaxa = length(Y)
Vy = phynetlm.Vy
Vyinv = inv(Vy)
logdetVy = logdet(Vy)
fittedValues =  zeros(length(Y))
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

@test_approx_eq nobs(phynetlm) ntaxa
@test_approx_eq residuals(phynetlm) resids
@test_approx_eq model_response(phynetlm) Y
@test_approx_eq predict(phynetlm) fittedValues
@test_approx_eq dof_residual(phynetlm) ntaxa
@test_approx_eq sigma2_estim(phynetlm) sigma2hat
@test_approx_eq loglikelihood(phynetlm) loglik
@test_approx_eq deviance(phynetlm) sigma2hat * ntaxa
@test_approx_eq nulldeviance(phynetlm) nullsigma2hat * ntaxa
@test_approx_eq nullloglikelihood(phynetlm) nullloglik
@test_approx_eq_eps r2(phynetlm) 1-sigma2hat / nullsigma2hat 1e-14
@test_approx_eq_eps adjr2(phynetlm) 1 - (1 - (1-sigma2hat/nullsigma2hat))*(ntaxa-1)/(ntaxa) 1e-14
@test_approx_eq aic(phynetlm) -2*loglik+2*(1)
@test_approx_eq aicc(phynetlm) -2*loglik+2*(1)+2(1)*((1)+1)/(ntaxa-(1)-1)
@test_approx_eq bic(phynetlm) -2*loglik+(1)*log(ntaxa)

# with data frames
dfr = DataFrame(trait = Y, tipNames = sim.M.tipNames)
fitbis = phyloNetworklm(trait ~ -1, dfr, net)
@show fitbis
#@test_approx_eq coef(phynetlm) coef(fitbis)
#@test_approx_eq vcov(phynetlm) vcov(fitbis)
@test_approx_eq nobs(phynetlm) nobs(fitbis)
@test_approx_eq residuals(phynetlm)[fitbis.model.ind] residuals(fitbis)
@test_approx_eq model_response(phynetlm)[fitbis.model.ind] model_response(fitbis)
@test_approx_eq predict(phynetlm)[fitbis.model.ind] predict(fitbis)
@test_approx_eq dof_residual(phynetlm) dof_residual(fitbis)
@test_approx_eq sigma2_estim(phynetlm) sigma2_estim(fitbis)
#@test_approx_eq stderr(phynetlm) stderr(fitbis)
#@test_approx_eq confint(phynetlm) confint(fitbis)
@test_approx_eq loglikelihood(phynetlm) loglikelihood(fitbis)
#@test_approx_eq dof(phynetlm)  dof(fitbis)
@test_approx_eq deviance(phynetlm)  deviance(fitbis)
@test_approx_eq nulldeviance(phynetlm)  nulldeviance(fitbis)
@test_approx_eq nullloglikelihood(phynetlm)  nullloglikelihood(fitbis)
@test_approx_eq_eps r2(phynetlm)  r2(fitbis) 1e-15
@test_approx_eq_eps adjr2(phynetlm)  adjr2(fitbis) 1e-15
@test_approx_eq aic(phynetlm)  aic(fitbis)
@test_approx_eq aicc(phynetlm)  aicc(fitbis)
@test_approx_eq bic(phynetlm)  bic(fitbis)
#@test_approx_eq mu_estim(phynetlm)  mu_estim(fitbis)
