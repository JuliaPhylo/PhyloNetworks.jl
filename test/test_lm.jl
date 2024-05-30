## tests of phylolm

@testset "phylolm on small network" begin
global net

tree_str= "(A:2.5,((B:1,#H1:0.5::0.1):1,(C:1,(D:0.5)#H1:0.5::0.9):1):0.5);"
net = readTopology(tree_str)
preorder!(net)

# Rk: Is there a way to check that the branch length are coherent with
# one another (Especialy for hybrids) ?
# see QuartetNetworkGoodnessFit.ultrametrize! which can detect if the network is
# time-consistent: all paths from the root to a given node have the same length
# https://github.com/cecileane/QuartetNetworkGoodnessFit.jl

# Ancestral state reconstruction with ready-made matrices
params = ParamsBM(10, 1)
Random.seed!(2468); # simulates the Y values below under julia v1.6
sim = simulate(net, params) # tests that the simulation runs, but results not used
Y = [11.239539657364706,8.600423079191044,10.559841251147608,9.965748423156297] # sim[:Tips]
X = ones(4, 1)
phynetlm = phylolm(X, Y, net; reml=false)
@test_logs show(devnull, phynetlm)
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


@test coef(phynetlm) ≈ betahat
@test nobs(phynetlm) ≈ ntaxa
@test residuals(phynetlm) ≈ resids
@test response(phynetlm) ≈ Y
@test predict(phynetlm) ≈ fittedValues
@test dof_residual(phynetlm) ≈ ntaxa-length(betahat)
@test sigma2_phylo(phynetlm) ≈ sigma2hat
@test loglikelihood(phynetlm) ≈ loglik
@test vcov(phynetlm) ≈ sigma2hat*ntaxa/(ntaxa-length(betahat))*inv(XtVyinv * X)
@test stderror(phynetlm) ≈ sqrt.(diag(sigma2hat*ntaxa/(ntaxa-length(betahat))*inv(XtVyinv * X)))
@test dof(phynetlm) ≈ length(betahat)+1
@test deviance(phynetlm, Val(true)) ≈ sigma2hat * ntaxa
@test nulldeviance(phynetlm) ≈ nullsigma2hat * ntaxa
@test nullloglikelihood(phynetlm) ≈ nullloglik
@test loglikelihood(phynetlm) ≈ nullloglikelihood(phynetlm)
@test deviance(phynetlm, Val(true)) ≈ nulldeviance(phynetlm)
@test r2(phynetlm) ≈ 1-sigma2hat / nullsigma2hat atol=1e-15
@test adjr2(phynetlm) ≈ 1 - (1 - (1-sigma2hat/nullsigma2hat))*(ntaxa-1)/(ntaxa-length(betahat)) atol=1e-15
@test aic(phynetlm) ≈ -2*loglik+2*(length(betahat)+1)
@test aicc(phynetlm) ≈ -2*loglik+2*(length(betahat)+1)+2(length(betahat)+1)*((length(betahat)+1)+1)/(ntaxa-(length(betahat)+1)-1)
@test bic(phynetlm) ≈ -2*loglik+(length(betahat)+1)*log(ntaxa)
@test hasintercept(phynetlm)

# with data frames
dfr = DataFrame(trait = Y, tipNames = ["A","B","C","D"]) # sim.M.tipNames
fitbis = phylolm(@formula(trait ~ 1), dfr, net; reml=false)
#@show fitbis

@test coef(phynetlm) ≈ coef(fitbis)
@test vcov(phynetlm) ≈ vcov(fitbis)
@test nobs(phynetlm) ≈ nobs(fitbis)
@test residuals(phynetlm)[fitbis.ind] ≈ residuals(fitbis)
@test response(phynetlm)[fitbis.ind] ≈ response(fitbis)
@test predict(phynetlm)[fitbis.ind] ≈ predict(fitbis)
@test dof_residual(phynetlm) ≈ dof_residual(fitbis)
@test sigma2_phylo(phynetlm) ≈ sigma2_phylo(fitbis)
@test stderror(phynetlm) ≈ stderror(fitbis)
@test confint(phynetlm) ≈ confint(fitbis)
@test loglikelihood(phynetlm) ≈ loglikelihood(fitbis)
@test dof(phynetlm) ≈ dof(fitbis)
@test deviance(phynetlm, Val(true)) ≈ deviance(fitbis, Val(true))
@test nulldeviance(phynetlm) ≈ nulldeviance(fitbis)
@test nullloglikelihood(phynetlm) ≈ nullloglikelihood(fitbis)
@test r2(phynetlm) ≈ r2(fitbis) atol=1e-15
@test adjr2(phynetlm) ≈ adjr2(fitbis) atol=1e-15
@test aic(phynetlm) ≈ aic(fitbis)
@test aicc(phynetlm) ≈ aicc(fitbis)
@test bic(phynetlm) ≈ bic(fitbis)
tmp = (@test_logs (:warn, r"^You fitted the data against a custom matrix") mu_phylo(phynetlm))
@test tmp ≈ mu_phylo(fitbis)
@test hasintercept(phynetlm)

## fixed values parameters
fitlam = phylolm(@formula(trait ~ 1), dfr, net, model = "lambda", fixedValue=1.0, reml=false)
@test_logs show(devnull, fitlam)

@test lambda_estim(fitlam) ≈ 1.0
@test coef(fitlam) ≈ coef(fitbis)
@test vcov(fitlam) ≈ vcov(fitbis)
@test nobs(fitlam) ≈ nobs(fitbis)
@test residuals(fitlam)[fitbis.ind] ≈ residuals(fitbis)
@test response(fitlam)[fitbis.ind] ≈ response(fitbis)
@test predict(fitlam)[fitbis.ind] ≈ predict(fitbis)
@test dof_residual(fitlam) ≈ dof_residual(fitbis)
@test sigma2_phylo(fitlam) ≈ sigma2_phylo(fitbis)
@test stderror(fitlam) ≈ stderror(fitbis)
@test confint(fitlam) ≈ confint(fitbis)
@test loglikelihood(fitlam) ≈ loglikelihood(fitbis)
@test dof(fitlam) ≈ dof(fitbis) + 1
@test deviance(fitlam, Val(true)) ≈ deviance(fitbis, Val(true))
@test nulldeviance(fitlam) ≈ nulldeviance(fitbis)
@test nullloglikelihood(fitlam) ≈ nullloglikelihood(fitbis)
@test r2(fitlam) ≈ r2(fitbis) atol=1e-15
@test adjr2(fitlam) ≈ adjr2(fitbis) - 0.5 atol=1e-15
@test aic(fitlam) ≈ aic(fitbis) + 2
#@test aicc(fitlam) ≈ aicc(fitbis)
@test bic(fitlam) ≈ bic(fitbis) + log(nobs(fitbis))
@test mu_phylo(fitlam) ≈ mu_phylo(fitbis)
@test hasintercept(fitlam)

fitSH = phylolm(@formula(trait ~ 1), dfr, net, model="scalingHybrid", fixedValue=1.0, reml=false)
@test loglikelihood(fitlam) ≈ loglikelihood(fitSH)
@test aic(fitlam) ≈ aic(fitSH)

@test modelmatrix(fitlam) == reshape(ones(4), (4,1))
s = IOBuffer(); show(s, formula(fitlam))
@test String(take!(s)) == "trait ~ 1"
PhyloNetworks.lambda!(fitlam, 0.5)
@test PhyloNetworks.lambda(fitlam) == 0.5

## Pagel's Lambda
fitlam = (@test_logs (:info, r"^Maximum lambda value") match_mode=:any phylolm(@formula(trait ~ 1), dfr, net, model="lambda", reml=false))
@test lambda_estim(fitlam) ≈ 1.24875

## Scaling Hybrid
fitSH = phylolm(@formula(trait ~ 1), dfr, net, model="scalingHybrid", reml=false)
@test lambda_estim(fitSH) ≈ 4.057891910001937 atol=1e-5

end

###############################################################################
### With shifts
###############################################################################
@testset "Shifts and Transgressive Evolution" begin
global net
net = readTopology("(((Ag:5,(#H1:1::0.056,((Ak:2,(E:1,#H2:1::0.004):1):1,(M:2)#H2:1::0.996):1):1):1,(((((Az:1,Ag2:1):1,As:2):1)#H1:1::0.944,Ap:4):1,Ar:5):1):1,(P:4,20:4):3,165:7);");
preorder!(net)

## Simulate
params = ParamsBM(10, 0.1, shiftHybrid([3.0, -3.0],  net))
Random.seed!(2468); # sets the seed for reproducibility, to debug potential error
sim = simulate(net, params) # checks for no error, but not used.
# values simulated using julia v1.6.4's RNG hardcoded below.
# Y = sim[:Tips]
Y = [11.640085037749985, 9.498284887480622, 9.568813792749083, 13.036916724865296, 6.873936265709946, 6.536647349405742, 5.95771939864956, 10.517318306450647, 9.34927049737206, 10.176238483133424, 10.760099940744308, 8.955543827353837]

## Construct regression matrix
dfr_shift = regressorShift(net.edge[[8,17]], net)
dfr_shift[!,:sum] = vec(sum(Matrix(dfr_shift[:,findall(DataFrames.propertynames(dfr_shift) .!= :tipNames)]), dims=2))
dfr_hybrid = regressorHybrid(net)

@test dfr_shift[!,:shift_8] ≈ dfr_hybrid[!,:shift_8]
@test dfr_shift[!,:shift_17] ≈ dfr_hybrid[!,:shift_17]
@test dfr_shift[!,:sum] ≈ dfr_hybrid[!,:sum]

## Data
dfr = DataFrame(trait = Y, tipNames = ["Ag","Ak","E","M","Az","Ag2","As","Ap","Ar","P","20","165"]) # sim.M.tipNames
dfr = innerjoin(dfr, dfr_hybrid, on=:tipNames)

## Simple BM
fitShift = phylolm(@formula(trait ~ shift_8 + shift_17), dfr, net; reml=false)
@test_logs show(devnull, fitShift)

## Test against fixed values lambda models
fitlam = phylolm(@formula(trait ~ shift_8 + shift_17), dfr, net, model="lambda", fixedValue=1.0,  reml=false)

@test lambda_estim(fitlam) ≈ 1.0
@test coef(fitlam) ≈ coef(fitShift)
@test vcov(fitlam) ≈ vcov(fitShift)
@test nobs(fitlam) ≈ nobs(fitShift)
@test residuals(fitlam) ≈ residuals(fitShift)
@test response(fitlam) ≈ response(fitShift)
@test predict(fitlam) ≈ predict(fitShift)
@test dof_residual(fitlam) ≈ dof_residual(fitShift)
@test sigma2_phylo(fitlam) ≈ sigma2_phylo(fitShift)
@test stderror(fitlam) ≈ stderror(fitShift)
@test confint(fitlam) ≈ confint(fitShift)
@test loglikelihood(fitlam) ≈ loglikelihood(fitShift)
@test dof(fitlam) ≈ dof(fitShift) + 1
@test deviance(fitlam, Val(true))  ≈ deviance(fitShift, Val(true))
@test nulldeviance(fitlam)  ≈ nulldeviance(fitShift)
@test nullloglikelihood(fitlam)  ≈ nullloglikelihood(fitShift)
@test r2(fitlam) ≈ r2(fitShift) atol=1e-15
#@test adjr2(fitlam) ≈ adjr2(fitShift) - 0.5 atol=1e-15
@test aic(fitlam) ≈ aic(fitShift) + 2
#@test aicc(fitlam)  ≈ aicc(fitShift)
@test bic(fitlam) ≈ bic(fitShift) + log(nobs(fitShift))
@test mu_phylo(fitlam)  ≈ mu_phylo(fitShift)
@test hasintercept(fitlam)

fitSH = phylolm(@formula(trait ~ shift_8 + shift_17), dfr, net, model="scalingHybrid", fixedValue=1.0, reml=false)
@test loglikelihood(fitlam) ≈ loglikelihood(fitSH)
@test aic(fitlam) ≈ aic(fitSH)

## ftest against own naive implementation
modnull = phylolm(@formula(trait ~ 1), dfr, net)
@test sigma2_phylo(modnull) ≈ 0.6517876326943942 atol=1e-6 # using REML
modhom = phylolm(@formula(trait ~ sum), dfr, net)
modhet = phylolm(@formula(trait ~ sum + shift_8), dfr, net)
#= 3 warnings thrown by ftest, one for each model, because after transforming the
   data to de-correlate the results, the intercept vector is not ∝ 1.
   Keep the warnings: because incorrect R² values in the ftest output
=#
table1 = redirect_stdio(stderr=devnull) do # to avoid seeing the warnings
    ftest(modhet, modhom, modnull)
end
table2 = PhyloNetworks.anova(modnull, modhom, modhet)

@test table1.fstat[2] ≈ table2[2,:F]
@test table1.fstat[3] ≈ table2[1,:F]
@test table1.pval[2] ≈ table2[2,Symbol("Pr(>F)")]
@test table1.pval[3] ≈ table2[1,Symbol("Pr(>F)")]
@test hasintercept(modnull) && hasintercept(modhom) && hasintercept(modhet)
# ## Replace next 4 lines with previous ones when GLM.ftest available
# @test table1[:F][2] ≈ table2[:F][2]
# @test table1[:F][1] ≈ table2[:F][1]
# @test table1[Symbol("Pr(>F)")][1] ≈ table2[Symbol("Pr(>F)")][1]
# @test table1[Symbol("Pr(>F)")][2] ≈ table2[Symbol("Pr(>F)")][2]

# Check that it is the same as doing shift_8 + shift_17
modhetbis = phylolm(@formula(trait ~ shift_8 + shift_17), dfr, net)
table2bis = PhyloNetworks.anova(modnull, modhom, modhetbis)
@test table2[!,:F] ≈ table2bis[!,:F]
@test table2[!,Symbol("Pr(>F)")] ≈ table2bis[!,Symbol("Pr(>F)")]
@test table2[!,:dof_res] ≈ table2bis[!,:dof_res]
@test table2[!,:RSS] ≈ table2bis[!,:RSS]
@test table2[!,:dof] ≈ table2bis[!,:dof]
@test table2[!,:SS] ≈ table2bis[!,:SS]

# re-fit with ML to do likelihood ratio test
modnull = phylolm(@formula(trait ~ 1), dfr, net; reml=false)
modhom = phylolm(@formula(trait ~ sum), dfr, net; reml=false)
modhet = phylolm(@formula(trait ~ sum + shift_8), dfr, net; reml=false)
table3 = (@test_logs lrtest(modhet, modhom, modnull))
@test all(isapprox.(table3.deviance, (25.10067039653046,47.00501928245542,47.0776339693065), atol=1e-6))
@test table3.dof == (4, 3, 2)
@test all(isapprox.(table3.pval[2:end], (2.865837220526082e-6,0.7875671600772386), atol=1e-6))

end

#################
### No intercept
#################
@testset "No Intercept" begin
global net
net = readTopology("(((Ag:5,(#H1:1::0.056,((Ak:2,(E:1,#H2:1::0.004):1):1,(M:2)#H2:1::0.996):1):1):1,(((((Az:1,Ag2:1):1,As:2):1)#H1:1::0.944,Ap:4):1,Ar:5):1):1,(P:4,20:4):3,165:7);");
preorder!(net)

## data
Y = [11.640085037749985, 9.498284887480622, 9.568813792749083, 13.036916724865296, 6.873936265709946, 6.536647349405742, 5.95771939864956, 10.517318306450647, 9.34927049737206, 10.176238483133424, 10.760099940744308, 8.955543827353837]
X = [9.199418112245104, 8.641506886650749, 8.827105915999073, 11.198420342332025, 5.8212242346434655, 6.130520100788492, 5.846098148463377, 9.125593652542882, 10.575371612483897, 9.198463833849347, 9.090317561636194, 9.603570747653789]
dfr = DataFrame(trait = Y, reg = X, tipNames = ["Ag","Ak","E","M","Az","Ag2","As","Ap","Ar","P","20","165"]) # sim.M.tipNames
phynetlm = phylolm(@formula(trait ~ -1 + reg), dfr, net; reml=false)
# Naive version (GLS): most of it hard-coded, but code shown below
ntaxa = length(Y)
X = phynetlm.X
# Vy = phynetlm.Vy; Vyinv = inv(Vy); XtVyinv = X' * Vyinv; logdetVy = logdet(Vy)
betahat = [1.073805579608655] # inv(XtVyinv * X) * XtVyinv * Y
fittedValues =  X * betahat
resids = Y - fittedValues
# sigma2hat = 1/ntaxa * (resids' * Vyinv * resids)
# loglik = - 1 / 2 * (ntaxa + ntaxa * log(2 * pi) + ntaxa * log(sigma2hat) + logdetVy)
#= null model: no X, and not even an intercept
nullX = zeros(ntaxa, 1); nullresids = Y; nullXtVyinv = nullX' * Vyinv
nullsigma2hat = 1/ntaxa * (nullresids' * Vyinv * nullresids) # 6.666261935713196
nullloglik = - 1 / 2 * (ntaxa + ntaxa * log(2 * pi) + ntaxa * log(nullsigma2hat) + logdetVy) # -38.01145980802529
=#
@test coef(phynetlm) ≈ betahat
@test nobs(phynetlm) ≈ 12 # ntaxa
@test residuals(phynetlm) ≈ resids
@test response(phynetlm) ≈ Y
@test predict(phynetlm) ≈ fittedValues
@test dof_residual(phynetlm) ≈ 11 # ntaxa-length(betahat)
@test sigma2_phylo(phynetlm) ≈ 0.1887449836519979 # sigma2hat
@test loglikelihood(phynetlm) ≈ -16.6249533603196 # loglik
@test vcov(phynetlm) ≈ [0.003054397019042955;;] # sigma2hat*ntaxa/(ntaxa-length(betahat))*inv(XtVyinv * X)
@test stderror(phynetlm) ≈ [0.05526659948868715] # sqrt.(diag(sigma2hat*ntaxa/(ntaxa-length(betahat))*inv(XtVyinv * X)))
@test dof(phynetlm) ≈ 2 # length(betahat)+1
@test deviance(phynetlm, Val(true)) ≈ 2.264939803823975 # sigma2hat * ntaxa
@test nulldeviance(phynetlm) ≈ 79.99514322855836 # nullsigma2hat * ntaxa
@test nullloglikelihood(phynetlm) ≈ -38.01145980802529 # nullloglik
@test r2(phynetlm) ≈ 0.9716865335517596 # 1-sigma2hat / nullsigma2hat
@test adjr2(phynetlm) ≈ 0.9716865335517596  atol=1e-15 # 1 - (1 - (1-sigma2hat/nullsigma2hat))*(ntaxa-1)/(ntaxa-length(betahat))
@test aic(phynetlm) ≈ 37.2499067206392 # -2*loglik+2*(length(betahat)+1)
@test aicc(phynetlm) ≈ 38.58324005397254 # -2*loglik+2*(length(betahat)+1)+2(length(betahat)+1)*((length(betahat)+1)+1)/(ntaxa-(length(betahat)+1)-1)
@test bic(phynetlm) ≈ 38.219720020215206 # -2*loglik+(length(betahat)+1)*log(ntaxa)
@test !hasintercept(phynetlm)

end

###############################################################################
#### Other Network
###############################################################################
@testset "phylolm and ancestralStateReconstruction" begin
global net
# originally: "(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);"
# followed by changes in net.edge[?].length values to make the network ultrametric
net = readTopology("(((Ag:5,(#H1:1::0.056,((Ak:2,(E:1,#H2:1::0.004):1):1,(M:2)#H2:1::0.996):1):1):1,(((((Az:1,Ag2:1):1,As:2):1)#H1:1::0.944,Ap:4):1,Ar:5):1):1,(P:4,20:4):3,165:7);");

#= Simulate correlated data in data frames
b0 = 1
b1 = 10
Random.seed!(5678)
sim = simulate(net, ParamsBM(1, 1))
A = sim[:Tips]
B = b0 .+ b1 * A .+ simulate(net,  ParamsBM(0, 0.1))[:Tips]
tipnam = sim.M.tipNames
=#
A = [2.626609842049044,0.6334773804400937,3.0577676668430476,0.8570052897626761,3.3415290038076875,2.7038939422417467,1.8694860778492748,3.354373836136418,7.436775409527188,2.6659127435884318,3.2298992674067417,-2.2323810599565013]
B = [27.60133970558981,8.228820310098914,32.42043423853238,10.249417359958978,33.52061781961048,27.008691929589997,19.11541648307886,35.38758567184537,75.04861071222199,27.68624399802581,33.03778377357321,-20.4001107607967]
tipnam = ["Ag","Ak","E","M","Az","Ag2","As","Ap","Ar","P","20","165"]

# With Matrices
X = hcat(ones(12), A)
fit_mat = phylolm(X, B, net; reml=false)
#@show fit_mat

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
@test coef(fit_mat) ≈ betahat
@test nobs(fit_mat) ≈ ntaxa
@test residuals(fit_mat) ≈ resids
@test response(fit_mat) ≈ B
@test predict(fit_mat) ≈ fittedValues
@test dof_residual(fit_mat) ≈ ntaxa-length(betahat)
@test sigma2_phylo(fit_mat) ≈ sigma2hat
@test loglikelihood(fit_mat) ≈ loglik
@test vcov(fit_mat) ≈ sigma2hat*ntaxa/(ntaxa-length(betahat)).*inv(XtVyinv * X)
@test stderror(fit_mat) ≈ sqrt.(diag(sigma2hat*ntaxa/(ntaxa-length(betahat)).*inv(XtVyinv * X)))
@test dof(fit_mat) ≈ length(betahat)+1
@test deviance(fit_mat, Val(true)) ≈ sigma2hat * ntaxa
@test nulldeviance(fit_mat) ≈ nullsigma2hat * ntaxa
@test nullloglikelihood(fit_mat) ≈ nullloglik
@test r2(fit_mat) ≈ 1-sigma2hat / nullsigma2hat atol=1e-15
@test adjr2(fit_mat) ≈ 1 - (1 - (1-sigma2hat/nullsigma2hat))*(ntaxa-1)/(ntaxa-length(betahat)) atol=1e-15
@test aic(fit_mat) ≈ -2*loglik+2*(length(betahat)+1)
@test aicc(fit_mat) ≈ -2*loglik+2*(length(betahat)+1)+2(length(betahat)+1)*((length(betahat)+1)+1)/(ntaxa-(length(betahat)+1)-1)
@test bic(fit_mat) ≈ -2*loglik+(length(betahat)+1)*log(ntaxa)

## perfect user using right format and formula
dfr = DataFrame(trait=B, pred=A, tipNames=tipnam)
phynetlm = phylolm(@formula(trait ~ pred), dfr, net; reml=false)
#@show phynetlm

@test coef(phynetlm) ≈ coef(fit_mat)
@test vcov(phynetlm) ≈ vcov(fit_mat)
@test nobs(phynetlm) ≈ nobs(fit_mat)
@test residuals(phynetlm) ≈ residuals(fit_mat)
@test response(phynetlm) ≈ response(fit_mat)
@test predict(phynetlm) ≈ predict(fit_mat)
@test dof_residual(phynetlm) ≈ dof_residual(fit_mat)
@test sigma2_phylo(phynetlm) ≈ sigma2_phylo(fit_mat)
@test stderror(phynetlm) ≈ stderror(fit_mat)
@test confint(phynetlm) ≈ confint(fit_mat)
@test loglikelihood(phynetlm) ≈ loglikelihood(fit_mat)
@test dof(phynetlm) ≈ dof(fit_mat)
@test deviance(phynetlm, Val(true)) ≈ deviance(fit_mat, Val(true))
@test nulldeviance(phynetlm) ≈ nulldeviance(fit_mat)
@test nullloglikelihood(phynetlm) ≈ nullloglikelihood(fit_mat)
@test r2(phynetlm) ≈ r2(fit_mat)
@test adjr2(phynetlm) ≈ adjr2(fit_mat)
@test aic(phynetlm) ≈ aic(fit_mat)
@test aicc(phynetlm) ≈ aicc(fit_mat)
@test bic(phynetlm) ≈ bic(fit_mat)

# Deprecated methods
@test (@test_logs (:warn,r"^accessing") phynetlm.model) === phynetlm
@test (@test_logs (:warn,r"^accessing") phynetlm.mf.f) == formula(phynetlm)
@test (@test_logs (:warn,r"^accessing") phynetlm.mm.m) == modelmatrix(phynetlm)

# unordered data
dfr = dfr[[2,6,10,5,12,7,4,11,1,8,3,9], :]
fitbis = phylolm(@formula(trait ~ pred), dfr, net; reml=false)

@test coef(phynetlm) ≈ coef(fitbis)
@test vcov(phynetlm) ≈ vcov(fitbis)
@test nobs(phynetlm) ≈ nobs(fitbis)
@test residuals(phynetlm)[fitbis.ind] ≈ residuals(fitbis)
@test response(phynetlm)[fitbis.ind] ≈ response(fitbis)
@test predict(phynetlm)[fitbis.ind] ≈ predict(fitbis)
@test dof_residual(phynetlm) ≈ dof_residual(fitbis)
@test sigma2_phylo(phynetlm) ≈ sigma2_phylo(fitbis)
@test stderror(phynetlm) ≈ stderror(fitbis)
@test confint(phynetlm) ≈ confint(fitbis)
@test loglikelihood(phynetlm) ≈ loglikelihood(fitbis)
@test dof(phynetlm) ≈ dof(fitbis)
@test deviance(phynetlm, Val(true)) ≈ deviance(fitbis, Val(true))
@test nulldeviance(phynetlm) ≈ nulldeviance(fitbis)
@test nullloglikelihood(phynetlm) ≈ nullloglikelihood(fitbis)
@test r2(phynetlm) ≈ r2(fitbis)
@test adjr2(phynetlm) ≈ adjr2(fitbis)
@test aic(phynetlm) ≈ aic(fitbis)
@test aicc(phynetlm) ≈ aicc(fitbis)
@test bic(phynetlm) ≈ bic(fitbis)
@test mu_phylo(phynetlm) ≈ mu_phylo(fitbis)

# unnamed ordered data
dfr = DataFrame(trait = B, pred = A)
fitter = (@test_logs (:info, r"^As requested \(no_names=true\)") match_mode=:any phylolm(@formula(trait ~ pred), dfr, net, no_names=true, reml=false))

@test coef(phynetlm) ≈ coef(fitter)
@test vcov(phynetlm) ≈ vcov(fitter)
@test nobs(phynetlm) ≈ nobs(fitter)
@test residuals(phynetlm) ≈ residuals(fitter)
@test response(phynetlm) ≈ response(fitter)
@test predict(phynetlm) ≈ predict(fitter)
@test dof_residual(phynetlm) ≈ dof_residual(fitter)
@test sigma2_phylo(phynetlm) ≈ sigma2_phylo(fitter)
@test stderror(phynetlm) ≈ stderror(fitter)
@test confint(phynetlm) ≈ confint(fitter)
@test loglikelihood(phynetlm) ≈ loglikelihood(fitter)
@test dof(phynetlm) ≈ dof(fitter)
@test deviance(phynetlm, Val(true)) ≈ deviance(fitter, Val(true))
@test nulldeviance(phynetlm) ≈ nulldeviance(fitter)
@test nullloglikelihood(phynetlm) ≈ nullloglikelihood(fitter)
@test r2(phynetlm) ≈ r2(fitter)
@test adjr2(phynetlm) ≈ adjr2(fitter)
@test aic(phynetlm) ≈ aic(fitter)
@test aicc(phynetlm) ≈ aicc(fitter)
@test bic(phynetlm) ≈ bic(fitter)

# unnamed un-ordered data
dfr = dfr[[9,6,5,10,1,11,12,7,2,3,8,4], :]
@test_throws ErrorException phylolm(@formula(trait ~ pred), dfr, net) # Wrong pred

### Add NAs
dfr = DataFrame(trait=B, pred=A, tipNames=tipnam)
allowmissing!(dfr, :pred)
dfr[[2, 8, 11], :pred] .= missing
fitna = phylolm(@formula(trait ~ pred), dfr, net)
#@show fitna

dfr = dfr[[8,2,3,4,6,10,9,1,11,12,7,5], :]
fitnabis = phylolm(@formula(trait ~ pred), dfr, net)
@test coef(fitna) ≈ coef(fitnabis)
@test vcov(fitna) ≈ vcov(fitnabis)
@test nobs(fitna) ≈ nobs(fitnabis)
@test sort(residuals(fitna)) ≈ sort(residuals(fitnabis))
@test sort(response(fitna)) ≈ sort(response(fitnabis))
@test sort(predict(fitna)) ≈ sort(predict(fitnabis))
@test dof_residual(fitna) ≈ dof_residual(fitnabis)
@test sigma2_phylo(fitna) ≈ sigma2_phylo(fitnabis)
@test stderror(fitna) ≈ stderror(fitnabis)
@test confint(fitna) ≈ confint(fitnabis)
@test loglikelihood(fitna) ≈ loglikelihood(fitnabis)
@test dof(fitna) ≈ dof(fitnabis)
@test deviance(fitna, Val(true)) ≈ deviance(fitnabis, Val(true))
@test nulldeviance(fitna) ≈ nulldeviance(fitnabis)
@test (@test_logs (:warn, r"^ML") nullloglikelihood(fitna)) ≈ (@test_logs (:warn, r"^ML") nullloglikelihood(fitnabis))
@test r2(fitna) ≈ r2(fitnabis)
@test adjr2(fitna) ≈ adjr2(fitnabis)
@test aic(fitna) ≈ aic(fitnabis)
@test aicc(fitna) ≈ aicc(fitnabis)
@test bic(fitna) ≈ bic(fitnabis)

## Tests against fixed values parameters
fitlam = phylolm(@formula(trait ~ pred), dfr, net, model="lambda", fixedValue=1.0)
#@show fitlam

@test lambda_estim(fitlam) ≈ 1.0
@test coef(fitlam) ≈ coef(fitnabis)
@test vcov(fitlam) ≈ vcov(fitnabis)
@test nobs(fitlam) ≈ nobs(fitnabis)
@test residuals(fitlam) ≈ residuals(fitnabis)
@test response(fitlam) ≈ response(fitnabis)
@test predict(fitlam) ≈ predict(fitnabis)
@test dof_residual(fitlam) ≈ dof_residual(fitnabis)
@test sigma2_phylo(fitlam) ≈ sigma2_phylo(fitnabis)
@test stderror(fitlam) ≈ stderror(fitnabis)
@test confint(fitlam) ≈ confint(fitnabis)
@test loglikelihood(fitlam) ≈ loglikelihood(fitnabis)
@test dof(fitlam) ≈ dof(fitnabis) + 1
@test deviance(fitlam, Val(true)) ≈ deviance(fitnabis, Val(true))
@test nulldeviance(fitlam) ≈ nulldeviance(fitnabis)
@test (@test_logs (:warn, r"^ML") nullloglikelihood(fitlam)) ≈ (@test_logs (:warn, r"^ML") nullloglikelihood(fitnabis))
@test r2(fitlam) ≈ r2(fitnabis) atol=1e-15
@test adjr2(fitlam)-1 ≈ (adjr2(fitnabis)-1)*(nobs(fitnabis)-dof(fitnabis)+1)/(nobs(fitnabis)-dof(fitlam)+1) atol=1e-15
@test aic(fitlam) ≈ aic(fitnabis) + 2
#@test aicc(fitlam) ≈ aicc(fitnabis)
@test bic(fitlam) ≈ bic(fitnabis) + log(nobs(fitnabis))
@test mu_phylo(fitlam) ≈ mu_phylo(fitnabis)

fitSH = phylolm(@formula(trait ~ pred), dfr, net, model="scalingHybrid", fixedValue=1.0)
@test loglikelihood(fitlam) ≈ loglikelihood(fitSH)
@test aic(fitlam) ≈ aic(fitSH)

## Pagel's Lambda
fitlam = (@test_logs (:info, r"^Maximum lambda value") match_mode=:any phylolm(@formula(trait ~ pred), dfr, net, model="lambda", reml=false))
#@show fitlam
@test lambda_estim(fitlam) ≈ 1.1135518305 atol=1e-6

## scaling Hybrid
fitSH = phylolm(@formula(trait ~ pred), dfr, net, model="scalingHybrid", reml=false)
@test_logs show(devnull, fitSH)
@test lambda_estim(fitSH) ≈ -52.81305448333567 atol=1e-6

### Ancestral State Reconstruction
params = ParamsBM(3, 1)
# sim = simulate(net, params); Y = sim[:Tips]; tipnam=tipLabels(sim)
Y = [7.49814057852738,7.713232061975018,7.4314117011628795,0.9850885689559203,4.970152778471174,5.384066549416034,4.326644522544125,0.6079385242666691,4.084254785718834,5.501648315448596,3.8732700346136597,4.790127215808698]
tipnam = ["Ag","Ak","E","M","Az","Ag2","As","Ap","Ar","P","20","165"]
# From known parameters
ancestral_traits = ancestralStateReconstruction(net, Y, params)
# BLUP
dfr = DataFrame(trait=Y, tipNames=tipnam)
phynetlm = phylolm(@formula(trait~1), dfr, net)
# prediction intervals larger with reml=true than with reml=false
blup = (@test_logs (:warn, r"^These prediction intervals show uncertainty in ancestral values") ancestralStateReconstruction(phynetlm));
# plot(net, blup)
@test_logs show(devnull, blup)

# BLUP same, using the function directly
blup_bis = (@test_logs (:warn, r"^These prediction intervals show uncertainty in ancestral values") match_mode=:any ancestralStateReconstruction(dfr, net));

@test expectations(blup)[!,:condExpectation] ≈ expectations(blup_bis)[!,:condExpectation]
@test expectations(blup)[!,:nodeNumber] ≈ expectations(blup_bis)[!,:nodeNumber]
@test blup.traits_tips ≈ blup_bis.traits_tips
@test blup.TipNumbers ≈ blup_bis.TipNumbers
@test predint(blup) ≈ predint(blup_bis)
@test predintPlot(blup)[!,:PredInt] == predintPlot(blup_bis)[!,:PredInt]
@test predintPlot(blup, withExp=true)[!,:PredInt] == predintPlot(blup_bis, withExp=true)[!,:PredInt]
@test expectationsPlot(blup)[!,:PredInt] == expectationsPlot(blup_bis)[!,:PredInt]

dfr = DataFrame(trait=Y, tipNames=tipnam, reg=Y)
@test_throws ErrorException ancestralStateReconstruction(dfr, net) # cannot handle a predictor

# Unordered
dfr2 = dfr[[5,4,9,2,6,12,8,11,7,1,3,10], :]
phynetlm = phylolm(@formula(trait~1), dfr2, net)
blup2 = (@test_logs (:warn, r"^These prediction intervals show uncertainty in ancestral values") ancestralStateReconstruction(phynetlm))

@test expectations(blup)[1:length(blup.NodeNumbers),:condExpectation] ≈ expectations(blup2)[1:length(blup.NodeNumbers),:condExpectation]
@test blup.traits_tips[phynetlm.ind] ≈ blup2.traits_tips
@test blup.TipNumbers[phynetlm.ind] ≈ blup2.TipNumbers
@test predint(blup)[1:length(blup.NodeNumbers), :] ≈ predint(blup2)[1:length(blup.NodeNumbers), :]

# With unknown tips
allowmissing!(dfr, :trait)
dfr[[2, 4], :trait] .= missing
phynetlm = phylolm(@formula(trait~1), dfr, net)
blup = (@test_logs (:warn, r"^These prediction intervals show uncertainty in ancestral values") ancestralStateReconstruction(phynetlm))
# plot(net, blup)

# Unordered
dfr2 = dfr[[1, 2, 5, 3, 4, 6, 7, 8, 9, 10, 11, 12], :]
phynetlm = phylolm(@formula(trait~1), dfr, net)
blup2 = (@test_logs (:warn, r"^These prediction intervals show uncertainty in ancestral values") ancestralStateReconstruction(phynetlm))

@test expectations(blup)[!,:condExpectation] ≈ expectations(blup2)[!,:condExpectation]
@test predint(blup) ≈ predint(blup2)
@test predintPlot(blup)[!,:PredInt] == predintPlot(blup2)[!,:PredInt]
@test predintPlot(blup, withExp=true)[!,:PredInt] == predintPlot(blup2, withExp=true)[!,:PredInt]

# Test mark on missing
ee = expectationsPlot(blup)
predMiss = ee[indexin([n.number for n in net.leaf][[2,4]], ee[!,:nodeNumber]),:PredInt]
for pp = predMiss
    @test pp[end] == '*'
end

end

#################
## Data with no phylogenetic signal
#################
@testset "lambda when no signal" begin
global net
net = readTopology("(((Ag:5,(#H1:1::0.056,((Ak:2,(E:1,#H2:1::0.004):1):1,(M:2)#H2:1::0.996):1):1):1,(((((Az:1,Ag2:1):1,As:2):1)#H1:1::0.944,Ap:4):1,Ar:5):1):1,(P:4,20:4):3,165:7);");

#= Simulate correlated data in data frames
b0 = 1
b1 = 10
Random.seed!(5678);
A = randn(size(tipLabels(net), 1))
B = b0 .+ (b1 .* A + randn(size(tipLabels(net), 1)))
=#
A = [-1.2217252038914663, 0.8431411538631137, 0.3847679754817904, 0.10277471357263539, 1.0944221266744778, 2.053347250198844, 1.4708882134841876, 1.1056475371071361, -0.94952153892202, -0.3477162381565148, -0.2742415177451819, 0.25034046948064764]
B = [-9.849415384443805, 10.765309004952346, 4.8269904926118565, 1.7279441642635127, 11.535570136728504, 20.16670120778599, 13.971404727143286, 13.019084912634444, -8.278125099304921, -4.784290010378141, -2.537139017477904, 2.9460706727827755]
dfr = DataFrame(trait = B, pred = A, tipNames = tipLabels(net))

## Network
phynetlm = (@test_logs (:info, r"^Maximum lambda value") match_mode=:any phylolm(@formula(trait ~ pred), dfr, net, model="lambda", reml=false))
@test lambda_estim(phynetlm) ≈ 0.5894200143 atol=1e-8
# using REML
phynetlm = (@test_logs (:info, r"^Max") phylolm(@formula(trait ~ pred), dfr, net, model="lambda"))
@test lambda_estim(phynetlm) ≈ 0.8356905283 atol=1e-8

## Major Tree
global tree
tree = majorTree(net)
phynetlm = (@test_logs (:info, r"^Maximum lambda value") match_mode=:any phylolm(@formula(trait ~ pred), dfr, tree, model="lambda", reml=false))

@test lambda_estim(phynetlm) ≈ 0.5903394415 atol=1e-6

## scaling Hybrid
lmtree = phylolm(@formula(trait ~ pred), dfr, tree, model = "BM")
lmnet = phylolm(@formula(trait ~ pred), dfr, net, model = "BM")
lmSHzero = phylolm(@formula(trait ~ pred), dfr, net, model = "scalingHybrid", fixedValue = 0.0)
lmSHone = phylolm(@formula(trait ~ pred), dfr, net, model = "scalingHybrid", fixedValue = 1.0)

@test loglikelihood(lmtree) ≈ loglikelihood(lmSHzero)
@test loglikelihood(lmnet) ≈ loglikelihood(lmSHone)

lmSH = phylolm(@formula(trait ~ pred), dfr, net, model="scalingHybrid", reml=false)
@test lambda_estim(lmSH) ≈ 23.46668204551696 atol=1e-5
lmSH = phylolm(@formula(trait ~ pred), dfr, net, model="scalingHybrid")
@test lambda_estim(lmSH) ≈ 24.61373831478016 atol=1e-5
# λ so large?? largest γ = 0.056, so λγ = 1.34 is > 1...
end

###############################################################################
### Undefined branch lengths
###############################################################################
@testset "Undefined branch length" begin
    ## No branch length
    net = readTopology("(A:2.5,((B,#H1:1::0.1):1,(C:1,(D:1)#H1:1::0.9):1):0.5);");
    dfr = DataFrame(trait = [11.6,8.1,10.3,9.1], tipNames = ["A","B","C","D"]);
    @test_throws ErrorException("""Branch(es) number 2 have no length.
        The variance-covariance matrix of the network is not defined.
        A phylogenetic regression cannot be done.""") phylolm(@formula(trait ~ 1), dfr, net);
    ## Negative branch length
    net.edge[2].length = -0.5;
    @test_throws ErrorException("""Branch(es) number 2 have negative length.
        The variance-covariance matrix of the network is not defined.
        A phylogenetic regression cannot be done.""") phylolm(@formula(trait ~ 1), dfr, net);
    ## Zero branch length: allowed
    net.edge[2].length = 0.0;
    fit = phylolm(@formula(trait ~ 1), dfr, net);
    @test loglikelihood(fit) ≈ -6.245746681512051
    net.edge[2].length = 0.1; # back to non-zero
    ## Illicit zero length for 1st edge: from root to single "outgroup" taxon
    net.edge[1].length = 0.0;
    @test_throws PosDefException phylolm(@formula(trait ~ 1), dfr, net);
end

############################
## Against no regressor
###########################
#= fixit: passes with GML up to v1.3, fails with GLM v1.4. `lm()` has by default
# `allowrankdeficient=false` in v1.3, but `dropcollinear=true` in v1.4
# We would need `dropcollinear=false` in this test. fixit later: pass kwargs... ?
# no predictors, so REML = ML in this case
@testset "phylolm with no regressor" begin
global net
net = readTopology("(((Ag:5,(#H1:1::0.056,((Ak:2,(E:1,#H2:1::0.004):1):1,(M:2)#H2:1::0.996):1):1):1,(((((Az:1,Ag2:1):1,As:2):1)#H1:1::0.944,Ap:4):1,Ar:5):1):1,(P:4,20:4):3,165:7);");

params = ParamsBM(10, 1)
Random.seed!(2468) # sets the seed for reproducibility, to debug potential error
sim = simulate(net, params)
Y = sim[:Tips]
phynetlm = phylolm(zeros(length(Y),0), Y, net)
#@show phynetlm
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

@test nobs(phynetlm) ≈ ntaxa
@test residuals(phynetlm) ≈ resids
@test response(phynetlm) ≈ Y
@test predict(phynetlm) ≈ fittedValues
@test dof_residual(phynetlm) ≈ ntaxa
@test sigma2_phylo(phynetlm) ≈ sigma2hat
@test loglikelihood(phynetlm) ≈ loglik
@test deviance(phynetlm, Val(true)) ≈ sigma2hat * ntaxa
@test nulldeviance(phynetlm) ≈ nullsigma2hat * ntaxa
@test nullloglikelihood(phynetlm) ≈ nullloglik
@test r2(phynetlm) ≈ 1-sigma2hat / nullsigma2hat atol=1e-14
@test adjr2(phynetlm) ≈ 1 - (1 - (1-sigma2hat/nullsigma2hat))*(ntaxa-1)/(ntaxa) atol=1e-14
@test aic(phynetlm) ≈ -2*loglik+2*(1)
@test aicc(phynetlm) ≈ -2*loglik+2*(1)+2(1)*((1)+1)/(ntaxa-(1)-1)
@test bic(phynetlm) ≈ -2*loglik+(1)*log(ntaxa)

# with data frames
dfr = DataFrame(trait = Y, tipNames = sim.M.tipNames)
fitbis = phylolm(@formula(trait ~ -1), dfr, net)
@test_logs show(devnull, fitbis)
#@test coef(phynetlm) ≈ coef(fitbis)
#@test vcov(phynetlm) ≈ vcov(fitbis)
@test nobs(phynetlm) ≈ nobs(fitbis)
@test residuals(phynetlm)[fitbis.ind] ≈ residuals(fitbis)
@test response(phynetlm)[fitbis.ind] ≈ response(fitbis)
@test predict(phynetlm)[fitbis.ind] ≈ predict(fitbis)
@test dof_residual(phynetlm) ≈ dof_residual(fitbis)
@test sigma2_phylo(phynetlm) ≈ sigma2_phylo(fitbis)
#@test stderror(phynetlm) ≈ stderror(fitbis)
#@test confint(phynetlm) ≈ confint(fitbis)
@test loglikelihood(phynetlm) ≈ loglikelihood(fitbis)
#@test dof(phynetlm) ≈ dof(fitbis)
@test deviance(phynetlm, Val(true)) ≈ deviance(fitbis, Val(true))
@test nulldeviance(phynetlm) ≈ nulldeviance(fitbis)
@test nullloglikelihood(phynetlm) ≈ nullloglikelihood(fitbis)
@test r2(phynetlm) ≈ r2(fitbis) atol=1e-15
@test adjr2(phynetlm) ≈ adjr2(fitbis) atol=1e-15
@test aic(phynetlm) ≈ aic(fitbis)
@test aicc(phynetlm) ≈ aicc(fitbis)
@test bic(phynetlm) ≈ bic(fitbis)
#@test mu_phylo(phynetlm)  mu_phylo(fitbis)
end
=#
