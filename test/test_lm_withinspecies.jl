## traits.jl : tools for within-species variation, continuous traits

#= to work on:
- add tests with reml=true with no within-species variation
- use some other network: on a tree (checked with pgls.SEy),
  and on a network with 1+ reticulation.
- check the estimation with independent implementation, like on a tree with
  pgls.SEy in the R package phytools. document here the R code (within comment block)
- check & update functions like: nobs, dof, dof_residual, deviance, residuals,
  sigma2 used to scale confidence intervals, etc.
- add documentation for within-species variation, including assumptions on X
  (no variation in X! but what if variation seen in practice, and what if
  missing in X but not Y for some individuals, or missing Y but not X in other individuals?)
- revise all docstrings: make sure signatures, options lists & jldoctest examples are correct
=#

@testset "phyloNetworklm: within-species variation, star" begin

#= simulation of data used below (test more reproducible when using fixed data)
# simulate traits at the species-level, then
# repeat identical trait values across m individuals
function simTraits(net, m, paramsprocess)
    sim = simulate(net, paramsprocess)
    trait = sim[:Tips] # simple vector now
    return repeat(trait, inner=m)
end

n = 4; m = 3
starnet = readTopology(PhyloNetworks.startree_newick(n))
Random.seed!(6591)
trait1 = simTraits(starnet, m, ParamsBM(2, 0.5)) # var: 0.5
trait2 = simTraits(starnet, m, ParamsBM(-2, 1))

phylo_noise_var = 2; meas_noise_var = 1
phylo_noise = simTraits(starnet, m, ParamsBM(0, phylo_noise_var))
meas_noise = randn(n*m) * sqrt(meas_noise_var)
trait3 = 10 .+ 2*trait1 + phylo_noise + meas_noise
print(round.(trait1, digits=4)) # etc
labels = repeat(starnet.names, inner=m)
df = DataFrame(trait1=trait1, trait2=trait2, trait3=trait3,tipNames=labels)
=#

n = 4; m = 3
starnet = readTopology(PhyloNetworks.startree_newick(n))
netnames = ["t1","t2","t3","t4"] # = tipLabels(starnet) # was generated to have length n
df = DataFrame(
    species = repeat(netnames, inner=m),
    trait1 = [2.8564, missing, 2.8564, 2.8457, 2.8457, 2.8457, 0.4197, 0.4197, 0.4197, 2.2359, 2.2359, 2.2359],
    trait2 = [-2.0935, -2.0935, -2.0935, 0.4955, 0.4955, 0.4955, -2.1977, -2.1977, -2.1977, -2.618, -2.618, -2.618],
    trait3 = [missing, 15.7869, 14.0615, 14.3547, 13.3721, 15.7062, 9.2764, 8.8553, 8.7627, 15.0298, 15.8258, 15.3248]
    # missing replaced 14.0467
)
Y = df[!,:trait3] # nm vector
X = fill(1.0, (n,2)); X[:,2] = df[1:m:n*m,:trait1] # nx2 matrix
#= reduced data: one row per species, extra columns for SD and n of trait 3
gdf = groupby(df, :species)
df_r = combine(gdf, :trait1 => (x -> mean(skipmissing(x))) => :trait1,
                    :trait2 => (x -> mean(skipmissing(x))) => :trait2,
                    :trait3 => (x -> mean(skipmissing(x))) => :trait3,
                    :trait3 => (x ->  std(skipmissing(x))) => :trait3_sd,
                    :trait3 => (x -> sum(.!ismissing.(x))) => :trait3_n)
=#
df_r = DataFrame(
  species = ["t1","t2","t3","t4"], trait1 = [2.8564,2.8457,.4197,2.2359],
  trait2 = [-2.0935,0.4955,-2.1977,-2.618],
  trait3 = [14.0615, 14.477666666666666, 8.9648, 15.393466666666667],
  #        [14.9242,... if we include ind 2
  trait3_sd = [0,1.1718985891848046,.2737966581242371,.4024181076111425],
  # [1.2200420402592675,... if we include ind 2
  trait3_n = [1,3,3,3]) # [2,3,3,3]

#= star tree: R code to check
df = data.frame(species = rep(c("t1","t2","t3","t4"), each=3),
  trait1 = c(2.8564,NA,2.8564,2.8457,2.8457,2.8457,.4197,.4197,.4197,2.2359,2.2359,2.2359),
  trait2 = c(-2.0935,-2.0935,-2.0935,.4955,.4955,.4955,-2.1977,-2.1977,-2.1977, -2.618,-2.618,-2.618),
  trait3 = c(NA,15.7869,14.0615,14.3547,13.3721,15.7062,9.2764,8.8553,8.7627,15.0298,15.8258,15.3248))
mR = lmer(trait3 ~ trait1 + (1|species), df)
fixef(mR) # 8.457020 2.296978
print(mR, digits=7, ranef.comp="Var") # variances: 2.1216068,0.5332865
logLik(mR) # -13.37294
logLik(mR, REML=FALSE) # -14.57265
vcov(mR) # matrix(c(3.111307,-1.219935,-1.219935,0.5913808), nrow=2)
mR = lmer(trait3 ~ trait1 + (1|species), df, REML=FALSE)
fixef(mR) # 8.439909 2.318488
print(mR, digits=7, ranef.comp="Var") # variances: 0.9470427,0.5299540
logLik(mR) # -14.30235
vcov(mR) # matrix(c(1.5237486,-0.6002803,-0.6002803,0.2941625), nrow=2)
=#

#= Alternatively: Julia code to check 
using MixedModels
mm1 = fit(MixedModel, @formula(trait3 ~ trait1 + (1|species)), df, REML=true) # compare with m1
mm3 = fit(MixedModel, @formula(trait3 ~ trait1 + (1|species)), df) # compare with m3
fixef(mm1) # fixed-effect param estimates: 8.45702,2.29698
VarCorr(mm1) # estimated variance-components: species-variance=2.121607, residual-variance=0.533287
# `objective(m)` returns -2 * (log-likelihood of model m)
objective(mm1)/(-2) # -13.3729434, `loglikelihood` not available for models fit by REML
loglikelihood(mm3) # -14.3023458
vcov(mm1) # var-cov matrix for fixef coeffs: [3.11131 -1.21993; -1.21993 0.591381]
=#
m1 = phyloNetworklm(@formula(trait3 ~ trait1), df, starnet; reml=true,
      tipnames=:species, msr_err=true)
m2 = phyloNetworklm(@formula(trait3 ~ trait1), df_r, starnet; reml=true,
      tipnames=:species, msr_err=true, y_mean_std=true)
@test m1.model.reml && m2.model.reml
@test coef(m1) ≈ [8.457020,2.296978] rtol=1e-5 # fixef(mR)
@test coef(m2) ≈ [8.457020,2.296978] rtol=1e-5
@test sigma2_phylo(m1) ≈ 2.1216068 rtol=1e-5 # print(mR, digits=7, ranef.comp="Var")
@test sigma2_phylo(m2) ≈ 2.1216068 rtol=1e-5
@test sigma2_within(m1) ≈ 0.5332865 rtol=1e-5
@test sigma2_within(m2) ≈ 0.5332865 rtol=1e-5
@test loglikelihood(m1) ≈ -13.37294 rtol=1e-5
@test loglikelihood(m2) ≈ -13.37294 rtol=1e-5
@test vcov(m1) ≈ [3.111307 -1.219935; -1.219935 0.5913808] rtol=1e-5
@test vcov(m2) ≈ [3.111307 -1.219935; -1.219935 0.5913808] rtol=1e-5
m3 = phyloNetworklm(@formula(trait3 ~ trait1), df_r, starnet; # reml=false
      tipnames=:species, msr_err=true, y_mean_std=true)
@test !m3.model.reml
@test coef(m3) ≈ [8.439909,2.318488] rtol=1e-5
@test sigma2_phylo(m3) ≈ 0.9470427 rtol=1e-5
@test sigma2_within(m3) ≈ 0.5299540 rtol=1e-5
@test loglikelihood(m3) ≈ -14.30235 rtol=1e-5
@test vcov(m3) ≈ [1.5237486 -0.6002803; -0.6002803 0.2941625] rtol=1e-5

end

@testset "phyloNetworklm: reml=true/false, msr_err=false, binary tree" begin

#= Rcode to generate the newick string and dataset:
library(phytools)
set.seed(1)

tree <- pbtree(n=5,scale=1) # topology and branch lengths
tree.newick <- write.tree(tree) # newick format
# tree.newick <- "((t1:0.8121974445,(t2:0.4806586387,t3:0.4806586387):0.3315388057):0.1878025555,(t4:0.1206907041,t5:0.1206907041):0.8793092959);"

X <- fastBM(tree,sig2=2,nsim=2); colnames(X) <- c("x1","x2") # non-intercept predictors
# X <- matrix(c(-0.5503706,-0.7503593,-0.7599224,-0.4217977,-0.9839226,-1.3102488,0.4544195,0.2151272,1.0384023,1.3138847),nrow=5)
# colnames(X) <- c("x1","x2"); rownames(X) <- c("t1","t2","t3","t4","t5")

bsperr <- fastBM(tree) # phylogenetic variation, true BM variance-rate is 1
y <- cbind(rep(1,5),X)%*%(1:3)+bsperr # response, true beta vector is c(1,2,3)
# y <- matrix(c(-4.1836578,0.6258014,-0.3070109,2.6120026,2.2391044),nrow=5)
# rownames(y) <- c("t1","t2","t3","t4","t5")

df <- data.frame(y,X,species=tree$tip.label)
=#
net = readTopology("((t1:0.8121974445,(t2:0.4806586387,t3:0.4806586387):0.3315388057):0.1878025555,(t4:0.1206907041,t5:0.1206907041):0.8793092959);")
df = DataFrame(
      y = [-4.18366,0.625801,-0.307011,2.612,2.2391],
      x1 = [-0.550371,-0.750359,-0.759922,-0.421798,-0.983923],
      x2 = [-1.31025,0.45442,0.215127,1.0384,1.31388],
      species = net.names # ["t1","t2","t3","t4","t5"]
)

#= Rcode to check model fit:
library(nlme)
m1 <- gls(y~x1+x2,data=df,
        correlation=corBrownian(1,tree,form=~species),method="REML")
coef(m1) # coefficients estimates: c(0.6469652,2.0420889,2.8285257)
sigma(m1)^2 # reml BM variance-rate estimate: 0.0438973
logLik(m1) # restricted log-likelihood: -0.07529961
vcov(m1) # cov mat for coef estimates
# matrix(c(0.030948901,0.024782246,0.003569513,0.024782246,0.039001092,0.009137043,0.003569513,0.009137043,0.013366014),nrow=3)
summary(m1)$tTable[,"t-value"] # t-values for coef estimates: c(3.677548,10.340374,24.465786)
summary(m1)$tTable[,"p-value"] # p-values for coef estimates: c(0.066635016,0.009223303,0.001666460)
coef(m1) + qt(p=0.025,df=5-3)*sqrt(diag(vcov(m1))) # lower limit 95%-CI: c(-0.1099703,1.1923712,2.3310897)
coef(m1) + qt(p=0.975,df=5-3)*sqrt(diag(vcov(m1))) # upper limit 95%-CI: c(1.403901,2.891807,3.325962)

m2 <- gls(y~x1+x2,data=df,
        correlation=corBrownian(1,tree,form=~species),method="ML")
coef(m2) # coefficients estimates: c(0.6469652,2.0420889,2.8285257)
sigma(m2)^2 # ml BM variance-rate estimate: 0.01755892
logLik(m2) # log-likelihood: 3.933531
vcov(m2) # cov mat for coef estimates, this evaluates to the same value as vcov(m1)
# Note: vcov(gls(...,method="REML)) == vcov(gls(...,method="ML"))
# The same holds for t/p-values and CIs for the coef estimates
summary(m2)$tTable[,"t-value"] # t-values for coef estimates: c(3.677548,10.340374,24.465786)
summary(m2)$tTable[,"p-value"] # p-values for coef estimates: c(0.066635016,0.009223303,0.001666460)
coef(m2) + qt(p=0.025,df=5-3)*sqrt(diag(vcov(m2))) # lower limit 95%-CI: c(-0.1099703,1.1923712,2.3310897)
coef(m2) + qt(p=0.975,df=5-3)*sqrt(diag(vcov(m2))) # upper limit 95%-CI: c(1.403901,2.891807,3.325962)
=#
m1 = phyloNetworklm(@formula(y~x1+x2),df,net;tipnames=:species,reml=true)
m2 = phyloNetworklm(@formula(y~x1+x2),df,net;tipnames=:species,reml=false)
@test m1.model.reml
@test coef(m1) ≈ [0.6469652,2.0420889,2.8285257] rtol=1e-5
@test sigma2_phylo(m1) ≈ 0.0438973 rtol=1e-4
@test isnothing(sigma2_within(m1))
@test loglikelihood(m1) ≈ -0.07529961 rtol=1e-3
@test vcov(m1) ≈ [0.030948901 0.024782246 0.003569513;0.024782246 0.039001092 0.009137043;0.003569513 0.009137043 0.013366014] rtol=1e-4
@test coeftable(m1).cols[coeftable(m1).teststatcol] ≈ [3.677548,10.340374,24.465786] rtol=1e-4
@test coeftable(m1).cols[coeftable(m1).pvalcol] ≈ [0.066635016,0.009223303,0.001666460] rtol=1e-4
@test coeftable(m1).cols[findall(coeftable(m1).colnms .== "Lower 95%")[1]] ≈ [-0.1099703,1.1923712,2.3310897] rtol=1e-5
@test coeftable(m1).cols[findall(coeftable(m1).colnms .== "Upper 95%")[1]] ≈ [1.403901,2.891807,3.325962] rtol=1e-5
@test !m2.model.reml
@test coef(m2) ≈ [0.6469652,2.0420889,2.8285257] rtol=1e-5
@test sigma2_phylo(m2) ≈ 0.01755892 rtol=1e-4
@test isnothing(sigma2_within(m2))
@test loglikelihood(m2) ≈ 3.933531 rtol=1e-4
@test vcov(m2) ≈ [0.030948901 0.024782246 0.003569513;0.024782246 0.039001092 0.009137043;0.003569513 0.009137043 0.013366014] rtol=1e-4
@test coeftable(m2).cols[coeftable(m2).teststatcol] ≈ [3.677548,10.340374,24.465786] rtol=1e-4
@test coeftable(m2).cols[coeftable(m2).pvalcol] ≈ [0.066635016,0.009223303,0.001666460] rtol=1e-4
@test coeftable(m2).cols[findall(coeftable(m2).colnms .== "Lower 95%")[1]] ≈ [-0.1099703,1.1923712,2.3310897] rtol=1e-5
@test coeftable(m2).cols[findall(coeftable(m2).colnms .== "Upper 95%")[1]] ≈ [1.403901,2.891807,3.325962] rtol=1e-5

end

@testset "phyloNetworklm: reml=true/false, msr_err=true, binary tree" begin

#= NOTE: 
This testset DOES NOT check for similarity of the variance-components estimates
from 'pgls.SEy' (phytools) and 'phyloNetworklm' (PhyloNetworks). In fact, 
'pgls.SEy' and 'phyloNetworklm' solve different problems! 
ML/REML variance-component estimates are obtained by running 'phyloNetworklm'.
The code for 'pgls.SEy' is adapted to calculate coef estimates, cov matrix for
coef estimates, t/p-values for coef estimates, and full/restricted-loglikelihood, 
given these variance-components estimates. These values are then compared against
those extracted from the output of 'phyloNetworklm'.
=#

#= Rcode to generate the newick string and dataset:
library(phytools)
set.seed(2)
m <- 5 # sample-size per taxa
n <- 5 # no. of taxa
p <- 3 # no. of predictors (including intercept)

tree <- pbtree(n=n,scale=1) # topology and branch lengths
tree.newick <- write.tree(tree) # newick format
# tree.newick <- "(((t4:0.2537636499,t5:0.2537636499):0.2103870459,t3:0.4641506959):0.5358493041,(t1:0.4807642475,t2:0.4807642475):0.5192357525);"

X <- fastBM(tree,sig2=2,nsim=2); colnames(X) <- c("x1","x2") # non-intercept predictors
# X <- matrix(c(0.7894387,0.3285286,2.1495131,-1.4935665,-2.3199935,2.8184268,0.4740687,2.5801004,1.9967379,-0.4121874),nrow=5)
# colnames(X) <- c("x1","x2"); rownames(X) <- c("t4","t5","t3","t1","t2")

bsperr <- fastBM(tree) # phylogenetic variation, true BM variance-rate is 1
wspvar <- 0.1 # measurement-error variance is 0.1
msrerr <- rnorm(n=n,sd=sqrt(wspvar/m)) # msr error in species-level mean responses
y_sd <- sqrt(wspvar*rchisq(n=n,df=m-1)/(m-1)) # sd in individual responses per species 
# y_sd <- c(0.2463003,0.3236629,0.2458547,0.4866844,0.3434582); 
y <- cbind(rep(1,5),X)%*%(1:3)+bsperr+msrerr # response, true beta vector is c(1,2,3)
# y <- matrix(c(11.399108,3.216645,13.648011,4.851454,-4.922803),nrow=5); rownames(y) <- c("t4","t5","t3","t1","t2")

df <- data.frame(y,y_sd,X,species=tree$tip.label)
=#
net = readTopology("(((t4:0.2537636499,t5:0.2537636499):0.2103870459,t3:0.4641506959):0.5358493041,(t1:0.4807642475,t2:0.4807642475):0.5192357525);")
df = DataFrame(
      y=[11.399108,3.216645,13.648011,4.851454,-4.922803],
      y_sd=[0.2463003,0.3236629,0.2458547,0.4866844,0.3434582],
      y_n=fill(5,5),
      x1=[0.7894387,0.3285286,2.1495131,-1.4935665,-2.3199935],
      x2=[2.8184268,0.4740687,2.5801004,1.9967379,-0.4121874],
      species=net.names
)

#= R/Julia code to check model fit:
(1) Julia code:
## To extract reml/ml msrerr variance estimates fitted by phyloNetworklm
m1 |> sigma2_within |> x -> round(x,sigdigits=6) # reml est msrerr var: 0.116721
m2 |> sigma2_within |> x -> round(x,sigdigits=6) # ml est msrerr var: 0.125628
m1 |> sigma2_phylo |> x -> round(x,sigdigits=6) # reml est BM var: 0.120746
m2 |> sigma2_phylo |> x -> round(x,sigdigits=6) # ml est BM var: 0.000399783

(2) R code:
## To check phyloNetworklm REML fit 
library(nlme)
T1 <- tree
wspvar1 <- 0.116721 # reml est msrerr var of indiv-lvl rsps 
bspvar1 <- 0.120746 # reml est BM var
se1 <- setNames(rep(sqrt(wspvar1/m),n),T1$tip.label) # msrerr sd of species-lvl mean rsp
T1$edge.length <- T1$edge.length*bspvar1 # scale all edges by est BM var
ii <- sapply(1:Ntip(T1),function(x,e) which(e==x),e=T1$edge[,2]) # indices of pendant edges
# extend pendant edges by msrerr sd  
T1$edge.length[ii] <- T1$edge.length[ii]+se1[T1$tip.label]^2
covmat1 <- vcv(T1) # extract est vars of species-lvl mean rsp
m1 <- gls(y~x1+x2,data=cbind(df,vf=diag(covmat1)),
          correlation=corBrownian(1,T1,form=~species),
          method="REML",
          weights=varFixed(~vf))

coef(m1) # reml coefficient estimates: c(1.079839,1.976719,3.217391)

Xp <- model.matrix(m1,df) # predictor matrix
vcov1 <- solve(t(Xp)%*%solve(covmat1)%*%Xp) # cov mat for coef estimates
# matrix(c(0.09386431,0.02273458,-0.02602937,0.02273458,0.02172123,-0.01133032,-0.02602937,-0.01133032,0.01584198),nrow=3)

teststat1 <- coef(m1)/sqrt(diag(vcov1)) # test stat for coef estimates: c(3.524591,13.412285,25.562252)
pval1 <- sapply(X=pt(teststat1,n-p),FUN=function(p) 2*min(p,1-p)) # pval for coef est:
# c(0.071921276,0.005513044,0.001526885)
lowerci1 <- coef(m1) + qt(p=0.025,df=n-p)*sqrt(diag(vcov1)) # lower limit 95%-CI: c(-0.2383769,1.3425890,2.6758379)
upperci1 <- coef(m1) + qt(p=0.975,df=n-p)*sqrt(diag(vcov1)) # upper limit 95%-CI: c(2.398055,2.610850,3.758944)

RSS <- sum((m-1)*(df$y_sd^2)) # residual sum-of-squares wrt to the species means
logLik(m1) # species-lvl cond restricted-ll: -3.26788
sigm1 <- sigma(m1) # this is not bspvar1!, but rather the best "scaling" for covmat1
# indiv-lvl cond restricted ll 
# "+ (n-p)*(2*log(sigm1)-sigm1^2+1)/2" un-scales the species-lvl rll returned by logLik
# "- n*(m-1)*(log(wspvar1)+log(2*pi))/2 - (n*log(m)+RSS/wspvar1)/2" corrects to indiv-lvl rll
rll.species <- logLik(m1)
rll.indiv <- (rll.species
              + (n-p)*(2*log(sigm1)-sigm1^2+1)/2
              - n*(m-1)*(log(wspvar1)+log(2*pi))/2 - (n*log(m)+RSS/wspvar1)/2
              )
round(rll.indiv,digits=6) # indiv-lvl rll of remles: -14.14184

## To check phyloNetworklm ML fit
T2 <- tree
wspvar2 <- 0.125628 # ml est msrerr var of indiv-lvl rsps 
bspvar2 <- 0.000399783 # ml est BM var
se2 <- setNames(rep(sqrt(wspvar2/m),n),T2$tip.label)
T2$edge.length <- T2$edge.length*bspvar2
ii <- sapply(1:Ntip(T2),function(x,e) which(e==x),e=T2$edge[,2])  
T2$edge.length[ii] <- T2$edge.length[ii]+se2[T2$tip.label]^2
covmat2 <- vcv(T2)
m2 <- gls(y~x1+x2,data=cbind(df,vf=diag(covmat2)),
          correlation=corBrownian(1,T2,form=~species),
          method="ML",
          weights=varFixed(~vf))

coef(m2) # ml coefficient estimates: c(0.9767352,1.9155142,3.2661862)

Xp <- model.matrix(m2,df) # predictor matrix
vcov2 <- solve(t(Xp)%*%solve(covmat2)%*%Xp) # cov mat for coef estimates
# matrix(c(0.019118171,0.004924625,-0.008977541,0.004924625,0.003577654,-0.003028413,-0.008977541,-0.003028413,0.005793475),nrow=3)

teststat2 <- coef(m2)/sqrt(diag(vcov2)) # test stat for coef est: c(7.064049,32.024784,42.911270)
pval2 <- sapply(X=pt(teststat2,n-p),FUN=function(p) 2*min(p,1-p)) # pval for coef est:
# c(0.0194568147,0.0009736278,0.0005426298)
lowerci2 <- coef(m2) + qt(p=0.025,df=n-p)*sqrt(diag(vcov2)) # lower limit 95%-CI: c(0.381814,1.658158,2.938690)
upperci2 <- coef(m2) + qt(p=0.975,df=n-p)*sqrt(diag(vcov2)) # upper limit 95%-CI: c(1.571656,2.172871,3.593682)

ll.species <- logLik(m2) # species-lvl cond ll: 1.415653
sigm2 <- sigma(m2) # this is not bspvar2!, but rather the best "scaling" for covmat2 
RSS <- sum((m-1)*(df$y_sd^2)) # residual sum-of-squares wrt to the species means
# indiv-lvl cond ll
# "+ n*(2*log(sigm2)-sigm2^2+1)/2" un-scales the species-lvl ll returned by logLik
# "- n*(m-1)*(log(wspvar2)+log(2*pi))/2 - (n*log(m)+RSS/wspvar2)/2" corrects to indiv-lvl ll
ll.indiv <- (ll.species 
             + n*(2*log(sigm2)-sigm2^2+1)/2 
             - n*(m-1)*(log(wspvar2)+log(2*pi))/2 - (n*log(m)+RSS/wspvar2)/2
             )
round(ll.indiv,digits=6) # indiv-lvl ll of mles: -9.582357
=#
m1 = phyloNetworklm(@formula(y~x1+x2),df,net;
                    tipnames=:species,reml=true,msr_err=true,y_mean_std=true)
m2 = phyloNetworklm(@formula(y~x1+x2),df,net;
                    tipnames=:species,reml=false,msr_err=true,y_mean_std=true)

@test coef(m1) ≈ [1.079839,1.976719,3.217391] rtol=1e-4
@test coef(m2) ≈ [0.9767352,1.9155142,3.2661862] rtol=1e-4
@test loglikelihood(m1) ≈ -14.14184 rtol=1e-4
@test loglikelihood(m2) ≈ -9.582357 rtol=1e-4
@test vcov(m1) ≈ [0.09386431 0.02273458 -0.02602937;0.02273458 0.02172123 -0.01133032;-0.02602937 -0.01133032 0.01584198] rtol=1e-4
@test vcov(m2) ≈ [0.019118171 0.004924625 -0.008977541;0.004924625 0.003577654 -0.003028413;-0.008977541 -0.003028413 0.005793475] rtol=1e-4
@test coeftable(m1).cols[coeftable(m1).teststatcol] ≈ [3.524591,13.412285,25.562252] rtol=1e-4
@test coeftable(m2).cols[coeftable(m2).teststatcol] ≈ [7.064049,32.024784,42.911270] rtol=1e-4
@test coeftable(m1).cols[coeftable(m1).pvalcol] ≈ [0.071921276,0.005513044,0.001526885] rtol=1e-4
@test coeftable(m2).cols[coeftable(m2).pvalcol] ≈ [0.0194568147,0.0009736278,0.0005426298] rtol=1e-4
@test coeftable(m1).cols[findall(coeftable(m1).colnms .== "Lower 95%")[1]] ≈ [-0.2383769,1.3425890,2.6758379] rtol=1e-4
@test coeftable(m2).cols[findall(coeftable(m2).colnms .== "Lower 95%")[1]] ≈ [0.381814,1.658158,2.938690] rtol=1e-4
@test coeftable(m1).cols[findall(coeftable(m1).colnms .== "Upper 95%")[1]] ≈ [2.398055,2.610850,3.758944] rtol=1e-4
@test coeftable(m2).cols[findall(coeftable(m2).colnms .== "Upper 95%")[1]] ≈ [1.571656,2.172871,3.593682] rtol=1e-4
end

@testset "phyloNetworklm: within-species variation, network (1 reticulation)" begin

#= Simulation of data used below:
(1) Individual-level data
function simTraits(net, m, paramsprocess)
      sim = simulate(net, paramsprocess)
      trait = sim[:Tips]
      return repeat(trait, inner=m)
end
n = 6; m = 3 # 6 species, 3 individuals per species 
net = readTopology("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
Random.seed!(18480224);
trait1 = simTraits(net, m, ParamsBM(2, 0.5)) # simulate a BM with mean 2 and variance 0.5 on net
trait2 = simTraits(net, m, ParamsBM(-2, 1.0)) # simulate a BM with mean -2 and variance 1.0 on net
phylo_noise_var = 0.1 # BM variance-rate
meas_noise_var = 0.01 # individual-level msrerr variance 
phylo_noise = simTraits(net, m, ParamsBM(0, phylo_noise_var))
meas_noise = randn(n*m)*sqrt(meas_noise_var)
trait3 = 10 .+ 2 * trait1 + phylo_noise + meas_noise
labels = repeat(names(vcv(net)), inner=m)
df = DataFrame(trait1=trait1, trait2=trait2, trait3=trait3, tipNames=labels)

(2) Species-level data (extra columns for SD and n of trait3)
gdf = groupby(df, :species)
df_r = combine(gdf, :trait1 => (x -> mean(x)) => :trait1,
                    :trait2 => (x -> mean(x)) => :trait2,
                    :trait3 => (x -> mean(x)) => :trait3,
                    :trait3 => (x -> std(x)) => :trait3_sd,
                    :trait3 => (x -> length(x)) => :trait3_n)
=#
n = 6; m = 3
net = readTopology("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
df = DataFrame(
      species = repeat(["D","C","A","B","E","O"],inner=m),
      trait1 = [4.08298,4.08298,4.08298,3.10782,3.10782,3.10782,2.17078,2.17078,2.17078,1.87333,1.87333,1.87333,2.8445,
                2.8445,2.8445,5.88204,5.88204,5.88204],
      trait2 = [-7.34186,-7.34186,-7.34186,-7.45085,-7.45085,-7.45085,-3.32538,-3.32538,-3.32538,-4.26472,-4.26472,
                -4.26472,-5.96857,-5.96857,-5.96857,-1.99388,-1.99388,-1.99388],
      trait3 = [18.8101,18.934,18.9438,17.0687,17.0639,17.0732,14.4818,14.1112,14.2817,13.0842,12.9562,12.9019,15.4373,
                15.4075,15.4317,24.2249,24.1449,24.1302]
)
df_r = DataFrame(
      species = ["D","C","A","B","E","O"],
      trait1 = [4.08298,3.10782,2.17078,1.87333,2.8445,5.88204],
      trait2 = [-7.34186,-7.45085,-3.32538,-4.26472,-5.96857,-1.99388],
      trait3 = [18.896,17.0686,14.2916,12.9808,15.4255,24.1667],
      trait3_sd = [0.074524,0.00465081,0.185497,0.0936,0.0158379,0.0509643],
      trait3_n = [3, 3, 3, 3, 3, 3]
)

# Check for agreement in fit when individual-level data is supplied vs when the
# corresponding species-level data is supplied.
m1 = phyloNetworklm(@formula(trait3 ~ trait1), df, net; reml=true, 
      tipnames=:species, msr_err=true)
m2 = phyloNetworklm(@formula(trait3 ~ trait1), df_r, net; reml=true,
      tipnames=:species, msr_err=true, y_mean_std=true)
@test coef(m1) ≈ [9.65347,2.30357] rtol=1e-4
@test coef(m2) ≈ [9.65347,2.30357] rtol=1e-4
@test sigma2_phylo(m1) ≈ 0.156188 rtol=1e-4
@test sigma2_phylo(m2) ≈ 0.156188 rtol=1e-4
@test sigma2_within(m1) ≈ 0.008634 rtol=1e-4
@test sigma2_within(m2) ≈ 0.008634 rtol=1e-4
@test loglikelihood(m1) ≈ 1.944626 rtol=1e-4
@test loglikelihood(m2) ≈ 1.944626 rtol=1e-4

m3 = phyloNetworklm(@formula(trait3 ~ trait1), df, net; reml=false, 
      tipnames=:species, msr_err=true)
m4 = phyloNetworklm(@formula(trait3 ~ trait1), df_r, net; reml=false,
      tipnames=:species, msr_err=true, y_mean_std=true)
@test coef(m3) ≈ [9.63523,2.30832] rtol=1e-4
@test coef(m4) ≈ [9.63523,2.30832] rtol=1e-4
@test sigma2_phylo(m3) ≈ 0.102255 rtol=1e-4
@test sigma2_phylo(m4) ≈ 0.102255 rtol=1e-4
@test sigma2_within(m3) ≈ 0.008677 rtol=1e-4
@test sigma2_within(m4) ≈ 0.008677 rtol=1e-4
@test loglikelihood(m3) ≈ 1.876606 rtol=1e-4
@test loglikelihood(m4) ≈ 1.876606 rtol=1e-4
end