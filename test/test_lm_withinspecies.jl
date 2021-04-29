## traits.jl : phylolm for within-species variation, continuous traits

@testset "phylolm: within-species variation, star" begin

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
m1 = phylolm(@formula(trait3 ~ trait1), df, starnet; reml=true,
      tipnames=:species, withinspecies_var=true)
m2 = phylolm(@formula(trait3 ~ trait1), df_r, starnet; reml=true,
      tipnames=:species, withinspecies_var=true, y_mean_std=true)
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
m3 = phylolm(@formula(trait3 ~ trait1), df_r, starnet; reml=false,
      tipnames=:species, withinspecies_var=true, y_mean_std=true)
@test !m3.model.reml
@test coef(m3) ≈ [8.439909,2.318488] rtol=1e-5
@test sigma2_phylo(m3) ≈ 0.9470427 rtol=1e-5
@test sigma2_within(m3) ≈ 0.5299540 rtol=1e-5
@test loglikelihood(m3) ≈ -14.30235 rtol=1e-5
@test vcov(m3) ≈ [1.5237486 -0.6002803; -0.6002803 0.2941625] rtol=1e-5

end

@testset "phylolm: binary tree, no withinspecies var" begin

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
m1 = phylolm(@formula(y~x1+x2),df,net; tipnames=:species, reml=true)
m2 = phylolm(@formula(y~x1+x2),df,net; tipnames=:species, reml=false)
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

@testset "phylolm: binary tree, withinspecies var" begin

#= NOTE:
This testset DOES NOT check for similarity of the variance-components estimates
from 'pgls.SEy' (phytools) and 'phylolm' (PhyloNetworks). In fact,
'pgls.SEy' and 'phylolm' solve different problems!
ML/REML variance-component estimates are obtained by running 'phylolm'.
The code for 'pgls.SEy' is adapted to calculate coef estimates, cov matrix for
coef estimates, t/p-values for coef estimates, and full/restricted-loglikelihood,
given these variance-components estimates. These values are then compared against
those extracted from the output of 'phylolm'.
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
wspvar <- 0.1 # within-species variance is 0.1
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
## To extract reml/ml msrerr variance estimates fitted by phylolm
m1 |> sigma2_within |> x -> round(x,sigdigits=6) # reml est msrerr var: 0.116721
m2 |> sigma2_within |> x -> round(x,sigdigits=6) # ml est msrerr var: 0.125628
m1 |> sigma2_phylo |> x -> round(x,sigdigits=6) # reml est BM var: 0.120746
m2 |> sigma2_phylo |> x -> round(x,sigdigits=6) # ml est BM var: 0.000399783

(2) R code:
## To check phylolm REML fit
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

## To check phylolm ML fit
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
m1 = phylolm(@formula(y~x1+x2),df,net;
      tipnames=:species, reml=true, withinspecies_var=true, y_mean_std=true)
m2 = phylolm(@formula(y~x1+x2),df,net;
      tipnames=:species, reml=false, withinspecies_var=true, y_mean_std=true)

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

# missing y (but y_n/SD present); missing/infinite sample size; missing/infinite y_SD
df = df[[1,2,3,4],:] # delete t2
allowmissing!(df, [:y,:y_n,:y_sd])
# m0 = phylolm(@formula(y~1),df,net;  tipnames=:species, withinspecies_var=true, y_mean_std=true)
df1 = DataFrame(df) # makes a copy
push!(df1, (missing,0.,0,0.,0.,"t2"))
m1 = phylolm(@formula(y~1),df1,net; tipnames=:species, withinspecies_var=true, y_mean_std=true)
@test coef(m1) ≈ [7.8116494901242195] atol=1e-6 # same as m0
@test sigma2_within(m1) ≈ 0.11568962074173368 atol=1e-6
@test nobs(m1) == 20        # 4 species x 5 indiv/species
@test dof_residual(m1) == 3
# sample sizes y_n: 1, 0, infinite
@test dof(m1) == 3 # intercept, s2phylo, s2within
df1[5,:] .= (0.,0.,1,0.,0.,"t2") # then almost same within-sp variation
m1 = phylolm(@formula(y~1),df1,net; tipnames=:species, withinspecies_var=true, y_mean_std=true)
@test sigma2_within(m1) ≈ 0.11569 atol=1e-4
@test nobs(m1) == 21
@test dof_residual(m1) == 4
df1[5,:] .= (0.,0.,0,0.,0.,"t2") # some species have 0 or <0 num_individuals, column y_n
@test_throws ErrorException phylolm(@formula(y~1),df1,net; tipnames=:species, withinspecies_var=true, y_mean_std=true)
df1[!,:y_n] = Float64.(df1[!,:y_n])
df1[5,:] .= (0.,0.,Inf,0.,0.,"t2")#some species have infinite num_individuals, column y_n
@test_throws ErrorException phylolm(@formula(y~1),df1,net; tipnames=:species, withinspecies_var=true, y_mean_std=true)
df1[5,:] .= (0.,missing,1,0.,0.,"t2") # some SD values are missing, column y_sd
@test_throws ErrorException phylolm(@formula(y~1),df1,net; tipnames=:species, withinspecies_var=true, y_mean_std=true)
df1[5,:] .= (0.0,Inf,1,0.,0.,"t2")   # some SD values are infinite, column y_sd
@test_throws ErrorException phylolm(@formula(y~1),df1,net; tipnames=:species, withinspecies_var=true, y_mean_std=true)

end

@testset "phylolm: within-species var, network h=1" begin

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
      trait3 = [18.895966666666666,17.0686,14.291566666666668,12.980766666666666,15.4255,24.166666666666668],
      trait3_sd = [0.07452397824414288,0.004650806381693211,0.18549690922851855,0.09360001780626576,0.01583792915756317,0.05096433393397276],
      trait3_n = [3, 3, 3, 3, 3, 3]
)

m1 = phylolm(@formula(trait3 ~ trait1), df, net; reml=true,
      tipnames=:species, withinspecies_var=true)
m2 = phylolm(@formula(trait3 ~ trait1), df_r, net; reml=true,
      tipnames=:species, withinspecies_var=true, y_mean_std=true)
m3 = phylolm(@formula(trait3 ~ trait1), df, net; reml=false,
      tipnames=:species, withinspecies_var=true)
m4 = phylolm(@formula(trait3 ~ trait1), df_r, net; reml=false,
      tipnames=:species, withinspecies_var=true, y_mean_std=true)

@testset "agreement across data input" begin # nested test set
#= - individual-level data vs when the species-level data: see m1, m2, m3, m4
- permuted data rows: see m1permute, m2permute
- w/o withinspecies var: when a species has missing entries vs when the species
  is completely absent from the data (see m3missingentry, m4missingrow). =#
m1permute = phylolm(@formula(trait3 ~ trait1), df[[1,6,11,17,16,18,8,5,9,3,12,7,13,10,2,14,4,15],:], # permute predictor rows
                  net; reml=true,tipnames=:species, withinspecies_var=true)
m2permute = phylolm(@formula(trait3 ~ trait1), df_r[[5,2,3,1,6,4],:], net; reml=true, # permute predictor rows
                  tipnames=:species, withinspecies_var=true, y_mean_std=true) 
@test coef(m1) ≈ [9.65347,2.30357] rtol=1e-4
@test coef(m1permute) ≈ [9.65347,2.30357] rtol=1e-4
@test coef(m2) ≈ [9.65347,2.30357] rtol=1e-4
@test coef(m2permute) ≈ [9.65347,2.30357] rtol=1e-4
@test sigma2_phylo(m1) ≈ 0.156188 rtol=1e-4
@test sigma2_phylo(m1permute) ≈ 0.156188 rtol=1e-4
@test sigma2_phylo(m2) ≈ 0.156188 rtol=1e-4
@test sigma2_phylo(m2permute) ≈ 0.156188 rtol=1e-4
@test sigma2_within(m1) ≈ 0.008634 rtol=1e-4
@test sigma2_within(m1permute) ≈ 0.008634 rtol=1e-4
@test sigma2_within(m2) ≈ 0.008634 rtol=1e-4
@test sigma2_within(m2permute) ≈ 0.008634 rtol=1e-4
@test loglikelihood(m1) ≈ 1.944626 rtol=1e-4
@test loglikelihood(m1permute) ≈ 1.944626 rtol=1e-4
@test loglikelihood(m2) ≈ 1.944626 rtol=1e-4
@test loglikelihood(m2permute) ≈ 1.944626 rtol=1e-4
@test stderror(m1) ≈ [1.3065987324433421,0.27616258597477233] atol=1e-6
@test stderror(m1permute) ≈ [1.3065987324433421,0.27616258597477233] atol=1e-6
@test stderror(m2) ≈ [1.3065987324433421,0.27616258597477233] atol=1e-6
@test stderror(m2permute) ≈ [1.3065987324433421,0.27616258597477233] atol=1e-6
@test nobs(m1) == 18 # 6 species, 18 ind
@test nobs(m1permute) == 18
@test residuals(m1) ≈ [-0.16295769587309603,0.2560302141911026,-0.36246083625543507,-0.9880623366773218,-0.78049231427113,0.9634719640663754] atol=1e-6
# indexin(df[[1,6,11,17,16,18,8,5,9,3,12,7,13,10,2,14,4,15],:].species |> unique, df.species |> unique) == [1,2,4,6,3,5]
@test residuals(m1permute) ≈ [-0.16295769587309603,0.2560302141911026,-0.36246083625543507,
                              -0.9880623366773218,-0.78049231427113,0.9634719640663754][[1,2,4,6,3,5]] atol=1e-6
@test residuals(m2) ≈ [-0.16295769587309603,0.2560302141911026,-0.36246083625543507,-0.9880623366773218,-0.78049231427113,0.9634719640663754] atol=1e-6
@test residuals(m2permute) ≈ [-0.16295769587309603,0.2560302141911026,-0.36246083625543507,
                              -0.9880623366773218,-0.78049231427113,0.9634719640663754][[5,2,3,1,6,4]] atol=1e-6
@test dof_residual(m1) == 4
@test dof_residual(m1permute) == 4
@test dof(m1) == 4

dfmissingentry = allowmissing(df); dfmissingentry.trait1[1:3] .= missing # remove entire trait1 col for species D
m3missingentry = phylolm(@formula(trait3 ~ trait1), dfmissingentry, net; reml=false, # not using data from species D
                        tipnames=:species, withinspecies_var=true)
df_rmissingrow = copy(df_r); df_rmissingrow = df_rmissingrow[2:6,:] # remove entire row for species D
m4missingrow = phylolm(@formula(trait3 ~ trait1), df_rmissingrow, net; reml=false, # not using data from species D
                        tipnames=:species, withinspecies_var=true, y_mean_std=true)
@test coef(m3) ≈ [9.63523,2.30832] rtol=1e-4
@test coef(m4) ≈ [9.63523,2.30832] rtol=1e-4
@test coef(m3missingentry) ≈ coef(m4missingrow) rtol=1e-4
@test sigma2_phylo(m3) ≈ 0.102255 rtol=1e-4
@test sigma2_phylo(m4) ≈ 0.102255 rtol=1e-4
@test sigma2_phylo(m3missingentry) ≈ sigma2_phylo(m4missingrow) rtol=1e-4
@test sigma2_within(m3) ≈ 0.008677 rtol=1e-4
@test sigma2_within(m4) ≈ 0.008677 rtol=1e-4
@test sigma2_within(m3missingentry) ≈ sigma2_within(m4missingrow) rtol=1e-4
@test loglikelihood(m3) ≈ 1.876606 rtol=1e-4
@test loglikelihood(m4) ≈ 1.876606 rtol=1e-4
@test loglikelihood(m3missingentry) ≈ loglikelihood(m4missingrow) rtol=1e-4
@test stderror(m3) ≈ [1.0619360781577734, 0.22496955609230126] atol=1e-6
@test stderror(m4) ≈ [1.0619360781577734, 0.22496955609230126] atol=1e-6
@test stderror(m3missingentry) ≈ stderror(m4missingrow) atol=1e-6
@test residuals(m3missingentry) ≈ residuals(m4missingrow) atol=1e-6
end # agreement test subset

@testset "model comparison & likelihood ratio test" begin
m3null = phylolm(@formula(trait3 ~ 1), df, net; reml=false, tipnames=:species, withinspecies_var=true)
m3full = phylolm(@formula(trait3 ~ trait1 + trait2), df, net; reml=false, tipnames=:species, withinspecies_var=true)
@test StatsModels.isnested(m3null, m3full)
tab = (@test_logs lrtest(m3null, m3, m3full))
@test all(isapprox.(tab.deviance, (13.720216379785523,-3.75321128763398,-4.243488375445074), atol=1e-6))
@test tab.dof == (3, 4, 5)
@test all(isapprox.(tab.pval[2:end], (2.9135155795020905e-5,0.4838037279625203), atol=1e-6))
m1null = phylolm(@formula(trait3 ~ 1), df, net; tipnames=:species, withinspecies_var=true) # REML
@test !(@test_logs (:error, r"fitted with ML") StatsModels.isnested(m1null, m1)) # REML, different FEs
@test_logs (:error, r"same criterion") (@test_throws ArgumentError lrtest(m3null, m1)) # ML and REML
m2w = phylolm(@formula(trait3 ~ trait1), df_r, net; tipnames=:species)
m5w = phylolm(@formula(trait3 ~ trait1), df[[1,2,4,7,13,18],:], net; tipnames=:species, withinspecies_var=true)
@test !(@test_logs (:error, r"same number of obs") StatsModels.isnested(m2, m2w))
@test !(@test_logs (:error, r"same response") StatsModels.isnested(m5w,m2w))
m6w = (@test_logs (:info,r"^Maximum lambda") phylolm(@formula(trait3 ~ trait1), df_r, net; tipnames=:species, model="lambda"))
tab = lrtest(m2w, m6w) # both REML, but same predictors
@test tab.pval[2] ≈ 0.03928341265297505 atol=1e-6
@test !PhyloNetworks.isnested(PhyloNetworks.PagelLambda(0.1),PhyloNetworks.BM())
@test !PhyloNetworks.isnested(PhyloNetworks.PagelLambda(),PhyloNetworks.ScalingHybrid())
@test !PhyloNetworks.isnested(PhyloNetworks.ScalingHybrid(2.0),PhyloNetworks.PagelLambda())
@test_throws ArgumentError ftest(m3null, m3, m3full) # not the same Y after transformation
end # lrt test subset

@testset "equivalence with Pagel's lambda on expanded net" begin
# find a given named tip or node
function getNode(name::String, net::HybridNetwork)
    i = findfirst(n -> n.name == name, net.node)
    isnothing(i) && error("Node $(name) was not found in the network.")
    return net.node[i]
end
# new column with tip ids, initialized to species names
dfbig = deepcopy(df)
insertcols!(dfbig, 1, :speciesIds => dfbig[!,:species])
# add tips with zero length branches on the network
netbig = deepcopy(net)
for i in 1:nrow(dfbig)
    sp_name = dfbig[i, :species]; tip_name = sp_name * string(i)
    PhyloNetworks.addleaf!(netbig, getNode(sp_name, netbig), tip_name, 0.0)
    dfbig[i, :speciesIds] = tip_name
end
# (((((D1:0.0,D2:0.0,D3:0.0)D:0.4,(C4:0.0,C5:0.0,C6:0.0)C:0.4):4.8,(((A7:0.0,A8:0.0,A9:0.0)A:0.8,(B10:0.0,B11:0.0,B12:0.0)B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0.0::0.3,(E13:0.0,E14:0.0,E15:0.0)E:3.0):6.2):2.0,(O16:0.0,O17:0.0,O18:0.0)O:11.2);
# fit lambda on bigger network
m1big = (@test_logs (:info, r"^M") phylolm(@formula(trait3 ~ trait1), dfbig, netbig, model="lambda"; tipnames=:speciesIds))
# equivalent phylo and within variances
net_height = maximum(PhyloNetworks.getHeights(net));
sig_phy_2 = sigma2_phylo(m1big) * lambda_estim(m1big);
sig_err_2 = sigma2_phylo(m1big) * (1 - lambda_estim(m1big)) * net_height;
@test sigma2_phylo(m1) ≈ sig_phy_2 atol=1e-7
@test sigma2_within(m1) ≈ sig_err_2 atol=1e-9
# equivalent lambda value; equivalent coefficients and other
w_var_lambda = sigma2_phylo(m1) / (sigma2_within(m1) / net_height + sigma2_phylo(m1));
@test w_var_lambda ≈ lambda_estim(m1big) atol=1e-8
@test coef(m1) ≈ coef(m1big) atol=1e-7
@test mu_phylo(m1) ≈ mu_phylo(m1big) atol=1e-7
@test vcov(m1) ≈ vcov(m1big) atol=1e-6
@test stderror(m1) ≈ stderror(m1big) atol=1e-6
@test nobs(m1) ≈ nobs(m1big) atol=1e-10
@test dof(m1) ≈ dof(m1big) atol=1e-10
@test deviance(m1) ≈ deviance(m1big) atol=1e-10
@test loglikelihood(m1) ≈ loglikelihood(m1big) atol=1e-10
# what differs: residuals, response, predict (not same length)
# also differs: dof_residuals (18-2=16 vs 6-2=4), and so predint
end

allowmissing!(df,  [:trait3]); df[4:6,:trait3] .= missing # species C missing
allowmissing!(df_r,[:trait3]); df_r[2,:trait3] = missing  # to check imputation

@testset "ancestral state prediction, intercept only" begin
m1 = phylolm(@formula(trait3 ~ 1), df_r, net; tipnames=:species, withinspecies_var=true, y_mean_std=true)
ar1 = (@test_logs (:warn, r"^T") ancestralStateReconstruction(m1))
# ar.NodeNumbers[8] == 2 (looking at node #2), m1.model.V.tipNames[indexin([2],m1.model.V.tipNumbers)[1]] == "C" (looking at tip "C")
@test ar1.traits_nodes[8] ≈ 18.74416393519304 rtol=1e-5 # masked sampled C_bar was 17.0686
@test predint(ar1)[8,:] ≈ [15.24005506417728,22.2482728062088] rtol=1e-5
# on dataframe with model passed as keyword args. must be individual data.
ar2 = (@test_logs (:warn, r"^T") ancestralStateReconstruction(df[!,[:species,:trait3]], net; tipnames=:species, withinspecies_var=true))
@test ar2.traits_nodes ≈ ar1.traits_nodes rtol=1e-5
@test predint(ar2) ≈ predint(ar1) rtol=1e-5
# When withinspecies_var=true, predicted values at the tips are part of
# "traits_nodes", not "traits_tips", and differ from the observed sample means.
@test length(ar1.traits_tips) == 0
@test length(ar2.traits_tips) == 0
@test length(ar1.traits_nodes) == 13
@test length(ar2.traits_nodes) == 13
@test ar1.traits_nodes[9] ≈ 18.895327175656757 # for tip D: observed 18.896
@test predint(ar1)[9,:] ≈ [18.73255168713768,19.058102664175834] # narrow
end

@testset "ancestral state prediction, more than intercept" begin
m3 = phylolm(@formula(trait3 ~ trait1 + trait2), df[[1,6,11,17,16,18,8,5,9,3,12,7,13,10,2,14,4,15],:], net; tipnames=:species, withinspecies_var=true)
X_n = [m3.model.X;m3.model.X[1:3,:]] # 8x3 Array
ar3 = (@test_logs (:warn, r"^T") ancestralStateReconstruction(m3, X_n))
m4 = phylolm(@formula(trait3 ~ trait1 + trait2), df_r[[1,4,6,3,2,5],:], net; tipnames=:species, withinspecies_var=true, y_mean_std=true)
ar4 = (@test_logs (:warn, r"^T") ancestralStateReconstruction(m4, X_n))
@test ar3.NodeNumbers == ar4.NodeNumbers
@test ar3.traits_nodes ≈ ar4.traits_nodes rtol=1e-5
@test ar3.variances_nodes ≈ ar4.variances_nodes rtol=1e-5
@test size(expectations(ar4)) == (13, 2)
@test expectationsPlot(ar4)[8,2] == "24.84*"
@test expectationsPlot(ar3)[8,2] == "24.84*"
@test predint(ar3)[13,:] ≈ [15.173280800793783,15.677786825808212] rtol=1e-5 # narrow at tip with data
@test predint(ar3)[7,:] ≈ [10.668220499837211,15.322037065478693] rtol=1e-5  # wide at a root
p3 = predintPlot(ar3)
p4 = predintPlot(ar4)
@test p3[!,:nodeNumber] == p4[!,:nodeNumber]
@test p3[13,2] == "[15.17, 15.68]"
@test p4[13,2] == "[15.17, 15.68]"
@test p3[7,2] == "[10.67, 15.32]"
@test p4[7,2] == "[10.67, 15.32]"
end # test subset

end # test set: withinspecies_var on network h=1
