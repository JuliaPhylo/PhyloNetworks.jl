# Test of PhyloNetworklm on trees

###############################################################################
## Caudata dataset - shared paths matrix
###############################################################################
@testset "phyloNetworklm: Caudata Dataset" begin
## Export "caudata" dataset (from geiger)
phy = readTopology(joinpath(@__DIR__, "..", "examples", "caudata_tree.txt"));

V = sharedPathMatrix(phy);
VR = CSV.read(joinpath(@__DIR__, "..", "examples", "caudata_shared_paths.txt");
              types=[Float64 for i in 1:393]); # to avoid bug in CSV
VR = convert(Matrix, VR);

# Tips
@test V[:Tips] ≈ VR[1:197, 1:197]

# Internal nodes
@test V[:InternalNodes] ≈ VR[-V.internalNodeNumbers .+ 196, -V.internalNodeNumbers .+ 196]

# Tips Nodes
@test V[:TipsNodes] ≈ VR[1:197, -V.internalNodeNumbers .+ 196]

### R Code to get those results
# library(geiger)
# ## Load data caudata (salamanders)
# data("caudata")
# ## Save tree
# write.tree(caudata$phy, file = "caudata_tree.txt", append = FALSE,
#            digits = 10, tree.names = FALSE)
# 
# ## Times shared
# V <- node.depth.edgelength(caudata$phy)
# prac <- mrca(caudata$phy, full = TRUE)
# V <- matrix(V[prac], dim(prac))
# write.table(V,
#             file = "caudata_shared_paths.txt",
#             sep = ",", row.names = FALSE,
#             col.names = TRUE)

###############################################################################
## Caudata dataset - BM
###############################################################################

## Export "caudata" dataset (from geiger)
phy = readTopology(joinpath(@__DIR__, "..", "examples", "caudata_tree.txt"));
dat = CSV.read(joinpath(@__DIR__, "..", "examples", "caudata_trait.txt"));

## Fit a BM
fitBM = phyloNetworklm(@formula(trait ~ 1), dat, phy)

# Tests against results obtained with geiger::fitContinuous or phylolm::phylolm
@test loglikelihood(fitBM) ≈ -78.9611507833 atol=1e-10
@test dof(fitBM) ≈ 2.0 atol=1e-10
@test aic(fitBM) ≈ 161.9223015666 atol=1e-10
@test aicc(fitBM) ≈ 161.9841572367 atol=1e-10
@test coef(fitBM) ≈ [4.6789989001] atol=1e-8
@test vcov(fitBM) ≈ [0.1093144100] atol=1e-10
@test nobs(fitBM) ≈ 197.0 atol=1e-10
@test sum(residuals(fitBM)) ≈ -115.5767321312 atol=1e-8
@test dof_residual(fitBM) ≈ 196.0 atol=1e-10
@test sigma2_estim(fitBM) ≈ 0.0029452097 atol=1e-10
@test stderror(fitBM) ≈ [0.3306272978] atol=1e-10
@test confint(fitBM)[1] ≈ 4.0269551772 atol=1e-10
@test confint(fitBM)[2] ≈ 5.3310426231 atol=1e-10
tmp = predict(fitBM);
# below: print for when the test was broken
#@show length(tmp) # fixit: tests fail here
#@show tmp[1:6], tmp[190:end] # 4.679 except for last 5: 9.31707, 9.33832, 9.33052, 9.34356, 9.31707
#println("are they all 4.6789989001?")
# next: looks random. sometimes passes, most times fails
@test predict(fitBM) ≈ [4.6789989001 for i in 1:197] atol=1e-8
# @test_approx_eq_eps deviance(fitBM)  deviance(fitbis)
# @test_approx_eq_eps nulldeviance(fitBM)  nulldeviance(fitbis)
# @test_approx_eq_eps nullloglikelihood(fitBM)  nullloglikelihood(fitbis)
# @test_approx_eq_eps r2(fitBM)  r2(fitbis)
# @test_approx_eq_eps adjr2(fitBM)  adjr2(fitbis)
# @test_approx_eq_eps bic(fitBM)  bic(fitbis)
# @test_approx_eq_eps mu_estim(fitBM)  mu_estim(fitbis)

### Ancestral state reconstruction (with Rphylopars)
anc = (@test_logs (:warn, r"^These prediction intervals show uncertainty in ancestral values") ancestralStateReconstruction(fitBM));
ancR = CSV.read(joinpath(@__DIR__, "..", "examples", "caudata_Rphylopars.txt"),
                types=[Float64,Float64]); # to avoid bug in CSV

## Expectations
expe = expectations(anc)
expeR = ancR[!,:trait]
# Matching tips ?
tipsR = expeR[expe[197:393, :nodeNumber]]
tipsJulia = expe[197:393, :condExpectation]
@test tipsR ≈ tipsJulia
# Matching nodes ?
nodesR = expeR[-expe[1:196, :nodeNumber] .+ 196]
nodesJulia = expe[1:196, :condExpectation]
# below: print for when the test was broken
#@show nodesR[1:6],    nodesR[190:end]
#@show nodesJulia[1:6],nodesJulia[190:end]
@test isapprox(nodesR, nodesJulia)

## Variances
vars = LinearAlgebra.diag(anc.variances_nodes)
# Rphylopars
varsR = ancR[!,:var]
# Matching nodes ?
nodesR = varsR[-expe[1:196, :nodeNumber] .+ 196]
@test nodesR ≈ vars atol=1e-3 ## RK: Small tol !!

### Ancestral state reconstruction (with Phytools)
ancRt = CSV.read(joinpath(@__DIR__, "..", "examples", "caudata_Phytools.txt"));

## Expectations
expe = expectations(anc)
expeRt = ancRt[!,:trait]
# Matching nodes ?
nodesRt = expeRt[-expe[1:196, :nodeNumber] .+ (196 - 197)]
nodesJulia = expe[1:196, :condExpectation]
# below: print for when the test was broken
#@show nodesRt[1:6],   nodesRt[190:end]
#@show nodesJulia[1:6],nodesJulia[190:end]
@test isapprox(nodesRt, nodesJulia)

## Variances
vars = LinearAlgebra.diag(anc.variances_nodes)
# Rphylopars
varsRt = ancRt[!,:var]
# Matching nodes ?
nodesRt = varsRt[-expe[1:196, :nodeNumber] .+ (196 - 197)]
@test nodesRt ≈ vars atol=2e-3 ## RK: Small tol !!

### Comparison between Rphylopars and Phytools:
@test nodesRt ≈ nodesR atol=0.003 ## RK: Small tol !!

### R script to get the above values:
# library(geiger)
# 
# ## Load data caudata (salamanders)
# data("caudata")
# 
# ## Save tree and data
# write.tree(caudata$phy, file = "caudata_tree.txt", append = FALSE,
#            digits = 10, tree.names = FALSE)
# 
# write.table(data.frame(tipsNames = names(caudata$dat),
#                        trait = unname(caudata$dat)),
#             file = "caudata_trait.txt",
#             sep = ",", row.names = FALSE)
# 
# ## Fit using Geiger
# fitgeiger <- fitContinuous(caudata$phy, caudata$dat, model = "BM")
# 
# ## Fit using phylolm
# library(phylolm)
# fitphylolm <- phylolm(trait ~ 1, 
#                       data.frame(trait = caudata$dat),
#                       caudata$phy, model = "BM")
# 
# ## Fit using Rphylopars
# library(Rphylopars)
# fitphylopars <- phylopars(data.frame(species = names(caudata$dat),
#                           trait = unname(caudata$dat)),
#                           caudata$phy,
#                           pheno_error = FALSE,
#                           pheno_correlated = FALSE,
#                           REML = FALSE)
# 
# # Save results of Rphylopars for ancestral trait reconstruction
# write.table(data.frame(trait = unname(fitphylopars$anc_recon),
#                        var = unname(fitphylopars$anc_var)),
#             file = "caudata_Rphylopars.txt",
#             sep = ",", row.names = FALSE)
# 
# ## Ancestral State reconstruction using phytools
# library(phytools)
# fitphytools <- fastAnc(caudata$phy, caudata$dat, vars = TRUE)
# 
# # Save results of Rphylopars for ancestral trait reconstruction
# write.table(data.frame(trait = unname(fitphytools$ace),
#                        var = unname(fitphytools$var)),
#             file = "caudata_Phytools.txt",
#             sep = ",", row.names = FALSE)
# 
# ## Quantities to compare
# sprintf("%.10f", fitgeiger$opt$ln) # log likelihood
# sprintf("%.10f", fitphylolm$logLik)
# sprintf("%.10f", fitphylopars$logLik)
# sprintf("%.10f", fitgeiger$opt$aic) # aic
# sprintf("%.10f", fitphylolm$aic)
# sprintf("%.10f", fitgeiger$opt$aicc) # aicc
# sprintf("%.10f", fitphylolm$coefficients) # coef
# sprintf("%.10f", fitgeiger$opt$z0)
# sprintf("%.10f", fitphylopars$mu)
# sprintf("%.10f", fitphylolm$vcov) # vcov
# sprintf("%.10f", fitphylolm$n) # nobs
# sprintf("%.10f", sum(fitphylolm$residuals)) # residuals (sum)
# sprintf("%.10f", fitphylolm$n - fitphylolm$d) # df residuals
# sprintf("%.10f", fitphylolm$sigma2) # sigma 2
# sprintf("%.10f", fitgeiger$opt$sigsq)
# sprintf("%.10f", fitphylopars$pars$phylocov)
# sprintf("%.10f", summary(fitphylolm)$coefficients[2]) # std error
# sprintf("%.10f", summary(fitphylolm)$df) # df
# sprintf("%.10f", fitgeiger$opt$k)
# sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[2] * qt(0.025, 196))
# sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[2] * qt(0.975, 196))
# sprintf("%.10f", predict(fitphylolm)) # df

###############################################################################
## Caudata dataset - Pagel's lambda
###############################################################################

## Fit Pagel's lambda
fitLambda = (@test_logs (:info, r"^Maximum lambda value") match_mode=:any phyloNetworklm(@formula(trait ~ 1), dat, phy, model = "lambda"));

@test lambda_estim(fitLambda) ≈ 0.9193 atol=1e-4 # Due to convergence issues, tolerance is lower.
@test loglikelihood(fitLambda) ≈ -51.684379 atol=1e-6
@test dof(fitLambda) ≈ 3.0 atol=1e-10
@test aic(fitLambda) ≈ 109.368759 atol=1e-6
@test aicc(fitLambda) ≈ 109.493111 atol=1e-6
@test coef(fitLambda) ≈ [4.66893] atol=1e-5
@test vcov(fitLambda) ≈ [0.05111] atol=1e-5
@test nobs(fitLambda) ≈ 197.0 atol=1e-10
# below: print for when the test was broken
#@show sum(residuals(fitLambda)) # locally: -115.91591040367894
#println("is this -113.594?")
@test isapprox(sum(residuals(fitLambda)), -113.594, atol=1e-2) ## Low Tolerance !!
@test dof_residual(fitLambda) ≈ 196.0 atol=1e-10 ## Correct Definition ?
@test sigma2_estim(fitLambda) ≈ 0.0014756 atol=1e-7
@test stderror(fitLambda) ≈ [0.22608] atol=1e-5
@test confint(fitLambda)[1] ≈ 4.2230 atol=1e-4
@test confint(fitLambda)[2] ≈ 5.114 atol=1e-3
tmp = predict(fitLambda);
# below: print for when the test was broken
#@show length(tmp)
#@show tmp[1:6], tmp[190:end] # all 4.66893 except for last 5: 8.42676, 8.44585, etc.
#println("are they all 4.66893?")
# next: looks random. sometimes passes, most times fails
@test predict(fitLambda) ≈ [4.66893 for i in 1:197] atol=1e-5 norm=x->LinearAlgebra.norm(x,Inf)

### R script to get the above values:
# library(geiger)
# 
# ## Load data caudata (salamanders)
# data("caudata")
# 
# ## Fit using Geiger
# fitgeiger <- fitContinuous(caudata$phy, caudata$dat, model = "lambda")
# 
# ## Fit using phylolm
# library(phylolm)
# fitphylolm <- phylolm(trait ~ 1, 
#                       data.frame(trait = caudata$dat),
#                       caudata$phy, model = "lambda",
#                       starting.value = 0.9)
# 
# ## Fit using Rphylopars
# library(Rphylopars)
# fitphylopars <- phylopars(data.frame(species = names(caudata$dat),
#                           trait = unname(caudata$dat)),
#                           caudata$phy,
#                           model = "lambda",
#                           model_par_start = 0.9,
#                           pheno_error = FALSE,
#                           pheno_correlated = FALSE,
#                           REML = FALSE)
# 
# ## Quantities to compare
# sprintf("%.10f", fitgeiger$opt$ln) # log likelihood
# sprintf("%.10f", fitphylolm$logLik)
# sprintf("%.10f", fitphylopars$logLik)
# sprintf("%.10f", fitgeiger$opt$aic) # aic
# sprintf("%.10f", fitphylolm$aic)
# sprintf("%.10f", fitgeiger$opt$aicc) # aicc
# sprintf("%.10f", fitphylolm$coefficients) # coef
# sprintf("%.10f", fitgeiger$opt$z0)
# sprintf("%.10f", fitphylopars$mu)
# sprintf("%.10f", fitphylolm$vcov) # vcov
# sprintf("%.10f", fitphylolm$n) # nobs
# sprintf("%.10f", sum(fitphylolm$residuals)) # residuals (sum)
# sprintf("%.10f", fitphylolm$n - fitphylolm$d) # df residuals
# sprintf("%.10f", fitphylolm$sigma2) # sigma 2
# sprintf("%.10f", fitgeiger$opt$sigsq)
# sprintf("%.10f", fitphylopars$pars$phylocov)
# sprintf("%.10f", summary(fitphylolm)$coefficients[2]) # std error
# sprintf("%.10f", summary(fitphylolm)$df) # df
# sprintf("%.10f", fitgeiger$opt$k)
# sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[2] * qt(0.025, 196))
# sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[2] * qt(0.975, 196))
# sprintf("%.10f", predict(fitphylolm)) # df

###############################################################################
## Caudata dataset - BM with shifts
###############################################################################

## Export "caudata" dataset (from geiger)
phy = readTopology(joinpath(@__DIR__, "..", "examples", "caudata_tree.txt"));
dat = CSV.read(joinpath(@__DIR__, "..", "examples", "caudata_trait.txt"));

## Add some shifts in the model
df_shift = regressorShift(phy.edge[[98, 326, 287]], phy)
dat = join(dat, df_shift, on=:tipNames)
## Fit a BM
fitBM = phyloNetworklm(@formula(trait ~ shift_98 + shift_326 + shift_287), dat, phy)

# Tests against results obtained with geiger::fitContinuous or phylolm::phylolm
@test loglikelihood(fitBM) ≈ -76.1541605207 atol=1e-10 
@test dof(fitBM) ≈ 5.0 atol=1e-10
@test aic(fitBM) ≈ 162.3083210414 atol=1e-10 
@test coef(fitBM) ≈ [4.8773804547 -0.6664678924 -1.0692950215 -0.0159349201]' atol=1e-10 
vcovR = [ 0.1165148753  -0.0431446679 -0.0305707092 0.0000000000
         -0.0431446679 0.3212796435  0.0340202694  -0.0000000000 
         -0.0305707092 0.0340202694 0.2472613252  -0.0236868278 
          0.0000000000  -0.0000000000 -0.0236868278 0.1263794810]
@test vcov(fitBM) ≈ vcovR atol=1e-8
@test nobs(fitBM) ≈ 197.0 atol=1e-10
@test sum(residuals(fitBM)) ≈ -10.1752334428 atol=1e-10
@test dof_residual(fitBM) ≈ 193.0 atol=1e-10
@test sigma2_estim(fitBM) ≈ 0.0028624636 atol=1e-10
@test stderror(fitBM) ≈ [0.3413427535 0.5668153522 0.4972537835 0.3554989184]' atol=1e-10
@test confint(fitBM)[:,1] ≈ [4.2041393297 -1.7844157659 -2.0500444097 -0.7170966977]' atol=1e-10
@test confint(fitBM)[:,2] ≈ [5.5506215797 0.4514799811 -0.0885456333 0.6852268574]' atol=1e-10
predictR = [4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 4.2109125624, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.7921505131, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 3.8080854332, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547, 4.8773804547]
tmp  = predict(fitBM)
tmp2 = predictR[fitBM.model.ind]
# below: print for when the test was broken
#@show tmp[1:6],  tmp[190:end]
#@show tmp2[1:6],tmp2[190:end] # the last 5 values are different
@test isapprox(predict(fitBM), predictR[fitBM.model.ind], atol=1e-8)

# ## R code to get those results
# library(geiger)
# ## Load data caudata (salamanders)
# data("caudata")
# ## Re-order with the tree
# dat <- caudata$dat[match(caudata$phy$tip.label, names(caudata$dat))]
# ## Incidence matrix
# Tm = PhylogeneticEM::incidence.matrix(caudata$phy) + 0
# ## Fit using phylolm
# library(phylolm)
# fitphylolm <- phylolm(trait ~ shift_37 + shift_105 + shift_222, 
#                       data.frame(trait = dat,
#                                  shift_37 = Tm[, 37],
#                                  shift_105 = Tm[, 105],
#                                  shift_222 = Tm[, 222]),
#                       caudata$phy, model = "BM")
# ## Quantities to compare
# sprintf("%.10f", fitphylolm$logLik)
# sprintf("%.10f", fitphylolm$aic)
# sprintf("%.10f", fitphylolm$coefficients) # coef
# sprintf("%.10f", fitphylolm$vcov) # vcov
# sprintf("%.10f", fitphylolm$n) # nobs
# sprintf("%.10f", sum(fitphylolm$residuals)) # residuals (sum)
# sprintf("%.10f", fitphylolm$n - fitphylolm$d) # df residuals
# sprintf("%.10f", fitphylolm$sigma2) # sigma 2
# sprintf("%.10f", summary(fitphylolm)$coefficients[,2]) # std error
# sprintf("%.10f", summary(fitphylolm)$df) # df
# sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[,2] * qt(0.025, 193))
# sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[,2] * qt(0.975, 193))
# sprintf("%.10f", predict(fitphylolm))

end

###############################################################################
## Lizard dataset - BM
###############################################################################

@testset "phyloNetworklm: Lizard Dataset" begin

## Export "lizard" dataset (Mahler et al 2013)
phy = readTopology(joinpath(@__DIR__, "..", "examples", "lizard_tree.txt"));
dat = CSV.read(joinpath(@__DIR__, "..", "examples", "lizard_trait.txt"));
categorical!(dat, :region)

## Fit a BM
fitBM = phyloNetworklm(@formula(AVG_SVL ~ AVG_ltoe_IV + AVG_lfing_IV * region), dat, phy)

# Tests against results obtained with geiger::fitContinuous or phylolm::phylolm
@test loglikelihood(fitBM) ≈ 105.17337853473711 atol=1e-10 
@test dof(fitBM) ≈ 10.0 atol=1e-10
@test aic(fitBM) ≈ -190.3467570695 atol=1e-10 
# @test aicc(fitBM)
@test coef(fitBM) ≈ [2.7925712673 -0.2010704391 0.9832555589 -0.1021226296 -0.3703658712 0.1557471731 0.0374549036 0.1805667675 -0.0495767233]' atol=1e-10 
vcovR =  [0.0200086273  -0.0136717540 0.0084815090  -0.0093192029 -0.0114417825 -0.0113346813 0.0041102304  0.0053787287  0.0050521693 
          -0.0136717540 0.0185396965  -0.0169114682 0.0020645005  0.0036352899  0.0026227856 -0.0012281620 -0.0018838231 -0.0014800242
          0.0084815090  -0.0169114682 0.0174647413  0.0016403871  0.0005624488  0.0015003301 -0.0005714235 -0.0003201257 -0.0005423354
          -0.0093192029 0.0020645005  0.0016403871  0.0167953394  0.0078012534  0.0086329399 -0.0080782771 -0.0037333495 -0.0039327836
          -0.0114417825 0.0036352899  0.0005624488  0.0078012534  0.0490482083  0.0092203882 -0.0033670465 -0.0191567265 -0.0040068947
          -0.0113346813 0.0026227856  0.0015003301  0.0086329399  0.0092203882  0.0331395502 -0.0037513830 -0.0041592743 -0.0146108207
          0.0041102304  -0.0012281620 -0.0005714235 -0.0080782771 -0.0033670465 -0.0037513830 0.0045172675  0.0018165174  0.0020857846
          0.0053787287  -0.0018838231 -0.0003201257 -0.0037333495 -0.0191567265 -0.0041592743 0.0018165174  0.0093292284  0.0020427637
          0.0050521693  -0.0014800242 -0.0005423354 -0.0039327836 -0.0040068947 -0.0146108207 0.0020857846  0.0020427637  0.0074817942]
@test vcov(fitBM) ≈ vcovR atol=1e-8
@test nobs(fitBM) ≈ 100.0 atol=1e-10
# below: print for when the test was broken (with @test_skip)
#@show sum(residuals(fitBM)) # looks random, e.g. 1.6091413520064477, or 0.8338189090359597
#println("is this equal to 0.6352899255?") # sometimes NO, yet the test passes below!!
@test sum(residuals(fitBM)) ≈ 0.6352899255 atol=1e-10
@test dof_residual(fitBM) ≈ 91.0 atol=1e-10
@test sigma2_estim(fitBM) ≈ 0.0003025014 atol=1e-10
@test stderror(fitBM) ≈ [0.1414518551,0.1361605540,0.1321542330,0.1295968341,0.2214683008,0.1820427154,0.0672106202,0.0965879311,0.0864973651] atol=1e-10
@test confint(fitBM)[:,1] ≈ [2.5115945339,-0.4715366529,0.7207474097,-0.3595508202,-0.8102854443,-0.2058583178,-0.0960507369,-0.0112932922,-0.2213931131] atol=1e-10 norm=x->LinearAlgebra.norm(x,Inf)
@test confint(fitBM)[:,2] ≈ [3.0735480006,0.0693957746,1.2457637082,0.1553055609,0.0695537019,0.5173526640,0.1709605441,0.3724268272,0.1222396666] atol=1e-10
# @test_approx_eq_eps predict(fitBM)[fitbis.model.ind] predict(fitbis)
# @test_approx_eq_eps deviance(fitBM)  deviance(fitbis)
# @test_approx_eq_eps nulldeviance(fitBM)  nulldeviance(fitbis)
# @test_approx_eq_eps nullloglikelihood(fitBM)  nullloglikelihood(fitbis)
# @test_approx_eq_eps r2(fitBM)  r2(fitbis)
# @test_approx_eq_eps adjr2(fitBM)  adjr2(fitbis)
# @test_approx_eq_eps bic(fitBM)  bic(fitbis)
# @test_approx_eq_eps mu_estim(fitBM)  mu_estim(fitbis)


### R script to get the above values
# ## Data
# dat <- read.csv(file = "GA_Anolis_traits.csv")
# geo <- read.csv(file = "GA_Anolis_biogeography.csv")
#
# dat <- merge(dat, geo, by = "species")
# # keep only toe and hand length
# dat <- dat[, c("species", "AVG.SVL", "AVG.ltoe.IV", "AVG.lfing.IV", "region")]
# colnames(dat)[1] <- "tipsNames"
#
# write.table(dat,
#             file = "lizard_trait.txt",
#             sep = ",", row.names = FALSE)
#
# ## Tree
# phy <- read.tree(file = "GA_Anolis_MCC.tre")
#
# write.tree(phy, file = "lizards_tree.txt", append = FALSE,
#            digits = 10, tree.names = FALSE)
#
# rownames(dat) <- dat$tipsNames
# dat <- dat[, -1]
# dat$region <- as.factor(dat$region)
#
# ## Fit
# fitphylolm <- phylolm(AVG.SVL ~ 1 + AVG.ltoe.IV + AVG.lfing.IV * region, dat, phy, model = "BM")
# ## Quantities to compare
# sprintf("%.10f", fitphylolm$logLik)
# sprintf("%.10f", summary(fitphylolm)$df) # df
# sprintf("%.10f", fitphylolm$aic)
# sprintf("%.10f", fitphylolm$coefficients) # coef
# matrix(sprintf("%.10f", fitphylolm$vcov), 9, 9) # vcov
# sprintf("%.10f", fitphylolm$n) # nobs
# sprintf("%.10f", sum(fitphylolm$residuals)) # residuals (sum)
# sprintf("%.10f", fitphylolm$n - fitphylolm$d) # df residuals
# sprintf("%.10f", fitphylolm$sigma2) # sigma 2
# sprintf("%.10f", summary(fitphylolm)$coefficients[,2]) # std error
# sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[, 2] * qt(0.025, 91))
# sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[, 2] * qt(0.975, 91))

###############################################################################
## Lizard dataset - lambda
###############################################################################

## Fit lambda
fitLambda = (@test_logs (:info, r"^Maximum lambda value") match_mode=:any phyloNetworklm(@formula(AVG_SVL ~ AVG_ltoe_IV + AVG_lfing_IV * region), dat, phy, model = "lambda"))

# Tests against results obtained with geiger::fitContinuous or phylolm::phylolm
@test lambda_estim(fitLambda) ≈ 0.9982715594 atol=1e-5
@test loglikelihood(fitLambda) ≈ 105.1769275564 atol=1e-8
@test dof(fitLambda) ≈ 11.0 atol=1e-10
@test aic(fitLambda) ≈ -188.3538551128 atol=1e-8 
# @test aicc(fitBM)
@test coef(fitLambda) ≈ [2.7940573420 -0.2066584606 0.9897083949 -0.1004840950 -0.3677991157 0.1576743022 0.0367633665 0.1792502383 -0.0505291142]' atol=1e-5 
vcovR =  [0.0200251600  -0.0137474015 0.0085637021  -0.0092973836 -0.0114259722 -0.0113056243 0.0041037877  0.0053740100  0.0050429112 
          -0.0137474015 0.0186885224  -0.0170645512 0.0020509207  0.0036334103 0.0026066694  -0.0012237488 -0.0018826836 -0.0014743137
          0.0085637021  -0.0170645512 0.0176200143  0.0016494733  0.0005604169 0.0015116125  -0.0005735573 -0.0003192320 -0.0005457420
          -0.0092973836 0.0020509207  0.0016494733  0.0167461876  0.0077885115 0.0086173037  -0.0080563819 -0.0037287856 -0.0039275469
          -0.0114259722 0.0036334103  0.0005604169  0.0077885115  0.0490092393 0.0092036032  -0.0033631662 -0.0191657329 -0.0040017905
          -0.0113056243 0.0026066694  0.0015116125  0.0086173037  0.0092036032 0.0330248707  -0.0037465110 -0.0041543671 -0.0145663751
          0.0041037877  -0.0012237488 -0.0005735573 -0.0080563819 -0.0033631662 -0.0037465110 0.0045042057  0.0018142470  0.0020823721 
          0.0053740100  -0.0018826836 -0.0003192320 -0.0037287856 -0.0191657329 -0.0041543671 0.0018142470  0.0093334212  0.0020404652 
          0.0050429112  -0.0014743137 -0.0005457420 -0.0039275469 -0.0040017905 -0.0145663751 0.0020823721  0.0020404652  0.0074600880]
@test vcov(fitLambda) ≈ vcovR atol=3e-7 
@test nobs(fitLambda) ≈ 100.0 atol=1e-10
# below: print for when the test was broken (with @test_skip)
#@show sum(residuals(fitLambda)) # looks random, eg 0.033126277561337916 or 0.644941961666333
#println("is this equal to 0.6369008979?")
@test sum(residuals(fitLambda)) ≈ 0.6369008979 atol=1e-5
@test dof_residual(fitLambda) ≈ 91.0 atol=1e-10
@test sigma2_estim(fitLambda) ≈ 0.0003009914 atol=1e-9
@test stderror(fitLambda) ≈ [0.1415102824,0.1367059706,0.1327404019,0.1294070617,0.2213803048,0.1817274626,0.0671133793,0.0966096332,0.0863718011] atol=1e-6
@test confint(fitLambda)[:,1] ≈ [2.5129645499,-0.4782080775,0.7260358930,-0.3575353260,-0.8075438955,-0.2033049779,-0.0965491169,-0.0126529301,-0.2220960868] atol=1e-5
@test confint(fitLambda)[:,2] ≈ [3.0751501341,0.0648911562,1.2533808968,0.1565671360,0.0719456641,0.5186535822,0.1700758500,0.3711534067,0.1210378584] atol=1e-5
# @test_approx_eq_eps predict(fitLambda)[fitbis.model.ind] predict(fitbis)
# @test_approx_eq_eps deviance(fitLambda)  deviance(fitbis)
# @test_approx_eq_eps nulldeviance(fitLambda)  nulldeviance(fitbis)
# @test_approx_eq_eps nullloglikelihood(fitLambda)  nullloglikelihood(fitbis)
# @test_approx_eq_eps r2(fitLambda)  r2(fitbis)
# @test_approx_eq_eps adjr2(fitLambda)  adjr2(fitbis)
# @test_approx_eq_eps bic(fitLambda)  bic(fitbis)
# @test_approx_eq_eps mu_estim(fitLambda)  mu_estim(fitbis)


### R script to get the above values
# ## Data
# dat <- read.csv(file = "GA_Anolis_traits.csv")
# geo <- read.csv(file = "GA_Anolis_biogeography.csv")
#
# dat <- merge(dat, geo, by = "species")
# # keep only toe and hand length
# dat <- dat[, c("species", "AVG.SVL", "AVG.ltoe.IV", "AVG.lfing.IV", "region")]
# colnames(dat)[1] <- "tipsNames"
#
# write.table(dat,
#             file = "lizard_trait.txt",
#             sep = ",", row.names = FALSE)
#
# ## Tree
# phy <- read.tree(file = "GA_Anolis_MCC.tre")
#
# write.tree(phy, file = "lizards_tree.txt", append = FALSE,
#            digits = 10, tree.names = FALSE)
#
# rownames(dat) <- dat$tipsNames
# dat <- dat[, -1]
# dat$region <- as.factor(dat$region)
#
# ## Fit
# fitphylolm <- phylolm(AVG.SVL ~ 1 + AVG.ltoe.IV + AVG.lfing.IV * region, dat, phy, model = "lambda")
# ## Quantities to compare
# sprintf("%.10f", fitphylolm$logLik)
# sprintf("%.10f", summary(fitphylolm)$df) # df
# sprintf("%.10f", fitphylolm$aic)
# sprintf("%.10f", fitphylolm$coefficients) # coef
# matrix(sprintf("%.10f", fitphylolm$vcov), 9, 9) # vcov
# sprintf("%.10f", fitphylolm$n) # nobs
# sprintf("%.10f", sum(fitphylolm$residuals)) # residuals (sum)
# sprintf("%.10f", fitphylolm$n - fitphylolm$d) # df residuals
# sprintf("%.10f", fitphylolm$sigma2) # sigma 2
# sprintf("%.10f", summary(fitphylolm)$coefficients[,2]) # std error
# sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[, 2] * qt(0.025, 91))
# sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[, 2] * qt(0.975, 91))

end
