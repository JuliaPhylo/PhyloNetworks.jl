# Test of PhyloNetworklm on trees

###############################################################################
## Caudata dataset - shared paths matrix
###############################################################################

## Export "caudata" dataset (from geiger)
phy = readTopology(joinpath(Pkg.dir("PhyloNetworks"), "examples", "caudata_tree.txt"));

V = sharedPathMatrix(phy);
VR = readtable(joinpath(Pkg.dir("PhyloNetworks"), "examples", "caudata_shared_paths.txt"));
VR = convert(Matrix, VR);

# Tips
@test_approx_eq V[:Tips] VR[1:197, 1:197]

# Internal nodes
@test_approx_eq V[:InternalNodes] VR[-V.internalNodeNumbers + 196, -V.internalNodeNumbers + 196]

# Tips Nodes
@test_approx_eq V[:TipsNodes] VR[1:197, -V.internalNodeNumbers + 196]

### R Code to get those results
# library(geiger)
# ## Load data caudata (salamanders)
# data("caudata")
# ## Save tree
# write.tree(caudata$phy, file = "caudata_tree.txt", append = FALSE,
# 					 digits = 10, tree.names = FALSE)
# 
# ## Times shared
# V <- node.depth.edgelength(caudata$phy)
# prac <- mrca(caudata$phy, full = TRUE)
# V <- matrix(V[prac], dim(prac))
# write.table(V,
# 						file = "caudata_shared_paths.txt",
# 						sep = ",", row.names = FALSE,
# 						col.names = TRUE)

###############################################################################
## Caudata dataset - BM
###############################################################################

## Export "caudata" dataset (from geiger)
phy = readTopology(joinpath(Pkg.dir("PhyloNetworks"), "examples", "caudata_tree.txt"));
dat = readtable(joinpath(Pkg.dir("PhyloNetworks"), "examples", "caudata_trait.txt"));

## Fit a BM
fitBM = phyloNetworklm(trait ~ 1, dat, phy)

# Tests against results obtained with geiger::fitContinuous or phylolm::phylolm
@test_approx_eq_eps loglikelihood(fitBM) -78.9611507833 1e-10 
@test_approx_eq_eps dof(fitBM)  2.0 1e-10
@test_approx_eq_eps aic(fitBM) 161.9223015666 1e-10 
@test_approx_eq_eps aicc(fitBM)  161.9841572367 1e-10
@test_approx_eq_eps coef(fitBM) 4.6789989001 1e-10 
@test_approx_eq_eps vcov(fitBM) 0.1093144100 1e-10 
@test_approx_eq_eps nobs(fitBM) 197.0 1e-10
@test_approx_eq_eps sum(residuals(fitBM)) -115.5767321312 1e-10
@test_approx_eq_eps dof_residual(fitBM) 196.0 1e-10
@test_approx_eq_eps sigma2_estim(fitBM) 0.0029452097 1e-10
@test_approx_eq_eps stderr(fitBM) 0.3306272978 1e-10
@test_approx_eq_eps confint(fitBM)[1] 4.0269551772 1e-10
@test_approx_eq_eps confint(fitBM)[2] 5.3310426231 1e-10
@test_approx_eq_eps predict(fitBM) [4.6789989001].*ones(197) 1e-10
# @test_approx_eq_eps model_response(fitBM)[fitbis.model.ind] model_response(fitbis)
# @test_approx_eq_eps deviance(fitBM)  deviance(fitbis)
# @test_approx_eq_eps nulldeviance(fitBM)  nulldeviance(fitbis)
# @test_approx_eq_eps nullloglikelihood(fitBM)  nullloglikelihood(fitbis)
# @test_approx_eq_eps r2(fitBM)  r2(fitbis)
# @test_approx_eq_eps adjr2(fitBM)  adjr2(fitbis)
# @test_approx_eq_eps bic(fitBM)  bic(fitbis)
# @test_approx_eq_eps mu_estim(fitBM)  mu_estim(fitbis)

## Ancestral state reconstruction (with Rphylopars)
anc = ancestralStateReconstruction(fitBM)
ancR = readtable(joinpath(Pkg.dir("PhyloNetworks"), "examples", "caudata_Rphylopars.txt"));

## Expectations
expe = expectations(anc)
expeR = ancR[:trait]
# Matching tips ?
tipsR = expeR[expe[197:393, :nodeNumber]]
tipsJulia = expe[197:393, :condExpectation]
for i in 1:197
    @test_approx_eq tipsR[i] tipsJulia[i]
end
# Matching nodes ?
nodesR = expeR[-expe[1:196, :nodeNumber] + 196]
nodesJulia = expe[1:196, :condExpectation]
for i in 1:196
    @test_approx_eq nodesR[i] nodesJulia[i]
end

## Variances
vars = diag(anc.variances_nodes)
# Rphylopars
varsR = ancR[:var]
# Matching nodes ?
nodesR = varsR[-expe[1:196, :nodeNumber] + 196]
for i in 1:196
    @test_approx_eq_eps nodesR[i] vars[i] 1e-3 ## RK: Small tol !!
end

### R script to get the above values:
# library(geiger)
# 
# ## Load data caudata (salamanders)
# data("caudata")
# 
# ## Save tree and data
# write.tree(caudata$phy, file = "caudata_tree.txt", append = FALSE,
# 					 digits = 10, tree.names = FALSE)
# 
# write.table(data.frame(tipsNames = names(caudata$dat),
# 											 trait = unname(caudata$dat)),
# 						file = "caudata_trait.txt",
# 						sep = ",", row.names = FALSE)
# 
# ## Fit using Geiger
# fitgeiger <- fitContinuous(caudata$phy, caudata$dat, model = "BM")
# 
# ## Fit using phylolm
# library(phylolm)
# fitphylolm <- phylolm(trait ~ 1, 
# 											data.frame(trait = caudata$dat),
# 											caudata$phy, model = "BM")
# 
# ## Fit using Rphylopars
# library(Rphylopars)
# fitphylopars <- phylopars(data.frame(species = names(caudata$dat),
# 													trait = unname(caudata$dat)),
# 													caudata$phy,
# 													pheno_error = FALSE,
# 													pheno_correlated = FALSE,
# 													REML = FALSE)
# 
# ## Save results of Rphylopars for ancestral trait reconstruction
# write.table(data.frame(trait = unname(fitphylopars$anc_recon),
# 											var = unname(fitphylopars$anc_var)),
# 					  file = "caudata_Rphylopars.txt",
# 						sep = ",", row.names = FALSE)
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
fitLambda = phyloNetworklm(trait ~ 1, dat, phy, model = "lambda")

@test_approx_eq_eps lambda_estim(fitLambda) 0.9193 1e-4 # Due to convergence issues, tolerance is lower.
@test_approx_eq_eps loglikelihood(fitLambda) -51.684379 1e-6
@test_approx_eq_eps dof(fitLambda)  3.0 1e-10
@test_approx_eq_eps aic(fitLambda) 109.368759 1e-6
@test_approx_eq_eps aicc(fitLambda)  109.493111 1e-6
@test_approx_eq_eps coef(fitLambda) 4.66893 1e-5
@test_approx_eq_eps vcov(fitLambda) 0.05111 1e-5
@test_approx_eq_eps nobs(fitLambda) 197.0 1e-10
@test_approx_eq_eps sum(residuals(fitLambda)) -113.59 1e-2 ## Low Tolerance !!
@test_approx_eq_eps dof_residual(fitLambda) 196.0 1e-10 ## Correct Definition ?
@test_approx_eq_eps sigma2_estim(fitLambda) 0.0014756 1e-7
@test_approx_eq_eps stderr(fitLambda) 0.22608 1e-5
@test_approx_eq_eps confint(fitLambda)[1] 4.2230 1e-4
@test_approx_eq_eps confint(fitLambda)[2] 5.114 1e-3
@test_approx_eq_eps predict(fitLambda) [4.66893].*ones(197) 1e-5

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
# 											data.frame(trait = caudata$dat),
# 											caudata$phy, model = "lambda",
#												starting.value = 0.9)
# 
# ## Fit using Rphylopars
# library(Rphylopars)
# fitphylopars <- phylopars(data.frame(species = names(caudata$dat),
# 													trait = unname(caudata$dat)),
# 													caudata$phy,
#                           model = "lambda",
#                           model_par_start = 0.9,
# 													pheno_error = FALSE,
# 													pheno_correlated = FALSE,
# 													REML = FALSE)
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
## Lizard dataset - BM
###############################################################################

## Export "lizard" dataset (Mahler et al 2013)
phy = readTopology(joinpath(Pkg.dir("PhyloNetworks"), "examples", "lizard_tree.txt"));
dat = readtable(joinpath(Pkg.dir("PhyloNetworks"), "examples", "lizard_trait.txt"));
dat[:region] = PooledDataArray(dat[:region]); # Pool by region

## Fit a BM
fitBM = phyloNetworklm(AVG_SVL ~ AVG_ltoe_IV + AVG_lfing_IV * region, dat, phy)

# Tests against results obtained with geiger::fitContinuous or phylolm::phylolm
@test_approx_eq_eps loglikelihood(fitBM) 105.17337853473711 1e-10 
@test_approx_eq_eps dof(fitBM) 10.0 1e-10
@test_approx_eq_eps aic(fitBM) -190.3467570695 1e-10 
# @test_approx_eq_eps aicc(fitBM)
@test_approx_eq_eps coef(fitBM) [2.7925712673 -0.2010704391 0.9832555589 -0.1021226296 -0.3703658712 0.1557471731 0.0374549036 0.1805667675 -0.0495767233]' 1e-10 
vcovR =  [0.0200086273  -0.0136717540 0.0084815090  -0.0093192029 -0.0114417825 -0.0113346813 0.0041102304  0.0053787287  0.0050521693 
 -0.0136717540 0.0185396965  -0.0169114682 0.0020645005  0.0036352899  0.0026227856 -0.0012281620 -0.0018838231 -0.0014800242
  0.0084815090  -0.0169114682 0.0174647413  0.0016403871  0.0005624488  0.0015003301 -0.0005714235 -0.0003201257 -0.0005423354
 -0.0093192029 0.0020645005  0.0016403871  0.0167953394  0.0078012534  0.0086329399 -0.0080782771 -0.0037333495 -0.0039327836
 -0.0114417825 0.0036352899  0.0005624488  0.0078012534  0.0490482083  0.0092203882 -0.0033670465 -0.0191567265 -0.0040068947
 -0.0113346813 0.0026227856  0.0015003301  0.0086329399  0.0092203882  0.0331395502 -0.0037513830 -0.0041592743 -0.0146108207
	0.0041102304  -0.0012281620 -0.0005714235 -0.0080782771 -0.0033670465 -0.0037513830 0.0045172675  0.0018165174  0.0020857846
  0.0053787287  -0.0018838231 -0.0003201257 -0.0037333495 -0.0191567265 -0.0041592743 0.0018165174  0.0093292284  0.0020427637
	0.0050521693  -0.0014800242 -0.0005423354 -0.0039327836 -0.0040068947 -0.0146108207 0.0020857846  0.0020427637  0.0074817942]
@test_approx_eq_eps vcov(fitBM) vcovR 1e-10 
@test_approx_eq_eps nobs(fitBM) 100.0 1e-10
@test_approx_eq_eps sum(residuals(fitBM)) 0.6352899255 1e-10
@test_approx_eq_eps dof_residual(fitBM) 91.0 1e-10
@test_approx_eq_eps sigma2_estim(fitBM) 0.0003025014 1e-10
@test_approx_eq_eps stderr(fitBM) [0.1414518551 0.1361605540 0.1321542330 0.1295968341 0.2214683008 0.1820427154 0.0672106202 0.0965879311 0.0864973651] 1e-10
@test_approx_eq_eps confint(fitBM)[:,1] [2.5115945339  -0.4715366529 0.7207474097  -0.3595508202 -0.8102854443 -0.2058583178 -0.0960507369 -0.0112932922 -0.2213931131] 1e-10
@test_approx_eq_eps confint(fitBM)[:,2] [3.0735480006 0.0693957746 1.2457637082 0.1553055609 0.0695537019 0.5173526640 0.1709605441 0.3724268272 0.1222396666] 1e-10
# @test_approx_eq_eps model_response(fitBM)[fitbis.model.ind] model_response(fitbis)
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
#	dat <- read.csv(file = "GA_Anolis_traits.csv")
#	geo <- read.csv(file = "GA_Anolis_biogeography.csv")
#
#	dat <- merge(dat, geo, by = "species")
#	# keep only toe and hand length
#	dat <- dat[, c("species", "AVG.SVL", "AVG.ltoe.IV", "AVG.lfing.IV", "region")]
#	colnames(dat)[1] <- "tipsNames"
#
#	write.table(dat,
#             file = "lizard_trait.txt",
#						  sep = ",", row.names = FALSE)
#
#	## Tree
#	phy <- read.tree(file = "GA_Anolis_MCC.tre")
#
#	write.tree(phy, file = "lizards_tree.txt", append = FALSE,
#            digits = 10, tree.names = FALSE)
#
#	rownames(dat) <- dat$tipsNames
#	dat <- dat[, -1]
#	dat$region <- as.factor(dat$region)
#
#	## Fit
#	fitphylolm <- phylolm(AVG.SVL ~ 1 + AVG.ltoe.IV + AVG.lfing.IV * region, dat, phy, model = "BM")
#	## Quantities to compare
#	sprintf("%.10f", fitphylolm$logLik)
#	sprintf("%.10f", summary(fitphylolm)$df) # df
#	sprintf("%.10f", fitphylolm$aic)
#	sprintf("%.10f", fitphylolm$coefficients) # coef
#	matrix(sprintf("%.10f", fitphylolm$vcov), 9, 9) # vcov
#	sprintf("%.10f", fitphylolm$n) # nobs
#	sprintf("%.10f", sum(fitphylolm$residuals)) # residuals (sum)
#	sprintf("%.10f", fitphylolm$n - fitphylolm$d) # df residuals
#	sprintf("%.10f", fitphylolm$sigma2) # sigma 2
#	sprintf("%.10f", summary(fitphylolm)$coefficients[,2]) # std error
#	sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[, 2] * qt(0.025, 91))
#	sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[, 2] * qt(0.975, 91))

###############################################################################
## Lizard dataset - lambda
###############################################################################

## Export "lizard" dataset (Mahler et al 2013)
phy = readTopology(joinpath(Pkg.dir("PhyloNetworks"), "examples", "lizard_tree.txt"));
dat = readtable(joinpath(Pkg.dir("PhyloNetworks"), "examples", "lizard_trait.txt"));
dat[:region] = PooledDataArray(dat[:region]); # Pool by region

## Fit lambda
fitLambda = phyloNetworklm(AVG_SVL ~ AVG_ltoe_IV + AVG_lfing_IV * region, dat, phy, model = "lambda")

# Tests against results obtained with geiger::fitContinuous or phylolm::phylolm
@test_approx_eq_eps lambda_estim(fitLambda) 0.9982715594 1e-5
@test_approx_eq_eps loglikelihood(fitLambda) 105.1769275564 1e-8
@test_approx_eq_eps dof(fitLambda) 11.0 1e-10
@test_approx_eq_eps aic(fitLambda) -188.3538551128 1e-8 
# @test_approx_eq_eps aicc(fitBM)
@test_approx_eq_eps coef(fitLambda) [2.7940573420 -0.2066584606 0.9897083949 -0.1004840950 -0.3677991157 0.1576743022 0.0367633665 0.1792502383 -0.0505291142]' 1e-5 
vcovR =  [0.0200251600  -0.0137474015 0.0085637021  -0.0092973836 -0.0114259722 -0.0113056243 0.0041037877  0.0053740100  0.0050429112 
					-0.0137474015 0.0186885224  -0.0170645512 0.0020509207  0.0036334103 0.0026066694  -0.0012237488 -0.0018826836 -0.0014743137
					0.0085637021  -0.0170645512 0.0176200143  0.0016494733  0.0005604169 0.0015116125  -0.0005735573 -0.0003192320 -0.0005457420
					-0.0092973836 0.0020509207  0.0016494733  0.0167461876  0.0077885115 0.0086173037  -0.0080563819 -0.0037287856 -0.0039275469
					-0.0114259722 0.0036334103  0.0005604169  0.0077885115  0.0490092393 0.0092036032  -0.0033631662 -0.0191657329 -0.0040017905
					-0.0113056243 0.0026066694  0.0015116125  0.0086173037  0.0092036032 0.0330248707  -0.0037465110 -0.0041543671 -0.0145663751
					0.0041037877  -0.0012237488 -0.0005735573 -0.0080563819 -0.0033631662 -0.0037465110 0.0045042057  0.0018142470  0.0020823721 
					0.0053740100  -0.0018826836 -0.0003192320 -0.0037287856 -0.0191657329 -0.0041543671 0.0018142470  0.0093334212  0.0020404652 
					0.0050429112  -0.0014743137 -0.0005457420 -0.0039275469 -0.0040017905 -0.0145663751 0.0020823721  0.0020404652  0.0074600880]
@test_approx_eq_eps vcov(fitLambda) vcovR 1e-7 
@test_approx_eq_eps nobs(fitLambda) 100.0 1e-10
@test_approx_eq_eps sum(residuals(fitLambda)) 0.6369008979 1e-6
@test_approx_eq_eps dof_residual(fitLambda) 91.0 1e-10
@test_approx_eq_eps sigma2_estim(fitLambda) 0.0003009914 1e-9
@test_approx_eq_eps stderr(fitLambda) [0.1415102824 0.1367059706 0.1327404019 0.1294070617 0.2213803048 0.1817274626 0.0671133793 0.0966096332 0.0863718011] 1e-6
@test_approx_eq_eps confint(fitLambda)[:,1] [2.5129645499 -0.4782080775 0.7260358930 -0.3575353260 -0.8075438955 -0.2033049779 -0.0965491169 -0.0126529301 -0.2220960868] 1e-5
@test_approx_eq_eps confint(fitLambda)[:,2] [3.0751501341 0.0648911562 1.2533808968 0.1565671360 0.0719456641 0.5186535822 0.1700758500 0.3711534067 0.1210378584] 1e-5
# @test_approx_eq_eps model_response(fitLambda)[fitbis.model.ind] model_response(fitbis)
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
#	dat <- read.csv(file = "GA_Anolis_traits.csv")
#	geo <- read.csv(file = "GA_Anolis_biogeography.csv")
#
#	dat <- merge(dat, geo, by = "species")
#	# keep only toe and hand length
#	dat <- dat[, c("species", "AVG.SVL", "AVG.ltoe.IV", "AVG.lfing.IV", "region")]
#	colnames(dat)[1] <- "tipsNames"
#
#	write.table(dat,
#             file = "lizard_trait.txt",
#						  sep = ",", row.names = FALSE)
#
#	## Tree
#	phy <- read.tree(file = "GA_Anolis_MCC.tre")
#
#	write.tree(phy, file = "lizards_tree.txt", append = FALSE,
#            digits = 10, tree.names = FALSE)
#
#	rownames(dat) <- dat$tipsNames
#	dat <- dat[, -1]
#	dat$region <- as.factor(dat$region)
#
#	## Fit
#	fitphylolm <- phylolm(AVG.SVL ~ 1 + AVG.ltoe.IV + AVG.lfing.IV * region, dat, phy, model = "lambda")
#	## Quantities to compare
#	sprintf("%.10f", fitphylolm$logLik)
#	sprintf("%.10f", summary(fitphylolm)$df) # df
#	sprintf("%.10f", fitphylolm$aic)
#	sprintf("%.10f", fitphylolm$coefficients) # coef
#	matrix(sprintf("%.10f", fitphylolm$vcov), 9, 9) # vcov
#	sprintf("%.10f", fitphylolm$n) # nobs
#	sprintf("%.10f", sum(fitphylolm$residuals)) # residuals (sum)
#	sprintf("%.10f", fitphylolm$n - fitphylolm$d) # df residuals
#	sprintf("%.10f", fitphylolm$sigma2) # sigma 2
#	sprintf("%.10f", summary(fitphylolm)$coefficients[,2]) # std error
#	sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[, 2] * qt(0.025, 91))
#	sprintf("%.10f", coef(fitphylolm) + summary(fitphylolm)$coefficients[, 2] * qt(0.975, 91))
