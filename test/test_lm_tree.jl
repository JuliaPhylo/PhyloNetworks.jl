# Test of PhyloNetworklm on trees

using PhyloNetworks
using GLM
using DataFrames
using Base.Test

###############################################################################
## Caudata dataset
###############################################################################

## Export "caudata" dataset (from geiger)
phy = readTopology("examples/caudata_tree.txt");
dat = readtable("examples/caudata_trait.txt");

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

###############################################################################
## Lizard dataset
###############################################################################

## Export "lizard" dataset (Mahler et al 2013)
phy = readTopology("examples/lizard_tree.txt");
dat = readtable("examples/lizard_trait.txt");
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
