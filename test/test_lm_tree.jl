# Test of PhyloNetworklm on trees

using PhyloNetworks
using GLM
using DataFrames
using Base.Test

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
# @test_approx_eq_eps model_response(fitBM)[fitbis.model.ind] model_response(fitbis)
# @test_approx_eq_eps predict(fitBM)[fitbis.model.ind] predict(fitbis)
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

