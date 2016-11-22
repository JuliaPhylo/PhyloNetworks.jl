# Test of PhyloNetworklm on trees

using PhyloNetworks
using GLM
using DataFrames
using Base.Test

## Export "caudata" dataset (from geiger)
phy = readTopology("examples/caudata_tree.txt")
dat = readtable("examples/caudata_trait.txt") 

## Fit a BM
fitBM = phyloNetworklm(trait ~ 1, dat, phy)

# Tests against results obtained with geiger::fitContinuous or phylolm::phylolm
@test_approx_eq_eps loglikelihood(fitBM) -78.9611507833 1e-10 
@test_approx_eq_eps aic(fitBM) 161.9223015666 1e-10 
@test_approx_eq_eps coef(fitBM) 4.6789989001 1e-10 
@test_approx_eq_eps vcov(fitBM) 0.1093144100 1e-10 
@test_approx_eq_eps nobs(fitBM) 197.0 1e-10
@test_approx_eq_eps sum(residuals(fitBM)) -115.5767321312 1e-10
@test_approx_eq_eps dof_residual(fitBM) 196.0 1e-10
@test_approx_eq_eps sigma2_estim(fitBM) 0.0029452097 1e-10
@test_approx_eq_eps stderr(fitBM) 0.3306272978 1e-10
@test_approx_eq_eps confint(fitBM)[1] 4.0269551772 1e-10
@test_approx_eq_eps confint(fitBM)[2] 5.3310426231 1e-10
@test_approx_eq_eps dof(fitBM)  2.0 1e-10
@test_approx_eq_eps aicc(fitBM)  161.9841572367 1e-10
# @test_approx_eq_eps model_response(fitBM)[fitbis.model.ind] model_response(fitbis)
# @test_approx_eq_eps predict(fitBM)[fitbis.model.ind] predict(fitbis)
# @test_approx_eq_eps deviance(fitBM)  deviance(fitbis)
# @test_approx_eq_eps nulldeviance(fitBM)  nulldeviance(fitbis)
# @test_approx_eq_eps nullloglikelihood(fitBM)  nullloglikelihood(fitbis)
# @test_approx_eq_eps r2(fitBM)  r2(fitbis)
# @test_approx_eq_eps adjr2(fitBM)  adjr2(fitbis)
# @test_approx_eq_eps bic(fitBM)  bic(fitbis)
# @test_approx_eq_eps mu_estim(fitBM)  mu_estim(fitbis)
