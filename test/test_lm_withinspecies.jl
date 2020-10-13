## tests of phyloNetworklm_wsp

#= to work on:
- hard-code the data for the tests: avoid randomness
- use some other network: star tree is an edge case with V=Identity,
  and has no reticulation. too easy.
- check the estimation with independent implementation, like on a tree with
  pgls.SEy in the R package phytools. document here the R code (within comment block)
- check & update functions like: nobs, dof, dof_residual, deviance, residuals, sigma2_estim, etc.
- finish handling of missing data: both within a species, and if a species missing entirely
- add documentation for within-species variation, including assumptions on X
  (no variation in X! but what if variation seen in practice, and what if
  missing in X but not Y for some individuals, or missing Y but not X in other individuals?)
=#

@testset "phyloNetworklm with measurement error on Star-topology" begin

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

gdf = groupby(df, :species)
# line below: does not handle missing values
# combine(gdf, nrow, valuecols(gdf) .=> mean, valuecols(gdf) .=> std)
# so doing more by hand:
df_r = combine(gdf, :trait1 => (x -> mean(skipmissing(x))) => :trait1,
                    :trait2 => (x -> mean(skipmissing(x))) => :trait2,
                    :trait3 => (x -> mean(skipmissing(x))) => :trait3,
                    :trait3 => (x ->  std(skipmissing(x))) => :trait3_sd,
                    :trait3 => (x -> sum(.!ismissing.(x))) => :trait3_n)
#=
4×6 DataFrame
│ Row │ species │ trait1  │ trait2  │ trait3  │ trait3_sd │ trait3_n │
│     │ String  │ Float64 │ Float64 │ Float64 │ Float64   │ Int64    │
├─────┼─────────┼─────────┼─────────┼─────────┼───────────┼──────────┤
│ 1   │ t1      │ 2.8564  │ -2.0935 │ 14.9242 │ 1.22004   │ 2        │
│ 2   │ t2      │ 2.8457  │ 0.4955  │ 14.4777 │ 1.1719    │ 3        │
│ 3   │ t3      │ 0.4197  │ -2.1977 │ 8.9648  │ 0.273797  │ 3        │
│ 4   │ t4      │ 2.2359  │ -2.618  │ 15.3935 │ 0.402418  │ 3        │
=#
fit = phyloNetworklm(@formula(trait3 ~ trait1), df, starnet;
                        tipnames=:species, msr_err=true)
fit_r = phyloNetworklm(@formula(trait3 ~ trait1), df_r, starnet;
                        tipnames=:species,
                        msr_err=true, response_std=true)
# relative tolerance for inexact equality comparison set to 1%
@test isapprox(coef(fit), coef(fit_r), rtol=0.01)
@test isapprox(sigma2_estim(fit), sigma2_estim(fit_r), rtol=0.01)
@test isapprox(wspvar_estim(fit), wspvar_estim(fit_r), rtol=0.01)

end
