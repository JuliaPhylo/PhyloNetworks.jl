## tests of phyloNetworkmem

@testset "phyloNetworklm with measurement error on Star-topology" begin
    global truenet

    # Generate Star-topology with n tips in Newick format
    function generateStar(n::Int64)
        topo_newick = reduce(*, ["t$i:1," for i=1:n])[1:(end-1)]
        topo_newick = "("*topo_newick*");"
    end
    
    # Simulate traits on a Star-topology
    function simTraits(net::PhyloNetworks.HybridNetwork,
                       m::Int64,
                       process::PhyloNetworks.ParamsProcess)
        # returns TraitSimulation object
        # M::MatrixTopologicalOrder, params::ParamsProcess,
        # model::AbstractString
        sim = simulate(net, process)
        # returns Array with one column
        trait = sim[:Tips]
        # replicate observations
        reduce((x,y)->vcat(x,y), [fill(obs, m) for obs in trait])
    end
    
    n = 30; m = 5
    truenet = readTopology(generateStar(n))
    trait1 = simTraits(truenet, m, ParamsBM(2, 0.5)) # var: 0.5
    trait2 = simTraits(truenet, m, ParamsBM(-2, 1))

    phylo_noise_var = 2; meas_noise_var = 1
    phylo_noise = simTraits(truenet, m, ParamsBM(0, phylo_noise_var)) # var: 2
    meas_noise = rand(Normal(0, sqrt(meas_noise_var)), n*m)
    trait3 = 10 .+ 2*trait1 + phylo_noise + meas_noise

    labels = reduce((x,y)->vcat(x,y), [fill(lab, m) for lab in truenet.names])
    Y = trait3
    X = fill(1.0, (30, 2)); X[:, 2] = [trait1[i] for i in collect(1:5:150)]

    df = DataFrame(trait1=trait1, trait2=trait2, trait3=trait3,
                   tipNames=labels)

    netnames = truenet.names
    a = length(netnames)
    trt1 = fill(0.0, n); trt1_sd = fill(0.0, n)
    trt2 = fill(0.0, n); trt2_sd = fill(0.0, n)
    trt3 = fill(0.0, n); trt3_sd = fill(0.0, n)
    for i in 1:a
        ind = (labels .== netnames[i])
        trt1[i] = mean(trait1[ind]); trt1_sd[i] = std(trait1[ind])
        trt2[i] = mean(trait2[ind]); trt2_sd[i] = std(trait2[ind])
        trt3[i] = mean(trait3[ind]); trt3_sd[i] = std(trait3[ind])
    end
    speciescts = repeat([m], n)
    df_r = DataFrame(trt1=trt1, trt1_sd=trt1_sd,
                     trt2=trt2, trt2_sd=trt2_sd,
                     trt3=trt3, trt3_sd=trt3_sd,
                     speciescts=speciescts, species=netnames)



    @test try
        fit = phyloNetworklm(@formula(trait3 ~ trait1), df, truenet;
                             msr_err=true)
        fit_r = phyloNetworklm(@formula(trt3 ~ trt1), df_r, truenet;
                               tipnames=:species,
                               msr_err=true, response_std=true)
        # relative tolerance for inexact equality comparison set to 1%
        @test isapprox(coef(fit), coef(fit_r), rtol=0.01) 
        @test isapprox(sigma2_estim(fit), sigma2_estim(fit_r), rtol=0.01)
        @test isapprox(wspvar_estim(fit), wspvar_estim(fit_r), rtol=0.01)
        true
    catch
        false
    end

end
