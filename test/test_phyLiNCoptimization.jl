@testset "phyLiNC" begin

fastasimple = joinpath(@__DIR__, "..", "examples", "simple.aln")
fasta8sites = joinpath(@__DIR__, "..", "examples", "Ae_bicornis_8sites.aln") # 8 sites only
fastacontig = joinpath(@__DIR__, "..", "examples", "Ae_bicornis_Tr406_Contig10132.aln")
fasta1missing = joinpath(@__DIR__, "..", "examples", "simple_missingone.aln")
fastaindiv    = joinpath(@__DIR__, "..", "examples", "individuals.aln")
mappingfile   = joinpath(@__DIR__, "..", "examples", "mappingIndividuals.csv")
#= for local trials only, not for automatic testing
  pkgpath = dirname(Base.find_package("PhyloNetworks"))
  fastasimple   = abspath(joinpath(pkgpath, "..", "examples", "simple.aln"))
  fasta8sites   = abspath(joinpath(pkgpath, "..", "examples", "Ae_bicornis_8sites.aln"))
  fastacontig   = abspath(joinpath(pkgpath, "..", "examples", "Ae_bicornis_Tr406_Contig10132.aln"))
  fasta1missing = abspath(joinpath(pkgpath, "..", "examples", "simple_missingone.aln"))
  fastaindiv    = abspath(joinpath(pkgpath, "..", "examples", "individuals.aln"))
  mappingfile   = abspath(joinpath(pkgpath, "..", "examples", "mappingIndividuals.csv"))
=#
emptyconstraint = PhyloNetworks.TopologyConstraint[]

@testset "optimize local BL & gammas, simple example" begin
Random.seed!(99)
net_simple = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
obj = PhyloNetworks.StatisticalSubstitutionModel(net_simple, fastasimple, :JC69)

## Local BL
obj.loglik = -Inf64
obj.net.edge[2].length= 0.0 # unzipping
e = PhyloNetworks.optimizelocalBL_LiNC!(obj, obj.net.edge[6],
        PhyloNetworks.CacheLengthLiNC(obj, 1e-6,1e-6,1e-2,1e-3, 10))
@test length(e) == 4     # not 5: edge below hybrid was excluded
@test obj.loglik > -40.0 # starting likelihood = -43.95468386633092

# ## Local Gamma
γcache = PhyloNetworks.CacheGammaLiNC(obj)
hybridmajorparent = PhyloNetworks.getMajorParentEdge(obj.net.hybrid[1])
@test_nowarn PhyloNetworks.optimizelocalgammas_LiNC!(obj, hybridmajorparent, 1e-6, γcache)
@test hybridmajorparent.gamma != 0.9
@test PhyloNetworks.getMinorParentEdge(obj.net.hybrid[1]).gamma != 0.1
end

@testset "optimize local BL & gammas, complex network and 8 sites" begin
Random.seed!(98)
dna_dat, dna_weights = readfastatodna(fasta8sites, true); # 22 species, 3 hybrid nodes, 103 edges
net = readTopology("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),(((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);");
PhyloNetworks.fuseedgesat!(93, net)
obj = (@test_logs (:warn, r"taxa with no data") PhyloNetworks.StatisticalSubstitutionModel(net, fasta8sites, :JC69))
@test length(obj.net.leaf) == 22
io = IOBuffer();
PhyloNetworks.showdata(io, obj)
@test String(take!(io)) == "data:\n  22 species\n  8 sites"
PhyloNetworks.showdata(io, obj, true)
@test String(take!(io)) ==
"data:
  22 species
  8 sites
  0 sites with no data (0.0%)
  2 invariant sites (25.0%)
  6 sites with 2 distinct states (75.0%)
  6 parsimony-informative sites (75.0%)
  6 sites with 1 or more missing values (75.0%)
  3.41% missing values overall"
close(io)

preorder!(obj.net)
PhyloNetworks.checknetwork_LiNC!(obj.net, 3, true, true, emptyconstraint)
# checknetwork removes degree-2 nodes (including root) and 2- and 3-cycles
# and requires that the network is preordered.
PhyloNetworks.updateSSM!(obj, true; constraints=emptyconstraint)
PhyloNetworks.startingBL!(obj.net, obj.trait, obj.siteweight)
PhyloNetworks.unzip_canonical!(obj.net)
## Local BL
lcache = PhyloNetworks.CacheLengthLiNC(obj, 1e-6,1e-6,1e-2,1e-3, 5)
obj.loglik = +Inf64
PhyloNetworks.optimizelocalBL_LiNC!(obj, obj.net.edge[27], lcache)
obj.loglik = -Inf64
@test_nowarn PhyloNetworks.optimizelocalBL_LiNC!(obj, obj.net.edge[27], lcache)
@test obj.net.edge[27].length ≈ 1e-8 atol=.0001
@test obj.net.edge[32].length ≈ 0.032 atol=.0005
@test obj.loglik > -80.4
# Local BL constrained edge
lengthe = obj.net.edge[44].length
@test_nowarn PhyloNetworks.optimizelocalBL_LiNC!(obj, obj.net.edge[44], lcache)
@test obj.net.edge[44].length == 0.0
@test obj.loglik > -78.95

# ## Local Gamma
# edge[4] = major parent edge of hybrid[1]
γcache = PhyloNetworks.CacheGammaLiNC(obj)
@test_nowarn PhyloNetworks.optimizelocalgammas_LiNC!(obj, obj.net.edge[4], 1e-6,γcache)
@test obj.net.edge[4].gamma == 0.0
@test PhyloNetworks.getMajorParentEdge(obj.net.hybrid[1]).gamma == 1.0

# gamma at a hybrid ladder: when some displayed trees don't have the focus edge
# 2 unzipped reticulations in a hybrid ladder, reasonable (small) branch lengths
net = readTopology("(#H2:0.0001::0.0244,((C:0.0262,((B:0.0)#H1:0.0::0.6)#H2:0.03::0.9756):0.4812,(#H1:0.0001::0.4,A:0.1274):0.0001):0.0274,D:0.151);")
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastasimple, :JC69)
obj.trait[1][3] = 3 # to create discordance across sites
obj.trait[4][4] = 2 #    for estimated γ to be within (0,1)
γcache = PhyloNetworks.CacheGammaLiNC(obj)
PhyloNetworks.discrete_corelikelihood!(obj) # -29.82754754416619
ll = PhyloNetworks.optimizegamma_LiNC!(obj, obj.net.edge[4], .001, γcache, 3)
@test ll ≈ -29.655467369763585
@test obj.net.edge[4].gamma ≈ 0.8090635871910823
@test γcache.hase == [true, true, false]
@test all(γcache.clikn .== 0.0)
ll = PhyloNetworks.optimizegamma_LiNC!(obj, obj.net.edge[1], .001, γcache, 3)
@test ll ≈ -29.468517891012983
@test obj.net.edge[1].gamma ≈ 0.19975558937688737
@test γcache.hase[[1,2]] == [false,true]
@test γcache.hase[3] === missing
@test γcache.clikn[[6,7]] ≈ [0.25, 0.0007] rtol=0.01
end #of local branch length and gamma optimization with localgamma! localBL! with 8 sites

@testset "global branch length and gamma optimization" begin
# to run locally on complex network:
# net = readTopology("(H_vulgare_HVens23:0.5,(((Ae_speltoides_Tr251:0.5):0.5,(Ae_mutica_Tr237:0.0)#H4:1.0::0.7):0.5,((((((Ae_caudata_Tr139:0.5,Ae_caudata_Tr275:0.5):0.0)#H1:1.0::0.7,#H2:1.0::0.3):0.5,#H1:1.0::0.3):0.5,((Ae_comosa_Tr271:0.5,Ae_comosa_Tr272:0.5):0.5,((Ae_uniaristata_Tr403:0.5,Ae_uniaristata_Tr357:0.5):0.5,Ae_uniaristata_Tr402:0.5):0.5):0.5):0.5,(((Ae_tauschii_Tr352:0.5,Ae_tauschii_Tr351:0.5):0.5,Ae_tauschii_Tr125:0.5):0.5,(((((((Ae_longissima_Tr241:0.5,Ae_longissima_Tr242:0.5):0.5,Ae_longissima_Tr355:0.5):0.5,Ae_sharonensis_Tr265:0.5):0.5,((Ae_bicornis_Tr408:0.5,Ae_bicornis_Tr407:0.5):0.5,Ae_bicornis_Tr406:0.5):0.5):0.5,(Ae_searsii_Tr164:0.5,Ae_searsii_Tr165:0.5):0.5):0.0)#H2:1.0::0.7,#H4:1.0::0.3):0.5):0.5):0.5):0.5);");
# obj = PhyloNetworks.StatisticalSubstitutionModel(net, fasta8sites, :JC69);
net = readTopology("(((A:0.5,(B:0.0)#H1:1.0::0.9):0.5,(C:0.5,#H1:1.0::0.1):0.5):0.5,D:0.5);")
# branch lengths set to 0.5, then unzipped -> some BL are 0, some 1, most 0.5
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastasimple, :JC69);

## optimizeBL
Random.seed!(5);
lcache = PhyloNetworks.CacheLengthLiNC(obj, 1e-6,1e-6,1e-2,1e-3, 5)
obj.loglik = -Inf64
@test_nowarn PhyloNetworks.optimizealllengths_LiNC!(obj, lcache);
@test all(e.length != 1.0 for e in obj.net.edge)
@test [e.length for e in obj.net.edge] ≈ [0.3727, 0.0, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8, 0.5201, 1.0e-8] rtol=.001

## optimizegammas -- and delete hybrid edges with γ=0
γcache = PhyloNetworks.CacheGammaLiNC(obj)
@test_nowarn PhyloNetworks.optimizeallgammas_LiNC!(obj,1e-6,γcache,100)
@test obj.net.numHybrids == 0
end

@testset "data to SSM pruning: complex network" begin
net = readTopology("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),(((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);");
PhyloNetworks.fuseedgesat!(93, net)
for edge in net.edge # reset network
    setLength!(edge,1.0)
end
for h in net.hybrid
    setGamma!(PhyloNetworks.getMajorParentEdge(h),0.6)
end
obj = (@test_logs (:warn, r"pruned") PhyloNetworks.StatisticalSubstitutionModel(net, fasta8sites, :JC69))
@test length(obj.net.leaf) == 22
@test length(obj.net.edge) == 52
@test length(obj.net.hybrid) == 3
@test !PhyloNetworks.hashybridladder(obj.net)
end

@testset "checknetwork LiNC" begin
tree = readTopology("(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);");
@test any(length(n.edge) == 2 for n in tree.node) # one node of degree 2
preorder!(tree)
PhyloNetworks.checknetwork_LiNC!(tree, 1, true, true)
@test all(length(n.edge) != 2 for n in tree.node) # no nodes of degree 2
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
@test_throws ErrorException PhyloNetworks.checknetwork_LiNC!(net, 0, true, true)

netstr = "(#H2:0.1::0.2,((C:0.2,((B:0.3)#H1:0.4)#H2:0.5::0.8):0.6,(#H1:0.7,((A1:0.8)#H3:0.01,(A2:0.9,#H3:0.02):0.03):1.0):1.1):1.2,O:1.3);"
net = readTopology(netstr)
# 2 reticulation in a hybrid ladder, and another isolated reticulation
undoinfo = PhyloNetworks.unzip_canonical!(net)
@test all(PhyloNetworks.getChildEdge(h).length == 0.0 for h in net.hybrid) # unzipped
@test writeTopology(net, round=true) == "(#H2:0.8::0.2,((C:0.2,((B:0.0)#H1:0.0)#H2:1.2::0.8):0.6,(#H1:1.0,((A1:0.0)#H3:0.81,(A2:0.9,#H3:0.82):0.03):1.0):1.1):1.2,O:1.3);"
PhyloNetworks.rezip_canonical!(undoinfo...)
@test writeTopology(net, round=true) == netstr
end

@testset "update root in SSM displayed trees" begin
# W structure
net = readTopology("(C:0.0262,(B:0.0)#H2:0.03::0.9756,(((D:0.1,A:0.1274):0.0)#H1:0.0::0.6,(#H2:0.0001::0.0244,#H1:0.151::0.4):0.0274):0.4812);")
# "((C:0.0262,(B:0.0)#H2:0.03::0.9756):0.4812,((D:0.1,A:0.1274):0.0)#H1:0.0::0.6,(#H2:0.0001::0.0244,#H1:0.151::0.4):0.0274);")
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastasimple, :JC69)
@test [t.node[t.root].number for t in obj.displayedtree] == [5,5,5,5]
# obj.displayedtree[1]: (C:0.026,(B:0.0):0.03,(((D:0.1,A:0.127):0.0):0.0):0.481);
# move the root to place the W structure at the root:
# the network's root node will be missing from some displayed trees.
rootatnode!(obj.net, 7) # node 7 = tree node whose 2 children are both hybrids
# "(#H2:0.0001::0.0244,((C:0.0262,(B:0.0)#H2:0.03::0.9756):0.4812,((D:0.1,A:0.1274):0.0)#H1:0.0::0.6):0.0274,#H1:0.151);"
PhyloNetworks.updateSSM_root!(obj) # re-root displayed trees in the same way
@test [t.node[t.root].number for t in obj.displayedtree] == [6,7,7,7]
@test writeTopology(obj.displayedtree[1]) == "(((D:0.1,A:0.1274):0.0)H1:0.0,(C:0.0262,(B:0.0)H2:0.03):0.4812);"
end

@testset "skip γ and lengths optimization when needed" begin
# W structure, with middle γs = 0
net = readTopology("(C:0.0262,(B:0.0)#H2:0.03::1,(((D:0.1,A:0.1274):0.0)#H1:0.004::1,(#H2:0.0001,#H1:0.151):0.0274):0.4812);")
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastasimple, :JC69)
for i in [8,9] setGamma!(obj.net.edge[i], 0.0); end
PhyloNetworks.updateSSM!(obj)
lcache = PhyloNetworks.CacheLengthLiNC(obj, 1e-6,1e-6,1e-2,1e-3, 5);
e = PhyloNetworks.optimizelocalBL_LiNC!(obj, obj.net.edge[10], lcache)
@test length(e) == 3 # not 5: the 2 edges with γ = 0 were excluded
@test obj.net.edge[10].length == 0.0274
# hybrid ladder, with lower γ = 0
net = readTopology("(#H2:0.0001::0.0244,((C:0.0262,((B:0.0)#H1:0.0::0.6)#H2:0.03::0.9756):0.4812,(#H1:0.0001::0.4,A:0.1274):0.0001):0.0274,D:0.151);")
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastasimple, :JC69)
setGamma!(obj.net.edge[4], 0.0); PhyloNetworks.updateSSM!(obj)
e = PhyloNetworks.optimizelocalBL_LiNC!(obj, obj.net.edge[5], lcache)
@test length(e) == 4
@test [obj.net.edge[i].length for i in [1,5]] == [0.0001, 0.03]
γcache = PhyloNetworks.CacheGammaLiNC(obj);
ll = PhyloNetworks.optimizegamma_LiNC!(obj, obj.net.edge[5], .001, γcache, 3)
@test ll ≈ -31.124547305074074 # same as before optimization
end

@testset "optimizestructure with simple example" begin
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastasimple, :JC69; maxhybrid=1)
PhyloNetworks.checknetwork_LiNC!(obj.net, 1, true, true)
PhyloNetworks.updateSSM!(obj, true; constraints=emptyconstraint)
PhyloNetworks.startingBL!(obj.net, obj.trait, obj.siteweight)
PhyloNetworks.unzip_canonical!(obj.net)
PhyloNetworks.discrete_corelikelihood!(obj)
@test obj.loglik ≈ -29.7762035
maxmoves = 2
Random.seed!(96)
γcache = PhyloNetworks.CacheGammaLiNC(obj)
lcache = PhyloNetworks.CacheLengthLiNC(obj, 1e-6,1e-6,1e-2,1e-3, 5)
PhyloNetworks.optimizestructure!(obj, maxmoves, 1, true, true, 0,100,
                                emptyconstraint, 1e-6, γcache, lcache)
@test obj.loglik > -27.42

# allow hybrid ladders
Random.seed!(110)
PhyloNetworks.optimizestructure!(obj, maxmoves, 1, true, false, 0,100,
                                emptyconstraint, 1e-6, γcache, lcache)
@test obj.loglik > -27.42
end # of optimizestructure with simple example

@testset "phyLiNCone with simple net, no constraints" begin
no3cycle = true
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");
seed = 102
for nohybridladder in [true, false]
    obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastasimple, :JC69)
    PhyloNetworks.checknetwork_LiNC!(obj.net, 1, no3cycle, nohybridladder)
    PhyloNetworks.updateSSM!(obj, true; constraints=emptyconstraint)
    PhyloNetworks.startingBL!(obj.net, obj.trait, obj.siteweight)
    PhyloNetworks.unzip_canonical!(obj.net)
    obj.loglik = -Inf # missing otherwise, which would cause an error below
    nullio = open("/dev/null", "w")
    γcache = PhyloNetworks.CacheGammaLiNC(obj)
    lcache = PhyloNetworks.CacheLengthLiNC(obj, 1e-6,1e-6,1e-2,1e-3, 5)
    @test_nowarn PhyloNetworks.phyLiNCone!(obj, 1, no3cycle,
            nohybridladder, 3, 2, false, false,
            nullio, seed, 0.5, emptyconstraint,
            1e-2, 1e-2, 1e-2, 1e-2, 0.0, 25.0, 0.01,0.9,
            γcache, lcache)
    @test obj.loglik > -27.45
end
end

@testset "phyLiNC no constraints: HKY, rate variation" begin
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");
obj = @test_nowarn PhyloNetworks.phyLiNC(net, fastasimple, :JC69, :G, 2; maxhybrid=2, # no missing BLs, so they're not re-estimated
                    no3cycle=true, nohybridladder=true, maxmoves=2,
                    nreject=1, nruns=1, filename="", verbose=false, seed=108)
@test obj.loglik > -27.27
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");
obj = @test_nowarn PhyloNetworks.phyLiNC(net, fastasimple, :HKY85; maxhybrid=2,
                    no3cycle=true, nohybridladder=true, maxmoves=2, probST=1.0, # not enough moves to get back to a good topology
                    nreject=1, nruns=1, filename="phyLiNC2", verbose=false, seed=0)
@test obj.loglik > -24.21
@test read("phyLiNC2.err", String) == ""
@test startswith(read("phyLiNC2.log", String), "PhyLiNC network estimation starting")
rm("phyLiNC2.log")
rm("phyLiNC2.err")

net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");
addprocs(1) # multiple cores
@everywhere using PhyloNetworks
#using Distributed; @everywhere begin; using Pkg; Pkg.activate("."); using PhyloNetworks; end
originalstdout = stdout  # verbose=true below
redirect_stdout(open("/dev/null", "w")) # not portable to Windows
obj = PhyloNetworks.phyLiNC(net, fastasimple, :JC69; maxhybrid=2, no3cycle=true,
                        nohybridladder=true, maxmoves=2, nreject=1, nruns=2,
                        filename="phyLiNCmult", verbose=true, seed=106)
redirect_stdout(originalstdout)
@test obj.loglik > -28.27
rmprocs(workers()) # remove extra processors
@test occursin("using 1 worker", read("phyLiNCmult.log", String))
rm("phyLiNCmult.log")
end

@testset "phyLiNC with simple net and one constraint" begin
net_level1_s = readTopology("(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));") # S1A S1B S1C go on leaf 1
# 3-cycle at degree-2 root -> 2-cycle after root deletion, removed within LiNC
# constraint
net_level1_i, c_species = PhyloNetworks.mapindividuals(net_level1_s, mappingfile)
PhyloNetworks.resetNodeNumbers!(net_level1_i)
net_level1_i.node[22].number = 100
PhyloNetworks.updateconstraints!(c_species, net_level1_i)
@test c_species[1].taxonnums == Set([8,9,100])
@test c_species[1].node.number == 21
@test PhyloNetworks.getParent(net_level1_i.node[22].edge[1]).number == 21

obj = PhyloNetworks.StatisticalSubstitutionModel(net_level1_i,fastaindiv,:JC69,:GI,2)
# obj.net = deepcopy of input net, so we need to rebuild the constraints
c_species[1] = PhyloNetworks.TopologyConstraint(0x01, c_species[1].taxonnames, obj.net)
# obj.net = deepcopy of input net, so we need to rebuild the constraints if done after
@test_logs (:warn, r"no 3-cycle") match_mode=:any PhyloNetworks.checknetwork_LiNC!(obj.net, 2,
                                                    true, true, c_species, true)
PhyloNetworks.updateSSM!(obj, true; constraints=emptyconstraint)

for e in obj.net.edge e.length = 0.1; end # was -1.0 for missing
PhyloNetworks.startingBL!(obj.net, obj.trait, obj.siteweight)
PhyloNetworks.unzip_canonical!(obj.net)
PhyloNetworks.setalpha!(obj.ratemodel, 0.48438)
obj.loglik = -Inf # actual likelihood -56.3068141288164. Need something non-missing
seed = 103
nullio = open("/dev/null", "w")
γcache = PhyloNetworks.CacheGammaLiNC(obj)
lcache = PhyloNetworks.CacheLengthLiNC(obj, 1e-2,1e-2,1e-2,1e-2, 5)
@test_nowarn PhyloNetworks.phyLiNCone!(obj, 2, true, true,
        3, 2, false, false, nullio,
        seed, 0.5, c_species, 1e-2, 1e-2,
        1e-2, 1e-2, 0.0,50.0, 0.01,.9, γcache, lcache)
@test obj.loglik > -65.0

obj = PhyloNetworks.phyLiNC(net_level1_s, # missing BLs, so BLs are re-estimated before starting
            fastaindiv, :JC69, :Inv; maxhybrid=2, no3cycle=true, nohybridladder=true,
            verbose=false, filename="", speciesfile=mappingfile, seed=138, nruns=1,
            maxmoves=10, nreject=2)
@test obj.loglik > -67.7 # -69.83824 with :noRV
@test obj.ratemodel.pinv[1] > 0.22 # 0.23753
# test that species stayed together after optimization, as the only polytomy
function polytomyS1(node)
    length(node.edge) > 3 || return false
    return Set(n.name for n in PhyloNetworks.getChildren(node)) == Set(["S1A", "S1B", "S1C"])
end
@test sum(polytomyS1(nod) for nod in obj.net.node) == 1
end

@testset "optimizeparameters" begin
net_simple = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
obj = PhyloNetworks.optimizeparameters(net_simple, fastasimple, :HKY85, :G, "", 123) # no file created
@test obj.net.edge[5].length ≈ 4.001604835835655 atol=.0001
@test obj.ratemodel.alpha[1] ≈ 0.2017 atol=.0001
@test obj.loglik ≈ -24.1758 atol=.0001

obj = PhyloNetworks.optimizeparameters(net_simple, fastasimple, :HKY85, :I, "", 123) # no file created
@test obj.loglik ≈ -23.04696 atol=.0001 # better because :I is a better fit for this data (50% invariant sites)
end

@testset "neighbornets" begin
seed = 1234
n6h1 = readTopology("(((1:0.013635011856564799,2:0.013337412436279021):0.02379142358426221,((3:0.014458424150016028,4:0.015899483647102967):0.0)#H1:0.01602497696107314::0.7876602143873201):0.01491008237330521,(#H1:0.014834901504447811::0.2123397856126799,5:0.005527884554698987):0.055340371373774774,6:0.0620387315019168);")
neighbors, distances = PhyloNetworks.neighbornets(n6h1, true, true)
# confirm that the arrays have the same order
@test hardwiredClusterDistance(n6h1, readTopology(neighbors[2]), true) == distances[2]
@test hardwiredClusterDistance(n6h1, readTopology(neighbors[4]), true) == distances[4]
@test length(neighbors) == 28
# confirm that the first neighbor is as expected
nni!(n6h1, n6h1.edge[3], 0x01, true, true)
@test hardwiredClusterDistance(n6h1, readTopology(neighbors[1]), true) == 0

n6h1 = readTopology("(((1:0.013635011856564799,2:0.013337412436279021):0.02379142358426221,((3:0.014458424150016028,4:0.015899483647102967):0.0)#H1:0.01602497696107314::0.7876602143873201):0.01491008237330521,(#H1:0.014834901504447811::0.2123397856126799,5:0.005527884554698987):0.055340371373774774,6:0.0620387315019168);")
neighbors, distances = PhyloNetworks.neighbornets(n6h1, true, false)
# confirm that the arrays have the same order
@test hardwiredClusterDistance(n6h1, readTopology(neighbors[5]), true) == distances[5]
@test hardwiredClusterDistance(n6h1, readTopology(neighbors[6]), true) == distances[6]
@test length(neighbors) == 40 # more neighbors because 3-cycles are allowed

n6h2 = readTopology("(((((1:1.4,#H1:0.7::0.2):1.3,((2:0.7)#H1:1.1,3:1.8):0.9):1.2,(4:3.4)#H2:0.5::0.7):2.0,(#H2:0.5,5:3.9):2.0):1.5,6:7.4);")
neighbors, distances = PhyloNetworks.neighbornets(n6h2, true, true)
# confirm that the arrays have the same order
@test hardwiredClusterDistance(n6h2, readTopology(neighbors[10]), true) == distances[10]
@test hardwiredClusterDistance(n6h2, readTopology(neighbors[11]), true) == distances[11]
@test length(neighbors) == 40
# confirm that the first neighbor is as expected
nni!(n6h2, n6h2.edge[2], 0x02, true, true)
@test hardwiredClusterDistance(n6h2, readTopology(neighbors[1]), true) == 0

n6h2 = readTopology("(((((1:1.4,#H1:0.7::0.2):1.3,((2:0.7)#H1:1.1,3:1.8):0.9):1.2,(4:3.4)#H2:0.5::0.7):2.0,(#H2:0.5,5:3.9):2.0):1.5,6:7.4);")
neighbors, distances = PhyloNetworks.neighbornets(n6h2, false, false)
# confirm that the arrays have the same order
@test hardwiredClusterDistance(n6h2, readTopology(neighbors[5]), true) == distances[5]
@test hardwiredClusterDistance(n6h2, readTopology(neighbors[6]), true) == distances[6]
@test length(neighbors) == 64 # multiple hybrids, so allowing hybrid ladders should increase the number of neighbors
end

@testset "uniqueneighbornets" begin
seed = 2345
n6h1 = readTopology("(((1:0.013635011856564799,2:0.013337412436279021):0.02379142358426221,((3:0.014458424150016028,4:0.015899483647102967):0.0)#H1:0.01602497696107314::0.7876602143873201):0.01491008237330521,(#H1:0.014834901504447811::0.2123397856126799,5:0.005527884554698987):0.055340371373774774,6:0.0620387315019168);")
neighbors, distances = PhyloNetworks.uniqueneighbornets(n6h1, true, true)
@test length(neighbors) == 12
end

@testset "nnistotruenet" begin
seed = 2345
truenet = readTopology("(((1:0.013635011856564799,2:0.013337412436279021):0.02379142358426221,((3:0.014458424150016028,4:0.015899483647102967):0.0)#H1:0.01602497696107314::0.7876602143873201):0.01491008237330521,(#H1:0.014834901504447811::0.2123397856126799,5:0.005527884554698987):0.055340371373774774,6:0.0620387315019168);")
startingnet = readTopology("(((1:0.009221630571539522,2:0.01090284156388857):0.03072521643460218,(3:1.0e-8,(4:0.0)#H1:0.01133912873252381::0.7292558128341733):0.030320601757657505):0.013615526953203836,6:0.07418479480620778,(5:0.016527926411687346,#H1:1.0000000050247593e-8::0.2707441871658267):0.026214884034880433);")

# test rooting check before rooting fix, HWCD = 4
@test hardwiredClusterDistance(startingnet, truenet, true) == 4
@test hardwiredClusterDistance(startingnet, truenet, false) == 2
rootatnode!(startingnet, "6") # confirm startingnet is rooted at outgroup
PhyloNetworks.removedegree2nodes!(startingnet)
@test length(PhyloNetworks.getChildren(startingnet.node[startingnet.root])) == 3

startingnet = readTopology("(((1:0.009221630571539522,2:0.01090284156388857):0.03072521643460218,(3:1.0e-8,(4:0.0)#H1:0.01133912873252381::0.7292558128341733):0.030320601757657505):0.013615526953203836,6:0.07418479480620778,(5:0.016527926411687346,#H1:1.0000000050247593e-8::0.2707441871658267):0.026214884034880433);")
nmoves = PhyloNetworks.nnistotruenet(startingnet, truenet, "6", true, true, 2)
@test nmoves == 1
end

@testset "hybridschangedbynnis" begin
startingnet = readTopology("(((1:0.009221630571539522,2:0.01090284156388857):0.03072521643460218,(3:1.0e-8,(4:0.0)#H1:0.01133912873252381::0.7292558128341733):0.030320601757657505):0.013615526953203836,6:0.07418479480620778,(5:0.016527926411687346,#H1:1.0000000050247593e-8::0.2707441871658267):0.026214884034880433);")
nflipped, nchanged = PhyloNetworks.hybridschangedbynnis(startingnet, true, true, 1)
@test nchanged == 15
@test nflipped == 0

nflipped, nchanged = PhyloNetworks.hybridschangedbynnis(startingnet, false, false, 2)
@test nchanged == 286
@test nflipped == 0
end

end # of overall phyLiNC test set
