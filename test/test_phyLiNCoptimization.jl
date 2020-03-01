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
lengthe = obj.net.edge[4].length
lengthep = obj.net.edge[4].node[1].edge[1].length
@test_nowarn PhyloNetworks.optimizelocalBL_LiNC!(obj, obj.net.edge[4], false, 1e-6,1e-6,1e-2,1e-3)
@test obj.net.edge[4].length != lengthe
@test obj.net.edge[4].node[1].edge[1].length != lengthep

# ## Local Gamma
hybridmajorparent = PhyloNetworks.getMajorParentEdge(obj.net.hybrid[1])
@test_nowarn PhyloNetworks.optimizelocalgammas_LiNC!(obj, hybridmajorparent, false, 1e-6,1e-6,1e-2,1e-3)
@test hybridmajorparent.gamma != 0.9
@test PhyloNetworks.getMinorParentEdge(obj.net.hybrid[1]).gamma != 0.1
end

@testset "optimize local BL & gammas, complex network and 8 sites" begin
Random.seed!(98)
dna_dat, dna_weights = readfastatodna(fasta8sites, true); # 22 species, 3 hybrid nodes, 103 edges
net = readTopology("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),(((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);");
PhyloNetworks.fuseedgesat!(93, net)
for edge in net.edge # reset network
    setLength!(edge,1.0)
end
for h in net.hybrid
    setGamma!(PhyloNetworks.getMajorParentEdge(h),0.7)
end
obj = (@test_logs (:warn, r"taxa with no data") PhyloNetworks.StatisticalSubstitutionModel(net, fasta8sites, :JC69))
@test length(obj.net.leaf) == 22

## Local BL
lengthe = obj.net.edge[48].length
lengthep = obj.net.edge[48].node[1].edge[1].length
@test_nowarn PhyloNetworks.optimizelocalBL_LiNC!(obj, obj.net.edge[48], false, 1e-6,1e-6,1e-2,1e-3)
@test obj.net.edge[48].length != lengthe
@test obj.net.edge[48].node[1].edge[1].length == 0.0 # below hybrid node

# ## Local Gamma
# edge[4] = major parent edge of hybrid[1]
originalstdout = stdout  # verbose=true below
redirect_stdout(open("/dev/null", "w")) # not portable to Windows
@test_nowarn PhyloNetworks.optimizelocalgammas_LiNC!(obj, obj.net.edge[4], true, 1e-6,1e-6,1e-2,1e-3)
redirect_stdout(originalstdout)
@test obj.net.edge[4].gamma != 0.7
@test PhyloNetworks.getMinorParentEdge(obj.net.hybrid[1]).gamma != 0.3
end #of local branch length and gamma optimization with localgamma! localBL! with 8 sites

@testset "global branch length and gamma optimization" begin
Random.seed!(97)
# to run locally on complex network:
# net = readTopology("(H_vulgare_HVens23:0.5,(((Ae_speltoides_Tr251:0.5):0.5,(Ae_mutica_Tr237:0.0)#H4:1.0::0.7):0.5,((((((Ae_caudata_Tr139:0.5,Ae_caudata_Tr275:0.5):0.0)#H1:1.0::0.7,#H2:1.0::0.3):0.5,#H1:1.0::0.3):0.5,((Ae_comosa_Tr271:0.5,Ae_comosa_Tr272:0.5):0.5,((Ae_uniaristata_Tr403:0.5,Ae_uniaristata_Tr357:0.5):0.5,Ae_uniaristata_Tr402:0.5):0.5):0.5):0.5,(((Ae_tauschii_Tr352:0.5,Ae_tauschii_Tr351:0.5):0.5,Ae_tauschii_Tr125:0.5):0.5,(((((((Ae_longissima_Tr241:0.5,Ae_longissima_Tr242:0.5):0.5,Ae_longissima_Tr355:0.5):0.5,Ae_sharonensis_Tr265:0.5):0.5,((Ae_bicornis_Tr408:0.5,Ae_bicornis_Tr407:0.5):0.5,Ae_bicornis_Tr406:0.5):0.5):0.5,(Ae_searsii_Tr164:0.5,Ae_searsii_Tr165:0.5):0.5):0.0)#H2:1.0::0.7,#H4:1.0::0.3):0.5):0.5):0.5):0.5);");
# obj = PhyloNetworks.StatisticalSubstitutionModel(net, fasta8sites, :JC69);
net = readTopology("(((A:0.5,(B:0.0)#H1:1.0::0.9):0.5,(C:0.5,#H1:1.0::0.1):0.5):0.5,D:0.5);")
# branch lengths set to 0.5, then unzipped -> some BL are 0, some 1, most 0.5
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastasimple, :JC69);

## optimizeBL
Random.seed!(5);
@test_nowarn PhyloNetworks.optimizeBL_LiNC!(obj, obj.net.edge, false,20,1e-2,1e-2,1e-2,1e-2);
@test all(e.length != 1.0 for e in obj.net.edge)

## optimizegammas
@test_nowarn PhyloNetworks.optimizeallgammas_LiNC!(obj,false,100,1e-6,1e-6,1e-2,1e-3)
@test PhyloNetworks.getMajorParentEdge(obj.net.hybrid[1]).gamma != 0.6
@test PhyloNetworks.getMinorParentEdge(obj.net.hybrid[1]).gamma != 0.4
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

@testset "checknetworkbeforeLiNC" begin
tree = readTopology("(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);");
@test any(length(n.edge) == 2 for n in tree.node) # one node of degree 2
preorder!(tree)
PhyloNetworks.checknetwork_LiNC!(tree, 1, true, true)
@test all(length(n.edge) != 2 for n in tree.node) # no nodes of degree 2

net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
@test any(length(n.edge) == 2 for n in net.node) # one node of degree 2
@test_throws ErrorException PhyloNetworks.unzip_canonical!(net)
preorder!(net)
PhyloNetworks.checknetwork_LiNC!(net, 1, true, true)
@test all(length(n.edge) != 2 for n in net.node) # no nodes of degree 2
# removed unzip from checknetwork_LiNC because its in startingBL! @test all(PhyloNetworks.getChildEdge(h).length == 0.0 for h in net.hybrid) # unzipped
@test_throws ErrorException PhyloNetworks.checknetwork_LiNC!(net, 0, true, true)
end

@testset "optimizestructure with simple example" begin
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastasimple, :JC69, 1)
PhyloNetworks.checknetwork_LiNC!(obj.net, 1, true, true)
PhyloNetworks.discrete_corelikelihood!(obj)
@test obj.loglik ≈ -32.934741316276956
maxmoves = 2
Random.seed!(92) # fixit: error with nan loglik if we do instead: Random.seed!(100)
PhyloNetworks.optimizestructure!(obj, maxmoves, 1, true, true, 0,100,false,
                                emptyconstraint, 1e-6,1e-6, 1e-2,1e-3)
@test obj.loglik > -32.934741316276956 # should be ≈ -28.596171545957727

# allow hybrid ladders
Random.seed!(110)
PhyloNetworks.optimizestructure!(obj, maxmoves, 1, true, false, 0,100,false,
                                emptyconstraint, 1e-6,1e-6, 1e-2,1e-3)
@test obj.loglik > -28.596171545957727
end # of optimizestructure with simple example

@testset "phyLiNCone with simple net, no constraints" begin
no3cycle = true
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");
seed = 102
for nohybridladder in [true, false]
    obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastasimple, :JC69, 1)
    PhyloNetworks.checknetwork_LiNC!(obj.net, 1, no3cycle, nohybridladder)
    PhyloNetworks.startingBL!(obj.net, true, obj.trait, obj.siteweight)
    obj.loglik = -Inf # missing otherwise, which would cause an error below
    @test_nowarn PhyloNetworks.phyLiNCone!(obj, 1, no3cycle,
                                           nohybridladder, 3, 2, false, seed,
                                           0.5, emptyconstraint, 1e-2, 1e-2,
                                           1e-2, 1e-2, 0.0, 25.0)
    @test obj.loglik > -21.209048958984734
end
end

@testset "phyLiNC multiple runs" begin
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");
@test_nowarn obj = PhyloNetworks.phyLiNC!(net, fastasimple, :JC69; maxhybrid=2,
                    no3cycle=true, nohybridladder=true, maxmoves=2,
                    nreject=1, nruns=1, filename="", verbose=false, seed=105)
@test obj.loglik > -21.209048958984734
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");
@test_nowarn obj = PhyloNetworks.phyLiNC!(net, fastasimple, :JC69; maxhybrid=2,
                    no3cycle=true, nohybridladder=true, maxmoves=2,
                    nreject=1, nruns=1, filename="phyLiNC2", verbose=false, seed=0)
@test obj.loglik > -21.209048958984734
@test read("phyLiNC2.err", String) == ""
@test startswith(read("phyLiNC2.log", String), "---------------------\nphyLiNC network estimation starting")
rm("phyLiNC2.log")
rm("phyLiNC2.err")

net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");
addprocs(1) # multiple cores
@everywhere using PhyloNetworks
#using Distributed; @everywhere begin; using Pkg; Pkg.activate("."); using PhyloNetworks; end
obj = PhyloNetworks.phyLiNC!(net, fastasimple, :JC69; maxhybrid=2, no3cycle=true,
                        nohybridladder=true, maxmoves=2, nreject=1, nruns=2,
                        filename="phyLiNCmult", verbose=false, seed=106)
@test obj.loglik > -21.209048958984734
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

obj = PhyloNetworks.StatisticalSubstitutionModel(net_level1_i,fastaindiv,:JC69,2)
# obj.net = deepcopy of input net, so we need to rebuild the constraints
c_species[1] = PhyloNetworks.TopologyConstraint(0x01, c_species[1].taxonnames, obj.net)
# obj.net = deepcopy of input net, so we need to rebuild the constraints if done after
@test_logs (:warn, r"no 3-cycle") PhyloNetworks.checknetwork_LiNC!(obj.net, 2,
                                                    true, true, c_species)
for e in obj.net.edge e.length = 0.1; end # was -1.0 for missing
PhyloNetworks.startingBL!(obj.net, true, obj.trait, obj.siteweight) # true: to unzip
obj.loglik = -Inf # actual likelihood -56.3068141288164. Need something non-missing
seed = 103
@test_nowarn PhyloNetworks.phyLiNCone!(obj, 2, true, true,
                                       3, 2, false, seed, 0.5, c_species,
                                       1e-2, 1e-2, 1e-2, 1e-2, 0.0, 50.0)

# # test phyLiNC only, not phyLiNCone
obj = (@test_logs (:warn, r"no 3-cycle") match_mode=:any phyLiNC!(net_level1_s,
            fastaindiv, :JC69; maxhybrid=2, no3cycle=true, nohybridladder=true,
            verbose=false, filename="", speciesfile=mappingfile, seed=106, nruns=1,
            maxmoves=10, nreject=2))
@test obj.loglik > -56.30702264264643
# test that species stayed together after optimization, as the only polytomy
function polytomyS1(node)
    length(node.edge) > 3 || return false
    return Set(n.name for n in PhyloNetworks.getChildren(node)) == Set(["S1A", "S1B", "S1C"])
end
@test sum(polytomyS1(nod) for nod in obj.net.node) == 1
end

end # of overall phyLiNC test set

@testset "phyLiNC with rate variation" begin
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");
@test_nowarn obj = PhyloNetworks.phyLiNC!(net, fastasimple, :JC69, 4; maxhybrid=2,
                    no3cycle=true, nohybridladder=true, maxmoves=2,
                    nreject=1, nruns=1, filename="phyLiNCratevariation", verbose=false, seed=105)
@test occursin("categories for Gamma discretization: 4", read("phyLiNCratevariation.log", String))
@test obj.loglik > -21.209048958984734
@test read("phyLiNCratevariation.err", String) == ""
rm("phyLiNCratevariation.log")
rm("phyLiNCratevariation.err")
end

@testset "phyLiNC" begin
@test !isempty("just to try")
end
