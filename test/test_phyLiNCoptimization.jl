@testset "optimizelocalBL! and optimizelocalgammas! with simple example" begin
net_simple = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple.aln"))
obj = PhyloNetworks.StatisticalSubstitutionModel(net_simple, fastafile, :JC69)

## Local BL: unzip = true
lengthe = obj.net.edge[4].length
lengthep = obj.net.edge[4].node[1].edge[1].length
@test typeof(PhyloNetworks.optimizelocalBL!(obj, obj.net, obj.net.edge[4], true)) == Vector{PhyloNetworks.Edge}
@test obj.net.edge[4].length != lengthe
@test obj.net.edge[4].node[1].edge[1].length != lengthep

## Local BL: unzip = false
lengthe = obj.net.edge[9].length
lengthep = obj.net.edge[9].node[1].edge[1].length
@test typeof(PhyloNetworks.optimizelocalBL!(obj, obj.net, obj.net.edge[9], false)) == Vector{PhyloNetworks.Edge}
@test obj.net.edge[9].length != lengthe
@test obj.net.edge[9].node[1].edge[1].length != lengthep

# ## Local Gamma: unzip = true
hybridmajorparent = PhyloNetworks.getMajorParentEdge(obj.net.hybrid[1])
@test typeof(PhyloNetworks.optimizelocalgammas!(obj, obj.net, hybridmajorparent,
    true)) == Vector{PhyloNetworks.Edge}
@test hybridmajorparent.gamma != 0.9
@test PhyloNetworks.getMinorParentEdge(obj.net.hybrid[1]).gamma != 0.1

# ## Local Gamma: unzip = false
obj = PhyloNetworks.StatisticalSubstitutionModel(net_simple, fastafile, :JC69)
hybridmajorparent = PhyloNetworks.getMajorParentEdge(obj.net.hybrid[1])
@test typeof(PhyloNetworks.optimizelocalgammas!(obj, obj.net, hybridmajorparent,
    false)) == Vector{PhyloNetworks.Edge}
@test hybridmajorparent.gamma != 0.9
@test PhyloNetworks.getMinorParentEdge(obj.net.hybrid[1]).gamma != 0.1
end

@testset "optimizelocalBL! optimizelocalgammas! with complex network and 8 sites" begin
fastafile = joinpath(@__DIR__, "..", "examples", "Ae_bicornis_8sites.aln") # 8 sites only
fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "Ae_bicornis_8sites.aln"))
dna_dat, dna_weights = readfastatodna(fastafile, true); # 22 species, 3 hybrid nodes, 103 edges
net = readTopology("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),(((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);");
PhyloNetworks.fuseedgesat!(93, net)
for edge in net.edge # reset network
    setLength!(edge,1.0)
end
for h in net.hybrid
    setGamma!(PhyloNetworks.getMajorParentEdge(h),0.7)
end
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69);
@test length(obj.net.leaf) == 22

## Local BL: unzip = true
lengthe = obj.net.edge[48].length
lengthep = obj.net.edge[48].node[1].edge[1].length
@test typeof(PhyloNetworks.optimizelocalBL!(obj, obj.net, obj.net.edge[48], true)) == Vector{PhyloNetworks.Edge}
@test obj.net.edge[48].length != lengthe
@test obj.net.edge[48].node[1].edge[1].length == 0.0 # below hybrid node

# ## Local Gamma: unzip = truE
@test typeof(PhyloNetworks.optimizelocalgammas!(obj, obj.net,
    PhyloNetworks.getMajorParentEdge(obj.net.hybrid[1]), true)) == Vector{PhyloNetworks.Edge}
@test PhyloNetworks.getMajorParentEdge(obj.net.hybrid[1]).gamma != 0.7
@test PhyloNetworks.getMinorParentEdge(obj.net.hybrid[1]).gamma != 0.3
end #of local branch length and gamma optimization with localgamma! localBL! with 8 sites

@testset "global branch length and gamma optimization with 8 sites" begin
fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "Ae_bicornis_Tr406_Contig10132.aln"))
dna_dat, dna_weights = readfastatodna(fastafile, true);
net = readTopology("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),(((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);");
PhyloNetworks.fuseedgesat!(93, net)
for edge in net.edge #adds branch lengths
    setLength!(edge,1.0)
end
for h in net.hybrid
    setGamma!(PhyloNetworks.getMajorParentEdge(h),0.6)
end
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69);

## optimizeBL: unzip = true
@test typeof(PhyloNetworks.optimizeBL!(obj, obj.net, obj.net.edge, true)) == Vector{PhyloNetworks.Edge}
@test obj.net.edge[10] != 1.0
@test obj.net.edge[40] != 1.0

## optimizegammas: unzip = true
@test typeof(PhyloNetworks.optimizeallgammas!(obj, obj.net, true)) == Vector{PhyloNetworks.Edge}
@test PhyloNetworks.getMajorParentEdge(obj.net.hybrid[1]).gamma != 0.6
@test PhyloNetworks.getMinorParentEdge(obj.net.hybrid[1]).gamma != 0.4
end

@testset "data to SSM pruning: simple example" begin
net_simple = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple_missingone.aln"))
obj = PhyloNetworks.StatisticalSubstitutionModel(net_simple, fastafile, :JC69)
@test length(obj.net.edge) == 7
@test length(obj.net.hybrid) == 1
@test length(obj.net.leaf) == 3
@test !PhyloNetworks.hashybridladder(obj.net)
end

@testset "data to SSM pruning: complex network" begin
fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "Ae_bicornis_8sites.aln"))
net = readTopology("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),(((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);");
PhyloNetworks.fuseedgesat!(93, net)
for edge in net.edge # reset network
    setLength!(edge,1.0)
end
for h in net.hybrid
    setGamma!(PhyloNetworks.getMajorParentEdge(h),0.6)
end
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69);
@test length(obj.net.leaf) == 22
@test length(obj.net.edge) == 52
@test length(obj.net.hybrid) == 3
@test !PhyloNetworks.hashybridladder(obj.net)
end

@testset "checknetworkbeforeLiNC" begin
tree = readTopology("(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);");
@test !all([!(length(n.edge) == 2) for n in tree.node]) # one node of degree 2
PhyloNetworks.checknetworkbeforeLiNC!(tree, 1, true, true, true)
@test all([!(length(n.edge) == 2) for n in tree.node]) # no nodes of degree 2

net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
@test !all([!(length(n.edge) == 2) for n in net.node]) # one node of degree 2
PhyloNetworks.checknetworkbeforeLiNC!(net, 1, true, true, true)
@test all([!(length(n.edge) == 2) for n in net.node]) # no nodes of degree 2
@test all([(PhyloNetworks.getChildEdge(h).length == 0.0) for h in net.hybrid])# edges below hybrid node are of length zero
@test_throws ErrorException PhyloNetworks.checknetworkbeforeLiNC!(net, 0, true, true, true)
end

@testset "optimizestructure with simple example" begin
maxmoves = 2
maxhybrid = 1
seed = 123
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "simple.aln"))
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69, maxhybrid)
PhyloNetworks.checknetworkbeforeLiNC!(obj.net, maxhybrid, true, true, true)
PhyloNetworks.discrete_corelikelihood!(obj)
@test typeof(PhyloNetworks.optimizestructure!(obj, maxmoves, maxhybrid, true, true, true)) == Bool
@test writeTopology(obj.net) != "(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);"

# unzip = false
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69, maxhybrid)
PhyloNetworks.checknetworkbeforeLiNC!(obj.net, maxhybrid, true, false, true)
PhyloNetworks.discrete_corelikelihood!(obj)
@test typeof(PhyloNetworks.optimizestructure!(obj, maxmoves, maxhybrid, true, false, true)) == Bool
@test writeTopology(obj.net) != "(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);"

# allow hybrid ladders
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
obj = PhyloNetworks.StatisticalSubstitutionModel(net, fastafile, :JC69, maxhybrid)
PhyloNetworks.checknetworkbeforeLiNC!(obj.net, maxhybrid, true, true, false)
PhyloNetworks.discrete_corelikelihood!(obj)
@test typeof(PhyloNetworks.optimizestructure!(obj, maxmoves, maxhybrid, true, true, false)) == Bool
@test writeTopology(obj.net) != "(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);"
end # of optimizestructure with simple example

@testset "phyLiNC with simple net, no constraints" begin
no3cycle = true
fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..",
            "examples", "simple.aln"));
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");
seed = 123
for unzip in [true, false]
    for nohybridladder in [true, false]
        maxhybrid = 1;
        @test typeof(PhyloNetworks.phyLiNC!(net, fastafile, :JC69, maxhybrid,
                                        no3cycle, unzip, nohybridladder, 5, 2, # maxmoves = 5, nreject = 2
                                        false, seed)) == PhyloNetworks.StatisticalSubstitutionModel

        maxhybrid = 0;
        @test_throws ErrorException PhyloNetworks.phyLiNC!(net, fastafile,
                                            :JC69, maxhybrid, no3cycle, unzip,
                                            nohybridladder, 5, 2, false, seed);
                                            # maxmoves = 5, nreject = 2
    end
end
end

@testset "multiphyLiNC" begin
fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..",
            "examples", "simple.aln"));
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");
@test typeof(PhyloNetworks.multiphyLiNC!(net, fastafile, :JC69; maxhybrid=2,
                    no3cycle=true, unzip=true, nohybridladder=true, maxmoves=5,
                    nreject=2, nruns=1, filename="phyLiNC", verbose=false,
                    seed=123)) == PhyloNetworks.StatisticalSubstitutionModel
end

# @testset "multiphyLiNC with simple net and one constraint" begin
# no3cycle = true
# unzip = true
# nohybridladder = true
# seed = 123
# fastafile = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..",
#             "examples", "individuals.aln"));
# maxhybrid = 1;
# # constraint
# str_level1_s = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));" # indviduals S1A S1B S1C go on leaf 1
# net_level1_s = readTopology(str_level1_s)
# filename = abspath(joinpath(dirname(Base.find_package("PhyloNetworks")), "..",
#             "examples", "mappingIndividuals.csv"));
# # for travis? filename = joinpath(@__DIR__, "..","examples", "mappingIndividuals.csv")
# net_level1_i, c_species = PhyloNetworks.mapindividuals(net_level1_s, filename)
# PhyloNetworks.deletehybridedge!(net_level1_i, net_level1_i.edge[18])
# # ^this hybrid creates a three cycle after root node of degree two is removed
#? how do we want to deal with this type of three cycle?

# @test typeof(phyLiNC!(net_level1_i, fastafile, :JC69, 2, no3cycle, unzip,
#                       nohybridladder, 5, 2, true, seed, 0.5, c_species)
#                       ) == PhyloNetworks.StatisticalSubstitutionModel
# #TODO ERROR: DomainError with -0.6984169736707944: log will only return a complex result if called with a complex argument. Try log(Complex(x)).
# end
