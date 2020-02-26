#= # for local testing, need this:
using Test
using PhyloNetworks
using PhyloPlots
using CSV
=#

@testset "setlengths and setgammas" begin
str_level1_s = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));" # indviduals S1A S1B S1C go on leaf 1
net = readTopology(str_level1_s)

PhyloNetworks.setlengths!([net.edge[1]], [1.1])
@test net.edge[1].length == 1.1
PhyloNetworks.setlengths!([net.edge[3], net.edge[4]], [3.3, 4.4])
@test net.edge[3].length == 3.3
@test net.edge[4].length == 4.4

PhyloNetworks.setmultiplegammas!([net.edge[18]], [0.25])
@test net.edge[18].gamma == 0.25
@test net.edge[16].gamma == 0.75

@test PhyloNetworks.getlengths([net.edge[1]]) == [net.edge[1].length]
@test PhyloNetworks.getlengths([net.edge[1], net.edge[5]]) == [net.edge[1].length, net.edge[5].length]

end

@testset "hashybridladder" begin
tree = readTopology("(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);");
@test !PhyloNetworks.hashybridladder(tree)
PhyloNetworks.addhybridedge!(tree, tree.edge[5], tree.edge[1], true)
PhyloNetworks.addhybridedge!(tree, tree.edge[2], tree.edge[1], true)
@test PhyloNetworks.hashybridladder(tree)
end # of testing hashybridladder

@testset "contain3cycles" begin
tree = readTopology("(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);");
@test !PhyloNetworks.contain3cycles(tree)
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
@test !PhyloNetworks.contain3cycles(net)
PhyloNetworks.addhybridedge!(net, net.edge[7], net.edge[4], true)
@test PhyloNetworks.contain3cycles(net)
net = readTopology("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),(((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);");
@test !PhyloNetworks.contain3cycles(net)
end

@testset "multiplygammas" begin
@test multiplygammas(0.25, 0.3) == 0.075
@test multiplygammas(-1.0, 0.3) == -1.0
@test multiplygammas(0.3, -1.0) == -1.0
@test multiplygammas(-1.0, -1.0) == -1.0
end
