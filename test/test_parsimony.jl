@testset "Testing parsimony score & reconstruction" begin

# on a tree:
net = readTopology("(A,(B,(C,D)));")
tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1)
score, states = parsimonyDiscrete(net, tips)
@test score==1
@test states==Dict(4=>Set([1]),-4=>Set([1]),-3=>Set([0]),
         2=>Set([0]),3=>Set([1]),-2=>Set([0]),1=>Set([0]))

# on a network:
net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
score, states = parsimonyDiscrete(net, tips)
@test score==1
@test states==Dict(4=>Set([1]),-4=>Set([0]),-3=>Set([1]),
         2=>Set([0]),-2=>Set([1]),5=>Set([1]),1=>Set([0]))

tips = Dict("A" => 0, "B" => 1, "C" => 0, "D" => 1)
score, states = parsimonyDiscrete(net, tips)
@test score==2
@test states==Dict(4=>Set([0]),-6=>Set([0]),-4=>Set([0]),
  -3=>Set([0]),2=>Set([1]),-2=>Set([0,1]),5=>Set([1]),1=>Set([0]))

# from a data frame and with missing data:
dat = DataFrame(taxon=["A","B","C","D"], trait=[0,0,1,1])
dat[1,2]=NA
score, states = parsimonyDiscrete(net, dat)
@test score==1
@test states==Dict(4=>Set([1]),-6=>Set([1]),-4=>Set([0]),
  -3=>Set([1]),2=>Set([0]),-2=>Set([1]),5=>Set([1]))

end

@testset "Testing Tarjan's biconnected components" begin

net = readTopology("(A,(B,(C,D)));")
a = biconnectedComponents(net);
@test [[e.number for e in b] for b in a] == [[2],[3],[4],[5],[6],[1]]
net = readTopology("(((A,(((C,(D)#H2),(E,#H2)))#H1),(B,#H1)),F);")
a = biconnectedComponents(net);
@test [[e.number for e in b] for b in a] == [[1],[2],[3],[6],
  [8, 7, 4, 5],[9],[12],[14, 13, 10, 11],[16],[15]]
net = readTopology("(((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2));")
a = biconnectedComponents(net);
@test [[e.number for e in b] for b in a] == [[1],
  [2],[5],[6],[12],[10, 14, 13, 7, 8, 9, 3, 4, 11]]
net = readTopology("((((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2)),(((F)#H3,G),(H,#H3)));")
a = biconnectedComponents(net);
@test [[e.number for e in b] for b in a] == [[1],[2],[5],[6],[12],
  [10, 14, 13, 7, 8, 9, 3, 4, 11],[16],[20],[18],
  [22, 21, 17, 19],[23],[15]]
a = biconnectedComponents(net, true);
@test [[e.number for e in b] for b in a] == [[10, 14, 13, 7, 8, 9, 3, 4, 11],
  [22, 21, 17, 19]]
# net = readTopology("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),((((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H3),#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))#H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);")
# above: parsing error because one hybrid has a hybrid child

a = blobRoots(net)
@test [n.number for n in a] == [-5, 3, -8, 6, -10, -3, 9, -14, -12, -11, -2, -2]
end
