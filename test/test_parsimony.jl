extrarun = false

@testset "Testing parsimony score & reconstruction" begin
global net
originalstdout = stdout

@testset "Fitch" begin
# on a tree:
net = readTopology("(A,(B,(C,D)));")
tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1)
redirect_stdout(open("/dev/null", "w")) # not portable to Windows
score, states = PhyloNetworks.parsimonyDiscreteFitch(net, tips)
redirect_stdout(originalstdout)
@test score==1
@test states==Dict(4=>Set([1]),-4=>Set([1]),-3=>Set([0]),
         2=>Set([0]),3=>Set([1]),-2=>Set([0]),1=>Set([0]))

# on a network:
net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
redirect_stdout(open("/dev/null", "w"))
score, states = PhyloNetworks.parsimonyDiscreteFitch(net, tips)
redirect_stdout(originalstdout)
@test score==1
@test states==Dict(4=>Set([1]),-4=>Set([0]),-3=>Set([1]),
         2=>Set([0]),-2=>Set([1]),5=>Set([1]),1=>Set([0]))

tips = Dict("A" => 0, "B" => 1, "C" => 0, "D" => 1)
redirect_stdout(open("/dev/null", "w"))
score, states = PhyloNetworks.parsimonyDiscreteFitch(net, tips)
redirect_stdout(originalstdout)
@test score==2
@test states==Dict(4=>Set([0]),-6=>Set([0]),-4=>Set([0]),
  -3=>Set([0]),2=>Set([1]),-2=>Set([0,1]),5=>Set([1]),1=>Set([0]))

# from a data frame and with missing data:
dat = DataFrame(taxon=["A","E","B","C","D"], trait=[missing,2,0,1,1])
redirect_stdout(open("/dev/null", "w"))
score, states = PhyloNetworks.parsimonyDiscreteFitch(net, dat)
redirect_stdout(originalstdout)
@test score==1
@test states==Dict(4=>Set([1]),-6=>Set([1]),-4=>Set([0]),
  -3=>Set([1]),2=>Set([0]),-2=>Set([1]),5=>Set([1]))
end

@testset "Testing Tarjan's biconnected components" begin

net = readTopology("(A,(B,(C,D)));");
a = biconnectedComponents(net);
@test [[e.number for e in b] for b in a] == [[2],[3],[4],[5],[6],[1]]
net = readTopology("(((A,(((C,(D)#H2),(E,#H2)))#H1),(B,#H1)),F);");
a = biconnectedComponents(net);
@test [[e.number for e in b] for b in a] == [[1],[2],[3],[6],
  [8, 7, 4, 5],[9],[12],[14, 13, 10, 11],[16],[15]]
net = readTopology("(((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2));");
a = biconnectedComponents(net);
@test [[e.number for e in b] for b in a] == [[1],
  [2],[5],[6],[12],[10, 14, 13, 7, 8, 9, 3, 4, 11]]
net = readTopology("((((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2)),(((F)#H3,G),(H,#H3)));");
a = biconnectedComponents(net);
@test [[e.number for e in b] for b in a] == [[1],[2],[5],[6],[12],
  [10, 14, 13, 7, 8, 9, 3, 4, 11],[16],[20],[18],
  [22, 21, 17, 19],[23],[15]]
a = biconnectedComponents(net, true);
@test [[e.number for e in b] for b in a] == [[10, 14, 13, 7, 8, 9, 3, 4, 11],
  [22, 21, 17, 19]]
# net = readTopology("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),((((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H3),#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))#H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);")
# above: parsing error because one hybrid has a hybrid child

r,major,minor = PhyloNetworks.blobInfo(net, false);
@test [n.number for n in r] == [-5, 3, -8, 6, -10, -3, 9, -14, -12, -11, -2, -2]
r,major,minor = PhyloNetworks.blobInfo(net);
@test [n.number for n in r] == [-3,-11,-2]
@test [[e.number for e in h] for h in major] == [[7, 3],[17],[]]
@test [[e.number for e in h] for h in minor] == [[13,9],[21],[]]
forest, blobs = blobDecomposition(net);
@test length(blobs)==3
@test writeTopology(forest) == "(dummy -3,dummy -11);"
s = IOBuffer()
writeSubTree!(s, blobs[1], nothing, false, true)
@test String(take!(s)) == "(((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2));"
writeSubTree!(s, blobs[2], nothing, false, true)
@test String(take!(s)) == "(((F)#H3,G),(H,#H3));"
writeSubTree!(s, blobs[3], nothing, false, true)
@test String(take!(s)) == "(dummy -3,dummy -11);"
end

@testset "Testing level-based softwired parsimony" begin

net = readTopology("(A,(B,(C,D)));");
tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1);
@test parsimonySoftwired(net, tips) == 1.0
net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
@test parsimonySoftwired(net, tips) == 1.0
net = readTopology("((((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2)),(((F)#H3,G),(H,#H3)));");
tips = Dict("A"=>0, "B"=>0, "C"=>0, "D"=>0, "E"=>0, "F"=>0, "G"=>0, "H"=>0);
@test parsimonySoftwired(net, tips) == 0.0
tips = Dict("A"=>0, "B"=>0, "C"=>0, "D"=>0, "E"=>0, "F"=>1, "G"=>1, "H"=>1);
@test parsimonySoftwired(net, tips) == 1.0
tips = Dict("A"=>1, "B"=>0, "C"=>0, "D"=>0, "E"=>0, "F"=>1, "G"=>1, "H"=>1);
@test parsimonySoftwired(net, tips) == 2.0
tips = Dict("A"=>"no", "B"=>"no", "C"=>"no", "D"=>"maybe", "E"=>"maybe", "F"=>"yes", "G"=>"yes", "H"=>"yes");
@test parsimonySoftwired(net, tips) == 2.0
tips = Dict("A"=>"notatall", "B"=>"no", "C"=>"no", "D"=>"maybe", "E"=>"maybe", "F"=>"yes", "G"=>"notatall", "H"=>"yes");
@test parsimonySoftwired(net, tips) == 3.0

net = readTopology("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),(((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);");
fastafile = joinpath(@__DIR__, "..", "examples", "Ae_bicornis_8sites.aln")
species, sequences = PhyloNetworks.readFastaToArray(fastafile);
@test parsimonySoftwired(net, species, sequences) == 11.0

end # of test set for softwired parsimony

if extrarun
  fastafile = "../examples/Ae_bicornis_Tr406_Contig10132.aln"
  @test parsimonySoftwired(net, species, sequences) == 209.0
  fastafile = "../examples/Ae_bicornis_Tr406_Contig10722.aln"
  species, sequences = PhyloNetworks.readFastaToArray(fastafile);
  @test parsimonySoftwired(net, species, sequences) == 583.0
  @time parsimonySoftwired(net, species, sequences)
  # 1.299817 seconds (16.93 M allocations: 414.510 MiB, 3.47% gc time)
end


@testset "Testing general framework parsimony" begin

net = readTopology("(A,(B,(C,D)));");
tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1);
@test parsimonyGF(net, tips) == 1.0
net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);");
@test parsimonyGF(net, tips) == 1.0
net = readTopology("((((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2)),(((F)#H3,G),(H,#H3)));");
tips = Dict("A"=>0, "B"=>0, "C"=>0, "D"=>0, "E"=>0, "F"=>0, "G"=>0, "H"=>0);
@test parsimonyGF(net, tips) == 0.0
tips = Dict("A"=>0, "B"=>0, "C"=>0, "D"=>0, "E"=>0, "F"=>1, "G"=>1, "H"=>1);
@test parsimonyGF(net, tips) == 1.0
tips = Dict("A"=>1, "B"=>0, "C"=>0, "D"=>0, "E"=>0, "F"=>1, "G"=>1, "H"=>1);
@test parsimonyGF(net, tips) == 2.0
tips = Dict("A"=>"no", "B"=>"no", "C"=>"no", "D"=>"maybe", "E"=>"maybe", "F"=>"yes", "G"=>"yes", "H"=>"yes");
@test parsimonyGF(net, tips) == 2.0
tips = Dict("A"=>"notatall", "B"=>"no", "C"=>"no", "D"=>"maybe", "E"=>"maybe", "F"=>"yes", "G"=>"notatall", "H"=>"yes");
@test parsimonyGF(net, tips) == 3.0
net = readTopology("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),(((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);");
fastafile = joinpath(@__DIR__, "..", "examples", "Ae_bicornis_8sites.aln")
species, sequences = PhyloNetworks.readFastaToArray(fastafile);
@test parsimonyGF(net, species, sequences) == 11.0

end # of test set for GF parsimony

if extrarun
  fastafile = "../examples/Ae_bicornis_Tr406_Contig10132.aln"
  species, sequences = PhyloNetworks.readFastaToArray(fastafile);
  @test parsimonyGF(net, species, sequences) == 209.0
  fastafile = "../examples/Ae_bicornis_Tr406_Contig10722.aln"
  species, sequences = PhyloNetworks.readFastaToArray(fastafile);
  @test parsimonyGF(net, species, sequences) == 583.0
  @time parsimonyGF(net, species, sequences)
  """
  8.890167 seconds (107.25 M allocations: 3.702 GiB, 7.70% gc time)
  after changing a = min(a, b) to if ... end in 2 places:
  reduced time and # allocations, but not reduced total GiB
  5.298533 seconds (84.46 M allocations: 3.360 GiB, 9.04% gc time)
  small time reduction by storing the # of detached parents:
  4.908907 seconds (82.30 M allocations: 3.193 GiB, 9.78% gc time)
  """
end
end


@testset "data from CSV, parsimony search" begin
originalstdout = stdout
dat = CSV.read(joinpath(@__DIR__, "..","examples","Swadesh.csv"));
# net = readTopology("(((English,German),Norwegian),(Spanish,Portuguese));")
# species, sequences = PhyloNetworks.readCSVtoArray(dat)
# @test parsimonyGF(net, species, sequences) == 17.0
# arguments chosen for a very short run
redirect_stdout(open("/dev/null", "w"))
@test_logs maxParsimonyNet(readTopology("(((English,Spanish),Norwegian),(German,Portuguese));"),
    dat, hmax=1, runs=1, Nfail=2, outgroup="Spanish", rootname="", seed=6)
redirect_stdout(originalstdout)
# best = PhyloNetworks.maxParsimonyNetRun1(net, dat, 100, 0.1, 1, 1234,stdout,false,0.3,"Spanish",:softwired)
end
