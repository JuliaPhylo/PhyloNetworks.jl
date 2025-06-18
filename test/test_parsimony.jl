extrarun = false

@testset "Testing parsimony score & reconstruction" begin
global net
originalstdout = stdout

@testset "utility: readfastatodna" begin
fasta8sites = joinpath(@__DIR__, "..", "examples", "Ae_bicornis_8sites.aln")
# fasta8sites = joinpath(dirname(pathof(PhyloNetworks)), "..", "examples", "Ae_bicornis_8sites.aln")
dna_dat, dna_weights = readfastatodna(fasta8sites, true) # 22 species
@test size(dna_dat) == (22, 9)
@test isa(dna_dat, DataFrame)
@test dna_weights == repeat([1], 8)
end

@testset "utility: readphylip" begin
  fastafile = "../examples/Ae_bicornis_truncated_names.aln"
  species, sequences = PhyloNetworks.readfastatoarray(fastafile);

  phylipfile = "../examples/Ae_bicornis_truncated_names.phylip"
  species2, sequences2 = PhyloNetworks.readphylip(phylipfile);

  dat1 = collect(zip(species , sequences ))
  dat2 = collect(zip(species2, sequences2))

  @test sort(dat1) == sort(dat2)
end

@testset "Fitch" begin
# on a tree:
net = readnewick("(A,(B,(C,D)));")
tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1)
redirect_stdout(devnull) # requires julia v1.6
score, states = PhyloNetworks.parsimonyfitch(net, tips)
redirect_stdout(originalstdout)
@test score==1
@test states==Dict(4=>Set([1]),-4=>Set([1]),-3=>Set([0]),
         2=>Set([0]),3=>Set([1]),-2=>Set([0]),1=>Set([0]))

# on a network:
net = readnewick("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
redirect_stdout(devnull)
score, states = PhyloNetworks.parsimonyfitch(net, tips)
redirect_stdout(originalstdout)
@test score==1
@test states==Dict(4=>Set([1]),-4=>Set([0]),-3=>Set([1]),
         2=>Set([0]),-2=>Set([1]),5=>Set([1]),1=>Set([0]))

tips = Dict("A" => 0, "B" => 1, "C" => 0, "D" => 1)
redirect_stdout(devnull)
score, states = PhyloNetworks.parsimonyfitch(net, tips)
redirect_stdout(originalstdout)
@test score==2
@test states==Dict(4=>Set([0]),-6=>Set([0]),-4=>Set([0]),
  -3=>Set([0]),2=>Set([1]),-2=>Set([0,1]),5=>Set([1]),1=>Set([0]))

# from a data frame and with missing data:
dat = DataFrame(taxon=["A","E","B","C","D"], trait=[missing,2,0,1,1])
redirect_stdout(devnull)
score, states = PhyloNetworks.parsimonyfitch(net, dat)
redirect_stdout(originalstdout)
@test score==1
@test states==Dict(4=>Set([1]),-6=>Set([1]),-4=>Set([0]),
  -3=>Set([1]),2=>Set([0]),-2=>Set([1]),5=>Set([1]))
end

@testset "Testing level-based softwired parsimony" begin

net = readnewick("(A,(B,(C,D)));");
tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1);
@test parsimonysoftwired(net, tips) == 1.0
net = readnewick("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
@test parsimonysoftwired(net, tips) == 1.0
net = readnewick("((((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2)),(((F)#H3,G),(H,#H3)));");
tips = Dict("A"=>0, "B"=>0, "C"=>0, "D"=>0, "E"=>0, "F"=>0, "G"=>0, "H"=>0);
@test parsimonysoftwired(net, tips) == 0.0
tips = Dict("A"=>0, "B"=>0, "C"=>0, "D"=>0, "E"=>0, "F"=>1, "G"=>1, "H"=>1);
@test parsimonysoftwired(net, tips) == 1.0
tips = Dict("A"=>1, "B"=>0, "C"=>0, "D"=>0, "E"=>0, "F"=>1, "G"=>1, "H"=>1);
@test parsimonysoftwired(net, tips) == 2.0
tips = Dict("A"=>"no", "B"=>"no", "C"=>"no", "D"=>"maybe", "E"=>"maybe", "F"=>"yes", "G"=>"yes", "H"=>"yes");
@test parsimonysoftwired(net, tips) == 2.0
tips = Dict("A"=>"notatall", "B"=>"no", "C"=>"no", "D"=>"maybe", "E"=>"maybe", "F"=>"yes", "G"=>"notatall", "H"=>"yes");
@test parsimonysoftwired(net, tips) == 3.0

net = readnewick("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),(((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);");
fastafile = joinpath(@__DIR__, "..", "examples", "Ae_bicornis_8sites.aln")
species, sequences = PhyloNetworks.readfastatoarray(fastafile);
@test parsimonysoftwired(net, species, sequences) == 11.0


end # of test set for softwired parsimony




if extrarun
  fastafile = "../examples/Ae_bicornis_Tr406_Contig10132.aln"
  species, sequences = PhyloNetworks.readfastatoarray(fastafile);
  @test parsimonysoftwired(net, species, sequences) == 209.0
  fastafile = "../examples/Ae_bicornis_Tr406_Contig10722.aln"
  species, sequences = PhyloNetworks.readfastatoarray(fastafile);
  @test parsimonysoftwired(net, species, sequences) == 583.0
  @time parsimonysoftwired(net, species, sequences)
  # 1.299817 seconds (16.93 M allocations: 414.510 MiB, 3.47% gc time)
end


@testset "Testing general framework parsimony" begin

net = readnewick("(A,(B,(C,D)));");
tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1);
@test parsimonyGF(net, tips) == 1.0
net = readnewick("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);");
@test parsimonyGF(net, tips) == 1.0
net = readnewick("((((A,(B)#H1),((C,(E)#H2),#H1)),(D,#H2)),(((F)#H3,G),(H,#H3)));");
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
net = readnewick("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),(((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);");
fastafile = joinpath(@__DIR__, "..", "examples", "Ae_bicornis_8sites.aln")
species, sequences = PhyloNetworks.readfastatoarray(fastafile);
@test parsimonyGF(net, species, sequences) == 11.0

end # of test set for GF parsimony

if extrarun
  fastafile = "../examples/Ae_bicornis_Tr406_Contig10132.aln"
  species, sequences = PhyloNetworks.readfastatoarray(fastafile);
  @test parsimonyGF(net, species, sequences) == 209.0
  fastafile = "../examples/Ae_bicornis_Tr406_Contig10722.aln"
  species, sequences = PhyloNetworks.readfastatoarray(fastafile);
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
dat = CSV.read(joinpath(@__DIR__,"..","examples","Swadesh.csv"), DataFrame)
# dat = CSV.read(joinpath(dirname(pathof(PhyloNetworks)),"..","examples","Swadesh.csv"), DataFrame)
# net = readnewick("(((English,German),Norwegian),(Spanish,Portuguese));")
# species, sequences = PhyloNetworks.readcsvtoarray(dat)
# @test parsimonyGF(net, species, sequences) == 17.0
# arguments chosen for a very short run
redirect_stdout(devnull)
@test_broken PhyloNetworks.maxParsimonyNet( # replace @test_broken by @test_logs when fixed
  readnewick("(((English,Spanish),Norwegian),(German,Portuguese));"),
  dat, hmax=1, runs=1, Nfail=2, outgroup="Spanish", rootname="", seed=6)
redirect_stdout(originalstdout)
# best = PhyloNetworks.maxParsimonyNetRun1(net, dat, 100, 0.1, 1, 1234,stdout,false,0.3,"Spanish",:softwired)
end
