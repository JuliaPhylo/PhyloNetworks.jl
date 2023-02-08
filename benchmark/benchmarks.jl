using BenchmarkTools, PhyloNetworks, DataFrames, Logging

#suppresses @warn and @info for benchmarks
logger = SimpleLogger(stderr, Logging.Error);
old_logger = global_logger(logger);

# Define a parent BenchmarkGroup to contain our SUITE
const SUITE = BenchmarkGroup()

SUITE["nasm"] = BenchmarkGroup(["JC69", "HKY85"])
SUITE["fitDiscreteFixed"] = BenchmarkGroup(["ERSM", "BTSM", "JC69", "HKY85"])
SUITE["fitdiscrete"] = BenchmarkGroup(["ERSM", "BTSM", "JC69", "HKY85"])

# Add benchmarks to nasm group
SUITE["nasm"]["JC69"] = @benchmarkable JC69([0.5])
m1 = HKY85([.5], [0.25, 0.25, 0.25, 0.25]);
SUITE["nasm"]["HKY85"] = @benchmarkable P!(P(m1, 1.0), m1, 3.0)

# fitDiscreteFixed benchmarks
net_dat = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
species_alone = ["C","A","B","D"]
dat_alone = DataFrame(trait=["hi","lo","lo","hi"])
SUITE["fitDiscreteFixed"]["ERSM"] = @benchmarkable fitdiscrete(net_dat, :ERSM, species_alone, dat_alone; optimizeQ=false, optimizeRVAS=false)
SUITE["fitDiscreteFixed"]["BTSM"] = @benchmarkable fitdiscrete(net_dat, :BTSM, species_alone, dat_alone; optimizeQ=false, optimizeRVAS=false)

fastafile = joinpath(@__DIR__, "..", "examples", "Ae_bicornis_Tr406_Contig10132.aln")
dna_dat, dna_weights = readfastatodna(fastafile, true);
net_dna = readTopology("((((((((((((((Ae_caudata_Tr275,Ae_caudata_Tr276),Ae_caudata_Tr139))#H1,#H2),(((Ae_umbellulata_Tr266,Ae_umbellulata_Tr257),Ae_umbellulata_Tr268),#H1)),((Ae_comosa_Tr271,Ae_comosa_Tr272),(((Ae_uniaristata_Tr403,Ae_uniaristata_Tr357),Ae_uniaristata_Tr402),Ae_uniaristata_Tr404))),(((Ae_tauschii_Tr352,Ae_tauschii_Tr351),(Ae_tauschii_Tr180,Ae_tauschii_Tr125)),(((((((Ae_longissima_Tr241,Ae_longissima_Tr242),Ae_longissima_Tr355),(Ae_sharonensis_Tr265,Ae_sharonensis_Tr264)),((Ae_bicornis_Tr408,Ae_bicornis_Tr407),Ae_bicornis_Tr406)),((Ae_searsii_Tr164,Ae_searsii_Tr165),Ae_searsii_Tr161)))#H2,#H4))),(((T_boeoticum_TS8,(T_boeoticum_TS10,T_boeoticum_TS3)),T_boeoticum_TS4),((T_urartu_Tr315,T_urartu_Tr232),(T_urartu_Tr317,T_urartu_Tr309)))),(((((Ae_speltoides_Tr320,Ae_speltoides_Tr323),Ae_speltoides_Tr223),Ae_speltoides_Tr251))H3,((((Ae_mutica_Tr237,Ae_mutica_Tr329),Ae_mutica_Tr244),Ae_mutica_Tr332))#H4))),Ta_caputMedusae_TB2),S_vavilovii_Tr279),Er_bonaepartis_TB1),H_vulgare_HVens23);");
for edge in net_dna.edge #adds branch lengths
    setLength!(edge,1.0)
    if edge.gamma < 0
        setGamma!(edge, 0.5)
    end
end
SUITE["fitDiscreteFixed"]["JC69"] = @benchmarkable fitdiscrete(net_dna, :JC69, dna_dat, dna_weights; optimizeQ=false, optimizeRVAS=false)
SUITE["fitDiscreteFixed"]["HKY85"] = @benchmarkable fitdiscrete(net_dna, :HKY85, dna_dat, dna_weights; optimizeQ=false, optimizeRVAS=false)

## fitdiscrete benchmarks
SUITE["fitdiscrete"]["ERSM"] = @benchmarkable fitdiscrete(net_dat, :ERSM, species_alone, dat_alone; optimizeQ=true, optimizeRVAS=true)
SUITE["fitdiscrete"]["BTSM"] = @benchmarkable fitdiscrete(net_dat, :BTSM, species_alone, dat_alone; optimizeQ=true, optimizeRVAS=true)
SUITE["fitdiscrete"]["JC69"] = @benchmarkable fitdiscrete(net_dna, :JC69, dna_dat, dna_weights; optimizeQ=true, optimizeRVAS=true)
SUITE["fitdiscrete"]["HKY85"] = @benchmarkable fitdiscrete(net_dna, :HKY85, dna_dat, dna_weights; optimizeQ=true, optimizeRVAS=true)

# If a cache of tuned parameters already exists, use it, otherwise, tune and cache
# the benchmark parameters. Reusing cached parameters is faster and more reliable
# than re-tuning `SUITE` every time the file is included.
paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
    loadparams!(SUITE, BenchmarkTools.load(paramspath)[1], :evals);
else
    tune!(SUITE)
    BenchmarkTools.save(paramspath, params(SUITE));
end

global_logger(old_logger) #restores typical logging at end of benchmarks
