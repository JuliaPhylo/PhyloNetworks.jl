```@setup concatdna
using PhyloNetworks
```

# fitting DNA on a network

The methods below model each DNA site as a trait, assuming that
sites are unlinked, that is, they evolve independently of each other.
In other words, this is a "concatenation" approach where
sites from the same locus do not share information about their
evolutionary path. This is appropriate if recombination is assumed
to have occurred within genes.

## DNA evolution: data and models

### reading in an alignment

As for trait evolution, [`fitdiscrete`](@ref) can be used. It can be
given data in a variety of ways. For DNA, this is one way:

```@repl concatdna
# read in network
dna_net = readTopology("((((((((((((((Ae_caudata_Tr275:1.0,Ae_caudata_Tr276:1.0):1.0,Ae_caudata_Tr139:1.0):1.0)#H1:1.0::0.6,((((((Ae_longissima_Tr241:1.0,Ae_longissima_Tr242:1.0):1.0,Ae_longissima_Tr355:1.0):1.0,(Ae_sharonensis_Tr265:1.0,Ae_sharonensis_Tr264:1.0):1.0):1.0,((Ae_bicornis_Tr408:1.0,Ae_bicornis_Tr407:1.0):1.0,Ae_bicornis_Tr406:1.0):1.0):1.0,((Ae_searsii_Tr164:1.0,Ae_searsii_Tr165:1.0):1.0,Ae_searsii_Tr161:1.0):1.0):1.0)#H2:1.0::0.6):1.0,(((Ae_umbellulata_Tr266:1.0,Ae_umbellulata_Tr257:1.0):1.0,Ae_umbellulata_Tr268:1.0):1.0,#H1:1.0::0.4):1.0):1.0,((Ae_comosa_Tr271:1.0,Ae_comosa_Tr272:1.0):1.0,(((Ae_uniaristata_Tr403:1.0,Ae_uniaristata_Tr357:1.0):1.0,Ae_uniaristata_Tr402:1.0):1.0,Ae_uniaristata_Tr404:1.0):1.0):1.0):1.0,(((Ae_tauschii_Tr352:1.0,Ae_tauschii_Tr351:1.0):1.0,(Ae_tauschii_Tr180:1.0,Ae_tauschii_Tr125:1.0):1.0):1.0,(#H2:1.0::0.4,((((Ae_mutica_Tr237:1.0,Ae_mutica_Tr329:1.0):1.0,Ae_mutica_Tr244:1.0):1.0,Ae_mutica_Tr332:1.0):1.0)#H4:1.0::0.6):1.0):1.0):1.0,(((T_boeoticum_TS8:1.0,(T_boeoticum_TS10:1.0,T_boeoticum_TS3:1.0):1.0):1.0,T_boeoticum_TS4:1.0):1.0,((T_urartu_Tr315:1.0,T_urartu_Tr232:1.0):1.0,(T_urartu_Tr317:1.0,T_urartu_Tr309:1.0):1.0):1.0):1.0):1.0,(((((Ae_speltoides_Tr320:1.0,Ae_speltoides_Tr323:1.0):1.0,Ae_speltoides_Tr223:1.0):1.0,Ae_speltoides_Tr251:1.0):1.0):1.0,#H4:1.0::0.4):1.0):1.0):1.0,Ta_caputMedusae_TB2:1.0):1.0,S_vavilovii_Tr279:1.0):1.0,Er_bonaepartis_TB1:1.0):1.0,H_vulgare_HVens23:1.0);");
# read in alignment in FASTA format
fastafile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","Ae_bicornis_Tr406_Contig10132.aln");
dna_dat, dna_weights = readfastatodna(fastafile, true);
dna_dat
dna_weights
```

Here, `dna_dat` is a single data frame containing both species names
and trait data (site patterns). The alignment was summarized by listing
each observed site pattern only once in `dna_dat`.
`dna_weights` is a vector of weights, containing
the number of times that each site pattern was observed.

### sequence substitution models

DNA-specific substitution models have 4 states: the 4 nucleotides from
[BioSymbols](https://github.com/BioJulia/BioSymbols.jl)
(listed [here](http://biojulia.net/BioSymbols.jl/stable/nucleicacids/)).
Each model has a relative and an absolute version.
- `:JC69` Jukes & Cantor 1969 model: one single rate for all transitions.
  The relative version has values -1 along the diagonal of the rate matrix
  (1 expected transition / unit of time). The absolute version has an extra
  parameter to scale the rate matrix.
- `:HKY85` Hasegawa, Kishino & Yano 1985: treats transitions differently
  from transversions. The relative is scaled to predict an average of
  1 transition / unit of time.

We may allow for rate variation across sites using the `:RV` option.

### likelihood of a fixed network

In the examples below, none of the rate parameters are optimized,
so we get to see the default starting values.

```@repl concatdna
d1 = fitdiscrete(dna_net, :JC69, dna_dat, dna_weights, :RV; optimizeQ=false, optimizeRVAS=false)
d2 = fitdiscrete(dna_net, :HKY85, dna_dat, dna_weights, :RV; optimizeQ=false, optimizeRVAS=false)
```
When allowing for rate variation across sites, the default α is 1.

In the more interesting examples below,
we optimize the evolutionary rates and the way rates vary across sites
(which is the default).
```@repl concatdna
d3 = fitdiscrete(dna_net, :JC69, dna_dat, dna_weights, :RV; ftolAbs=0.1, xtolAbs=0.01)
```
Lenient tolerance parameters `ftolAbs` etc. have been chosen here to
make this example faster.
Note that the fitted object contains a separate version of the input network,
where any taxon without data has been pruned, and where branch length numbers
may have been modified.

```@repl concatdna
d3.net
```
