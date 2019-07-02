```@setup traitevol_fixednet
using PhyloNetworks
using DataFrames
mkpath("../assets/figures")
```

# Discrete Trait Evolution

With a phylogenetic network structure inferred, we can now estimate how quickly traits
have evolved over time using a likelihood model. These traits should be discrete
characteristics of a species such as feather color, diet type,
or DNA in aligned genetic sequences.

## Discrete Trait Data

As with continuous trait evolution, we assume a fixed network, correctly rooted,
with branch lengths proportional to calendar time. We start with a network, then
add data about the tips of this network. We allow data of two types.

1. A vector of species names with a data frame of traits:

   ```@example traitevol_fixednet
   # read in network
   net = readTopology("(A:3,((B:0.4)#H1:1.6::0.92,((C:0.4,#H1:0::0.08):0.6,D:1):1):1);");
   # read in trait data
   species = ["C","A","B","D"]
   dat = DataFrame(trait=["hi","lo","lo","hi"])
   ```

   If your species names and trait data are in the same data frame,
   read in your data frame then subset the data like this:
   ```@example traitevol_fixednet
   dat = DataFrame(species=["C","A","B","D"], trait=["hi","lo","lo","hi"])
   species = dat[:species]
   dat = DataFrame(trait = dat[:trait])
   ```

2. To use dna data, read in the network structure then start with a fasta
   file. Reading the data from this file using the `readfastatodna` function.
   This creates a data frame of dna data and a vector of dna pattern weights.

   ```@example traitevol_fixednet
   # read in network
   dna_net = readTopology("((((((((((((((Ae_caudata_Tr275:1.0,Ae_caudata_Tr276:1.0):1.0,Ae_caudata_Tr139:1.0):1.0)#H1:1.0::0.6,((((((Ae_longissima_Tr241:1.0,Ae_longissima_Tr242:1.0):1.0,Ae_longissima_Tr355:1.0):1.0,(Ae_sharonensis_Tr265:1.0,Ae_sharonensis_Tr264:1.0):1.0):1.0,((Ae_bicornis_Tr408:1.0,Ae_bicornis_Tr407:1.0):1.0,Ae_bicornis_Tr406:1.0):1.0):1.0,((Ae_searsii_Tr164:1.0,Ae_searsii_Tr165:1.0):1.0,Ae_searsii_Tr161:1.0):1.0):1.0)#H2:1.0::0.6):1.0,(((Ae_umbellulata_Tr266:1.0,Ae_umbellulata_Tr257:1.0):1.0,Ae_umbellulata_Tr268:1.0):1.0,#H1:1.0::0.4):1.0):1.0,((Ae_comosa_Tr271:1.0,Ae_comosa_Tr272:1.0):1.0,(((Ae_uniaristata_Tr403:1.0,Ae_uniaristata_Tr357:1.0):1.0,Ae_uniaristata_Tr402:1.0):1.0,Ae_uniaristata_Tr404:1.0):1.0):1.0):1.0,(((Ae_tauschii_Tr352:1.0,Ae_tauschii_Tr351:1.0):1.0,(Ae_tauschii_Tr180:1.0,Ae_tauschii_Tr125:1.0):1.0):1.0,(#H2:1.0::0.4,((((Ae_mutica_Tr237:1.0,Ae_mutica_Tr329:1.0):1.0,Ae_mutica_Tr244:1.0):1.0,Ae_mutica_Tr332:1.0):1.0)#H4:1.0::0.6):1.0):1.0):1.0,(((T_boeoticum_TS8:1.0,(T_boeoticum_TS10:1.0,T_boeoticum_TS3:1.0):1.0):1.0,T_boeoticum_TS4:1.0):1.0,((T_urartu_Tr315:1.0,T_urartu_Tr232:1.0):1.0,(T_urartu_Tr317:1.0,T_urartu_Tr309:1.0):1.0):1.0):1.0):1.0,(((((Ae_speltoides_Tr320:1.0,Ae_speltoides_Tr323:1.0):1.0,Ae_speltoides_Tr223:1.0):1.0,Ae_speltoides_Tr251:1.0):1.0):1.0,#H4:1.0::0.4):1.0):1.0):1.0,Ta_caputMedusae_TB2:1.0):1.0,S_vavilovii_Tr279:1.0):1.0,Er_bonaepartis_TB1:1.0):1.0,H_vulgare_HVens23:1.0);");
   # read in dna data
   fastafile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","Ae_bicornis_Tr406_Contig10132.aln")
   dna_dat, dna_weights = readfastatodna(fastafile, true);
   dna_dat
   dna_weights
   ```

## Choosing a Substitution Model

After reading in your data, choose a model to describe how evolutionary changes
(or substitutions, in the case of DNA) happened over time.
Available Markov substitution models are described below.

### Generic Trait Models

These models works well for any type of trait we may want to model. For general
trait types, use one of these three models:
- `:BTSM` Binary Trait Substitution Model (2 states, rates unconstrained)
- `:ERSM` Equal Rates Substitution Model
  (`k` states, all transitions possible with equal rates)
- `:TBTSM` Two Binary Trait Substituion Model (though not fully implemented yet)

### DNA-Specific Models

The DNA-specific models are optimized for aligned sequence data.
The 4 nucleotide states are from
[BioSymbols](https://github.com/BioJulia/BioSymbols.jl)
(listed [here](http://biojulia.net/BioSymbols.jl/stable/nucleicacids/)).
Each model has a relative and an absolute version.
- `:JC69` Jukes & Cantor 1969 model: one single rate for all transitions.
  The relative version has values -1 along the diagonal of the rate matrix
  (1 expected transition / unit of time). The absolute version has an extra
  parameter to scale the rate matrix.
- `:HKY85` Hasegawa, Kishino & Yano 1985: treats transitions differently
  from transversions.

## Fitting the model

To infer evolutionary rates, run the `fitdiscrete` function on the network and data.
It will calculate the maximum likelihood score of a fixed network
given one or more discrete trait characters at the tips.
Along each edge, evolutionary changes
are modeled with a continous time Markov model, with parameters estimated by
maximizing the likelihood. At each hybrid node, the trait is assumed to be
inherited from the immediate parent (or parents, in the case of a hybrid node).
At a hybrid node, the trait is assumed to be inherited from one or the other
parent, with probabilities equal to the inheritance γ of each parent edge
(which is given by the network).
The model ignores incomplete lineage sorting (e.g. hemiplasy).

### General Trait Data

```@repl traitevol_fixednet
s1 = fitdiscrete(net, :ERSM, species, dat; optimizeQ=false)
s2 = fitdiscrete(net, :BTSM, species, dat; optimizeQ=false)
```
In this `fitdiscrete` call, we do not optimize rates or allow for rate variation
across sites. The default rates (which act as starting value if rates
were to be optimized) are chosen equal to the inverse of the total edge lengths
in the network (or 1/ntax if all branch lengths are missing).

If `optimizeQ = true` (which is the default), the `fitdiscrete`
function estimates the parameters of the rate matrix.
Because we didn't allow for rate variation across sites in these models,
there is nothing to optimize in the way rates may vary across traits (sites).

```@repl traitevol_fixednet
s3 = fitdiscrete(net, :ERSM, species, dat)
s4 = fitdiscrete(net, :BTSM, species, dat)
```

### DNA Data

For DNA data, use one of `:JC69` or `:HKY85`.
To allow for rate variation across sites, use the `:RV` option.

```@example traitevol_fixednet
d1 = fitdiscrete(dna_net, :JC69, dna_dat, dna_weights, :RV; optimizeQ=false, optimizeRVAS=false)
d2 = fitdiscrete(dna_net, :HKY85, dna_dat, dna_weights, :RV; optimizeQ=false, optimizeRVAS=false)
```
In these `fitdiscrete` models, we do not optimize rates (`optimizeQ=false`), but
we do allow for rate variation across sites, with a default α of 1.

### Rate Variation Across Sites

In its default version, `fitdiscrete` does not allow for rate variation across sites.
To allow for rate variation across sites in your estimate of evolutionary rates
(or rate variation across traits in the case of general traits),
include `:RV`. If you include `:RV` and `optimizeRVAS = true`,
the model will allow for rate variation and
also optimize the parameter α of the distribution of rates across sites.

We optimize the evolutionary rates and the way rates vary across sites for the
DNA data here:
```@repl traitevol_fixednet
d3 = fitdiscrete(dna_net, :JC69, dna_dat, dna_weights, :RV; optimizeRVAS=false)
d4 = fitdiscrete(dna_net, :HKY85, dna_dat, dna_weights, :RV)
```
