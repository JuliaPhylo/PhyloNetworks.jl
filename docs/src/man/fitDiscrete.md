# Discrete Trait Evolution

With a phylogenetic network structure inferred, we can now estimate how quickly traits
have evolved over time using a likelihood model. These traits should be discrete 
characteristics of a species such as feather color, diet type, or dna in aligned 
genetic sequences.

## Discrete Trait Data

As with continuous trait evolution, we assume a fixed network, correctly rooted, 
with branch lengths proportional to calendar time. We start with a network, then
add data about the tips of this network. We allow data of two types.
1. A vector of species names with a data frame of traits:

```@setup fitDiscrete
using PhyloNetworks, DataFrames
mkpath("../assets/figures")
```

```@example fitDiscrete
#read in network
simple_net = readTopology("(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);");

#read in trait data
simple_species = ["C","A","B","D"]
simple_dat = DataFrame(trait=["hi","lo","lo","hi"])
```

If your species names and trait data are in the same data frame, read in your data 
frame then subset the data like this:
```@example fitDiscrete
dat = DataFrame(species=["C","A","B","D"], trait=["hi","lo","lo","hi"])
simple_species = dat[:species]
simple_dat = DataFrame(trait = dat[:trait])
```

2. To use dna data, read in the network structure then start with a fasta 
file. Reading the data from this file using the `readfastatodna` function. This 
creates a data frame of dna data and a vector of dna pattern weights.

```@example fitDiscrete
#read in network
dna_net = readTopology("((((((((((((((Ae_caudata_Tr275:1.0,Ae_caudata_Tr276:1.0):1.0,Ae_caudata_Tr139:1.0):1.0)#H1:1.0::0.6,((((((Ae_longissima_Tr241:1.0,Ae_longissima_Tr242:1.0):1.0,Ae_longissima_Tr355:1.0):1.0,(Ae_sharonensis_Tr265:1.0,Ae_sharonensis_Tr264:1.0):1.0):1.0,((Ae_bicornis_Tr408:1.0,Ae_bicornis_Tr407:1.0):1.0,Ae_bicornis_Tr406:1.0):1.0):1.0,((Ae_searsii_Tr164:1.0,Ae_searsii_Tr165:1.0):1.0,Ae_searsii_Tr161:1.0):1.0):1.0)#H2:1.0::0.6):1.0,(((Ae_umbellulata_Tr266:1.0,Ae_umbellulata_Tr257:1.0):1.0,Ae_umbellulata_Tr268:1.0):1.0,#H1:1.0::0.4):1.0):1.0,((Ae_comosa_Tr271:1.0,Ae_comosa_Tr272:1.0):1.0,(((Ae_uniaristata_Tr403:1.0,Ae_uniaristata_Tr357:1.0):1.0,Ae_uniaristata_Tr402:1.0):1.0,Ae_uniaristata_Tr404:1.0):1.0):1.0):1.0,(((Ae_tauschii_Tr352:1.0,Ae_tauschii_Tr351:1.0):1.0,(Ae_tauschii_Tr180:1.0,Ae_tauschii_Tr125:1.0):1.0):1.0,(#H2:1.0::0.4,((((Ae_mutica_Tr237:1.0,Ae_mutica_Tr329:1.0):1.0,Ae_mutica_Tr244:1.0):1.0,Ae_mutica_Tr332:1.0):1.0)#H4:1.0::0.6):1.0):1.0):1.0,(((T_boeoticum_TS8:1.0,(T_boeoticum_TS10:1.0,T_boeoticum_TS3:1.0):1.0):1.0,T_boeoticum_TS4:1.0):1.0,((T_urartu_Tr315:1.0,T_urartu_Tr232:1.0):1.0,(T_urartu_Tr317:1.0,T_urartu_Tr309:1.0):1.0):1.0):1.0):1.0,(((((Ae_speltoides_Tr320:1.0,Ae_speltoides_Tr323:1.0):1.0,Ae_speltoides_Tr223:1.0):1.0,Ae_speltoides_Tr251:1.0):1.0):1.0,#H4:1.0::0.4):1.0):1.0):1.0,Ta_caputMedusae_TB2:1.0):1.0,S_vavilovii_Tr279:1.0):1.0,Er_bonaepartis_TB1:1.0):1.0,H_vulgare_HVens23:1.0);");

#read in dna data
fastafile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples",
"Ae_bicornis_Tr406_Contig10132.aln")
dna_dat, dna_weights = readfastatodna(fastafile, true);
```

## Choosing a Substitution Model

After reading in your data, choose a model to describe how evolutionary changes 
(or substitutions, in the case of DNA) happened over time. We offer a selection 
of Markov substitution models to describe the evolutionary process.

### Generic Trait Models  

These models works well for any type of trait we may want to model. For general
trait types, use one of these three models:    
`:ERSM` Equal Rates Substitution Model  
`:BTSM` Binary Trait Substitution Model    
`:TBTSM` Two Binary Trait Substituion Model    

### DNA-Specific Models  

The DNA-specific models are optimized for aligned sequence data. We offer JC69 
and HKY85 models in both relative and absolute versions. The JC69 model was 
developed by Jukes and Cantor in 1969 and uses one rate for all type of substitutions. 
The HKY85 model was developed in 1985 by Hasegawa, Kishino, & Yano. It treats 
transitions differently from transversions.    
`:JC69` Jukes Cantor 69 Model    
`:HKY85` Hasegawa, Kishino and Yano 1985     

## Running FitDiscrete

To infer evolutionary rates, run the `fitDiscrete` function on the network and data. 
It will calculate the maximum likelihood score of a network given one or more 
discrete trait characters at the tips. Along each edge, evolutionary changes
are modeled with a continous time Markov model, with parameters estimated by 
maximizing the likelihood. At each hybrid node, the trait is assumed to be 
inherited from the immediate parent (or parents, in the case of a hybrid edge).
If there is a hybrid edge, the trait is modeled according to the parents' weighted 
average genetic contributions, as measured by inheritance gamma γ. The model 
currently ignores incomplete lineage sorting.

### General Trait Data

```@example fitDiscrete
s1 = fitDiscrete(simple_net, :ERSM, simple_species, simple_dat; optimizeQ=false, 
optimizeRVAS=false)
s2 = fitDiscrete(simple_net, :BTSM, simple_species, simple_dat; optimizeQ=false, 
optimizeRVAS=false)
```
In this `fitDiscrete` call, we do not optimize rates or allow for rate variation
across sites.

If `optimizeQ = true`, the `fitDiscrete` function estimates the evolutionary rate
or rates. Because we didn't allow for rate variation across sites in these models, 
we do not optimize the way rates may vary across trait types.

```@example fitDiscrete
s3 = fitDiscrete(simple_net, :ERSM, simple_species, simple_dat; optimizeQ=true, 
optimizeRVAS=true)
s4 = fitDiscrete(simple_net, :BTSM, simple_species, simple_dat; optimizeQ=true, 
optimizeRVAS=true)
```

### DNA Data

For DNA data, use `:JC69` or `:HKY85` models. 
```@example fitDiscrete
d1 = fitDiscrete(dna_net, :JC69, dna_dat, dna_weights, :RV; optimizeQ=false, 
optimizeRVAS=false)
d2 = fitDiscrete(dna_net, :HKY85, dna_dat, dna_weights, :RV; optimizeQ=false, 
optimizeRVAS=false)
```
In these `fitDiscrete` models, we do not optimize rates (`optimizeQ=false`), but
we do allow for rate variation across sites.

### Rate Variation Across Sites

In its default version, `fitDiscrete` does not allow for rate variation across sites.
To allow for rate variation across sites in your estimate of evolutionary rates 
(or rate variation across trait types, in the case of general trait types), 
include `:RV`. If you include `:RV` and `optimizeRVAS = true`, the model will 
not only allow for rate variation, but it will also optimize how rates vary across 
sites.

We optimize the evolutionary rates and the way rates vary across sites for the
DNA data here:    
```@example fitDiscrete
d3 = fitDiscrete(dna_net, :JC69, dna_dat, dna_weights, :RV; optimizeQ=true, 
optimizeRVAS=false)
d4 = fitDiscrete(dna_net, :HKY85, dna_dat, dna_weights, :RV; optimizeQ=true, 
optimizeRVAS=false)
```