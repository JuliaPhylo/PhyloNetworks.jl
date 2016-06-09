# Multiple alleles

The usual settings for SNaQ consider each allele in the gene trees to
be its own tip in the network. If instead each allele can be mapped confidently
to a species, and if only the species-level network needs
to be estimated, this can be done with the following functions:
```julia
new_df = mapAllelesCFtable(mappingFile, CFtable);
new_d = readTableCF(new_df);
```
where the mapping file can be a text file (or csv) with two columns
named *allele* and *species*, mapping each allele name to a species
name. The CF table is the original table with allele names for each
4-taxon subset. This function will create a new CF data frame with the
species names instead of allele names, and will modify *new_df* by
removing rows like *sp1,sp1,sp1,sp1*, which contain no information about
between-species relationships.

Estimation will work the same way:
```julia
new_net = snaq!(new_T,new_d);
```
where *new_T* should be a starting topology with one tip per species, labelled with the species names.
<!--
WARNING: the current function works best if all alleles from the same
individual are given the same name (the individual's 'name') across
all genes for which that individual was sequenced.
-->

