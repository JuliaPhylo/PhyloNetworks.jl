# Multiple alleles per species

The default setting for SNaQ considers that each allele in a gene tree corresponds
to a taxon (a tip) in the network. If instead each allele/individual can be mapped confidently
to a species, and if only the species-level network needs to be estimated,
then the following functions should be used:

```julia
df_sp = mapAllelesCFtable(mappingFile, CFtable_ind);
d_sp = readTableCF!(df_sp);
```
where the mapping file can be a text (or `csv`) file with two columns
named `allele` and `species`, mapping each allele name to a species name.
The CF table `CFtable_ind` should be a table of concordance factors at the level of individuals.
In other words, it should list CFs using one row for each set of 4 alleles/individuals.
The first command creates a new data frame `df_sp` of quartet concordance factors at the
species level: with the allele names replaced by the appropriate species names.

The second command modifies this data frame `df_sp` by deleting rows that are uninformative
about between-species relationships, such as rows corresponding to 4 individuals from the
same species. The output `d_sp` of this second command is an object of type `DataCF` at the
species level, which can be used as input for networks estimation with `snaq!`:

```julia
net = snaq!(T_sp, d_sp);
```
where `T_sp` should be a starting topology with one tip per species,
labelled with the same species names as the names used in the mapping file.

Warnings:

- This feature has not been fully tested
- This procedure is slow and should be made faster, when working with gene trees as input.
- If input data are gene trees, the CF table at the individual level should be created first,
  like this:
  ```julia
  CFtable_ind = readTrees2CF(gene tree file);
  ```
  before applying the two commands above.
  At this time, however, this procedure requires that all alleles from the same
  individual are given the same name (the individual's 'name') across
  all genes for which that individual was sequenced.