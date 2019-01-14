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
species level, which can be used as input for networks estimation with `snaq!`.
But before, it is safe to save the concordance factor of quartets of species,
which can be calculated by averaging the CFs of quartets of individuals
from the associated species:

```julia
df_sp = writeTableCF(d_sp) # data frame, quartet CFs averaged across individuals of same species
CSV.write("CFtable_species.csv", df_sp) # save to file
```

Some quartets have the same species repeated twice,
representing cases when 2 of the 4 individuals came from the same species.
These quartets, with repeated species, are informative about the population
size of extant populations, i.e. about the lengths of external branches in
coalescent units.

now we can run snaq:

```julia
net = snaq!(T_sp, d_sp);
```
where `T_sp` should be a starting topology with one tip per species,
labelled with the same species names as the names used in the mapping file.

If `snaq!` takes too long that way, we can try a less ambitious estimation
that does not estimate the external branch lengths, that is,
*without* using quartets that have 2 individuals from the same species.
To do so, we can use the quartet concordance factors at the species level,
but filter out the quartets with one (or more) species repeated:

```julia
df_sp = writeTableCF(d_sp) # some quartets have the same species twice
function hasrep(row) # see if a row (4-taxon set) has a species name ending with "__2": repeated species
  occursin(r"__2$", row[:tx1]) || occursin(r"__2$", row[:tx2]) ||
    occursin(r"__2$", row[:tx3]) || occursin(r"__2$", row[:tx4])
end
df_sp_reduced = filter(!hasrep, df_sp) # removes rows with repeated species
df_sp_reduced # should have fewer rows than df_sp
CSV.write("CFtable_species_norep.csv", df_sp_reduced) # to save to file
d_sp_reduced = readTableCF(df_sp_reduced) # DataCF object, for input to snaq!
```

and now we can run `snaq!` on the reduced set of quartets without repeats,
which should be faster:

```julia
net = snaq!(T_sp, d_sp_reduced);
```


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