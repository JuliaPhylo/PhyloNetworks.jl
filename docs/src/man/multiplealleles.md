# Multiple alleles per species

## between-species 4-taxon sets

The default setting for SNaQ considers that each allele in a gene tree corresponds
to a taxon (a tip) in the network. If instead each allele/individual can be mapped confidently
to a species, and if only the species-level network needs to be estimated,
then the following functions can be used:

```julia
tm = CSV.read(mappingFile) # taxon map as a data frame
taxonmap = Dict(tm[i,:allele] => tm[i,:species] for i in 1:110) # taxon map as a dictionary
```

The mapping file can be a text (or `csv`) file with two columns (at least):
one column named `allele` and one named `species`,
mapping each allele name to a species name. Next, read in the gene trees
and calculate the quartet CFs at the species level:


```julia
genetrees = readMultiTopology(genetreeFile)
df_sp = writeTableCF(countquartetsintrees(genetrees, taxonmap)...)
```

Now `df_sp` is a data frame containing the quartet concordance factors
at the species level only, that is, considering sets made of 4 distinct species,
even if the gene trees may have multiple alleles from the same species.
For 4 distinct species `A,B,C,D`, all alleles from each species (`A` etc.)
will be used to calculate the quartet CF. If a given gene tree has
`n_a` alleles from `a`, `n_b` alleles from `b` etc., then
each set of 4 alleles is given a weight of `1/(n_a n_b b_c n_c)`
to calculated of the CF for `A,B,C,D` (such that the total weight from
this particular gene trees is 1).
It is save to save this data frame, then use it for `snaq!` like this:

```julia
CSV.write("tableCF_species.csv", df)      # to save the data frame to a file
d_sp = readTableCF("tableCF_species.csv") # to get a "DataCF" object for use in snaq!.
```

## within-species 4-taxon sets

Four-taxon sets involving 2 individuals per species can provide more
information about the underlying network, including external branch
length in coalescent units. However, `snaq!` runs more slowly when
using this extra information. To get quartet CFs from sets of 4 individuals
in which 2 individuals are from the same species, the following functions
should be used:

```julia
df_ind = writeTableCF(countquartetsintrees(genetrees)...) # no mapping here: so quartet CFs across individuals
CSV.write("tableCF_individuals.csv", df)               # to save to a file
df_sp = mapAllelesCFtable(mappingFile, "tableCF_individuals.csv");
d_sp = readTableCF!(df_sp);
```
where the mapping file can be a text (or `csv`) file with two columns
named `allele` and `species`, mapping each allele name to a species name.
The data in `df_ind` is the table of concordance factors at the level of individuals.
In other words, it list CFs using one row for each set of 4 alleles/individuals.

`mapAllelesCFtable` creates a new data frame `df_sp` of quartet concordance factors at the
species level: with the allele names replaced by the appropriate species names.

**Warnings**:
- This procedure requires that all alleles from the same
  individual are given the same name (the individual's 'name') across
  all genes for which that individual was sequenced.
- For a four-taxon set `A,B,C,D`, all the individuals from `A`, `B`, `C` and `D`
  are considered, say `a1,b1,c1,d1`, `a2,b1,c1,d1`, `a1,b2,c1,d1`, `a2,b2,c1,d1`
  and so on. The CFs of these 4-taxon sets are averaged together to obtain the
  CFs at the species level. This procedures gives more weight to genes that have
  many alleles (because they contribute to more sets of 4 individuals) and less
  weight to genes that have few alleles.

The last command modifies this data frame `df_sp` by deleting rows that are uninformative
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
