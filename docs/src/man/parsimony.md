# Parsimony on networks

## Parsimony score of a given network

We can calculate the parsimony score of a given
network topology and a given set of characters.
The characters could be in a CSV file in this format:

| taxon | trait1 | trait2 | trait3 | trait4 | trait5 | trait6 ...
|:-------|:-------|:-------|:-------|:--------|:--------|:-------
| English      | 1| 2 | 1|   1 |       1 |       3
| ...    |  |   |  |         |              |       ...

The trait values can be integer numbers, or strings.
The data table may have missing data, and may contain extra taxa
that we might want to exclude.

An example file comes with the package, available
[here](https://github.com/crsl4/PhyloNetworks/blob/master/examples/Swadesh.csv)
or
[here](https://raw.githubusercontent.com/crsl4/PhyloNetworks/master/examples/Swadesh.csv).

```@setup parsimony
using PhyloNetworks
mkpath("../assets/figures")
using RCall
R"name <- function(x) file.path('..', 'assets', 'figures', x)"
# net1 = readTopology(joinpath(dirname(pathof(PhyloNetworks)), "..","examples","swadesh.out"))
# we would get net1 from analyzing the complete data, but not available with the package
```

First, we need to read the trait table as a DataFrame object:

```@repl parsimony
using CSV, DataFrames
csvfile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","Swadesh.csv");
dat = CSV.read(csvfile);
first(dat, 6) # to see the first 6 rows
```

Then, we need to convert the DataFrame object `dat`
into a vector of species and traits.
The species names are in column 1 named `taxon`,
and the traits are in columns 2-11. The trait data need to be
converted to a list of vectors, with one vector for each species.
An internal function is provided for this:

```@repl parsimony
species, traits = PhyloNetworks.readCSVtoArray(dat);
species
traits
```

Then, we read the network as usual:

```@repl parsimony
net = readTopology("(Spanish,((English)#H1,(Norwegian,(German,#H1))));");
```

```@example parsimony
using PhyloPlots, RCall
R"svg(name('parsimony-fixed-net.svg'), width=4, height=4)"; # hide
R"par"(mar = [0,0,0,0]);
plot(net, :R, xlim=[0.8,7.5]);
R"dev.off"(); # hide
```
![parsimony-fixed-net](../assets/figures/parsimony-fixed-net.svg)

There are different types of parsimony scores on networks.
Currently, we have implemented the **softwired** criterion only,
with two different functions:
[`parsimonySoftwired`](@ref) and
[`parsimonyGF`](@ref).

The function `parsimonySoftwired` uses a faster algorithm than
`parsimonyGF`, but can solve the softwired criterion only.

```@repl parsimony
score = parsimonySoftwired(net, species, traits)
score = parsimonyGF(net,species,traits,:softwired)
```


## Finding the most parsimonious network

The function [`maxParsimonyNet`](@ref) searches for the most parsimonious
level-1 network. It uses the `parsimonyGF` function, with softwired criterion
as default, which will be extended to other criteria later.

Just like [`snaq!`](@ref), [`maxParsimonyNet`](@ref) requires a
starting topology, which can be a tree or a level-1 network,
and returns a level-1 network.
Taxa present in the data but absent from the starting topology
will be ignored during the search.

```julia
starttree = readTopology("(((English,German),Norwegian),Spanish);");
net1 = maxParsimonyNet(starttree, dat, hmax=1, outgroup="Spanish", rootname="swadesh")
```

The example data is very small: only 1 of the 11 traits is parsimony informative,
on the 4 taxa specified by the starting topology. So these data happen
to be compatible with a tree, and that tree is returned despite allowing
for up to 1 reticulation:
`(Spanish,((English,Norwegian),German));`.
