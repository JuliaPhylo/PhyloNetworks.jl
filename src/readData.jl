# function to write a csv table from the expCF of an
# array of quartets
# warning: does not check if the expCF have been calculated
function writeExpCF(quartets::Array{Quartet,1})
    df = DataFrames.DataFrame(t1=String[],t2=String[],t3=String[],t4=String[],
                              CF12_34=Float64[],CF13_24=Float64[],CF14_23=Float64[])
    for q in quartets
        length(q.taxon) == 4 || error("quartet $(q.number) does not have 4 taxa")
        length(q.qnet.expCF) == 3 || error("quartet $(q.number) does have qnet with 3 expCF")
        push!(df, [q.taxon[1],q.taxon[2],q.taxon[3],q.taxon[4],q.qnet.expCF[1],q.qnet.expCF[2],q.qnet.expCF[3]])
    end
    return df
end

writeExpCF(d::DataCF) = writeExpCF(d.quartet)

"""
    writeTableCF(vector of Quartet objects)
    writeTableCF(DataCF)
    writeTableCF(vector of QuartetT objects [, taxonames])

Build a DataFrame containing observed quartet concordance factors,
with columns named:
- `:tx1`, `:tx2`, `:tx3`, `:tx4` for the four taxon names in each quartet
-  `:CF12_34`, `:CF13_24`, `:CF14_23` for the 3 quartets of a given four-taxon set
- `:ngenes` if this information is available for some quartets

If the input are [`QuartetT`](@ref) objects, their `data` field needs to
contain vectors of 4 values.
"""
function writeTableCF(quartets::Array{Quartet,1})
    df = DataFrames.DataFrame(t1=String[],t2=String[],t3=String[],t4=String[],
                              CF12_34=Float64[],CF13_24=Float64[],CF14_23=Float64[],
                              ngenes=Union{Missing, Float64}[])
    for q in quartets
        length(q.taxon) == 4 || error("quartet $(q.number) does not have 4 taxa")
        length(q.obsCF) == 3 || error("quartet $(q.number) does have qnet with 3 expCF")
        push!(df, [q.taxon[1],q.taxon[2],q.taxon[3],q.taxon[4],q.obsCF[1],q.obsCF[2],q.obsCF[3],
                   (q.ngenes==-1.0 ? missing : q.ngenes)])
    end
    if all(ismissing, df[!,:ngenes])
        select!(df, Not(:ngenes))
    end
    return df
end

writeTableCF(d::DataCF) = writeTableCF(d.quartet)

function writeTableCF(quartets::Vector{QuartetT{T}},
                      taxa=Vector{String}()::AbstractVector{String}) where T<:AbstractVector
    V = eltype(T) # would expect Float64, but Int would be reasonable if counts are not normalized
    V <: Real || error("CFs need to take real values")
    df = DataFrames.DataFrame(t1=String[],t2=String[],t3=String[],t4=String[],
                              CF12_34=V[],CF13_24=V[],CF14_23=V[], ngenes=V[])
    ntaxa = length(taxa)
    taxstring = x -> ( x>ntaxa ? string(x) : taxa[x] )
    for q in quartets
        length(q.data) == 4 || error("quartet $(q.data) does not have 4 data points: CF12, CF13, CF14, ngenes")
        qn = taxstring.(q.taxonnumber)
        push!(df, [qn[1],qn[2],qn[3],qn[4],q.data[1],q.data[2],q.data[3],q.data[4]])
    end
    return df
end

"""
    readTableCF(file)
    readTableCF(data frame)
    readTableCF!(data frame)

Read a file or DataFrame object containing a table of concordance factors (CF),
with one row per 4-taxon set. The first 4 columns are assumed to give the labels
of the 4 taxa in each set (tx1, tx2, tx3, tx4).
Columns containing the CFs are assumed to be named
`CF12_34`, `CF13_24` and `CF14_23`;
or `CF12.34`, `CF13.24` and `CF14.23`;
or else are assumed to be columns 5,6,7.
If present, a column named 'ngenes' will be used to get the number of loci
used to estimate the CFs for each 4-taxon set.

Output: [`DataCF`](@ref) object

Optional arguments:

- summaryfile: if specified, a summary file will be created with that name.
- delim (for the first form only): to specify how columns are delimited,
  with single quotes: delim=';'. Default is a `csv` file, i.e. `delim=','`.

The last version modifies the input data frame, if species are represented by multiple alleles
for instance (see [`readTableCF!`](@ref)(data frame, columns)).
"""
function readTableCF(file::AbstractString; delim=','::Char, summaryfile=""::AbstractString)
    df = CSV.read(file, delim=delim)
    readTableCF!(df, summaryfile=summaryfile)
end

function readTableCF(df0::DataFrames.DataFrame; summaryfile=""::AbstractString)
    df = deepcopy(df0)
    readTableCF!(df, summaryfile=summaryfile)
end

function readTableCF!(df::DataFrames.DataFrame; summaryfile=""::AbstractString)
    @debug "assume the numbers for the taxon read from the observed CF table match the numbers given to the taxon when creating the object network"
    alternativecolnames = [ # obsCF12 is as exported by fittedQuartetCF()
        [:CF12_34, Symbol("CF12.34"), :obsCF12],
        [:CF13_24, Symbol("CF13.24"), :obsCF13],
        [:CF14_23, Symbol("CF14.23"), :obsCF14]
    ]
    obsCFcol = [findfirst(x-> x ∈ alternativecolnames[1], DataFrames.names(df)),
                findfirst(x-> x ∈ alternativecolnames[2], DataFrames.names(df)),
                findfirst(x-> x ∈ alternativecolnames[3], DataFrames.names(df))]
    ngenecol =  findfirst(isequal(:ngenes), DataFrames.names(df))
    withngenes = ngenecol !== nothing
    if nothing in obsCFcol # one or more col names for CFs were not found
        size(df,2) == (withngenes ? 8 : 7) ||
          @warn """Column names for quartet concordance factors (CFs) were not recognized.
          Was expecting CF12_34, CF13_24 and CF14_23 for the columns with CF values,
          or CF12.34 or obsCF12, etc.
          Will assume that the first 4 columns give the taxon names, and that columns 5-7 give the CFs."""
        obsCFcol = [5,6,7] # assuming CFs are in columns 5,6,7, with colname mismatch
    end
    minimum(obsCFcol) > 4 ||
        error("CFs found in columns $obsCFcol, but taxon labels expected in columns 1-4")
    # fixit: what about columns giving the taxon names: always assumed to be columns 1-4? No warning if not?
    columns = [[1,2,3,4]; obsCFcol]
    if withngenes  push!(columns, ngenecol)  end

    d = readTableCF!(df, columns)

    if withngenes # && d.numTrees == -1
        m1 = minimum([q.ngenes for q in d.quartet])
        m2 = maximum([q.ngenes for q in d.quartet])
        if m1<m2 print("between $m1 and ") end
        println("$m2 gene trees per 4-taxon set")
        # other info printed by show() on a DataCF object: num quartets and num gene trees
    end
    if(summaryfile != "")
        descData(d,summaryfile)
    end
    return d
end

# see docstring below, for readTableCF!
# takes in df and 7 or 8 column numbers (4 labels + 3 CFs + ngenes possibly)
function readTableCF!(df::DataFrames.DataFrame, co::Vector{Int})
    withngenes = (length(co)==8) # true if column :ngenes exists, false ow
    repSpecies = cleanAlleleDF!(df,co) # removes uninformative rows from df (not df0)
    # fixit: cleanAlleleDF! is time consuming but many times not needed
    # add option to skip it, if the user knows that each tip appears once only?
    if !isempty(repSpecies)
        df = mergeRows(df,co)   # warning: this 'df' is *not* changed externally
    end                         # we cannot move to mapAllelesCFtable because we need repSpecies in here
    quartets = Quartet[]
    for i in 1:size(df,1)
        push!(quartets,Quartet(i,string(df[i,co[1]]),string(df[i,co[2]]),string(df[i,co[3]]),string(df[i,co[4]]),
                               [df[i,co[5]],df[i,co[6]],df[i,co[7]]]))
        if withngenes
            quartets[end].ngenes = df[i,co[8]]
        end
    end
    d = DataCF(quartets)
    if(!isempty(repSpecies))
        d.repSpecies = repSpecies
    end
    return d  # return d, df ## to save memory & gc with readTableCF! for bootstrapping?
end

"""
    readTableCF!(data frame, columns)

Read in quartet CFs from data frame, assuming information is in columns numbered `columns`,
of length **7 or 8**: 4 taxon labels then 3 CFs then ngenes possibly.

If some species appears more than once in the same 4-taxon set (e.g. t1,t1,t2,t3),
then the data frame is modified to remove rows (4-taxon sets) that are uninformative about
between-species relationships. This situation may occur if multiple individuals are sampled from
the same species. A 4-taxon set is uninformative (and its row is removed)
if one taxon is repeated 3 or 4 times (like t1,t1,t1,t1 or t1,t2,t2,t2).
The list of species appearing twice in some 4-taxon sets is stored in the output DataCF object.
For these species, the length of their external edge is identifiable (in coalescent units).
If multiple rows correspond to the same 4-taxon set, these rows are merged and their CF values
(and number of genes) are averaged.

    readTableCF!(DataCF, data frame, columns)

Modify the `.quartet.obsCF` values in the `DataCF` object with those read from the data frame
in columns numbered `columns`.
`columns` should have **3** columns numbers for the 3 CFs in this order:
`12_34`, `13_24` and `14_23`.

Assumptions:
- same 4-taxon sets in `DataCF` and in the data frame, and in the same order,
  but this assumption is *not checked* (for speed, e.g. during bootstrapping).
- one single row per 4-taxon set (multiple individuals representatives
  of the same 4-taxon set should have been already merged);
  basically: the DataCF should have been created from the data frame by `readTableCF!(df, colums)`
"""
function readTableCF!(datcf::DataCF, df::DataFrame, cols::Vector{Int})
    for i in 1:size(df,1)
        for j in 1:3
            datcf.quartet[i].obsCF[j] = df[i,cols[j]]
        end
    end
end


# ---------------- read input gene trees and calculate obsCF ----------------------

"""
    readInputTrees(file)

Read a text file with a list of trees/networks in parenthetical format
(one tree per line) and transform them like [`readTopologyLevel1`](@ref)
does: to be unrooted, with resolved polytomies, missing branch lengths
set to 1.0, etc. See [`readMultiTopology`](@ref) to read multiple
trees or networks with no modification.

Output: array of HybridNetwork objects.

Each line starting with "(" will be considered as describing one topology.
The file can have extra lines that are ignored.
"""
function readInputTrees(file::AbstractString)
    try
        s = open(file)
    catch
        error("Could not find or open $(file) file");
    end
    vnet = HybridNetwork[];
    s = open(file)
    numl = 1
    for line in eachline(s)
        line = strip(line) # remove spaces
        @debug "$(line)"
        c = isempty(line) ? "" : line[1]
        if(c == '(')
           try
               push!(vnet, readTopologyUpdate(line,false))
           catch(err)
               error("could not read tree in line $(numl). The error is $(err)")
           end
        end
        numl += 1
    end
    close(s)
    return vnet # consistent output type: HybridNetwork vector. might be of length 0.
end

# function to list all quartets for a set of taxa names
# return a vector of quartet objects, and if writeFile=true, writes a file
# warning: taxon has to be vector of String, vector of AbstractString do not work
function allQuartets(taxon::Union{Vector{String},Vector{Int}}, writeFile::Bool)
    quartets = combinations(taxon,4)
    vquartet = Quartet[];
    if(writeFile)
        #allName = "allQuartets$(string(integer(time()/1000))).txt"
        allName = "allQuartets.txt"
        f = open(allName,"w")
    end
    i = 1
    for q in quartets
        if(writeFile)
            write(f,"$(q[1]),$(q[2]),$(q[3]),$(q[4])\n")
        end
        push!(vquartet,Quartet(i,string(q[1]),string(q[2]),string(q[3]),chomp(string(q[4])),[1.0,0.0,0.0]))
        i += 1 # overflow error if # quartets > typemax(Int), i.e. if 121,978+ taxa with Int64, 478+ taxa with Int32
    end
    if(writeFile)
        close(f)
    end
    return vquartet
end

allQuartets(numTaxa::Integer, writeFile::Bool) = allQuartets(1:numTaxa, writeFile)

# function to list num randomly selected quartets for a vector of all quartets
# return a vector of Quartet randomly chosen (and a file if writeFile=true)
function randQuartets(allquartets::Vector{Quartet},num::Integer, writeFile::Bool)
    randquartets = Quartet[]
    n = length(allquartets)
    if num == 0
        num = Integer(floor(0.1*n))
    end
    num <= n || error("you cannot choose a sample of $(num) quartets when there are $(n) in total")
    indx = [ones(Bool,num); zeros(Bool,n-num)]
    indx = indx[sortperm(randn(n))]
    if writeFile
        randName = "rand$(num)Quartets.txt"
        println("list of randomly selected quartets in file $(randName)")
        out = open(randName,"w")
    end
    for i in 1:n
        if indx[i]
            if writeFile
                q = allquartets[i].taxon
                write(out,"$(q[1]),$(q[2]),$(q[3]),$(q[4])\n")
            end
            push!(randquartets,allquartets[i])
        end
    end
    if writeFile
        close(out)
    end
    return randquartets
end

# function that will not use randQuartets(list of quartets,...)
# this function uses whichQuartet to avoid making the list of all quartets
# fixit: i think we should write always the file, but not sure
function randQuartets(taxon::Union{Vector{String},Vector{Int}},num::Integer, writeFile::Bool)
    randquartets = Quartet[]
    n = length(taxon)
    ntotal = binom(n,4)
    num <= ntotal || error("you cannot choose a sample of $(num) quartets when there are $(ntotal) in total")
    # indx = [rep(1,num);rep(0,ntotal-num)] # requires much more memory than necessary:
    # indx = indx[sortperm(randn(ntotal))]  # several arrays of size ntotal !!
    # rq = findall(x -> x==1, indx)
    rq = sample(1:ntotal, num, replace=false, ordered=true)
    randName = "rand$(num)Quartets.txt"
    println("list of randomly selected quartets in file $(randName)")
    out = open(randName,"w")
    i = 1
    for q in rq
        qind = whichQuartet(n,q) # vector of int
        quartet = createQuartet(taxon,qind,i)
        write(out,"$(quartet.taxon[1]), $(quartet.taxon[2]), $(quartet.taxon[3]), $(quartet.taxon[4])\n")
        push!(randquartets,quartet)
        i += 1
    end
    close(out)
    return randquartets
end

randQuartets(numTaxa::Integer,num::Integer, writeFile::Bool) = randQuartets(1:numTaxa,num, writeFile)

# function to read list of quartets from a file
# and create Quartet type objects
function readListQuartets(file::AbstractString)
    try
        f = open(file)
    catch
        error("could not open file $(file)")
    end
    f = open(file)
    quartets = Quartet[];
    i = 1
    for line in eachline(f)
        l = split(line,",")
        length(l) == 4 || error("quartet with $(length(l)) elements, should be 4: $(line)")
        push!(quartets,Quartet(i,string(l[1]),string(l[2]),string(l[3]),chomp(string(l[4])),[1.0,0.0,0.0]))
        i += 1
    end
    return quartets
end

"""
    sameTaxa(Quartet, HybridNetwork)

Return `true` if all taxa in the quartet are represented in the network,
`false` if one or more taxa in the quartet does not appear in the network.

warning: the name can cause confusion. A more appropriate name might be
"in", or "taxain", or "taxonsubset", or etc.
"""
function sameTaxa(q::Quartet, t::HybridNetwork)
    for name in q.taxon
        in(name,t.names) || return false
    end
    return true
end

"""
    taxadiff(Vector{Quartet}, network; multiplealleles=true)
    taxadiff(DataCF, network; multiplealleles=true)

Return 2 vectors:

- taxa in at least 1 of the quartets but not in the network, and
- taxa in the network but in none of the quartets.

When `multiplealleles` is true, the taxon names that end with "__2"
are ignored in the quartets: they are not expected to appear in the
networks that users give as input, or get as output.
"""
function taxadiff(quartets::Vector{Quartet}, t::HybridNetwork;
                  multiplealleles=true::Bool)
    tq = tipLabels(quartets)
    secondallele = occursin.(Ref(r"__2$"), tq)
    for i in length(secondallele):-1:1
        secondallele[i] || continue
        basetax = match(r"(.*)__2$", tq[i]).captures[1]
        # if tq[i] = "mouse__2" for instance, then basetax = "mouse"
        if basetax in tq     # some other taxon is "mouse"
            deleteat!(tq, i) # delete "mouse__2" from tq IF "mouse" is present
        end
    end
    tn = tipLabels(t)
    return (setdiff(tq,tn), setdiff(tn,tq))
end

taxadiff(d::DataCF, t::HybridNetwork; multiplealleles=true::Bool) =
    taxadiff(d.quartet, t; multiplealleles=multiplealleles)


# extract & sort the union of taxa of list of gene trees
function unionTaxa(trees::Vector{HybridNetwork})
    taxa = reduce(union, tipLabels(t) for t in trees)
    return sort_stringasinteger!(taxa)
end

"""
    sort_stringasinteger!(taxa)

Sort a vector of strings `taxa`, numerically if
elements can be parsed as an integer, alphabetically otherwise.
"""
function sort_stringasinteger!(taxa)
    sortby = x->parse(Int,x)
    try
        parse.(Int,taxa)
    catch
        sortby = identity
    end
    sort!(taxa, by=sortby)
    return taxa
end

# extract & sort the union of taxa of list of quartets
function unionTaxa(quartets::Vector{Quartet})
    taxa = reduce(union, q.taxon for q in quartets)
    return sort_stringasinteger!(taxa)
end

unionTaxaTree(file::AbstractString) = unionTaxa(readInputTrees(file))

tipLabels(t::Vector{HybridNetwork}) = unionTaxa(t)
tipLabels(q::Vector{Quartet}) = unionTaxa(q)
tipLabels(d::DataCF) = unionTaxa(d.quartet)

"""
    calculateObsCFAll!(DataCF, taxa::Union{Vector{String}, Vector{Int}})

Calculate observed concordance factors:
update the `.quartet[i].obsCF` values of the `DataCF` object based on its .tree vector.

    calculateObsCFAll!(vector of quartets, vector of trees, taxa)

Calculate observed concordance factors:
update the `.obsCF` values of the quartets, based on the trees, and returns a new `DataCF` object
with these updated quartets and trees.

    calculateObsCFAll_noDataCF!(vector of quartets, vector of trees, taxa)

update the `.obsCF` values of the quartets based on the trees, but returns nothing.

Warning: all these functions need input trees (without any reticulations: h=0).

See also: [`countquartetsintrees`](@ref), which uses a faster algorithm,
processing each input tree only once.
`calculateObsCFAll_noDataCF!` processes each input tree `# quartet` times.
"""
function calculateObsCFAll!(dat::DataCF, taxa::Union{Vector{String}, Vector{Int}})
    calculateObsCFAll_noDataCF!(dat.quartet, dat.tree, taxa)
end

function calculateObsCFAll!(quartets::Vector{Quartet}, trees::Vector{HybridNetwork}, taxa::Union{Vector{String}, Vector{Int}})
    calculateObsCFAll_noDataCF!(quartets, trees, taxa)
    d = DataCF(quartets,trees)
    return d
end

function calculateObsCFAll_noDataCF!(quartets::Vector{Quartet}, trees::Vector{HybridNetwork}, taxa::Union{Vector{String}, Vector{Int}})
    println("calculating obsCF from $(length(trees)) gene trees and for $(length(quartets)) quartets")
    index = 1
    totalq = length(quartets)
    println("Reading in quartets...")
    r = round(1/totalq, digits=2)
    numq = (r > 0.02 ? totalq : 50)
    print("0+")
    for i in 1:numq
        print("-")
    end
    print("+100%")
    println("  ")
    print("  ")
    for q in quartets
        if round(index/totalq, digits=2) > 0.02
            print("*")
            index = 1
        end
        suma = 0
        sum12 = 0
        sum13 = 0
        sum14 = 0
        for t in trees
            isTree(t) || error("gene tree found in file that is a network $(writeTopology(t))")
            if sameTaxa(q,t)
                M = tree2Matrix(t,taxa) #fixit: way to reuse M? length(t.edge) will be different across trees
                res = extractQuartetTree(q,M,taxa)
                @debug "res is $(res)"
                if(res == 1)
                    sum12 += 1
                elseif(res == 2)
                    sum13 += 1
                elseif(res == 3)
                    sum14 += 1
                end
                suma += (res == 0) ? 0 : 1
            end
        end
        q.obsCF = [sum12/suma, sum13/suma, sum14/suma]
        q.ngenes = suma
        index += 1
    end
    println("  ")
    return nothing
end

"""
    countquartetsintrees(trees [, taxonmap=Dict{String,String}]; which=:all, weight_byallele=true)

Calculate the quartet concordance factors (CF) observed in the `trees` vector.
If present, `taxonmap` should map each allele name to it's species name.
To save to a file, first convert to a data frame using [`writeTableCF`](@ref).
When `which=:all`, quartet CFs are calculated for all 4-taxon sets.
(Other options are not implemented yet.)

The algorithm runs in O(mn⁴) where m is the number of trees and n is the number
of tips in the trees.

CFs are calculated at the species level only, that is, considering 4-taxon sets
made of 4 distinct species, even if the gene trees may have multiple alleles
from the same species. For 4 distinct species `a,b,c,d`, all alleles from
each species (`a` etc.) will be considered to calculate the quartet CF.

By default, each gene has a weight of 1. So if there are `n_a` alleles from `a`,
`n_b` alleles from `b` etc. in a given gene, then each set of 4 alleles has a
weight of `1/(n_a n_b b_c n_c)` in the calculation of the CF for `a,b,c,d`.
With option `weight_byallele=true`, then each set of 4 alleles is given a
weight of 1 instead. This inflates the total number of sets used to calculate
the quartet CFs (to something larger than the number of genes). This may also
affect the CF values if the number of alleles varies across genes: genes with
more alleles will be given more weight.

# examples
```jldoctest
julia> tree1 = readTopology("(E,(A,B),(C,D),O);"); tree2 = readTopology("(((A,B),(C,D)),E);");

julia> q,t = countquartetsintrees([tree1, tree2]);
Reading in trees, looking at 15 quartets in each...
0+--+100%
  **

julia> t # taxon order: t[i] = name of taxon number i
6-element Array{String,1}:
 "A"
 "B"
 "C"
 "D"
 "E"
 "O"

julia> length(q) # 15 four-taxon sets on 6 taxa
15

julia> q[1] # both trees agree on AB|CD: resolution 1
4-taxon set number 1; taxon numbers: 1,2,3,4
data: [1.0, 0.0, 0.0, 2.0]

julia> q[8] # tree 2 is missing O (taxon 6), tree 1 wants resolution 3: AO|CD
4-taxon set number 8; taxon numbers: 1,3,4,6
data: [0.0, 0.0, 1.0, 1.0]

julia> q[11] # tree 1 has ACEO unresolved, and tree 2 is missing O: no data for this quartet
4-taxon set number 11; taxon numbers: 1,3,5,6
data: [0.0, 0.0, 0.0, 0.0]

julia> tree1 = readTopology("(E,(a1,B),(a2,D),O);"); tree2 = readTopology("(((a1,a2),(B,D)),E);");

julia> q,t = countquartetsintrees([tree1, tree2], Dict("a1"=>"A", "a2"=>"A"); showprogressbar=false);

julia> t
5-element Array{String,1}:
 "A"
 "B"
 "D"
 "E"
 "O"

julia> q[1] # tree 1 has discordance: a1B|DE and a2D|BE. tree 2 has AE|BD for both alleles of A
4-taxon set number 1; taxon numbers: 1,2,3,4
data: [0.25, 0.25, 0.5, 2.0]

julia> q[3] # tree 2 is missing O (taxon 5), and a2 is unresolved in tree 1. There's only a1B|EO
4-taxon set number 3; taxon numbers: 1,2,4,5
data: [1.0, 0.0, 0.0, 0.5]

julia> df = writeTableCF(q,t); # to get a DataFrame that can be saved to a file later

julia> show(df, allcols=true, splitcols=false)
5×8 DataFrames.DataFrame
│ Row │ t1     │ t2     │ t3     │ t4     │ CF12_34 │ CF13_24 │ CF14_23 │ ngenes  │
│     │ String │ String │ String │ String │ Float64 │ Float64 │ Float64 │ Float64 │
├─────┼────────┼────────┼────────┼────────┼─────────┼─────────┼─────────┼─────────┤
│ 1   │ A      │ B      │ D      │ E      │ 0.25    │ 0.25    │ 0.5     │ 2.0     │
│ 2   │ A      │ B      │ D      │ O      │ 0.5     │ 0.5     │ 0.0     │ 1.0     │
│ 3   │ A      │ B      │ E      │ O      │ 1.0     │ 0.0     │ 0.0     │ 0.5     │
│ 4   │ A      │ D      │ E      │ O      │ 1.0     │ 0.0     │ 0.0     │ 0.5     │
│ 5   │ B      │ D      │ E      │ O      │ 0.0     │ 0.0     │ 0.0     │ 0.0     │

julia> # using CSV; CSV.write(df, "filename.csv");

julia> tree2 = readTopology("((A,(B,D)),E);");

julia> q,t = countquartetsintrees([tree1, tree2], Dict("a1"=>"A", "a2"=>"A"); weight_byallele=true);
Reading in trees, looking at 5 quartets in each...
0+--+100%
  **

julia> show(writeTableCF(q,t), allcols=true, splitcols=false)
5×8 DataFrames.DataFrame
│ Row │ t1     │ t2     │ t3     │ t4     │ CF12_34  │ CF13_24  │ CF14_23  │ ngenes  │
│     │ String │ String │ String │ String │ Float64  │ Float64  │ Float64  │ Float64 │
├─────┼────────┼────────┼────────┼────────┼──────────┼──────────┼──────────┼─────────┤
│ 1   │ A      │ B      │ D      │ E      │ 0.333333 │ 0.333333 │ 0.333333 │ 3.0     │
│ 2   │ A      │ B      │ D      │ O      │ 0.5      │ 0.5      │ 0.0      │ 2.0     │
│ 3   │ A      │ B      │ E      │ O      │ 1.0      │ 0.0      │ 0.0      │ 1.0     │
│ 4   │ A      │ D      │ E      │ O      │ 1.0      │ 0.0      │ 0.0      │ 1.0     │
│ 5   │ B      │ D      │ E      │ O      │ 0.0      │ 0.0      │ 0.0      │ 0.0     │
```
"""
function countquartetsintrees(tree::Vector{HybridNetwork},
                           taxonmap=Dict{String,String}()::Dict{String,String};
                           whichQ=:all::Symbol, weight_byallele=false::Bool,
                           showprogressbar=true::Bool)
    whichQ in [:all, :intrees] || error("whichQ must be either :all or :intrees, but got $whichQ")
    if isempty(taxonmap)
        taxa = unionTaxa(tree)
    else
        taxa = sort!(collect(Set(haskey(taxonmap, l.name) ? taxonmap[l.name] : l.name
                                 for t in tree for l in t.leaf)))
    end
    taxonnumber = Dict(taxa[i] => i for i in eachindex(taxa))
    ntax = length(taxa)
    nCk = nchoose1234(ntax) # matrix used to ranks 4-taxon sets
    qtype = MVector{4,Float64} # 4 floats: CF12_34, CF13_24, CF14_23, ngenes; initialized at 0.0
    if whichQ == :all
        numq = nCk[ntax+1,4]
        quartet = Vector{QuartetT{qtype}}(undef, numq)
        ts = [1,2,3,4]
        for qi in 1:numq
            quartet[qi] = QuartetT(qi, SVector{4}(ts), MVector(0.,0.,0.,0.))
            # next: find the 4-taxon set with the next rank,
            #       faster than using the direct mapping function
            ind = findfirst(x -> x>1, diff(ts))
            if ind === nothing ind = 4; end
            ts[ind] += 1
            for j in 1:(ind-1)
                ts[j] = j
            end
        end
    else
        error("whichQ = :intrees not implemented yet")
        # fixit: read all the trees, go through each edge, check if the current quartet list covers the edge, if not: add a quartet
    end
    totalt = length(tree)
    if showprogressbar
        nstars = (totalt < 50 ? totalt : 50)
        ntrees_perstar = (totalt/nstars)
        println("Reading in trees, looking at $numq quartets in each...")
        print("0+" * "-"^nstars * "+100%\n  ")
        stars = 0
        nextstar = Integer(ceil(ntrees_perstar))
    end
    for i in 1:totalt # number of times each quartet resolution is seen in each tree
        countquartetsintrees!(quartet, tree[i], whichQ, weight_byallele, nCk, taxonnumber, taxonmap)
        if showprogressbar && i >= nextstar
            print("*")
            stars += 1
            nextstar = Integer(ceil((stars+1) * ntrees_perstar))
        end
    end
    showprogressbar && print("\n")
    # normalize counts to frequencies & number of genes
    for q in quartet
        d = q.data
        d[4] = d[1]+d[2]+d[3] # number of genes
        if d[4] > 0.0
            d[1:3] /= d[4]
        end # otherwise: no genes with data on this quartet (missing taxa or polytomy): leave all zeros (NaN if /0)
    end
    return quartet, taxa
end
function countquartetsintrees!(quartet::Vector{QuartetT{MVector{4,Float64}}},
            tree::HybridNetwork, whichQ::Symbol, weight_byallele::Bool, nCk::Matrix,
            taxonnumber::Dict{String,Int64}, taxonmap::Dict{String,String})
    tree.numHybrids == 0 || error("input phylogenies must be trees")
    # next: reset node & edge numbers so that they can be used as indices: 1,2,3,...
    resetNodeNumbers!(tree; checkPreorder=true, ape=false) # leaves first & post-order
    resetEdgeNumbers!(tree)
    # next: build list leaf number -> species ID, using the node name then taxon map
    nleaves = length(tree.leaf)
    nnodes  = length(tree.node)
    taxID = Vector{Int}(undef, nleaves)
    for n in tree.leaf
        taxID[n.number] = haskey(taxonmap, n.name) ? taxonnumber[taxonmap[n.name]] : taxonnumber[n.name]
    end
    # number of individuals from each species: needed to weigh the quartets at the individual level
    # weight of t1,t2,t3,t4: 1/(taxcount[t1]*taxcount[t2]*taxcount[t3]*taxcount[t4])
    taxcount = zeros(Int, length(taxonnumber))
    for ti in taxID taxcount[ti] += 1; end
    # next: build data structure to get descendant / ancestor clades
    below,above = ladderpartition(tree) # re-checks that node numbers can be used as indices, with leaves first
    # below[n][1:2 ]: left & clades below node number n
    # above[n][1:end]: grade of clades above node n
    for n in (nleaves+1):nnodes # loop through internal nodes indices only
        bn = below[n]
        an = above[n]
        for c1 in 2:length(bn) # c = child clade, loop over all pairs of child clades
          for t1 in bn[c1]     # pick 1 tip from each child clade
            s1 = taxID[t1]
            for c2 in 1:(c1-1)
              for t2 in bn[c2]
                s2 = taxID[t2]
                s1 != s2 || continue # skip quartets that have repeated species
                t12max = max(t1,t2)
                leftweight = 1/(taxcount[s1]*taxcount[s2])
                for p1 in 1:length(an) # p = parent clade
                  for t3i in 1:length(an[p1])
                    t3 = an[p1][t3i]
                    s3 = taxID[t3]
                    (s3 != s1 && s3 != s2) || continue
                    for t4i in 1:(t3i-1) # pick 2 distinct tips from the same parent clade
                        t4 = an[p1][t4i]
                        t3 > t12max || t4 > t12max || continue   # skip: would be counted twice otherwise
                        s4 = taxID[t4]
                        (s4 != s1 && s4 != s2 && s4 != s3) || continue
                        rank,res = quartetRankResolution(s1, s2, s3, s4, nCk)
                        weight = ( weight_byallele ? 1.0 : leftweight / (taxcount[s3]*taxcount[s4]) )
                        quartet[rank].data[res] += weight
                    end
                    for p2 in 1:(p1-1) # distinct parent clade: no risk of counting twice
                      for t4 in an[p2]
                        s4 = taxID[t4]
                        (s4 != s1 && s4 != s2 && s4 != s3) || continue
                        rank,res = quartetRankResolution(s1, s2, s3, s4, nCk)
                        weight = ( weight_byallele ? 1.0 : leftweight / (taxcount[s3]*taxcount[s4]) )
                        quartet[rank].data[res] += weight
                      end
                    end
                  end
                end
              end
            end
          end
        end
    end
end

function quartetRankResolution(t1::Int, t2::Int, t3::Int, t4::Int, nCk::Matrix)
    # quartet: t1 t2 | t3 t4, but indices have not yet been ordered: t1<t2, t3<t4, t1<min(t3,t4)
    if t3 > t4 # make t3 smallest of t3, t4
        (t3,t4) = (t4,t3)
    end
    if t1 > t2 # make t1 smallest of t1, t2
        (t1,t2) = (t2,t1)
    end
    if t1 > t3 # swap t1 with t3, t2 with t4 - makes t1 smallest
        (t1,t3) = (t3,t1)
        (t2,t4) = (t4,t2)
    end
    if t2 < t3 # t2 2nd smallest: order t1 < t2 < t3 < t4
        resolution = 1; # 12|34 after ordering indices
        rank = quartetrank(t1, t2, t3, t4, nCk)
    else # t3 2nd smallest
        if t2 < t4 # order t1 < t3 < t2 < t4
            resolution = 2; # 13|24 after ordering
            rank = quartetrank(t1, t3, t2, t4, nCk);
        else # order t1 < t3 < t4 < t2
            resolution = 3; # 14|23 after ordering
            rank = quartetrank(t1, t3, t4, t2, nCk);
        end
    end
    return rank, resolution
end

"""
    readInputData(trees, quartetfile, whichQuartets, numQuartets, writetable, tablename, writeQfile, writesummary)
    readInputData(trees, whichQuartets, numQuartets, taxonlist,   writetable, tablename, writeQfile, writesummary)

Read gene trees and calculate the observed quartet concordance factors (CF),
that is, the proportion of genes (and the number of genes) that display each
quartet for a given list of four-taxon sets.

Input:

- `trees`: name of a file containing a list of input gene trees,
  or vector of trees (`HybridNetwork` objects)

Optional arguments (defaults):

- `quartetfile`: name of a file containing a list of quartets, or more precisely,
  a list of four-taxon sets
- `whichQuartets` (`:all`): which quartets to sample.
  `:all` for all of them, `:rand` for a random sample.
- `numQuartets`: number of quartets in the sample.
  default: total number of quartets if `whichQuartets=:all`
  and 10% of total if `whichQuartets=:rand`
- `taxonlist` (all in the input gene trees):
  If `taxonlist` is used, `whichQuartets` will consist of *all* sets of 4 taxa in the `taxonlist`. 
- `writetable` (true): write the table of observed CF?
- `tablename` ("tableCF.txt"): if `writetable` is true, the table of observed CFs is write to file `tablename`
- `writeQfile` (false): write intermediate file with sampled quartets?
- `writesummary` (true): write a summary file?
  if so, the summary will go in file "summaryTreesQuartets.txt".

Uses [`calculateObsCFAll!`](@ref), which implements a slow algorithm.

See also:
[`countquartetsintrees`](@ref), which uses a much faster algorithm;
[`readTrees2CF`](@ref), which is basically a re-naming of `readInputData`, and
[`readTableCF`](@ref) to read a table of quartet CFs directly.
"""
function readInputData(treefile::AbstractString, quartetfile::AbstractString, whichQ::Symbol, numQ::Integer, writetab::Bool, filename::AbstractString, writeFile::Bool, writeSummary::Bool)
    if writetab
        if(filename == "none")
            filename = "tableCF.txt" # "tableCF$(string(integer(time()/1000))).txt"
        end
        if (isfile(filename) && filesize(filename) > 0)
           error("""file $(filename) already exists and is non-empty. Cannot risk to erase data.
                    Choose a different CFfile name, use writeTab=false, or read the existing file
                    with readTableCF(\"$(filename)\")""")
        end
    end
    println("read input trees from file $(treefile)\nand quartetfile $(quartetfile)")
    trees = readInputTrees(treefile)
    readInputData(trees, quartetfile, whichQ, numQ, writetab, filename, writeFile, writeSummary)
end

readInputData(treefile::AbstractString, quartetfile::AbstractString, whichQ::Symbol, numQ::Integer, writetab::Bool) = readInputData(treefile, quartetfile, whichQ, numQ, writetab, "none", false, true)
readInputData(treefile::AbstractString, quartetfile::AbstractString, whichQ::Symbol, numQ::Integer) = readInputData(treefile, quartetfile, whichQ, numQ, true, "none", false, true)
readInputData(treefile::AbstractString, quartetfile::AbstractString, writetab::Bool, filename::AbstractString) = readInputData(treefile, quartetfile, :all, 0, writetab, filename, false, true)

function readInputData(trees::Vector{HybridNetwork}, quartetfile::AbstractString, whichQ::Symbol, numQ::Integer, writetab::Bool, filename::AbstractString, writeFile::Bool, writeSummary::Bool)
    if(whichQ == :all)
        numQ == 0 || @warn "set numQ=$(numQ) but whichQ is not rand, so all quartets will be used and numQ will be ignored. If you want a specific number of 4-taxon subsets not random, you can input with the quartetfile option"
        println("will use all quartets in file $(quartetfile)")
        quartets = readListQuartets(quartetfile)
    elseif(whichQ == :rand)
        if(numQ == 0)
            @warn "not specified numQ but whichQ=rand, so 10% of quartets will be sampled" #handled inside randQuartets
        else
            println("will take a random sample of $(numQ) 4-taxon sets from file $(quartetfile)")
        end
        allquartets = readListQuartets(quartetfile)
        quartets = randQuartets(allquartets,numQ,writeFile)
    else
        error("unknown symbol for whichQ $(whichQ), should be either all or rand")
    end
    d = calculateObsCFAll!(quartets,trees, unionTaxa(trees))
    if(writetab)
        if(filename == "none")
            filename = "tableCF.txt" # "tableCF$(string(integer(time()/1000))).txt"
        end
        if (isfile(filename) && filesize(filename) > 0)
           error("""file $(filename) already exists and is non-empty. Cannot risk to erase data.
                    Choose a different CFfile name, use writeTab=false, or read the existing file
                    with readTableCF(\"$(filename)\")""")
        end
        println("\ntable of obsCF printed to file $(filename)")
        df = writeTableCF(d)
        CSV.write(filename,df)
    end
    #descData(d,"summaryTreesQuartets$(string(integer(time()/1000))).txt")
    writeSummary && descData(d,"summaryTreesQuartets.txt")
    return d
end


function readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer, taxa::Union{Vector{String}, Vector{Int}}, writetab::Bool, filename::AbstractString, writeFile::Bool, writeSummary::Bool)
    if writetab
        if(filename == "none")
            filename = "tableCF.txt" # "tableCF$(string(integer(time()/1000))).txt"
        end
        if (isfile(filename) && filesize(filename) > 0)
           error("""file $(filename) already exists and is non-empty. Cannot risk to erase data.
                    Choose a different CFfile name, use writeTab=false, or read the existing file
                    with readTableCF(\"$(filename)\")""")
        end
    end
    println("read input trees from file $(treefile). no quartet file given.")
    trees = readInputTrees(treefile)
    readInputData(trees, whichQ, numQ, taxa, writetab, filename, writeFile, writeSummary)
end

readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer, taxa::Union{Vector{String}, Vector{Int}}, writetab::Bool) = readInputData(treefile, whichQ, numQ, taxa, writetab, "none", false, true)
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer, taxa::Union{Vector{String}, Vector{Int}}) = readInputData(treefile, whichQ, numQ, taxa, true, "none",false, true)
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer, writetab::Bool, filename::AbstractString) = readInputData(treefile, whichQ, numQ, unionTaxaTree(treefile), writetab, filename,false, true)
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer, writetab::Bool) = readInputData(treefile, whichQ, numQ, unionTaxaTree(treefile), writetab, "none",false, true)
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer) = readInputData(treefile, whichQ, numQ, unionTaxaTree(treefile), true, "none",false, true)
readInputData(treefile::AbstractString) = readInputData(treefile, :all, 0, unionTaxaTree(treefile), true, "none",false, true)
readInputData(treefile::AbstractString,taxa::Union{Vector{String}, Vector{Int}}) = readInputData(treefile, :all, 0, taxa, true, "none",false, true)
# above: the use of unionTaxaTree to set the taxon set
#        is not good: need to read the tree file twice: get the taxa, then get the trees
#        this inefficiency was fixed in readTrees2CF

function readInputData(trees::Vector{HybridNetwork}, whichQ::Symbol, numQ::Integer, taxa::Union{Vector{String}, Vector{Int}}, writetab::Bool, filename::AbstractString, writeFile::Bool, writeSummary::Bool)
    if(whichQ == :all)
        numQ == 0 || @warn "set numQ=$(numQ) but whichQ=all, so all quartets will be used and numQ will be ignored. If you want a specific number of 4-taxon subsets not random, you can input with the quartetfile option"
        quartets = allQuartets(taxa,writeFile)
        println("will use all quartets on $(length(taxa)) taxa")
    elseif(whichQ == :rand)
        if(numQ == 0)
            @warn "not specified numQ with whichQ=rand, so 10% of quartets will be sampled" #handled inside randQuartets
        else
            println("will use a random sample of $(numQ) 4-taxon sets ($(round((100*numQ)/binomial(length(taxa),4), digits=2)) percent) on $(length(taxa)) taxa")
        end
        quartets = randQuartets(taxa,numQ, writeFile)
    else
        error("unknown symbol for whichQ $(whichQ), should be either all or rand")
    end
    d = calculateObsCFAll!(quartets,trees,taxa)
    if writetab
        if(filename == "none")
            filename = "tableCF.txt"
        end
        println("table of obsCF printed to file $(filename)")
        df = writeTableCF(d)
        CSV.write(filename,df)
    end
    #descData(d,"summaryTreesQuartets$(string(integer(time()/1000))).txt")
    writeSummary && descData(d,"summaryTreesQuartets.txt")
    return d
end



# rename the function readInputData to make it more user-friendly
"""
    readTrees2CF(treefile)
    readTrees2CF(vector of trees)

Read trees in parenthetical format from a file, or take a vector of trees already read,
and calculate the proportion of these trees having a given quartet (concordance factor: CF),
for all quartets or for a sample of quartets.
Optional arguments include:

- quartetfile: name of text file with list of 4-taxon subsets to be analyzed. If none is specified, the function will list all possible 4-taxon subsets.
- whichQ="rand": to choose a random sample of 4-taxon subsets
- numQ: size of random sample (ignored if whichQ is not set to "rand")
- writeTab=false: does not write the observedCF to a table (default true)
- CFfile: name of file to save the observedCF (default tableCF.txt)
- writeQ=true: save intermediate files with the list of all 4-taxon subsets and chosen random sample (default false).
- writeSummary: write descriptive stats of input data (default: true)
- nexus: if true, it assumes the gene trees are written in nexus file (default: false)

See also:
[`countquartetsintrees`](@ref), which uses a much faster algorithm;
[`readTableCF`](@ref) to read a table of quartet CFs directly.
"""
function readTrees2CF(treefile::AbstractString; quartetfile="none"::AbstractString, whichQ="all"::AbstractString, numQ=0::Integer,
                      writeTab=true::Bool, CFfile="none"::AbstractString,
                      taxa=Vector{String}()::Union{Vector{String},Vector{Int}},
                      writeQ=false::Bool, writeSummary=true::Bool, nexus=false::Bool)
    trees = (nexus ? readNexusTrees(treefile, readTopologyUpdate, false, false) : readInputTrees(treefile))
    if length(taxa)==0        # unionTaxa(trees) NOT default argument:
      taxa = unionTaxa(trees) # otherwise: tree file is read twice
    end
    readTrees2CF(trees, quartetfile=quartetfile, whichQ=whichQ, numQ=numQ, writeTab=writeTab,
                 CFfile=CFfile, taxa=taxa, writeQ=writeQ, writeSummary=writeSummary)
end

# same as before, but with input vector of HybridNetworks
function readTrees2CF(trees::Vector{HybridNetwork}; quartetfile="none"::AbstractString, whichQ="all"::AbstractString, numQ=0::Integer, writeTab=true::Bool, CFfile="none"::AbstractString, taxa=unionTaxa(trees)::Union{Vector{String},Vector{Int}}, writeQ=false::Bool, writeSummary=true::Bool)
    whichQ == "all" || whichQ == "rand" ||
        error("whichQ should be all or rand, not $(whichQ)")
    if(quartetfile == "none")
        readInputData(trees, Symbol(whichQ), numQ, taxa, writeTab, CFfile, writeQ, writeSummary)
    else
        readInputData(trees, quartetfile, Symbol(whichQ), numQ, writeTab, CFfile, writeQ, writeSummary)
    end
end

# ---------------------- descriptive stat for input data ----------------------------------

# function to check how taxa is represented in the input trees
function taxaTreesQuartets(trees::Vector{HybridNetwork}, quartets::Vector{Quartet},s::IO)
    taxaT = unionTaxa(trees)
    taxaQ = unionTaxa(quartets)
    dif = symdiff(taxaT,taxaQ)
    isempty(dif) ? write(s,"\n same taxa in gene trees and quartets: $(taxaT)\n") :
                   write(s,"\n $(length(dif)) different taxa found in gene trees and quartets. \n Taxa $(intersect(taxaT,dif)) in trees, not in quartets; and taxa $(intersect(taxaQ,dif)) in quartets, not in trees\n")
    u = union(taxaT,taxaQ)
    for taxon in u
        numT = taxonTrees(taxon,trees)
        #numQ = taxonQuartets(taxon,quartets)
        write(s,"Taxon $(taxon) appears in $(numT) input trees ($(round(100*numT/length(trees), digits=2)) %)\n")  #and $(numQ) quartets ($(round(100*numQ/length(quartets), digits=2)) %)\n")
    end
end

taxaTreesQuartets(trees::Vector{HybridNetwork}, quartets::Vector{Quartet}) = taxaTreesQuartets(trees, quartets, stdout)

# function that counts the number of trees in which taxon appears
function taxonTrees(taxon::AbstractString, trees::Vector{HybridNetwork})
    suma = 0
    for t in trees
        suma += in(taxon,t.names) ? 1 : 0
    end
    return suma
end

# function that counts the number of quartets in which taxon appears
function taxonQuartets(taxon::AbstractString, quartets::Vector{Quartet})
    suma = 0
    for q in quartets
        suma += in(taxon,q.taxon) ? 1 : 0
    end
    return suma
end


# function to create descriptive stat from input data, will save in stream sout
# which can be a file or stdout
# default: send to stdout
# pc: only 4-taxon subsets with percentage of gene trees less than pc will be printed (default 70%)
function descData(d::DataCF, sout::IO, pc::Float64)
    0<=pc<=1 || error("percentage of missing genes should be between 0,1, not: $(pc)")
    if !isempty(d.tree)
        print(sout,"data consists of $(d.numTrees) gene trees and $(d.numQuartets) 4-taxon subsets\n")
        taxaTreesQuartets(d.tree,d.quartet,sout)
        print(sout,"----------------------------\n\n")
        print(sout,"will print below only the 4-taxon subsets with data from <= $(round((pc)*100, digits=2))% genes\n")
        for q in d.quartet
            percent  = q.ngenes == -1.0 ? 0.0 : round(q.ngenes/d.numTrees*100, digits=2)
            if percent < pc
                print(sout,"4-taxon subset $(q.taxon) obsCF constructed with $(round(q.ngenes)) gene trees ($(percent)%)\n")
            end
        end
        print(sout,"----------------------------\n\n")
    else
        if !isempty(d.quartet)
            print(sout,"data consists of $(d.numQuartets) 4-taxon subsets")
            taxa=unionTaxa(d.quartet)
            print(sout,"\nTaxa: $(taxa)\n")
            print(sout,"Number of Taxa: $(length(taxa))\n")
            numQ = binomial(length(taxa),4);
            print(sout,"Maximum number of 4-taxon subsets: $(numQ). Thus, $(round(100*d.numQuartets/numQ, digits=2)) percent of 4-taxon subsets sampled\n")
        end
    end
end

function descData(d::DataCF, filename::AbstractString,pc::Float64)
    println("descriptive stat of input data printed to file $(filename)")
    s = open(filename, "w")
    descData(d,s,pc)
    close(s)
end

descData(d::DataCF, sout::IO) = descData(d, sout,0.7)
descData(d::DataCF) = descData(d, stdout,0.7)
descData(d::DataCF,pc::Float64) = descData(d, stdout,pc)
descData(d::DataCF, filename::AbstractString) = descData(d, filename,0.7)

"""
`summarizeDataCF(d::DataCF)`

function to summarize the information contained in a DataCF object. It has the following optional arguments:
- filename: if provided, the summary will be saved in the filename, not to screen
- pc (number between (0,1)): threshold of percentage of missing genes to identify 4-taxon subsets with fewer genes than the threshold
"""
function summarizeDataCF(d::DataCF; filename="none"::AbstractString, pc=0.7::Float64)
    0<=pc<=1 || error("percentage of missing genes should be between 0,1, not: $(pc)")
    if filename == "none"
        descData(d,stdout,pc)
    else
        descData(d,filename,pc)
    end
end

# -------- branch length estimate in coalescent units on species tree ------

"""
    updateBL!(net::HybridNetwork, d::DataCF)

Update internal branch lengths of `net` based on the average quartet concordance
factor (CF) across all quartets that exactly correspond to a given branch:
new branch length = `-log(3/2(1-mean(CF observed in d)))`.
`net` is assumed to be a tree, such that the above equation holds.
"""
function updateBL!(net::HybridNetwork,d::DataCF)
    if !isTree(net)
        @error "updateBL! was created for a tree, and net here is not a tree, so no branch lengths updated"
    end
    parts = edgesParts(net)
    df = makeTable(net,parts,d)
    x = by(df, :edge, Nquartets= :CF => length,
                      edgeL = :CF => x -> -log(3/2*(1. - mean(x))))
    # ommitting columns: meanCF= :CF => mean, sdCF= :CF => std
    edges = x[!,:edge]
    lengths = x[!,:edgeL]
    for i in 1:length(edges)
        ind = getIndexEdge(edges[i],net) # helpful error if not found
        if net.edge[ind].length < 0.0 || net.edge[ind].length==1.0
            # readTopologyLevel1 changes missing branch length to 1.0
            setLength!(net.edge[ind], (lengths[i] > 0 ? lengths[i] : 0.0))
        end
    end
    for e in net.edge
        if e.length < 0.0 # some edges might have *no* quartet in the data
            setLength!(e, 1.0)
        end
    end
    return x
end


# function to get part1,part2,part3,part4 for each edge in net.edge
# returns a EdgeParts object
function edgesParts(net::HybridNetwork)
    parts = EdgeParts[] #vector to hold part1,...,part4 for each edge
    for e in net.edge
        if(isInternalEdge(e))
            length(e.node) == 2 || error("strange edge with $(length(e.node)) nodes instead of 2")
            n1 = e.node[1]
            n2 = e.node[2]
            e11,e12 = hybridEdges(n1,e)
            e21,e22 = hybridEdges(n2,e)
            part1 = Node[]
            part2 = Node[]
            part3 = Node[]
            part4 = Node[]
            getDescendants!(getOtherNode(e11,n1),e11,part1)
            getDescendants!(getOtherNode(e12,n1),e12,part2)
            getDescendants!(getOtherNode(e21,n2),e21,part3)
            getDescendants!(getOtherNode(e22,n2),e22,part4)
            push!(parts, EdgeParts(e.number,part1,part2,part3,part4))
        end
    end
    return parts
end

# aux function to traverse the network from a node and an edge
# based on traverseContainRoot
# warning: it does not go accross hybrid node, minor hybrid edge
# there is another getDescendants in update.jl for updatePartition
function getDescendants!(node::Node, edge::Edge, descendants::Array{Node,1})
    if(node.leaf)
        push!(descendants, node)
    else
        for e in node.edge
            if(!isEqual(edge,e) && e.isMajor)
                other = getOtherNode(e,node);
                getDescendants!(other,e, descendants);
            end
        end
    end
end

# function to make table to later use in updateBL
# uses vector parts obtained from edgeParts function
function makeTable(net::HybridNetwork, parts::Vector{EdgeParts},d::DataCF)
    df = DataFrame(edge=Int[],t1=AbstractString[],t2=AbstractString[],t3=AbstractString[],t4=AbstractString[],resolution=AbstractString[],CF=Float64[])
    sortedDataQ = [sort(q.taxon) for q in d.quartet]
    for p in parts #go over internal edges too
        for t1 in p.part1
            for t2 in p.part2
                for t3 in p.part3
                    for t4 in p.part4
                        tx1 = net.names[t1.number]
                        tx2 = net.names[t2.number]
                        tx3 = net.names[t3.number]
                        tx4 = net.names[t4.number]
                        nam = [tx1,tx2,tx3,tx4]
                        snam = sort(nam)
                        row = findall(isequal(snam), sortedDataQ)
                        for r in row # nothing if tax set not found: length(row)=0
                          col,res = resolution(nam,d.quartet[r].taxon)
                          push!(df, [p.edgenum,tx1,tx2,tx3,tx4,res,d.quartet[r].obsCF[col]])
                        end
                    end
                end
            end
        end
    end
    return df
end

# function to determine the resolution of taxa picked from part1,2,3,4 and DataCF
# names: taxa from part1,2,3,4
# rownames: taxa from table of obsCF
function resolution(names::Vector{String},rownames::Vector{String})
    length(names) == length(rownames) || error("names and rownames should have the same length")
    length(names) == 4 || error("names should have 4 entries, not $(length(names))")
    bin = [n == names[1] || n == names[2] ? 1 : 0 for n in rownames]
    if(bin == [1,1,0,0] || bin == [0,0,1,1])
        return 1,"12|34"
    elseif(bin == [1,0,1,0] || bin == [0,1,0,1])
        return 2,"13|24"
    elseif(bin == [1,0,0,1] || bin == [0,1,1,0])
        return 3,"14|23"
    else
        error("strange resolution $(bin)")
    end
end


# function to extract a quartet from a matrix M
# obtained from tree2Matrix (defined in file compareNetworks.jl)
# this function is meant to replace extractQuartet! in calculateObsCFAll
# input: Quartet, Matrix, vector of taxa names
# returns 1 if quartet found is 12|34, 2 if 13|24, 3 if 14|23, and 0 if not found
function extractQuartetTree(q::Quartet, M::Matrix{Int},S::Union{Vector{String},Vector{Int}})
    @debug "extractQuartet: $(q.taxon)"
    @debug "matrix: $(M)"
    inds = indexin(q.taxon, S)
    if any(isnothing, inds)
        error("some taxon in quartet $(q.taxon) not found in list of all species $(S)")
    end
    subM = M[:, inds.+1]
    @debug "subM: $(subM)"
    for r in 1:size(subM,1) #rows in subM
        @debug "subM[r,:]: $(subM[r,:])"
        if subM[r,:] == [0,0,1,1] || subM[r,:] == [1,1,0,0]
            return 1
        elseif subM[r,:] == [0,1,0,1] || subM[r,:] == [1,0,1,0]
            return 2
        elseif subM[r,:] == [0,1,1,0] || subM[r,:] == [1,0,0,1]
            return 3
        end
    end
    return 0
end


# function that will give the qth quartet without making a list of all quartets
# input: n number of taxa, q desired index of quartet
# returns vector of int, e.g. 1234
function whichQuartet(n::Int, q::Int)
    p = 4
    q <= binom(n,p) || error("the index for the quartet $(q) needs to be less than choose(n,4)=$(binom(n,p))")
    n > 4 || error("there must be at least 5 taxa, not $(n)")
    quartet = Int[]
    while(n > 1)
        abs = binom(n-1,p) #fixit: we don't want to compute this, we want to look for it in a table
        if(q > abs)
            push!(quartet,n)
            n -= 1
            p -= 1
            q = q-abs
        else
            n -= 1
        end
    end
    if(length(quartet) == 3)
        push!(quartet,1)
    end
    quartet = quartet[[4,3,2,1]] #sort
    return quartet
end

# function to write a quartet on integer to taxon names
# it creates a Quartet type
# input: list of taxa names, vector of integers (each integer corresponds to a taxon name),
# num is the number of the Quartet
# assumes taxa is sorted already
function createQuartet(taxa::Union{Vector{String},Vector{Int}},qvec::Vector{Int}, num::Int)
    length(qvec) == 4 || error("a quartet should have only 4 taxa, not $(length(qvec))")
    names = String[]
    for i in qvec
        i <= length(taxa) || error("want taxon number $(i) in list of taxon names $(taxa) which has only $(length(taxa)) names")
        push!(names,string(taxa[i]))
    end
    return Quartet(num,names,[1.0,0.0,0.0])
end

"""
    readNexusTrees(filename::AbstractString, treereader=readTopology::Function [, args...])

Read trees in nexus-formatted file and return a vector of `HybridNetwork`s.
For the nexus format, see Maddison, Swofford & Maddison (1997)
https://doi.org/10.1093/sysbio/46.4.590.
The optional arguments are passed onto the individual tree reader.

Warnings:
- "translate" tables are not supported yet
- only the first tree block is read
"""
function readNexusTrees(file::AbstractString, treereader=readTopology::Function, args...)
    vnet = HybridNetwork[]
    rx_start = r"^\s*begin\s+trees\s*;"i
    rx_end = r"^\s*end\s*;"i
    rx_tree = r"^\s*tree\s+[^(]+(\([^;]*;)"i
    # spaces,"Tree",spaces,any_symbols_other_than_(, then we capture:
    # ( any_symbols_other_than_; ;
    treeblock = false # whether we are currently reading the TREE block or not
    open(file) do s
        numl = 0
        for line in eachline(s)
            numl += 1
            if treeblock # currently reading trees, check for END signal
                occursin(rx_end, line) && break # break if end of tree block
            else # not reading trees: check for the BEGIN signal
                if occursin(rx_start, line) treeblock=true; end
                continue # to next line, either way
            end
            # if we get there, it's that we are inside the treeblock (true) and no END signal yet
            m = match(rx_tree, line)
            m != nothing || continue # continue to next line if no match
            phy = m.captures[1] # string
            try
                push!(vnet, treereader(phy, args...)) # readTopologyUpdate(phy,false)
            catch err
                print("skipped phylogeny on line $(numl) of file $file: ")
                if :msg in fieldnames(typeof(err)) println(err.msg); else println(typeof(err)); end
            end
        end
    end
    return vnet # consistent output type: HybridNetwork vector. might be of length 0.
end

