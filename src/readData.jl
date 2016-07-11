# functions to read data
# originally in functions.jl
# Claudia March 2015
# based on Julia 0.4
# Claudia October 2015

# ----- read data --------

# function to write a csv table from the expCF of an
# array of quartets
# warning: does not check if the expCF have been calculated
function writeExpCF(quartets::Array{Quartet,1})
    df = DataFrames.DataFrame(t1=ASCIIString[],t2=ASCIIString[],t3=ASCIIString[],t4=ASCIIString[],
                              CF12_34=Float64[],CF13_24=Float64[],CF14_23=Float64[])
    for(q in quartets)
        length(q.taxon) == 4 || error("quartet $(q.number) does not have 4 taxa")
        length(q.qnet.expCF) == 3 || error("quartet $(q.number) does have qnet with 3 expCF")
        push!(df, [q.taxon[1],q.taxon[2],q.taxon[3],q.taxon[4],q.qnet.expCF[1],q.qnet.expCF[2],q.qnet.expCF[3]])
    end
    return df
end

writeExpCF(d::DataCF) = writeExpCF(d.quartet)

# function to write a csv table from the obsCF of an
# array of quartets
function writeObsCF(quartets::Array{Quartet,1})
    df = DataFrames.DataFrame(t1=ASCIIString[],t2=ASCIIString[],t3=ASCIIString[],t4=ASCIIString[],
                              CF12_34=Float64[],CF13_24=Float64[],CF14_23=Float64[],ngenes=Int[])
    for(q in quartets)
        length(q.taxon) == 4 || error("quartet $(q.number) does not have 4 taxa")
        length(q.obsCF) == 3 || error("quartet $(q.number) does have qnet with 3 expCF")
        push!(df, [q.taxon[1],q.taxon[2],q.taxon[3],q.taxon[4],q.obsCF[1],q.obsCF[2],q.obsCF[3],q.ngenes])
    end
    return df
end

writeObsCF(d::DataCF) = writeObsCF(d.quartet)

# function that takes a dataframe and creates a DataCF object
# has ! because it can modify the dataframe inside
"""
    readTableCF(file)
    readTableCF(data frame)
    readTableCF!(data frame)

Read a file or DataFrame object containing a table of concordance factors (CF),
with one row per 4-taxon set. The first 4 columns are assumed to give the labels
of the 4 taxa in each set (tx1, tx2, tx3, tx4).
Columns containing the CFs are assumed to be named
'CF12_34', 'CF13_24' and 'CF14_23', or else are assumed to be columns 5,6,7.
If present, a column named 'ngenes' will be used to get the number of loci
used to estimate the CFs for each 4-taxon set.

Optional arguments:

- summaryfile: if specified, a summary file will be created with that name.
- sep (for the first form only): to specify the type of separator in the file,
  with single quotes: sep=';'. Default is a `csv` file, i.e. `sep=','`.

The last version modifies the input data frame, if species are represented by multiple alleles
for instance (see `readTableCF!`).
"""
# warning: file AbstractString bc it can be read as UTF8String
function readTableCF(file::AbstractString; sep=','::Char, summaryfile=""::AbstractString)
    df = readtable(file,separator=sep)
    readTableCF!(df, summaryfile=summaryfile)
end

function readTableCF(df0::DataFrames.DataFrame; summaryfile=""::AbstractString)
    df = deepcopy(df0)
    readTableCF!(df, summaryfile=summaryfile)
end

function readTableCF!(df::DataFrames.DataFrame; summaryfile=""::AbstractString)
    DEBUG && println("assume the numbers for the taxon read from the observed CF table match the numbers given to the taxon when creating the object network")
    obsCFcol = [findfirst(DataFrames.names(df), :CF12_34),
                findfirst(DataFrames.names(df), :CF13_24),
                findfirst(DataFrames.names(df), :CF14_23)] # use DataFrames.index( ).names ow
    ngenecol =  findfirst(DataFrames.names(df), :ngenes)
    withngenes = ngenecol>0
    if findfirst(obsCFcol, 0) > 0 # one or more col names for CFs were not found
        size(df,2) == (withngenes ? 8 : 7) ||
          warn("""Column names for quartet concordance factors (CFs) were not recognized.
          Was expecting CF12_34, CF13_24 and CF14_23 for the columns with CF values.
          Will assume that the first 4 columns give the taxon names, and that columns 5-7 give the CFs.""")
        obsCFcol = [5,6,7] # assuming CFs are in columns 5,6,7, with colname mismatch
    end
    if minimum(obsCFcol) <= 4
        error("CFs found in columns $obsCFcol, but taxon labels expected in columns 1-4")
    end
    # fixit: what about columns giving the taxon names: always assumed to be columns 1-4? No warning if not?
    columns = [[1,2,3,4]; obsCFcol]
    if withngenes  push!(columns, ngenecol)  end

    d = readTableCF!(df, columns)

    if withngenes && d.numTrees == -1
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
    if(!isempty(repSpecies))    # fixit: cleanAlleleDF! and mergeRows! are time consuming but many times not needed.
        mergeRows!(df,co)       # add option to skip them, for the user to say that each tip appears once only?
    end
    quartets = Quartet[]
    for(i in 1:size(df,1))
        push!(quartets,Quartet(i,string(df[i,co[1]]),string(df[i,co[2]]),string(df[i,co[3]]),string(df[i,co[4]]),
                               [df[i,co[5]],df[i,co[6]],df[i,co[7]]]))
        if(withngenes)
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
`columns` should have **3** columns numbers for the 3 CFs in this order: 12_34, 13_24 and 14_23.
Assume the same 4-taxon sets in `DataCF` and in the data frame, and in the same order,
but this assumption is *not checked* (for speed, e.g. during bootstrapping).
"""
# WARNING: assumes df has a single row per 4-taxon sets (multiple ind already merged)
#                  dcf created earlier from df.
function readTableCF!(datcf::DataCF, df::DataFrame, cols::Vector{Int})
    for(i in 1:size(df,1))
        for (j in 1:3)
            datcf.quartet[i].obsCF[j] = df[i,cols[j]]
        end
    end
end


# ---------------- read input gene trees and calculate obsCF ----------------------

# function to read a file and create one object per line read
# (each line starting with "(" will be considered a topology)
# the file can have extra lines that are ignored
# returns an array of HybridNetwork objects (that can be trees)
"""
`readInputTrees(file)`

function to read a text file with a list of trees in parenthetical format (one tree per line), it returns an array of HybridNetwork object.
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
        DEBUG && println("$(line)")
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
# warning: taxon has to be vector of ASCIIString, vector of AbstractString do not work
function allQuartets(taxon::Union{Vector{ASCIIString},Vector{Int}}, writeFile::Bool)
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
    if(num == 0)
        num = integer(floor(0.1*n))
    end
    num <= n || error("you cannot choose a sample of $(num) quartets when there are $(n) in total")
    indx = [ones(Int,num); zeros(Int,n-num)]
    indx = indx[sortperm(randn(n))]
    if(writeFile)
        #randName = "rand$(numQ)Quartets$(string(integer(time()/1000))).txt"
        randName = "rand$(numQ)Quartets.txt"
        println("list of randomly selected quartets in file $(randName)")
        out = open(randName,"w")
    end
    for i in 1:n
        if(indx[i] == 1)
            if(writeFile)
                write(out,"$(lines[i])")
            end
            push!(randquartets,allquartets[i])
        end
    end
    if(writeFile)
        close(out)
    end
    return randquartets
end

# function that will not use randQuartets(list of quartets,...)
# this function uses whichQuartet to avoid making the list of all quartets
# fixit: i think we should write always the file, but not sure
function randQuartets(taxon::Union{Vector{ASCIIString},Vector{Int}},num::Integer, writeFile::Bool)
    randquartets = Quartet[]
    n = length(taxon)
    ntotal = binom(n,4)
    num <= ntotal || error("you cannot choose a sample of $(num) quartets when there are $(ntotal) in total")
    # indx = [rep(1,num);rep(0,ntotal-num)] # requires much more memory than necessary:
    # indx = indx[sortperm(randn(ntotal))]  # several arrays of size ntotal !!
    # rq = find(indx .== 1)
    rq = sample(1:ntotal, num, replace=false, ordered=true)
    randName = "rand$(num)Quartets.txt"
    println("list of randomly selected quartets in file $(randName)")
    out = open(randName,"w")
    i = 1
    for(q in rq)
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


# function to check if taxa in quartet is in tree t
function sameTaxa(q::Quartet, t::HybridNetwork)
    for name in q.taxon
        in(name,t.names) || return false
    end
    return true
end

# function to check if taxa in all quartets is in tree t
function sameTaxa(quartets::Vector{Quartet}, t::HybridNetwork)
    suc = true
    for(q in quartets)
        for name in q.taxon
            suc = in(name,t.names) ? true : false
        end
    end
    return suc
end

sameTaxa(d::DataCF, t::HybridNetwork) = sameTaxa(d.quartet, t)


# function to extract the union of taxa of list of gene trees
function unionTaxa(trees::Vector{HybridNetwork})
    taxa = trees[1].names
    for t in trees
        taxa = union(taxa,t.names)
    end
    integers = true
    try
        parse(Int,taxa[1])
    catch
        integers = false
    end
    # to sort taxa
    if(integers)
        taxa = map(x->string(x),sort(map(x->parse(Int,x),taxa)))
    else
        taxa = sort(taxa)
    end
    return taxa
end

# function to extract the union of taxa of list of quartets
function unionTaxa(quartets::Vector{Quartet})
    taxa = quartets[1].taxon
    for q in quartets
        taxa = union(taxa,q.taxon)
    end
    integers = true
    try
        parse(Int,taxa[1])
    catch
        integers = false
    end
    # to sort taxa
    if(integers)
        taxa = map(x->string(x),sort(map(x->parse(Int,x),taxa)))
    else
        taxa = sort(taxa)
    end
    return taxa
end

unionTaxaTree(file::AbstractString) = unionTaxa(readInputTrees(file))


"""
    calculateObsCFAll!(DataCF, taxa::Union{Vector{ASCIIString}, Vector{Int}})

update the .quartet[i].obsCF values of the DataCF object, based on its .tree vector.

    calculateObsCFAll!(vector of quartets, vector of trees, taxa)

update the .obsCF values of the quartets, based on the trees, and returns a new DataCF object
with these updated quartets and trees.

    calculateObsCFAll_noDataCF!(vector of quartets, vector of trees, taxa)

update the .obsCF values of the quartets based on the trees, but returns nothing.
"""
# function to calculate the obsCF from a file with a set of gene trees
# returns a DataCF object (*no* longer writes a csv table with the obsCF)
# warning: it needs trees (not networks) as input
function calculateObsCFAll!(dat::DataCF, taxa::Union{Vector{ASCIIString}, Vector{Int}})
    calculateObsCFAll_noDataCF!(dat.quartet, dat.tree, taxa)
end

function calculateObsCFAll!(quartets::Vector{Quartet}, trees::Vector{HybridNetwork}, taxa::Union{Vector{ASCIIString}, Vector{Int}})
    calculateObsCFAll_noDataCF!(quartets, trees, taxa)
    d = DataCF(quartets,trees)
    return d
end

function calculateObsCFAll_noDataCF!(quartets::Vector{Quartet}, trees::Vector{HybridNetwork}, taxa::Union{Vector{ASCIIString}, Vector{Int}})
    println("calculating obsCF from $(length(trees)) gene trees and for $(length(quartets)) quartets")
    index = 1
    totalq = length(quartets)
    println("Reading in quartets...")
    r = round(1/totalq,2)
    if(r > 0.02)
        numq = totalq
    else
        numq = 50
    end
    print("0+")
    for(i in 1:numq)
        print("-")
    end
    print("+100%")
    println("  ")
    print("  ")
    for q in quartets
        if(round(index/totalq,2)>0.02)
            print("*")
            index = 1
        end
        suma = 0
        sum12 = 0
        sum13 = 0
        sum14 = 0
        for t in trees
            isTree(t) || error("gene tree found in file that is a network $(writeTopology(t))")
            if(sameTaxa(q,t))
                M = tree2Matrix(t,taxa) #fixit: way to reuse M? length(t.edge) will be different across trees
                res = extractQuartetTree(q,M,taxa)
                DEBUG && println("res is $(res)")
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

# function to read input list of gene trees/quartets and calculates obsCF
# as opposed to readTableCF that read the table of obsCF directly
# input: treefile (with gene trees), quartetfile (with list of quartets),
# whichQ (:add/:rand to decide if all or random sample of quartets, default all)
# numQ: number of quartets in random sample
# writetab = true to write the table of obsCF as file with name filename
# does it by default
# writeFile=true writes file with sampled quartets, default false
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
readInputData(treefile::AbstractString, quartetfile::AbstractString) = readInputData(treefile, quartetfile, :all, 0, true, "none", false, true)
readInputData(treefile::AbstractString, quartetfile::AbstractString, writetab::Bool, filename::AbstractString) = readInputData(treefile, quartetfile, :all, 0, writetab, filename, false, true)

# function to read input list of gene trees/quartets and calculates obsCF
# as opposed to readTableCF that read the table of obsCF directly
# input: trees Vector of HybridNetwork, quartetfile (with list of quartets),
# whichQ (:add/:rand to decide if all or random sample of quartets, default all)
# numQ: number of quartets in random sample
# writetab = true to write the table of obsCF as file with name filename
# does it by default
# writeFile=true writes file with sampled quartets, default false
function readInputData(trees::Vector{HybridNetwork}, quartetfile::AbstractString, whichQ::Symbol, numQ::Integer, writetab::Bool, filename::AbstractString, writeFile::Bool, writeSummary::Bool)
    if(whichQ == :all)
        numQ == 0 || warn("set numQ=$(numQ) but whichQ is not rand, so all quartets will be used and numQ will be ignored. If you want a specific number of 4-taxon subsets not random, you can input with the quartetfile option")
        println("will use all quartets in file $(quartetfile)")
        quartets = readListQuartets(quartetfile)
    elseif(whichQ == :rand)
        if(numQ == 0)
            warn("not specified numQ but whichQ=rand, so 10% of quartets will be sampled") #handled inside randQuartets
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
        df = writeObsCF(d)
        writetable(filename,df)
    end
    #descData(d,"summaryTreesQuartets$(string(integer(time()/1000))).txt")
    writeSummary && descData(d,"summaryTreesQuartets.txt")
    return d
end



# function to read input list of gene trees, and not the list of quartets
# so it creates the list of quartets inside and calculates obsCF
# as opposed to readTableCF that read the table of obsCF directly
# input: treefile (with gene trees), whichQ (:add/:rand to decide if all or random sample of quartets, default all)
# numQ: number of quartets in random sample
# taxa: list of taxa, if not given, all taxa in gene trees used
# writetab = true to write the table of obsCF as file with name filename
# does it by default
# writeFile= true, writes intermediate files with the quartets info (default false)
function readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer, taxa::Union{Vector{ASCIIString}, Vector{Int}}, writetab::Bool, filename::AbstractString, writeFile::Bool, writeSummary::Bool)
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

readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer, taxa::Union{Vector{ASCIIString}, Vector{Int}}, writetab::Bool) = readInputData(treefile, whichQ, numQ, taxa, writetab, "none", false, true)
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer, taxa::Union{Vector{ASCIIString}, Vector{Int}}) = readInputData(treefile, whichQ, numQ, taxa, true, "none",false, true)
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer, writetab::Bool, filename::AbstractString) = readInputData(treefile, whichQ, numQ, unionTaxaTree(treefile), writetab, filename,false, true)
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer, writetab::Bool) = readInputData(treefile, whichQ, numQ, unionTaxaTree(treefile), writetab, "none",false, true)
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer) = readInputData(treefile, whichQ, numQ, unionTaxaTree(treefile), true, "none",false, true)
readInputData(treefile::AbstractString) = readInputData(treefile, :all, 0, unionTaxaTree(treefile), true, "none",false, true)
readInputData(treefile::AbstractString,taxa::Union{Vector{ASCIIString}, Vector{Int}}) = readInputData(treefile, :all, 0, taxa, true, "none",false, true)
readInputData(treefile::AbstractString, filename::AbstractString) = readInputData(treefile, :all, 0, unionTaxaTree(treefile), true, filename,false, true)

# function to read input vector of HybridNetworks, and not the list of quartets
# so it creates the list of quartets inside and calculates obsCF
# as opposed to readTableCF that read the table of obsCF directly
# input: trees (Vector of HybridNetwork), whichQ (:add/:rand to decide if all or random sample of quartets, default all)
# numQ: number of quartets in random sample
# taxa: list of taxa, if not given, all taxa in gene trees used
# writetab = true to write the table of obsCF as file with name filename
# does it by default
# writeFile= true, writes intermediate files with the quartets info (default false)
function readInputData(trees::Vector{HybridNetwork}, whichQ::Symbol, numQ::Integer, taxa::Union{Vector{ASCIIString}, Vector{Int}}, writetab::Bool, filename::AbstractString, writeFile::Bool, writeSummary::Bool)
    if(whichQ == :all)
        numQ == 0 || warn("set numQ=$(numQ) but whichQ=all, so all quartets will be used and numQ will be ignored. If you want a specific number of 4-taxon subsets not random, you can input with the quartetfile option")
        quartets = allQuartets(taxa,writeFile)
        println("will use all quartets on $(length(taxa)) taxa")
    elseif(whichQ == :rand)
        if(numQ == 0)
            warn("not specified numQ with whichQ=rand, so 10% of quartets will be sampled") #handled inside randQuartets
        else
            println("will use a random sample of $(numQ) 4-taxon sets ($(round((100*numQ)/binomial(length(taxa),4),2)) percent) on $(length(taxa)) taxa")
        end
        quartets = randQuartets(taxa,numQ, writeFile)
    else
        error("unknown symbol for whichQ $(whichQ), should be either all or rand")
    end
    d = calculateObsCFAll!(quartets,trees,taxa)
    if writetab
        println("table of obsCF printed to file $(filename)")
        df = writeObsCF(d)
        writetable(filename,df)
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
"""
function readTrees2CF(treefile::AbstractString; quartetfile="none"::AbstractString, whichQ="all"::AbstractString, numQ=0::Integer, writeTab=true::Bool, CFfile="none"::AbstractString, taxa=unionTaxaTree(treefile)::Union{Vector{ASCIIString},Vector{Int}}, writeQ=false::Bool, writeSummary=true::Bool)
    whichQ == "all" || whichQ == "rand" ||
        error("whichQ should be all or rand, not $(whichQ)")
    if(quartetfile == "none")
        readInputData(treefile, symbol(whichQ), numQ, taxa, writeTab, CFfile, writeQ, writeSummary)
    else
        readInputData(treefile, quartetfile, symbol(whichQ), numQ, writeTab, CFfile, writeQ, writeSummary)
    end
end

# same as before, but with input vector of HybridNetworks
function readTrees2CF(trees::Vector{HybridNetwork}; quartetfile="none"::AbstractString, whichQ="all"::AbstractString, numQ=0::Integer, writeTab=true::Bool, CFfile="none"::AbstractString, taxa=unionTaxa(trees)::Union{Vector{ASCIIString},Vector{Int}}, writeQ=false::Bool, writeSummary=true::Bool)
    whichQ == "all" || whichQ == "rand" ||
        error("whichQ should be all or rand, not $(whichQ)")
    if(quartetfile == "none")
        readInputData(trees, symbol(whichQ), numQ, taxa, writeTab, CFfile, writeQ, writeSummary)
    else
        readInputData(trees, quartetfile, symbol(whichQ), numQ, writeTab, CFfile, writeQ, writeSummary)
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
        write(s,"Taxon $(taxon) appears in $(numT) input trees ($(round(100*numT/length(trees),2)) %)\n")  #and $(numQ) quartets ($(round(100*numQ/length(quartets),2)) %)\n")
    end
end

taxaTreesQuartets(trees::Vector{HybridNetwork}, quartets::Vector{Quartet}) = taxaTreesQuartets(trees, quartets, STDOUT)

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
# which can be a file or STDOUT
# default: send to STDOUT
# pc: only 4-taxon subsets with percentage of gene trees less than pc will be printed (default 70%)
function descData(d::DataCF, sout::IO, pc::Float64)
    0<=pc<=1 || error("percentage of missing genes should be between 0,1, not: $(pc)")
    if(!isempty(d.tree))
        print(sout,"data consists of $(d.numTrees) gene trees and $(d.numQuartets) 4-taxon subsets\n")
        taxaTreesQuartets(d.tree,d.quartet,sout)
        print(sout,"----------------------------\n\n")
        print(sout,"will print below only the 4-taxon subsets with data from <= $(round((pc)*100,2))% genes\n")
        for q in d.quartet
            percent  = q.ngenes == -1 ? 0.0 : round(q.ngenes/d.numTrees*100,2)
            if(percent < pc)
                print(sout,"4-taxon subset $(q.taxon) obsCF constructed with $(q.ngenes) gene trees ($(percent)%)\n")
            end
        end
        print(sout,"----------------------------\n\n")
    else
        if(!isempty(d.quartet))
            print(sout,"data consists of $(d.numQuartets) 4-taxon subsets")
            taxa=unionTaxa(d.quartet)
            print(sout,"\nTaxa: $(taxa)\n")
            print(sout,"Number of Taxa: $(length(taxa))\n")
            numQ = binomial(length(taxa),4);
            print(sout,"Maximum number of 4-taxon subsets: $(numQ). Thus, $(round(100*d.numQuartets/numQ,2)) percent of 4-taxon subsets sampled\n")
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
descData(d::DataCF) = descData(d, STDOUT,0.7)
descData(d::DataCF,pc::Float64) = descData(d, STDOUT,pc)
descData(d::DataCF, filename::AbstractString) = descData(d, filename,0.7)

"""
`summarizeDataCF(d::DataCF)`

function to summarize the information contained in a DataCF object. It has the following optional arguments:
- filename: if provided, the summary will be saved in the filename, not to screen
- pc (number between (0,1)): threshold of percentage of missing genes to identify 4-taxon subsets with fewer genes than the threshold
"""
function summarizeDataCF(d::DataCF; filename="none"::AbstractString, pc=0.7::Float64)
    0<=pc<=1 || error("percentage of missing genes should be between 0,1, not: $(pc)")
    if(filename == "none")
        descData(d,STDOUT,pc)
    else
        descData(d,filename,pc)
    end
end

# ------------------ read starting tree from astral and put branch lengths ------------------------------

# function to read the starting topology (can be tree/network)
# if updateBL=true, updates the branch lengths with the obsCF in d
# by default, updateBL=true
function readStartTop(file::AbstractString,d::DataCF,updateBL::Bool)
    net = readTopologyUpdate(file)
    if(updateBL)
        updateBL!(net,d)
    end
    return net
end

"""
`readStartTop(treefile,d::DataCF)`

function to read a tree in parenthetical format from text file treefile and a DataCF object to update the branch lengths according to the average observed CF for any given edge.
"""
readStartTop(file::AbstractString,d::DataCF) = readStartTop(file,d,true)
#readStartTop(file::AbstractString) = readStartTop(file,DataCF(),false) #not sure why we need this one

# function to update starting branch lengths for starting tree read from ASTRAL
# BL are updated as -log(3/2(1-mean(obsCF)))
# based on Cecile's chr4-species-tree-units.r
# input: starting tree, data after read table of obsCF
function updateBL!(net::HybridNetwork,d::DataCF)
    if(isTree(net))
        parts = edgesParts(net)
        df = makeTable(net,parts,d)
        x=by(df,[:edge],df->DataFrame(Nquartets=length(df[:CF]),edgeL=-log(3/2*(1-mean(df[:CF])))))
        # ommitting columns: meanCF=mean(df[:CF]), sdCF=std(df[:CF])
        edges = x[:edge]
        lengths = x[:edgeL]
        for(i in 1:length(edges))
            try
                ind = getIndexEdge(edges[i],net)
            catch
                error("edge $(edges[i]) not in net")
            end
            ind = getIndexEdge(edges[i],net)
            if (net.edge[ind].length < 0.0 || net.edge[ind].length==1.0)
                # readTopologyLevel1 changes missing branch length to 1.0
                setLength!(net.edge[ind], (lengths[i] > 0 ? lengths[i] : 0.0))
            end
        end
        for (e in net.edge)
            if e.length < 0.0 # some edges might have *no* quartet in the data
                setLength!(e, 1.0)
            end
        end
        return x
    else
        warn("updateStartBL was created for a tree, and net here is not a tree, so no branch lengths updated")
    end
end


# function to get part1,part2,part3,part4 for each edge in net.edge
# returns a EdgeParts object
function edgesParts(net::HybridNetwork)
    parts = EdgeParts[] #vector to hold part1,...,part4 for each edge
    for(e in net.edge)
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
        for(e in node.edge)
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
    for(p in parts) #go over internal edges too
        for(t1 in p.part1)
            for(t2 in p.part2)
                for(t3 in p.part3)
                    for(t4 in p.part4)
                        tx1 = net.names[t1.number]
                        tx2 = net.names[t2.number]
                        tx3 = net.names[t3.number]
                        tx4 = net.names[t4.number]
                        nam = [tx1,tx2,tx3,tx4]
                        snam = sort(nam)
                        # row = getIndex(true,[sort(nam) == sort(q.taxon) for q in d.quartet])
                        # getIndex was getting the first index only: like findfirst
                        row = findin([dnam==snam for dnam in sortedDataQ], true)
                        for (r in row) # nothing if tax set not found: length(row)=0
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
function resolution(names::Vector{ASCIIString},rownames::Vector{ASCIIString})
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
function extractQuartetTree(q::Quartet, M::Matrix{Int},S::Union{Vector{ASCIIString},Vector{Int}})
    DEBUG && println("extractQuartet: $(q.taxon)")
    try
        ind1 = getIndex(q.taxon[1],S)
        ind2 = getIndex(q.taxon[2],S)
        ind3 = getIndex(q.taxon[3],S)
        ind4 = getIndex(q.taxon[4],S)
    catch
        error("some taxon in quartet $(q.taxon) not found in list of all species $(S)")
    end
    ind1 = getIndex(q.taxon[1],S)
    ind2 = getIndex(q.taxon[2],S)
    ind3 = getIndex(q.taxon[3],S)
    ind4 = getIndex(q.taxon[4],S)
    subM = M[:,[ind1+1,ind2+1,ind3+1,ind4+1]]
    DEBUG && println("subM: $(subM)")
    for(r in 1:size(subM,1)) #rows in subM
        if(subM[r,:] == [0 0 1 1] || subM[r,:] == [1 1 0 0])
            return 1
        elseif(subM[r,:] == [0 1 0 1] || subM[r,:] == [1 0 1 0])
            return 2
        elseif(subM[r,:] == [0 1 1 0] || subM[r,:] == [1 0 0 1])
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
function createQuartet(taxa::Union{Vector{ASCIIString},Vector{Int}},qvec::Vector{Int}, num::Int)
    length(qvec) == 4 || error("a quartet should have only 4 taxa, not $(length(qvec))")
    names = ASCIIString[]
    for(i in qvec)
        i <= length(taxa) || error("want taxon number $(i) in list of taxon names $(taxa) which has only $(length(taxa)) names")
        push!(names,string(taxa[i]))
    end
    return Quartet(num,names,[1.0,0.0,0.0])
end



