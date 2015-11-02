# functions to read data
# originally in functions.jl
# Claudia March 2015

# same as readData.jl, but modified the Union(...) to Union{...}
# based on Julia 0.4
# Claudia October 2015

# ----- read data --------

# function to write a csv table from the expCF of an
# array of quartets
# warning: does not check if the expCF have been calculated
function writeExpCF(quartets::Array{Quartet,1})
    df = DataFrames.DataFrame(t1="",t2="",t3="",t4="",CF1234=0.,CF1324=0.,CF1423=0.)
    for(q in quartets)
        length(q.taxon) == 4 || error("quartet $(q.number) does not have 4 taxa")
        length(q.qnet.expCF) == 3 || error("quartet $(q.number) does have qnet with 3 expCF")
        append!(df,DataFrames.DataFrame(t1=q.taxon[1],t2=q.taxon[2],t3=q.taxon[3],t4=q.taxon[4],CF1234=q.qnet.expCF[1],CF1324=q.qnet.expCF[2],CF1423=q.qnet.expCF[3]))
    end
    df = df[2:size(df,1),1:size(df,2)]
    return df
end

writeExpCF(d::DataCF) = writeExpCF(d.quartet)

# function to write a csv table from the obsCF of an
# array of quartets
function writeObsCF(quartets::Array{Quartet,1})
    df = DataFrames.DataFrame(t1="",t2="",t3="",t4="",CF1234=0.,CF1324=0.,CF1423=0.,numGT=0)
    for(q in quartets)
        length(q.taxon) == 4 || error("quartet $(q.number) does not have 4 taxa")
        length(q.obsCF) == 3 || error("quartet $(q.number) does have qnet with 3 expCF")
        append!(df,DataFrames.DataFrame(t1=q.taxon[1],t2=q.taxon[2],t3=q.taxon[3],t4=q.taxon[4],CF1234=q.obsCF[1],CF1324=q.obsCF[2],CF1423=q.obsCF[3],numGT=q.numGT))
    end
    df = df[2:size(df,1),1:size(df,2)]
    return df
end

writeObsCF(d::DataCF) = writeObsCF(d.quartet)

# function that takes a dataframe and creates a DataCF object
function readTableCF(df::DataFrames.DataFrame,writeTab::Bool)
    DEBUG && println("assume the numbers for the taxon read from the observed CF table match the numbers given to the taxon when creating the object network")
    size(df,2) == 7 || warn("Dataframe should have 7 columns: 4taxa, 3CF, will ignore columns from 8th on")
    quartets = Quartet[]
    for(i in 1:size(df,1))
        push!(quartets,Quartet(i,string(df[i,1]),string(df[i,2]),string(df[i,3]),string(df[i,4]),[df[i,5],df[i,6],df[i,7]]))
    end
    d = DataCF(quartets)
    println("DATA: data consists of $(d.numTrees) gene trees and $(d.numQuartets) quartets")
    #descData(d,"summaryCFtable$(string(integer(time()/1000))).txt")
    if(writeTab)
        descData(d,"summaryCFtable.txt")
    end
    return d
end

readTableCF(df::DataFrames.DataFrame) = readTableCF(df,true)

# warning: when file needs to be AbstractString bc it can be read as UTF8String
"""
`readTableCF(file)`

read a file with a table of CF. It has one optional argument: sep to specify the type of separator in the table with single quotes: sep=';'
"""
readTableCF(file::AbstractString;sep=','::Char) = readTableCF(readtable(file,separator=sep))

# ---------------- read input gene trees and calculate obsCF ----------------------

# function to read a file and create one object per line read
# (each line starting with "(" will be considered a topology)
# the file can have extra lines that are ignored
# returns an array of HybridNetwork objects (that can be trees)
"""
`readInputTrees(file)`

function to read a text file with a list of trees in parenthetical format (one tree per line), it returns an array of HybridNetwork object.
Careful to put ; after to avoid output written to screen
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
        DEBUG && println("$(line)")
        c = line[1]
        if(c == '(')
           try
               net = readTopologyUpdate(line,false)
               push!(vnet,deepcopy(net))
           catch(err)
               error("could not read tree in line $(numl). The error is $(err)")
           end
        end
        numl += 1
    end
    close(s)
    !isempty(vnet) || return nothing
    return vnet
end

# function to list all quartets for a set of taxa names
# return a vector of quartet objects, and if writeFile=true, writes a file
# warning: taxon has to be vector of ASCIIString, vector of AbstractString do not work
function allQuartets(taxon::Union{Vector{ASCIIString},Vector{Int64}}, writeFile::Bool)
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
        i += 1
    end
    if(writeFile)
        close(f)
    end
    return vquartet
end

allQuartets(numTaxa::Int64, writeFile::Bool) = allQuartets(1:numTaxa, writeFile)

# function to list num randomly selected quartets for a vector of all quartets
# return a vector of Quartet randomly chosen (and a file if writeFile=true)
function randQuartets(allquartets::Vector{Quartet},num::Int64, writeFile::Bool)
    randquartets = Quartet[]
    n = length(allquartets)
    if(num == 0)
        num = integer(floor(0.1*n))
    end
    num <= n || error("you cannot choose a sample of $(num) quartets when there are $(n) in total")
    indx = [rep(1,num);rep(0,n-num)]
    indx = indx[sortperm(randn(n))]
    if(writeFile)
        #randName = "rand$(numQ)Quartets$(string(integer(time()/1000))).txt"
        randName = "rand$(numQ)Quartets.txt"
        println("DATA: chosen list of random quartets in file $(randName)")
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
function randQuartets(taxon::Union{Vector{ASCIIString},Vector{Int64}},num::Int64, writeFile::Bool)
    randquartets = Quartet[]
    n = length(taxon)
    ntotal = binom(n,4)
    num <= ntotal || error("you cannot choose a sample of $(num) quartets when there are $(ntotal) in total")
    indx = [rep(1,num);rep(0,ntotal-num)]
    indx = indx[sortperm(randn(ntotal))]
    rq = find(indx .== 1)
    randName = "rand$(num)Quartets.txt"
    println("DATA: chosen list of random quartets in file $(randName)")
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

randQuartets(numTaxa::Int64,num::Int64, writeFile::Bool) = randQuartets(1:numTaxa,num, writeFile)

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
    return taxa
end

# function to extract the union of taxa of list of quartets
function unionTaxa(quartets::Vector{Quartet})
    taxa = quartets[1].taxon
    for q in quartets
        taxa = union(taxa,q.taxon)
    end
    return taxa
end

unionTaxaTree(file::AbstractString) = unionTaxa(readInputTrees(file))


# function to calculate the obsCF from a file with a set of gene trees
# returns a DataCF object and write a csv table with the obsCF
# warning: it needs trees (not networks) as input
function calculateObsCFAll!(quartets::Vector{Quartet}, trees::Vector{HybridNetwork}, taxa::Union{Vector{ASCIIString}, Vector{Int64}})
    println("DATA: calculating obsCF from set of $(length(trees)) gene trees and list of $(length(quartets)) quartets")
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
#    println("0    10   20   30   40   50   60   70   80   90   100")
#    println("+----+----+----+----+----+----+----+----+----+----+")
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
        q.numGT = suma
        index += 1
    end
    println("  ")
    d = DataCF(quartets,trees)
    return d
end

# function to read input list of gene trees/quartets and calculates obsCF
# as opposed to readTableCF that read the table of obsCF directly
# input: treefile (with gene trees), quartetfile (with list of quartets),
# whichQ (:add/:rand to decide if all or random sample of quartets, default all)
# numQ: number of quartets in random sample
# writetab = true to write the table of obsCF as file with name filename
# does it by default
# writeFile=true writes file with sampled quartets, default false
function readInputData(treefile::AbstractString, quartetfile::AbstractString, whichQ::Symbol, numQ::Int64, writetab::Bool, filename::AbstractString, writeFile::Bool)
    println("DATA: reading input data for treefile $(treefile) \nand quartetfile $(quartetfile)")
    trees = readInputTrees(treefile)
    if(whichQ == :all)
        numQ == 0 || warn("set numQ=$(numQ) but whichQ is not rand, so all quartets will be used and numQ will be ignored. If you want a specific number of 4-taxon subsets not random, you can input with the quartetfile option")
        println("DATA: will use all quartets in quartetfile $(quartetfile)")
        quartets = readListQuartets(quartetfile)
    elseif(whichQ == :rand)
        if(numQ == 0)
            warn("not specified numQ but whichQ=rand, so 10% of quartets will be sampled") #handled inside randQuartets
        else
            println("DATA: will use a random sample of $(numQ) from quartetfile $(quartetfile)")
        end
        allquartets = readListQuartets(quartetfile)
        quartets = randQuartets(allquartets,numQ,writeFile)
    else
        error("unknown symbol for whichQ $(whichQ), should be either all or rand")
    end
    d = calculateObsCFAll!(quartets,trees, unionTaxa(trees))
    if(writetab)
        if(filename == "none")
            filename = "tableCF$(string(integer(time()/1000))).txt"
        end
        println("\nDATA: printing table of obsCF in file $(filename)")
        df = writeObsCF(d)
        writetable(filename,df)
    end
    #descData(d,"summaryTreesQuartets$(string(integer(time()/1000))).txt")
    descData(d,"summaryTreesQuartets.txt")
    return d
end

readInputData(treefile::AbstractString, quartetfile::AbstractString, whichQ::Symbol, numQ::Int64, writetab::Bool) = readInputData(treefile, quartetfile, whichQ, numQ, writetab, "none", false)
readInputData(treefile::AbstractString, quartetfile::AbstractString, whichQ::Symbol, numQ::Int64) = readInputData(treefile, quartetfile, whichQ, numQ, true, "none", false)
readInputData(treefile::AbstractString, quartetfile::AbstractString) = readInputData(treefile, quartetfile, :all, 0, true, "none", false)
readInputData(treefile::AbstractString, quartetfile::AbstractString, writetab::Bool, filename::AbstractString) = readInputData(treefile, quartetfile, :all, 0, writetab, filename, false)

# function to read input list of gene trees, and not the list of quartets
# so it creates the list of quartets inside and calculates obsCF
# as opposed to readTableCF that read the table of obsCF directly
# input: treefile (with gene trees), whichQ (:add/:rand to decide if all or random sample of quartets, default all)
# numQ: number of quartets in random sample
# taxa: list of taxa, if not given, all taxa in gene trees used
# writetab = true to write the table of obsCF as file with name filename
# does it by default
# writeFile= true, writes intermediate files with the quartets info (default false)
function readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Int64, taxa::Union{Vector{ASCIIString}, Vector{Int64}}, writetab::Bool, filename::AbstractString, writeFile::Bool)
    println("DATA: reading input data for treefile $(treefile) and no quartetfile given: will get quartets here")
    trees = readInputTrees(treefile)
    if(whichQ == :all)
        numQ == 0 || warn("set numQ=$(numQ) but whichQ=all, so all quartets will be used and numQ will be ignored. If you want a specific number of 4-taxon subsets not random, you can input with the quartetfile option")
        quartets = allQuartets(taxa,writeFile)
        println("DATA: will use all quartets based on $(length(taxa)) taxa")
    elseif(whichQ == :rand)
        if(numQ == 0)
            warn("not specified numQ with whichQ=rand, so 10% of quartets will be sampled") #handled inside randQuartets
        else
            println("DATA: will use a random sample of $(numQ) quartets ($(round((100*numQ)/binomial(length(taxa),4),2)) percent) based on $(length(taxa)) taxa")
        end
        quartets = randQuartets(taxa,numQ, writeFile)
    else
        error("unknown symbol for whichQ $(whichQ), should be either all or rand")
    end
    d = calculateObsCFAll!(quartets,trees,taxa)
    if(writetab)
        if(filename == "none")
            #filename = "tableCF$(string(integer(time()/1000))).txt"
            filename = "tableCF.txt"
        end
        println("DATA: printing table of obsCF in file $(filename)")
        df = writeObsCF(d)
        writetable(filename,df)
    end
    #descData(d,"summaryTreesQuartets$(string(integer(time()/1000))).txt")
    descData(d,"summaryTreesQuartets.txt")
    return d
end

readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Int64, taxa::Union{Vector{ASCIIString}, Vector{Int64}}, writetab::Bool) = readInputData(treefile, whichQ, numQ, taxa, writetab, "none", false)
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Int64, taxa::Union{Vector{ASCIIString}, Vector{Int64}}) = readInputData(treefile, whichQ, numQ, taxa, true, "none",false)
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Int64, writetab::Bool, filename::AbstractString) = readInputData(treefile, whichQ, numQ, unionTaxaTree(treefile), writetab, filename,false)
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Int64, writetab::Bool) = readInputData(treefile, whichQ, numQ, unionTaxaTree(treefile), writetab, "none",false)
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Int64) = readInputData(treefile, whichQ, numQ, unionTaxaTree(treefile), true, "none",false)
readInputData(treefile::AbstractString) = readInputData(treefile, :all, 0, unionTaxaTree(treefile), true, "none",false)
readInputData(treefile::AbstractString,taxa::Union{Vector{ASCIIString}, Vector{Int64}}) = readInputData(treefile, :all, 0, taxa, true, "none",false)
readInputData(treefile::AbstractString, filename::AbstractString) = readInputData(treefile, :all, 0, unionTaxaTree(treefile), true, filename,false)


# rename the function readInputData to make it more user-friendly
"""
`readTrees2CF(treefile)`

function to read the trees in parenthetical format from treefile (text file) and calculate the observed CF. It has many optional arguments:

- quartetfile: name of text file with list of 4-taxon subsets to be analyzed. If none is specified, the function will list all possible 4-taxon subsets.
- whichQ="rand": to choose a random sample of 4-taxon subsets
- numQ: size of random sample (ignored if whichQ is not set to "rand")
- writeTab=false: does not write the observedCF to a table (default true)
- CFfile: name of file to save the observedCF (default tableCF.txt)
- writeFile=true: save intermediate files with the list of all 4-taxon subsets and chosen random sample (default false).
"""
function readTrees2CF(treefile::AbstractString; quartetfile="none"::AbstractString, whichQ="all"::AbstractString, numQ=0::Int64, writetab=true::Bool, CFfile="none"::AbstractString, taxa=unionTaxaTree(treefile)::Union{Vector{ASCIIString},Vector{Int64}}, writeFile=false::Bool)
    if(quartetfile == "none")
        if(whichQ == "all")
            readInputData(treefile, :all, numQ, taxa, writetab, CFfile, writeFile)
        elseif(whichQ == "rand")
            readInputData(treefile, :rand, numQ, taxa, writetab, CFfile, writeFile)
        else
            error("whichQ should be all or rand, not $(whichQ)")
        end
    else
        if(whichQ == "all")
            readInputData(treefile, quartetfile, :all, numQ, writetab, CFfile, writeFile)
        elseif(whichQ == "rand")
            readInputData(treefile, quartetfile, :rand, numQ, writetab, CFfile, writeFile)
        else
            error("whichQ should be all or rand, not $(whichQ)")
        end
    end
end

# ---------------------- descriptive stat for input data ----------------------------------

# function to check how taxa is represented in the input trees
function taxaTreesQuartets(trees::Vector{HybridNetwork}, quartets::Vector{Quartet},s::IO)
    taxaT = unionTaxa(trees)
    taxaQ = unionTaxa(quartets)
    dif = symdiff(taxaT,taxaQ)
    isempty(dif) ? write(s,"\nDATA: same taxa in gene trees and quartets: $(taxaT)\n") : write(s,"\nDATA: $(length(dif)) different taxa found in gene trees and quartets. \n Taxa $(intersect(taxaT,dif)) in trees, not in quartets; and taxa $(intersect(taxaQ,dif)) in quartets, not in trees\n")
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
        print(sout,"DATA: data consists of $(d.numTrees) gene trees and $(d.numQuartets) 4-taxon subsets\n")
        taxaTreesQuartets(d.tree,d.quartet,sout)
        print(sout,"----------------------------\n\n")
        print(sout,"will print below only the 4-taxon subsets with data from <= $(round((pc)*100,2))% genes\n")
        for q in d.quartet
            percent  = round(q.numGT/d.numTrees*100,2)
            if(percent < pc)
                print(sout,"4-taxon subset $(q.taxon) obsCF constructed with $(q.numGT) gene trees ($(percent)%)\n")
            end
        end
        print(sout,"----------------------------\n\n")
    else
        if(!isempty(d.quartet))
            print(sout,"DATA: data consists of $(d.numQuartets) 4-taxon subsets")
            taxa=unionTaxa(d.quartet)
            print(sout,"\nTaxa: $(taxa)\n")
            print(sout,"Number of Taxa: $(length(taxa))\n")
            numQ = binomial(length(taxa),4);
            print(sout,"Maximum number of 4-taxon subsets: $(numQ). Thus, $(round(100*d.numQuartets/numQ,2)) percent of 4-taxon subsets sampled\n")
        end
    end
end

function descData(d::DataCF, filename::AbstractString,pc::Float64)
    println("DATA: printing descriptive stat of input data in file $(filename)")
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
    isTree(net) || warn("updateStartBL was created for a tree, and net here is not a tree")
    parts = edgesParts(net)
    df = makeTable(net,parts,d)
    x=by(df,[:edge],df->DataFrame(meanCF=mean(df[:CF]),sdCF=std(df[:CF]),Nquartets=length(df[:CF]),edgeL=-log(3/2*(1-mean(df[:CF])))))
    edges = x[1]
    lengths = x[5]
    for(i in 1:length(edges))
        try
            ind = getIndexEdge(edges[i],net)
        catch
            error("edge $(edges[i]) not in net")
        end
        ind = getIndexEdge(edges[i],net)
        if(lengths[i] > 0)
            setLength!(net.edge[ind],lengths[i])
        else
            setLength!(net.edge[ind],0.0)
        end
    end
    return x
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
    df = DataFrames.DataFrame(edge=1,t1="",t2="",t3="",t4="",resolution="",CF=0.)
    for(p in parts) #go over internal edges too
        for(t1 in p.part1)
            for(t2 in p.part2)
                for(t3 in p.part3)
                    for(t4 in p.part4)
                        tx1 = net.names[t1.number]
                        tx2 = net.names[t2.number]
                        tx3 = net.names[t3.number]
                        tx4 = net.names[t4.number]
                        names = [tx1,tx2,tx3,tx4]
                        row = getIndex(true,[sort(names) == sort(q.taxon) for q in d.quartet])
                        col,res = resolution(names,d.quartet[row].taxon)
                        append!(df,DataFrames.DataFrame(edge=p.edgenum,t1=tx1,t2=tx2,t3=tx3,t4=tx4,resolution=res,CF=d.quartet[row].obsCF[col]))
                    end
                end
            end
        end
    end
    df = df[2:size(df,1),1:size(df,2)]
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

# function to traverse a tree postorder and modify the matrix M
# with edges as rows and species as columns (see tree2Matrix)
# S should be sorted
function traverseTree2Matrix!(node::Node, edge::Edge, M::Matrix{Int}, S::Union{Vector{ASCIIString},Vector{Int64}})
    if(!node.leaf)
        for(e in node.edge) #postorder traversal
            if(!isEqual(e,edge))
                traverseTree2Matrix!(getOtherNode(e,node),e,M,S)
            end
        end
    end
    try
        indedge = getIndex(edge.number,M[:,1])
    catch
        error("either edge $(edge.number) not in tree")
    end
    indedge = getIndex(edge.number,M[:,1])
    if(node.leaf)
        try
            indsp = getIndex(node.name,S)
        catch
            error("leaf $(node.name) not in species list $(S)")
        end
        indsp = getIndex(node.name,S)
        M[indedge,indsp+1] = 1 #indsp+1 bc first column is edge numbers
    else
        for(e in node.edge)
            if(!isEqual(e,edge))
                inde = getIndex(e.number,M[:,1])
                M[indedge,2:size(M,2)] |= M[inde,2:size(M,2)]
            end
        end
    end
end

# function to take a tree and a list of species as input,
# and produce a matrix with edges as rows and species as columns, and in cells:
# Mij=1 if species j is descendant of edge i, 0 ow.
# Mij=-1 if species not present in tree
# if tree has missing taxa, they will keep 0, but this is handled in calculateObsCFAll with sameTaxa function
function tree2Matrix(T::HybridNetwork, S::Union{Vector{ASCIIString},Vector{Int64}})
    sort!(S)
    M = zeros(Int,length(T.edge),length(S)+1)
    M[:,1] = sort!([e.number for e in T.edge])
    for(e in T.node[T.root].edge)
        traverseTree2Matrix!(getOtherNode(e,T.node[T.root]),e,M,S)
    end
    return M
end


# function to extract a quartet from a matrix M
# obtained from tree2Matrix
# this function is meant to replace extractQuartet! in calculateObsCFAll
# input: Quartet, Matrix, vector of taxa names
# returns 1 if quartet found is 12|34, 2 if 13|24, 3 if 14|23, and 0 if not found
function extractQuartetTree(q::Quartet, M::Matrix{Int},S::Union{Vector{ASCIIString},Vector{Int64}})
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
function createQuartet(taxa::Union{Vector{ASCIIString},Vector{Int64}},qvec::Vector{Int}, num::Int)
    length(qvec) == 4 || error("a quartet should have only 4 taxa, not $(length(qvec))")
    names = ASCIIString[]
    for(i in qvec)
        i <= length(taxa) || error("want taxon number $(i) in list of taxon names $(taxa) which has only $(length(taxa)) names")
        push!(names,string(taxa[i]))
    end
    return Quartet(num,names,[1.0,0.0,0.0])
end



