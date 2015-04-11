# functions to read data
# originally in functions.jl
# Claudia March 2015


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
    df = DataFrames.DataFrame(t1="",t2="",t3="",t4="",CF1234=0.,CF1324=0.,CF1423=0.)
    for(q in quartets)
        length(q.taxon) == 4 || error("quartet $(q.number) does not have 4 taxa")
        length(q.obsCF) == 3 || error("quartet $(q.number) does have qnet with 3 expCF")
        append!(df,DataFrames.DataFrame(t1=q.taxon[1],t2=q.taxon[2],t3=q.taxon[3],t4=q.taxon[4],CF1234=q.obsCF[1],CF1324=q.obsCF[2],CF1423=q.obsCF[3]))
    end
    df = df[2:size(df,1),1:size(df,2)]
    return df
end

writeObsCF(d::DataCF) = writeObsCF(d.quartet)

# function that takes a dataframe and creates a DataCF object
function readTableCF(df::DataFrame)
    warn("assume the numbers for the taxon read from the observed CF table match the numbers given to the taxon when creating the object network")
    size(df,2) == 7 || error("Dataframe should have 7 columns: 4taxa, 3CF")
    quartets = Quartet[]
    for(i in 1:size(df,1))
        push!(quartets,Quartet(i,string(df[i,1]),string(df[i,2]),string(df[i,3]),string(df[i,4]),[df[i,5],df[i,6],df[i,7]]))
    end
    return DataCF(quartets)
end

# ---------------- read input gene trees and calculate obsCF ----------------------

# function to read a file and create one object per line read
# (each line starting with "(" will be considered a topology)
# the file can have extra lines that are ignored
# returns an array of HybridNetwork objects (that can be trees)
function readInputTrees(file::String)
    try
        s = open(file)
    catch
        error("Could not find or open $(file) file");
    end
    vnet = HybridNetwork[];
    s = open(file)
    for line in eachline(s)
        c = line[1]
        if(c == '(')
           net = readTopology(line)
           push!(vnet,net)
        end
    end
    close(s)
    return vnet
end

# function to list all quartets for a set of taxa names
# return a text file with the list of quartets, one per line
function allQuartets(taxon::Union(Vector{ASCIIString},Vector{Int64}))
    quartets = combinations(taxon,4)
    f = open("allQuartets.txt","w")
    for q in quartets
        write(f,"$(q[1]),$(q[2]),$(q[3]),$(q[4])\n")
    end
    close(f)
end

# function to list all quartets for a set of taxa names
# return a text file with the list of quartets, one per line
# instead of taxon names, only needs numTaxa
function allQuartets(numTaxa::Int64)
    quartets = combinations(1:numTaxa,4)
    f = open("allQuartets.txt","w")
    for q in quartets
        write(f,"$(q[1]),$(q[2]),$(q[3]),$(q[4])\n")
    end
    close(f)
end

# function to list num randomly selected quartets for a set of taxa names
# return a text file with the list of quartets, one per line
# input: file with the list of all quartets
function randQuartets(allQ::ASCIIString,num::Int64)
    try
        f = open(allQ)
    catch
        error("could not open file $(allQ)")
    end
    f = open(allQ)
    lines = readlines(f)
    n = length(lines)
    num <= n || error("you cannot choose a sample of $(num) quartets when there are $(n) in total")
    indx = [rep(1,num),rep(0,n-num)]
    indx = indx[sortperm(randn(n))]
    out = open("rand$(num)Quartets.txt","w")
    for i in 1:n
        if(indx[i] == 1)
            write(out,"$(lines[i])")
        end
    end
    close(out)
    close(f)
end

# function to list num randomly selected quartets for a set of taxa names
# return a text file with the list of quartets, one per line
# no input file with all quartets because it will create it
function randQuartets(numTaxa::Int64,num::Int64)
    allQuartets(numTaxa)
    f = open("allQuartets.txt")
    lines = readlines(f)
    n = length(lines)
    num <= n || error("you cannot choose a sample of $(num) quartets when there are $(n) in total")
    indx = [rep(1,num),rep(0,n-num)]
    indx = indx[sortperm(randn(n))]
    out = open("rand$(num)Quartets.txt","w")
    for i in 1:n
        if(indx[i] == 1)
            write(out,"$(lines[i])")
        end
    end
    close(out)
    close(f)
end
# function to list num randomly selected quartets for a set of taxa names
# return a text file with the list of quartets, one per line
# no input file with all quartets because it will create it
function randQuartets(taxon::Union(Vector{ASCIIString},Vector{Int64}),num::Int64)
    allQuartets(taxon)
    f = open("allQuartets.txt")
    lines = readlines(f)
    n = length(lines)
    num <= n || error("you cannot choose a sample of $(num) quartets when there are $(n) in total")
    indx = [rep(1,num),rep(0,n-num)]
    indx = indx[sortperm(randn(n))]
    out = open("rand$(num)Quartets.txt","w")
    for i in 1:n
        if(indx[i] == 1)
            write(out,"$(lines[i])")
        end
    end
    close(out)
    close(f)
end


# function to read list of quartets from a file
# and create Quartet type objects
function readListQuartets(file::ASCIIString)
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
