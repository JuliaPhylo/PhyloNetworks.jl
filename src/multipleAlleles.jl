# functions to extend snaq to multiple alleles case
# Claudia November 2015


# function to read a table of allele-species matchs (dataframe with 2 columns)
# and a table of CF in the allele names, and replace all the allele names
# to the species names
# this will create a new CF table, will not rewrite on the original one
# filename is the name to give to the new table, if write=true
function mapAllelesCFtable!(alleleDF::DataFrame, cfDF::DataFrame,write::Bool,filename::AbstractString)
    compareTaxaNames(alleleDF,cfDF)
    newt1 = map(x->replace(string(x),string(alleleDF[1,:allele]),alleleDF[1,:species]),cfDF[1])
    newt2 = map(x->replace(string(x),string(alleleDF[1,:allele]),alleleDF[1,:species]),cfDF[2])
    newt3 = map(x->replace(string(x),string(alleleDF[1,:allele]),alleleDF[1,:species]),cfDF[3])
    newt4 = map(x->replace(string(x),string(alleleDF[1,:allele]),alleleDF[1,:species]),cfDF[4])
    if(size(alleleDF,1) > 1)
        for(j in 2:size(alleleDF,1)) #for all the allele matchings
            newt1 = map(x->replace(string(x),string(alleleDF[j,:allele]),alleleDF[j,:species]),newt1)
            newt2 = map(x->replace(string(x),string(alleleDF[j,:allele]),alleleDF[j,:species]),newt2)
            newt3 = map(x->replace(string(x),string(alleleDF[j,:allele]),alleleDF[j,:species]),newt3)
            newt4 = map(x->replace(string(x),string(alleleDF[j,:allele]),alleleDF[j,:species]),newt4)
        end
    end
    newdf = DataFrames.DataFrame(t1=newt1,t2=newt2,t3=newt3,t4=newt4,CF12_34=cfDF[5],CF13_24=cfDF[6],CF14_23=cfDF[7])
    if(write)
        filename != "" || error("want to write new table of CF with alleles mapped but filename is empty")
        writetable(filename,newdf)
    end
    return newdf
end

# function to clean a df after changing allele names to species names
# inside mapAllelesCFtable
# by deleting rows that are not informative like sp1 sp1 sp1 sp2
function cleanAlleleDF!(newdf::DataFrame)
    global DEBUG
    keeprows = ones(Bool,size(newdf,1)) # repeats 'true'
    ##keeprows =  Bool[]
    repSpecies = ASCIIString[]
    if(isa(newdf[1,1],Int)) #taxon names as integers: we need this to be able to add _2
        newdf[:1] = map(x->string(x),newdf[:1])
        newdf[:2] = map(x->string(x),newdf[:2])
        newdf[:3] = map(x->string(x),newdf[:3])
        newdf[:4] = map(x->string(x),newdf[:4])
    end
    for(i in 1:size(newdf,1)) #check all rows
        DEBUG && println("row number: $i")
        row = convert(Array,DataArray(newdf[i,1:4]))
        DEBUG && println("row $(row)")
        uniq = unique(row)
        DEBUG && println("unique $(uniq)")
        if(length(uniq) == 1)
            ##push!(keeprows,false)
            keeprows[i]= false
        elseif(length(uniq) == 4)
            ##push!(keeprows,true)
            continue
        elseif(length(uniq) == 3) #sp1 sp1 sp2 sp3
            for(u in uniq)
                DEBUG && println("u $(u), typeof $(typeof(u))")
                ind = row .== u #taxon names matching u
                DEBUG && println("taxon names matching u $(ind)")
                if(sum(ind) == 2)
                    push!(repSpecies,string(u))
                    found = false
                    for(k in 1:4)
                        if(ind[k])
                            if(found)
                                DEBUG && println("found the second one in k $(k), will change newdf[i,k] $(newdf[i,k]), typeof $(typeof(newdf[i,k]))")
                                newdf[i,k] = string(u,"_2")
                                break
                            else
                                found = true
                            end
                        end
                    end
                    break
                end
            end
        elseif(length(uniq) == 2)
            keep = true
            for(u in uniq)
                DEBUG && println("length uniq is 2, u $(u)")
                ind = row .== u
                if(sum(ind) == 1 || sum(ind) == 3)
                    DEBUG && println("ind $(ind) is 1 or 3, should not keep")
                    keep = false
                    break
                elseif(sum(ind) == 2)
                    DEBUG && println("ind $(ind) is 2, should keep")
                    keep = true
                    found = false
                    push!(repSpecies,string(u))
                    for(k in 1:4)
                        if(ind[k])
                            if(found)
                                newdf[i,k] = string(u,"_2")
                                break
                            else
                                found = true
#                                newdf[i,k] = u
                            end
                        end
                    end
                end
            end
            DEBUG && println("after if, keep is $(keep)")
            ##push!(keeprows,keep)
            keeprows[i] = keep
        end
        DEBUG && println("$(keeprows[i])")
    end
    DEBUG && println("keeprows is $(keeprows)")
    DEBUG && println("repSpecies is $(repSpecies)")
    if(!all(keeprows))
        warn("found $(length(keeprows)-sum(keeprows)) troublesome 4-taxon subsets out of $(size(newdf,1)) 4-taxon subsets. for the moment, we will ignore these 4-taxon subsets: will use $(sum(keeprows)) 4-taxon subsets.")
        size(newdf,1) > (length(keeprows)-sum(keeprows)) || warn("4-taxon subsets with repeated taxon names are all the 4-taxon subsets, so resulting new dataframe is empty")
    end
    return unique(repSpecies),keeprows
end


# function to merge rows that have repeated taxon names by using the weigthed average of CF
# (if info on number of genes is provided) or simple average
function mergeRows!(df::DataFrame)
    keeprow = ones(Bool,size(df,1)) # repeats 'true'
    addeddf = df[1,1:size(df,2)] #df that will have the averaged rows, need to delete the first row at the end
    for(i in 1:size(df,1)) #rows
        if(keeprow[i])
            partialdf = df[i,1:size(df,2)]
            for(j in (i+1):size(df,1)) #other rows
                if(df[i,1:4] == df[j,1:4])
                    keeprow[i] = false
                    append!(partialdf,df[j,1:size(df,2)])
                    keeprow[j] = false
                end
            end
            averagePartialDF!(partialdf) # transform all rows in partialdf to only one row with average
            append!(addeddf,partialdf)
        end
    end
    if(!all(keeprow))
        warn("found repeated 4-taxon subsets in CF table, the CF values were averaged between rows")
        df = df[keeprow,:]
        addeddf = addeddf[2:size(addeddf,1),1:size(addeddf,2)]
        append!(df,addeddf) #append averaged rows
    end
end

# function to average all rows in a dataframe
# assumes all rows have repeated taxon names, and does not
# verify this
function averagePartialDF!(df::DataFrame)
    global DEBUG
    df2 = deepcopy(df)
    deleterows!(df,2:size(df,1)) #only keep one row to hold the average
    numfound = true
    try
        df[:ngenes] #column with number of genes
    catch
        numfound = false
    end
    if(numfound)
        DEBUG && println("knows it has ngenes column")
        num = df2[:ngenes]
        suma = sum(num)
        delete!(df2,1:4)
        delete!(df2,:ngenes) #delete all columns except CF
        for(j in 1:size(df2,2)) #columns in df2
            df[1,j+4] = sum(df2[:,j].*num)/suma
        end
        df[1,:ngenes] = suma
    else
        DEBUG && println("knows it does not have ngenes column")
        delete!(df2,1:4)
        for(j in 1:size(df2,2)) #columns in df2
            df[1,j+4] = mean(df2[:,j])
        end
    end
end



# function to expand leaves in tree to two individuals
# based on cf table with alleles mapped to species names
function expandLeaves!(repSpecies::Union{Vector{ASCIIString},Vector{Int64}},tree::HybridNetwork)
    for(sp in repSpecies)
        for(n in tree.node)
            if(n.name == sp) #found leaf with sp name
                n.leaf || error("name $(sp) should correspond to a leaf, but it corresponds to an internal node")
                length(n.edge) == 1 || error("leaf $(sp) should have only one edge attached and it has $(length(n.edge))")
                if(n.edge[1].length == -1.0)
                    setLength!(n.edge[1],1.0)
                end
                removeLeaf!(tree,n)
                n.leaf = false
                n.edge[1].istIdentifiable = true
                n.name = ""
                max_node = maximum([e.number for e in tree.node]);
                max_edge = maximum([e.number for e in tree.edge]);
                e1 = Edge(max_edge+1,0.0)
                e2 = Edge(max_edge+2,0.0)
                n1 = Node(max_node+1,true,false,[e1])
                n2 = Node(max_node+2,true,false,[e2])
                setNode!(e1,n1)
                setNode!(e1,n)
                setNode!(e2,n2)
                setNode!(e2,n)
                setEdge!(n,e1)
                setEdge!(n,e2)
                pushNode!(tree,n1)
                pushNode!(tree,n2)
                pushEdge!(tree,e1)
                pushEdge!(tree,e2)
                n1.name = string(sp)
                n2.name = string(sp,"_2")
                break
            end
        end
    end
end


"""
`mapAllelesCFtable(mapping file, cf file)`

function that change the allele names in the CF table to species names.
The new DataFrame object is returned.
Optional argument: filename for the resulting CF table. If not specified, then no CF is saved as file.
"""
function mapAllelesCFtable(alleleDF::AbstractString, cfDF::AbstractString; filename=""::AbstractString)
    d = readtable(alleleDF)
    d2 = readtable(cfDF)
    if(filename=="")
        mapAllelesCFtable!(d,d2,false,filename)
    else
        mapAllelesCFtable!(d,d2,true,filename)
    end
end

# function to compare the taxon names in the allele-species matching table
# and the CF table
function compareTaxaNames(alleleDF::DataFrame, cfDF::DataFrame)
    checkMapDF(alleleDF)
    size(cfDF,2) == 7 || warn("CF Dataframe should have 7 columns: 4taxa, 3CF, will ignore columns from 8th on")
    d = readTableCF(cfDF)
    println("ALLELE MAP: there are $(length(alleleDF[1])) allele-species matches")
    CFtaxa = unionTaxa(d.quartet)
    CFtaxa = map(x->string(x),CFtaxa) #treat as string
    alleleTaxa = map(x->string(x),alleleDF[:allele]) #treat as string
    sizeCF = length(CFtaxa)
    sizeAllele = length(alleleTaxa)
    if(sizeAllele > sizeCF)
        println("ALLELE MAP: there are more taxa in the allele-species mapping file: $(sizeAllele) than in the CF table: $(sizeCF), so extra allele names will be ignored")
        alleleTaxa = intersect(alleleTaxa,CFtaxa)
    elseif(sizeAllele < sizeCF)
        println("ALLELE MAP: there are fewer taxa in the allele-species map: $(sizeAllele) than in the CF table: $(sizeCF), so some names in the CF table will remained unchanged")
    end
    unchanged = setdiff(CFtaxa,alleleTaxa)
    if(length(unchanged) == length(CFtaxa))
        warn("no allele names in CF table match with the mapping file")
    end
    if(isempty(unchanged))
        println("TAXON NAMES MATCH: taxon names in the CF table were changed according to the allele-species mapping file")
    else
        warn("not all alleles mapped")
        println("the following taxa in the CF table were not modified by the allele-species map (since they are absent in the mapping file):\n $(unchanged)")
    end
end

# function to check that the allele df has one column labelled alleles and one column labelled species
function checkMapDF(alleleDF::DataFrame)
    size(alleleDF,2) <= 2 || error("Allele-Species matching Dataframe should have at least 2 columns")
    size(alleleDF,2) >= 2 || warn("allele mapping file contains more than two columns: will ignore all columns not labelled allele or species")
    try
        alleleDF[:allele]
    catch
        error("In allele mapping file there is no column named allele")
    end
    try
        alleleDF[:species]
    catch
        error("In allele mapping file there is no column named species")
    end
end

