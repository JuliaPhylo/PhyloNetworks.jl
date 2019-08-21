# functions to extend snaq to multiple alleles case
# Claudia November 2015

global repeatAlleleSuffix = "__2"
repeatAlleleSuffix_re = Regex("$repeatAlleleSuffix\$")


"""
    mapAllelesCFtable(mapping file, CF file; filename, columns, delim)

Create a new DataFrame containing the same concordance factors as in the input CF file,
but with modified taxon names. Each allele name in the input CF table is replaced by the
species name that the allele maps onto, based on the mapping file. The mapping file should have column names: allele and species.

Optional arguments:

- file name to write/save resulting CF table. If not specified, then the output
  data frame is not saved to a file.
- column numbers for the taxon names. 1-4 by default.
- any keyword arguments that `CSV.read` would accept.
  For example, delim=',' by default: columns are delimited by commas.
  Unless specified otherwise by the user, `categorical`=false
  (to read taxon names as Strings, not levels of a categorical factor,
  for combining the 4 columns with taxon names more easily).
  The same CSV arguments are used to read both input file (mapping file and quartet file)

See also [`mapAllelesCFtable!`](@ref) to input DataFrames instead of file names.

If a `filename` is specified, such as "quartetCF_speciesNames.csv"
in the example below, this file is best read later with the option
`categorical=false`. example:

```julia
mapAllelesCFtable("allele-species-map.csv", "allele-quartet-CF.csv";
                  filename = "quartetCF_speciesNames.csv")
df_sp = CSV.read("quartetCF_speciesNames.csv", categorical=false); # DataFrame object
dataCF_specieslevel = readTableCF!(df_sp); # DataCF object
```
"""
function mapAllelesCFtable(alleleDF::AbstractString, cfDF::AbstractString;
        filename=""::AbstractString, columns=Int[]::Vector{Int}, CSVargs...)
    # force categorical=false unless the user wants otherwise
    if :categorical ∉ [pair[1] for pair in CSVargs]
        CSVargs = (CSVargs..., (:categorical, false))
    end
    d = CSV.read(alleleDF; CSVargs...)
    d2 = CSV.read(cfDF; CSVargs...)
    mapAllelesCFtable!(d2,d, columns, filename != "", filename)
end

"""
    mapAllelesCFtable!(quartet CF DataFrame, mapping DataFrame, columns, write?, filename)

Modify (and return) the quartet concordance factor (CF) DataFrame:
replace each allele name by the species name that the allele maps onto
based on the mapping data frame. This mapping data frame should have columns
named "allele" and "species" (see `rename!` to change column names if need be).

If `write?` is `true`, the modified data frame is written to a file named "filename".

Warning: [`mapAllelesCFtable`](@ref) takes the quartet data file as its second
argument, while `mapAllelesCFtable!` takes the quartet data (which it modifies)
as its first argument.
"""
function mapAllelesCFtable!(cfDF::DataFrame, alleleDF::DataFrame, co::Vector{Int},write::Bool,filename::AbstractString)
    size(cfDF,2) >= 7 || error("CF DataFrame should have 7+ columns: 4taxa, 3CF, and possibly ngenes")
    if length(co)==0 co=[1,2,3,4]; end
    compareTaxaNames(alleleDF,cfDF,co)
    for j in 1:4
        for ia in 1:size(alleleDF,1) # for all alleles
            cfDF[!,co[j]] = map(x->replace(string(x),
                                         Regex("^$(string(alleleDF[ia,:allele]))\$") =>
                                         alleleDF[ia,:species]),
                                cfDF[!,co[j]])
        end
    end
    if write
        filename != "" || error("cannot write quartet CF with mapped alleles: empty filename")
        CSV.write(filename, cfDF)
    end
    return cfDF
end

# function to clean a df after changing allele names to species names
# inside readTableCF!
# by deleting rows that are not informative like sp1 sp1 sp1 sp2
# keepOne=true: we only keep one allele per species
function cleanAlleleDF!(newdf::DataFrame, cols::Vector{Int};keepOne=false::Bool)
    withngenes = (length(cols)==8)
    delrows = Int[] # indices of rows to delete
    repSpecies = String[]
    if(isa(newdf[1,cols[1]],Integer)) #taxon names as integers: we need this to be able to add __2
        newdf[!,cols[1]] = map(string, newdf[!,cols[1]])
        newdf[!,cols[2]] = map(string, newdf[!,cols[2]])
        newdf[!,cols[3]] = map(string, newdf[!,cols[3]])
        newdf[!,cols[4]] = map(string, newdf[!,cols[4]])
    end
    row = Vector{String}(undef, 4)
    for i in 1:size(newdf,1) #check all rows
        @debug "row number: $i"
        # fixit: check for no missing value, or error below
        map!(j -> newdf[i,cols[j]], row, 1:4)
        @debug "row $(row)"
        uniq = unique(row)
        @debug "unique $(uniq)"

        keep = false # default: used if 1 unique name, or 2 in some cases
        if(length(uniq) == 4)
            keep = true
        else
            if(!keepOne)
                if(length(uniq) == 3) #sp1 sp1 sp2 sp3
                    keep = true
                    for u in uniq
                        @debug "u $(u), typeof $(typeof(u))"
                        ind = row .== u #taxon names matching u
                        @debug "taxon names matching u $(ind)"
                        if(sum(ind) == 2)
                            push!(repSpecies,string(u))
                            found = false
                            for k in 1:4
                                if(ind[k])
                                    if(found)
                                        @debug "found the second one in k $(k), will change newdf[i,cols[k]] $(newdf[i,cols[k]]), typeof $(typeof(newdf[i,cols[k]]))"
                                        newdf[i,cols[k]] = string(u, repeatAlleleSuffix)
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
                    # keep was initialized to false
                    for u in uniq
                        @debug "length uniq is 2, u $(u)"
                        ind = row .== u
                        if(sum(ind) == 1 || sum(ind) == 3)
                            @debug "ind $(ind) is 1 or 3, should not keep"
                            break
                        elseif(sum(ind) == 2)
                            @debug "ind $(ind) is 2, should keep"
                            keep = true
                            found = false
                            push!(repSpecies,string(u))
                            for k in 1:4
                                if(ind[k])
                                    if(found)
                                        newdf[i,cols[k]] = string(u, repeatAlleleSuffix)
                                        break
                                    else
                                        found = true
                                    end
                                end
                            end
                        end
                    end
                end
                @debug "after if, keep is $(keep)"
            end
        end
        keep || push!(delrows, i)
        @debug "" keep
    end
    @debug "" delrows
    @debug "" repSpecies
    nrows = size(newdf,1)
    nkeep = nrows - length(delrows)
    if nkeep < nrows
        print("""found $(length(delrows)) 4-taxon sets uninformative about between-species relationships, out of $(nrows).
              These 4-taxon sets will be deleted from the data frame. $nkeep informative 4-taxon sets will be used.
              """)
        nkeep > 0 || @warn "All 4-taxon subsets are uninformative, so the dataframe will be left empty"
        deleterows!(newdf, delrows)
    end
    # @show size(newdf)
    return unique(repSpecies)
end


# function to merge rows that have repeated taxon names by using the weigthed average of CF
# (if info on number of genes is provided) or simple average
function mergeRows(df::DataFrame, cols::Vector{Int})
    sorttaxa!(df, cols) # sort taxa alphabetically within each row
    colnam = names(df)[cols[5:end]]
    df = aggregate(df, names(df)[cols[1:4]], mean);
    rename!(df, Dict((Symbol(n, "_mean"), n) for n in colnam) )
    n4tax = size(df,1) # total number of 4-taxon sets
    print("$n4tax unique 4-taxon sets were found. CF values of repeated 4-taxon sets will be averaged")
    println((length(cols)>7 ? " (ngenes too)." : "."))
    return df
end


# function to expand leaves in tree to two individuals
# based on cf table with alleles mapped to species names
function expandLeaves!(repSpecies::Union{Vector{String},Vector{Int}},tree::HybridNetwork)
    for sp in repSpecies
        for n in tree.node
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
                n2.name = string(sp,repeatAlleleSuffix)
                push!(tree.names,n2.name)
                break
            end
        end
    end
end


# function to compare the taxon names in the allele-species matching table
# and the CF table
function compareTaxaNames(alleleDF::DataFrame, cfDF::DataFrame, co::Vector{Int})
    checkMapDF(alleleDF)
    #println("found $(length(alleleDF[1])) allele-species matches")
    CFtaxa = string.(mapreduce(x -> unique(skipmissing(x)), union, eachcol(cfDF[!,co[1:4]])))
    alleleTaxa = map(string, alleleDF[!,:allele]) # as string, too
    sizeCF = length(CFtaxa)
    sizeAllele = length(alleleTaxa)
    if sizeAllele > sizeCF
        @warn "some alleles in the mapping file do not occur in the quartet CF data. Extra allele names will be ignored"
        alleleTaxa = intersect(alleleTaxa,CFtaxa)
    end
    unchanged = setdiff(CFtaxa,alleleTaxa)
    if length(unchanged) == length(CFtaxa)
        @warn "None of the taxon names in CF table match with allele names in the mapping file"
    end
    if !isempty(unchanged)
        warnmsg = "not all alleles were mapped\n"
        warnmsg *= "alleles not mapped to a species name:"
        for n in unchanged warnmsg *= " $n"; end
        @warn warnmsg
    end
    return nothing
end

# function to check that the allele df has one column labelled alleles and one column labelled species
function checkMapDF(alleleDF::DataFrame)
    size(alleleDF,2) <= 2 || error("Allele-Species matching Dataframe should have at least 2 columns")
    size(alleleDF,2) >= 2 || @warn "allele mapping file contains more than two columns: will ignore all columns not labelled allele or species"
    try
        alleleDF[!,:allele]
    catch
        error("In allele mapping file there is no column named allele")
    end
    try
        alleleDF[!,:species]
    catch
        error("In allele mapping file there is no column named species")
    end
end



## function to check if a new proposed topology satisfies the condition for
## multiple alleles: no gene flow to either allele, and both alleles as sister
## returns false if the network is not ok
function checkTop4multAllele(net::HybridNetwork)
    for n in net.leaf
        if occursin(repeatAlleleSuffix_re, n.name)
            n.leaf || error("weird node $(n.number) not leaf in net.leaf list")
            length(n.edge) == 1 || error("weird leaf with $(length(n.edge)) edges")
            par = getOtherNode(n.edge[1],n)
            if(par.hybrid) ## there is gene flow into n
                return false
            end
            nameOther = replace(n.name, repeatAlleleSuffix_re => "")
            foundOther = false
            for i in 1:3
                other = getOtherNode(par.edge[i],par)
                if(other.leaf && other.name == nameOther)
                    foundOther = true
                end
            end
            foundOther || return false
        end
    end
    return true
end



## function to merge the two alleles into one
function mergeLeaves!(net::HybridNetwork)
    leaves = copy(net.leaf) # bc we change this list
    for n in leaves
        if occursin(repeatAlleleSuffix_re, n.name)
            n.leaf || error("weird node $(n.number) not leaf in net.leaf list")
            length(n.edge) == 1 || error("weird leaf with $(length(n.edge)) edges")
            par = getOtherNode(n.edge[1],n)
            foundOther = false
            other = Node()
            nameOther = replace(n.name, repeatAlleleSuffix_re => "")
            for i in 1:3
                other = getOtherNode(par.edge[i],par)
                if(other.leaf && other.name == nameOther)
                    foundOther = true
                    break
                end
            end
            if(!foundOther)
                checkTop4multAllele(net) || error("current network does not comply with multiple allele condition")
                error("strange network that passes checkTop4multAllele, but cannot find the other allele for $(n.name)")
            end
            removeEdge!(par,n.edge[1])
            removeEdge!(par,other.edge[1])
            deleteNode!(net,n)
            deleteNode!(net,other)
            deleteEdge!(net,n.edge[1])
            deleteEdge!(net,other.edge[1])
            par.name = other.name
            par.leaf = true
            push!(net.leaf,par)
            net.numTaxa += 1
        end
    end
end
