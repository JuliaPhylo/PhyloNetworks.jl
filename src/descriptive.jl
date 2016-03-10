# functions to describe a HybridNetwork object to avoid accessing the attributes directly
# Claudia August 2015

"""
`tipLabels(net::HybridNetwork)`

returns a vector of taxon names (at the leaves) from a HybridNetwork object
"""
function tipLabels(net::HybridNetwork)
    return ASCIIString[l.name for l in net.leaf] # AbstractString does not work for use by tree2Matrix
end



"""
`dfObsExpCF(d::DataCF)`

function that will create a table with the observed and expected CF after estimation of a network with snaq(T,d).
"""
function dfObsExpCF(d::DataCF)
    df=DataFrame(obsCF1=[q.obsCF[1] for q in d.quartet],obsCF2=[q.obsCF[2] for q in d.quartet],obsCF3=[q.obsCF[3] for q in d.quartet], expCF1=[q.qnet.expCF[1] for q in d.quartet],expCF2=[q.qnet.expCF[2] for q in d.quartet],expCF3=[q.qnet.expCF[3] for q in d.quartet])
    return df
end


# function to set nonidentifiable edges BL to -1.0
# used at the end of optTopRuns
function setNonIdBL!(net::HybridNetwork)
    for(e in net.edge)
        if(!e.istIdentifiable)
            e.length = -1.0 #do not use setLength because it does not allow BL too negative
        end
    end
end

# function that we need to overwrite to avoid printing useless scary
# output for HybridNetworks

function Base.show(io::IO, obj::HybridNetwork)
    disp = "$(typeof(obj)), "
    if obj.isRooted
        disp = disp * "Rooted Network"
    else
        disp = disp * "Un-rooted Network"
    end
    disp = disp * "\n$(obj.numEdges) edges\n"
    disp = disp * "$(obj.numNodes) nodes: $(obj.numTaxa) tips, "
    disp = disp * "$(obj.numHybrids) hybrid nodes, "
    disp = disp * "$(obj.numNodes - obj.numTaxa - obj.numHybrids) internal tree nodes.\n"
    tipslabels = [n.name for n in obj.leaf]
    if length(tipslabels) > 1 || !all(tipslabels .== "")
        disptipslabels = "$(tipslabels[1])"
        for i in 2:min(obj.numTaxa, 4)
            disptipslabels = disptipslabels * ", $(tipslabels[i])"
        end
        if obj.numTaxa > 4 disptipslabels = disptipslabels * ", ..." end
        disp = disp * "\nTips labels: \n " * disptipslabels
    end
    disp = disp * "\n\n $(writeTopology(obj))"
    println(io, disp)
end

# function show(io::IO, net::HybridNetwork)
#     print(io,"$(writeTopology(net))")
# end

# and QuartetNetworks (which cannot be just written because they do not have root)
function Base.show(io::IO, net::QuartetNetwork)
    print(io,"taxa: $(net.quartetTaxon)\n")
    print(io,"number of hybrid nodes: $(net.numHybrids)\n")
    if(net.split != [-1,-1,-1,-1])
        print(io,"split: $(net.split)\n")
    end
end

function Based.show(io::IO,d::DataCF)
    print(io,"-----------------------\n")
    print(io,"Object DataCF\n")
    print(io,"number of quartets: $(d.numQuartets)\n")
    if(d.numTrees != -1)
        print(io,"number of trees: $(d.numTrees)\n")
        print(io,"For your object DataCF, you can access the list of Quartet types with the attribute .quartet. \nFor example, if your DataCF object is named d, d.quartet[1] will print the first quartet.")
    else
        print(io, "CF not computed from input gene trees, information on gene trees could be present in the quartets if the column ngenes was in the CF table.\n")
        print(io,"For your object DataCF, you can access the list of Quartet types with the attribute .quartet, and the list of HybridNetwork types (if the input was a list of trees) with the attribute .tree. \nFor example, if your DataCF object is named d, d.quartet[1] will print the first quartet, and d.tree[1] will print the first tree.")
    end
end

function Base.show(io::IO,q::Quartet)
    print(io,"number: $(q.number)\n")
    print(io,"taxon names: $(q.taxon)\n")
    print(io,"observed CF: $(q.obsCF)\n")
    print(io,"-logPseudo-dev under best estimated network $(q.logPseudoLik) (meaningless before estimation)\n")
    print(io,"expected CF under best estimated network: $(q.qnet.expCF) (meaningless before estimation)\n")
    if(q.ngenes != -1)
        print(io,"number of genes used to compute observed CF: $(q.ngenes)\n")
    end
end

function Base.show(io::IO, obj::Node)
    disp = "$(typeof(obj)):"
    disp = disp * "\n Node Number:$(obj.number)"
    if (obj.name != "") disp = disp * "\n Node Name:$(obj.name)" end
    if (obj.hybrid) disp = disp * "\n This node is a hybrid" end
    if (obj.leaf) disp = disp * "\n This node is a leaf" end
    println(io, disp)
end

function Base.show(io::IO, obj::Edge)
    disp = "$(typeof(obj)):"
    disp = disp * "\n Edge Number:$(obj.number)"
    disp = disp * "\n Edge length:$(obj.length)"
    if (obj.hybrid) disp = disp * "\n This edge is a hybrid edge with gamma=$(obj.gamma)" end
    println(io, disp)
end

# function to check (and change) all attributes in a set of trees
    function updateTrees!(trees::Vector{HybridNetwork})
        for(t in trees)
            if(isTree(t))
                for(e in t.edge)
                    e.containRoot = true
                    e.inCycle = -1
                end
                for(n in t.node)
                    n.inCycle = -1
                end
            end
        end
    end


