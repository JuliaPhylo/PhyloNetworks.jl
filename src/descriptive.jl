# functions to describe a HybridNetwork object to avoid accessing the attributes directly
# Claudia August 2015

function tipLabels(net::HybridNetwork)
    return [l.name for l in net.leaf]
end

# function to re root on a node
# resolve=true, a branch of length=0 is
# added if node is internal if node is leaf, root(net,outgroup) is
# called, so a new node is created on external edge
function root!(net::HybridNetwork, node::Node, resolve::Bool)
    node.hybrid && error("node $(node.number) is hybrid, cannot root network on a hybrid node")
    if(node.leaf)
        warn("node $(node.number) is a leaf, so we will root as an outgroup if possible")
        root!(net,node.name)
    else
        if(!isTree(net))
            if(!net.cleaned)
                DEBUG && println("net not cleaned inside root, need to run updateCR")
                for(n in graph.hybrid)
                    flag,edges = updateContainRoot!(graph,n)
                    flag || error("hybrid node $(n.hybrid) has conflicting containRoot")
                end
            end
        end
        if(canBeRoot(node))
            try
                ind = getIndex(node,net)
            catch
                error("cannot set node $(node.number) as root because it is not part of net")
            end
            ind = getIndex(node,net)
            net.root = ind
            if(resolve)
                resolve!(net,node)
            end
        else
            warn("node $(node.number) cannot be root, will leave root as is")
        end
    end
end

function root!(net::HybridNetwork, nodeNum::Int64, resolve::Bool)
    try
        ind = getIndexNode(nodeNum,net)
    catch
        error("cannot set node $(nodeNum) as root because it is not part of net")
    end
    ind = getIndexNode(nodeNum,net)
    root!(net,net.node[ind],resolve)
end

root!(net::HybridNetwork, nodeNum::Int64) = root!(net, nodeNum, false)

# function to resolve an internal node for root function
function resolve!(net::HybridNetwork, node::Node)
    length(node.edge) == 3 || error("node $(node.number) has $(length(node.edge)) edges instead of 3")
    !node.hybrid || error("node $(node.number) is hybrid and cannot root/resolve on hybrid node")
    canBeRoot(node) || error("cannot resolve root in node $(node.number) because it cannot be root to begin with")
    if(net.cleaned)
        if(node.inCycle == -1) #node not in cycle
            done=false
            for(e in node.edge)
                if(e.containRoot)
                    done = true
                    newNodeResolve!(net,e,node)
                    net.root = length(net.node) #last node is root
                    break
                end
            end
            done || error("not found edge to contain root for node $(node.number)")
        else
            done = false
            for(e in node.edge)
                if(e.inCycle == -1 && e.containRoot)
                    done = true
                    newNodeResolve!(net,e,node)
                    net.root = length(net.node) #last node is root
                    break
                end
            end
            done || error("strange: not found edge to contain root not in cycle for node $(node.number) even when node can be root")
        end
    else # net not cleaned
        if(!node.hasHybEdge)
            done=false
            for(e in node.edge)
                if(!e.hybrid)
                    done = true
                    newNodeResolve!(net,e,node)
                    net.root = length(net.node) #last node is root
                    break
                end
            end
            done || error("not found edge to contain root for node $(node.number)")
        else
            done = false
            for(e in node.edge)
                if(e.inCycle == -1 && e.containRoot)
                    done = true
                    newNodeResolve!(net,e,node)
                    net.root = length(net.node) #last node is root
                    break
                end
            end
            done || error("strange: not found edge to contain root not in cycle for node $(node.number) even when node can be root")
        end
    end
end

# function to create new node in edge for resolve
function newNodeResolve!(net::HybridNetwork,e::Edge, node::Node)
    removeEdge!(node,e)
    removeNode!(node,e)
    max_edge = maximum([e.number for e in net.edge]);
    max_node = maximum([e.number for e in net.node]);
    newedge = Edge(max_edge+1)
    newnode = Node(max_node+1,false,false,[e,newedge])
    if(net.cleaned && !isTree(net))
        part = whichPartition(net,e)
        push!(net.partition[part].edges,newedge)
    end
    setNode!(e,newnode)
    setNode!(newedge,newnode)
    setEdge!(node,newedge)
    setNode!(newedge,node)
    pushEdge!(net,newedge)
    pushNode!(net,newnode)
    t = e.length
    setLength!(e,t/2)
    setLength!(newedge,t/2)
end

# function to root a network on an outgroup
# (single taxon)
function root!(net::HybridNetwork, outgroup::String)
    if(!isTree(net))
        if(!net.cleaned)
            DEBUG && println("net not cleaned inside root, need to run updateCR")
            for(n in graph.hybrid)
                flag,edges = updateContainRoot!(graph,n)
                flag || error("hybrid node $(n.hybrid) has conflicting containRoot")
            end
        end
    end
    updateRoot!(net,outgroup)
end


# function to read a table of allele-species matchs (dataframe with 2 columns)
# and a table of CF in the allele names, and replace all the allele names
# to the species names
# this will create a new CF table, will not rewrite on the original one
# filename is the name to give to the new table, by default "newCFtable"
function renameTaxaCFtable(alleleDF::DataFrame, cfDF::DataFrame,filename::String)
    size(alleleDF,2) == 2 || error("Allele-Species matching Dataframe should have 2 columns only: allele names, species names")
    size(cfDF,2) == 7 || warn("CF Dataframe should have 7 columns: 4taxa, 3CF, will ignore columns from 8th on")
    compareTaxaNames(alleleDF,cfDF)
    newt1 = map(x->replace(string(x),string(alleleDF[1,1]),alleleDF[1,2]),cfDF[1])
    newt2 = map(x->replace(string(x),string(alleleDF[1,1]),alleleDF[1,2]),cfDF[2])
    newt3 = map(x->replace(string(x),string(alleleDF[1,1]),alleleDF[1,2]),cfDF[3])
    newt4 = map(x->replace(string(x),string(alleleDF[1,1]),alleleDF[1,2]),cfDF[4])
    if(size(alleleDF,1) > 1)
        for(j in 2:size(alleleDF,1)) #for all the allele matchings
            newt1 = map(x->replace(string(x),string(alleleDF[j,1]),alleleDF[j,2]),newt1)
            newt2 = map(x->replace(string(x),string(alleleDF[j,1]),alleleDF[j,2]),newt2)
            newt3 = map(x->replace(string(x),string(alleleDF[j,1]),alleleDF[j,2]),newt3)
            newt4 = map(x->replace(string(x),string(alleleDF[j,1]),alleleDF[j,2]),newt4)
        end
    end
    newdf = DataFrames.DataFrame(t1=newt1,t2=newt2,t3=newt3,t4=newt4,CF1234=cfDF[5],CF1324=cfDF[6],CF1423=cfDF[7])
    writetable(filename,newdf)
    return newdf
end

function renameTaxaCFtable(alleleDF::String, cfDF::String; filename="newCFtable.txt"::String)
    d = readtable(alleleDF)
    d2 = readtable(cfDF)
    renameTaxaCFtable(d,d2,filename)
end

# function to compare the taxon names in the allele-species matching table
# and the CF table
function compareTaxaNames(alleleDF::DataFrame, cfDF::DataFrame)
    size(alleleDF,2) == 2 || error("Allele-Species matching Dataframe should have 2 columns only: allele names, species names")
    size(cfDF,2) == 7 || warn("CF Dataframe should have 7 columns: 4taxa, 3CF, will ignore columns from 8th on")
    d = readTableCF(cfDF)
    println("there are $(length(alleleCF[1])) allele-species matches")
    CFtaxa = unionTaxa(d.quartet)
    CFtaxa = map(x->string(x),CFtaxa) #treat as string
    alleleTaxa = map(x->string(x),alleleDF[1]) #treat as string
    sizeCF = length(CFtaxa)
    sizeAllele = length(alleleTaxa)
    if(sizeAllele > sizeCF)
        println("there are more taxa in the allele-species map: $(sizeAllele) than in the CF table: $(sizeCF), so extra allele names will be ignored")
        alleleTaxa = intersect(alleleTaxa,CFtaxa)
    elseif(sizeAllele < sizeCF)
        println("there are fewer taxa in the allele-species map: $(sizeAllele) than in the CF table: $(sizeCF), so names in the CF table will remained unchanged")
    else
        println("there are the same number of taxa in the allele-species map: $(sizeAllele) as in the CF table: $(sizeCF), but we still need to verify that the names match")
    end
    unchanged = setdiff(CFtaxa,alleleTaxa)
    if(isempty(unchanged))
        println("TAXON NAMES MATCH: all taxa in the CF table was changed according to the allele-species map")
    else
        println("the following taxa in the CF table were not modified to the allele-species map (since they are not present in the allele-species table):\n $(unchanged)")
    end
end


