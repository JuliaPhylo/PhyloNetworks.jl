# functions to describe a HybridNetwork object to avoid accessing the attributes directly
# Claudia August 2015

"""
`tipLabels(net::HybridNetwork)`

returns a vector of taxon names (at the leaves) from a HybridNetwork object
"""
function tipLabels(net::HybridNetwork)
    return ASCIIString[l.name for l in net.leaf] # AbstractString does not work for use by tree2Matrix
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
                for(n in net.hybrid)
                    flag,edges = updateContainRoot!(net,n)
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

"""
`root!(net::HybridNetwork, nodeNumber)`

root the network/tree object at the node with nodeNumber
"""
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
"""
`root!(net::HybridNetwork, outgroup)`

root the network/tree object with the outgroup taxon name given as argument.
"""
function root!(net::HybridNetwork, outgroup::AbstractString)
    if(!isTree(net))
        if(!net.cleaned)
            DEBUG && println("net not cleaned inside root, need to run updateCR")
            for(n in net.hybrid)
                flag,edges = updateContainRoot!(net,n)
                flag || error("hybrid node $(n.hybrid) has conflicting containRoot")
            end
        end
    end
    updateRoot!(net,outgroup)
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
function show(io::IO, net::HybridNetwork)
    print(io,"$(writeTopology(net))")
end

# and QuartetNetworks (which cannot be just written because they do not have root)
function show(io::IO, net::QuartetNetwork)
    print(io,"taxa: $(net.quartetTaxon)\n")
    print(io,"number of hybrid nodes: $(net.numHybrids)\n")
    if(net.split != [-1,-1,-1,-1])
        print(io,"split: $(net.split)\n")
    end
end

function show(io::IO,d::DataCF)
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

function show(io::IO,q::Quartet)
    print(io,"number: $(q.number)\n")
    print(io,"taxon names: $(q.taxon)\n")
    print(io,"observed CF: $(q.obsCF)\n")
    print(io,"-logPseudo-dev under best estimated network $(q.logPseudoLik) (meaningless before estimation)\n")
    print(io,"expected CF under best estimated network: $(q.qnet.expCF) (meaningless before estimation)\n")
    if(q.ngenes != -1)
        print(io,"number of genes used to compute observed CF: $(q.ngenes)\n")
    end
end


