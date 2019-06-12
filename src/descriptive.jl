# functions to describe a HybridNetwork object to avoid accessing the attributes directly
# Claudia August 2015

"""
    tipLabels(x)


Return a vector of taxon names at the leaves, for objects of various types:
`HybridNetwork`,
Vector of `HybridNetwork`s (in which case the union is taken then sorted),
Vector of `Quartet`s, `DataCF`,
`TraitSimulation`, `MatrixTopologicalOrder`.

For a network, the taxon names are coerced to strings.
"""
function tipLabels(net::HybridNetwork)
    return String[l.name for l in net.leaf] # AbstractString does not work for use by tree2Matrix
end



"""
`fittedQuartetCF(d::DataCF, format::Symbol)`

return a data frame with the observed and expected quartet concordance factors
after estimation of a network with snaq(T,d).
The format can be :wide (default) or :long.

- if wide, the output has one row per 4-taxon set, and each row has 10 columns: 4 columns
  for the taxon names, 3 columns for the observed CFs and 3 columns for the expected CF.
- if long, the output has one row per quartet, i.e. 3 rows per 4-taxon sets, and 7 columns:
  4 columns for the taxon names, one column to give the quartet resolution, one column for
  the observed CF and the last column for the expected CF.

see also: `topologyQPseudolik!` and `topologyMaxQPseudolik!` to update the fitted CF expected
under a specific network, inside the DataCF object `d`.
"""
function fittedQuartetCF(d::DataCF, format=:wide::Symbol)
    if format == :wide
        df=DataFrame(
                 tx1 = [q.taxon[1] for q in d.quartet],
                 tx2 = [q.taxon[2] for q in d.quartet],
                 tx3 = [q.taxon[3] for q in d.quartet],
                 tx4 = [q.taxon[4] for q in d.quartet],
                obsCF12=[q.obsCF[1] for q in d.quartet],
                obsCF13=[q.obsCF[2] for q in d.quartet],
                obsCF14=[q.obsCF[3] for q in d.quartet],
                expCF12=[q.qnet.expCF[1] for q in d.quartet],
                expCF13=[q.qnet.expCF[2] for q in d.quartet],
                expCF14=[q.qnet.expCF[3] for q in d.quartet])
    elseif format == :long
        nQ = length(d.quartet)
        df = DataFrame(
            tx1 = repeat([q.taxon[1] for q in d.quartet], inner=[3]),
            tx2 = repeat([q.taxon[2] for q in d.quartet], inner=[3]),
            tx3 = repeat([q.taxon[3] for q in d.quartet], inner=[3]),
            tx4 = repeat([q.taxon[4] for q in d.quartet], inner=[3]),
            quartet = repeat(["12_34","13_24","14_23"], outer=[nQ]),
            obsCF = Array{Float64}(undef, 3*nQ),
            expCF = Array{Float64}(undef, 3*nQ)  )
        row = 1
        for i in 1:nQ
            for j in 1:3
                df[row, 6] = d.quartet[i].obsCF[j]
                df[row, 7] = d.quartet[i].qnet.expCF[j]
                row += 1
            end
        end
    else
        error("format $(format) was not recognized. Should be :wide or :long.")
    end
    return df
end

"""
`setNonIdBL!(net)`

Set non-identifiable edge branch lengths to -1.0 (i.e. missing) for a level-1 network `net`,
except for edges in

- a good triangle: the edge below the hybrid is constrained to 0.
- a bad diamond II: the edge below the hybrid is constrained to 0
- a bad diamond I: the edges across from the hybrid node have non identifiable lengths
  but are kept, because the two γ*(1-exp(-t)) values are identifiable.

will break if `inCycle` attributes are not initialized (at -1) or giving a correct node number.

see [`Node`](@ref) for the meaning of boolean attributes
`isBadTriangle` (which corresponds to a "good" triangle above),
`isBadDiamondI` and `isBadDiamondII`.
"""
function setNonIdBL!(net::HybridNetwork)
    for e in net.edge
        if !e.istIdentifiable
            keeplength = any(n -> (n.isBadDiamondII || n.isBadTriangle), e.node)
            # if below 'bad' hybrid node, length=0 by constraint. If above, length estimated.
            # next: keep length if edge across from hybrid node in bad diamond I.
            if  !keeplength && e.inCycle != -1
                hyb = net.node[getIndexNode(e.inCycle, net)]
                if hyb.isBadDiamondI
                    keeplength |= !any(n -> (n==hyb), e.node) # only if e across hyb: no touching it
                end
            end
            if !keeplength
                e.length = -1.0 # not setLength, which does not allow too negative BLs
            end
        end
    end
end

"""
setBLGammaParsimony!(net::HybridNetwork)

Maximum parsimony function does not provide estimates for branch lengths, or gamma.
But since the `maxParsimonyNet` function is using snaq move functions, branch lengths
and gamma values are set randomly (to initialize optimization).
We need to remove these random values before returning the maximum parsimony network.
"""
function setBLGammaParsimony!(net::HybridNetwork)
    for e in net.edge
        e.length = -1.0
        if e.hybrid
            e.gamma = -1.0
        end
    end
end


# function that we need to overwrite to avoid printing useless scary
# output for HybridNetworks
# PROBLEM: writeTopology changes the network and thus show changes the network
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
        disp *= "tip labels: " * disptipslabels
    end
    par = ""
    try
        # par = writeTopology(obj,round=true) # but writeTopology changes the network, not good
        s = IOBuffer()
        writeSubTree!(s, obj, false,true, true,3,false)
        par = String(take!(s))
    catch err
        println("ERROR with writeSubTree!:")
        showerror(stdout, err)
        println("Trying writeTopologyLevel1")
        par = writeTopologyLevel1(obj)
    end
    disp *= "\n$par"
    println(io, disp)
end

# and QuartetNetworks (which cannot be just written because they do not have root)
function Base.show(io::IO, net::QuartetNetwork)
    print(io,"taxa: $(net.quartetTaxon)\n")
    print(io,"number of hybrid nodes: $(net.numHybrids)\n")
    if(net.split != [-1,-1,-1,-1])
        print(io,"split: $(net.split)\n")
    end
end

function Base.show(io::IO,d::DataCF)
    print(io,"Object DataCF\n")
    print(io,"number of quartets: $(d.numQuartets)\n")
    if(d.numTrees != -1)
        print(io,"number of trees: $(d.numTrees)\n")
    end
end

function Base.show(io::IO,q::Quartet)
    print(io,"number: $(q.number)\n")
    print(io,"taxon names: $(q.taxon)\n")
    print(io,"observed CF: $(q.obsCF)\n")
    print(io,"pseudo-deviance under last used network: $(q.logPseudoLik) (meaningless before estimation)\n")
    print(io,"expected CF under last used network: $(q.qnet.expCF) (meaningless before estimation)\n")
    if(q.ngenes != -1)
        print(io,"number of genes used to compute observed CF: $(q.ngenes)\n")
    end
end

function Base.show(io::IO, obj::Node)
    disp = "$(typeof(obj)):"
    disp = disp * "\n number:$(obj.number)"
    if (obj.name != "") disp *= "\n name:$(obj.name)" end
    if (obj.hybrid)     disp *= "\n hybrid node" end
    if (obj.leaf)       disp *= "\n leaf node" end
    disp *= "\n attached to $(length(obj.edge)) edges, numbered:"
    for e in obj.edge disp *= " $(e.number)"; end
    println(io, disp)
end

function Base.show(io::IO, obj::Edge)
    disp = "$(typeof(obj)):"
    disp *= "\n number:$(obj.number)"
    disp *= "\n length:$(obj.length)"
    if (obj.hybrid)
        disp *= "\n " * (obj.isMajor ? "major" : "minor")
        disp *= " hybrid edge with gamma=$(obj.gamma)"
    elseif (!obj.isMajor)
        disp *= "\n minor tree edge"
    end
    disp *= "\n attached to $(length(obj.node)) node(s) (parent first):"
    if (length(obj.node)==1) disp *= " $(obj.node[1].number)";
    elseif (length(obj.node)==2)
        disp *= " $(obj.node[obj.isChild1 ? 2 : 1].number)"
        disp *= " $(obj.node[obj.isChild1 ? 1 : 2].number)"
    end
    println(io, disp)
end

