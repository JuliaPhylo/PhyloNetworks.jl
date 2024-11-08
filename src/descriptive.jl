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

# function that we need to overwrite to avoid printing useless scary
# output for HybridNetworks
# PROBLEM: writenewick changes the network and thus show changes the network
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
        # par = writenewick(obj,round=true) # but writenewick changes the network, not good
        s = IOBuffer()
        writeSubTree!(s, obj, false,true, true,3,true)
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

