# functions to describe a HybridNetwork object to avoid accessing the attributes directly
# Claudia August 2015

"""
    tiplabels(x)

Vector of taxon names at the leaves, defined for objects of various types:
`HybridNetwork`,
`MatrixTopologicalOrder`.

For a network, the taxon names are coerced to strings.
"""
function tiplabels(net::HybridNetwork)
    return String[l.name for l in net.leaf] # AbstractString does not work for use by tree2Matrix
end

# extract & sort the union of taxa of list of gene trees
function uniontaxa(trees::Vector{HybridNetwork})
    taxa = reduce(union, tiplabels(t) for t in trees)
    return sort_stringasinteger!(taxa)
end

"""
    sort_stringasinteger!(taxa)

Sort a vector of strings `taxa`, numerically if
elements can be parsed as an integer, alphabetically otherwise.
"""
function sort_stringasinteger!(taxa)
    sortby = x->parse(Int,x)
    try
        parse.(Int,taxa)
    catch
        sortby = identity
    end
    sort!(taxa, by=sortby)
    return taxa
end

# 'show' should not (and does not) modify the object
function Base.show(io::IO, obj::HybridNetwork)
    disp = "$(typeof(obj)), "
    if obj.isrooted
        disp = disp * "Rooted Network"
    else
        disp = disp * "Semidirected Network"
    end
    disp = disp * "\n$(obj.numedges) edges\n"
    disp = disp * "$(obj.numnodes) nodes: $(obj.numtaxa) tips, "
    disp = disp * "$(obj.numhybrids) hybrid nodes, "
    disp = disp * "$(obj.numnodes - obj.numtaxa - obj.numhybrids) internal tree nodes.\n"
    tipslabels = [n.name for n in obj.leaf]
    if length(tipslabels) > 1 || !all(tipslabels .== "")
        disptipslabels = "$(tipslabels[1])"
        for i in 2:min(obj.numtaxa, 4)
            disptipslabels = disptipslabels * ", $(tipslabels[i])"
        end
        if obj.numtaxa > 4 disptipslabels = disptipslabels * ", ..." end
        disp *= "tip labels: " * disptipslabels
    end
    par = ""
    try
        #= *not* 'writenewick(obj,round=true)' because writenewick changes
        the network (runs directedges! and if it fails, tries other rootings).
        writesubtree! does *not* modify the network.
        =#
        s = IOBuffer()
        writesubtree!(s, obj, false,true, true,3,true)
        par = String(take!(s))
    catch err
        println("ERROR with writesubtree!:")
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
        disp *= "\n " * (obj.ismajor ? "major" : "minor")
        disp *= " hybrid edge with gamma=$(obj.gamma)"
    elseif (!obj.ismajor)
        disp *= "\n minor tree edge"
    end
    disp *= "\n attached to $(length(obj.node)) node(s) (parent first):"
    if (length(obj.node)==1) disp *= " $(obj.node[1].number)";
    elseif (length(obj.node)==2)
        disp *= " $(obj.node[obj.ischild1 ? 2 : 1].number)"
        disp *= " $(obj.node[obj.ischild1 ? 1 : 2].number)"
    end
    println(io, disp)
end

