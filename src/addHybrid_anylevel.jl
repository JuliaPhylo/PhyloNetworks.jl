# functions to add hybridizations to a semi-directed network
# subject to topological constraints; no level restriction

"""
    addhybridedge!(net::HybridNetwork, no3cycle::Bool, nohybridladder::Bool,
                   constraints=TopologyConstraint[]::Vector{TopologyConstraint};
                   maxattempts=10::Int)

Randomly choose two edges in `net` then: add hybrid edge from edge 1 to edge 2
and randomly decide which "half" of edge 2 will serve as partner hybrid edge.
This partner edge will point toward the newly-created node
on the middle of the original edge 2.

If the resulting network is a DAG, satisfies the constraint(s),
does not contain any 3-cycle (if `no3cycle=true`), and does not have
a hybrid ladder (if `nohybridladder=true`) then the proposal is successful:
`net` is modified, and the function returns the newly created hybrid node.

If the resulting network is not acceptable, then a new set of edges
is proposed (using a blacklist), until one is found acceptable, or until
a maximum number of attempts have been made (`maxattempts`).
If none of the attempted proposals are successful, `nothing` is returned
(without causing an error).

After a pair of edges is picked, the "top" half of edge2 is proposed
as the partner hybrid edge with probability 0.8
(to avoid changing the direction of edge2 with more than 50% chance).
If this choice does not work, the "bottom" half of edge2 is proposed as
the partner hybrid edge (which would require to flip the direction of edge 2).
"""
function addhybridedge!(net::HybridNetwork, no3cycle::Bool, nohybridladder::Bool,
        constraints=TopologyConstraint[]::Vector{TopologyConstraint};
        maxattempts=10::Int)
    all(con.type == 2 for con in constraints) || error("only type-2 constraints implemented so far")
    numedges = length(net.edge)
    blacklist = Set{Tuple{Int,Int}}()
    nmax_blacklist = numedges * (numedges-1) # all sets of edge1 -> edge2
    nattempts = 0
    while nattempts < maxattempts && length(blacklist) < nmax_blacklist
        e1 = Random.rand(1:numedges) # presumably faster than Random.randperm or Random.shuffle
        edge1 = net.edge[e1]
        e2 = Random.rand(1:(numedges-1)) # e2 must be different from e1: only numedges-1 options
        edge2 = net.edge[(e2<e1 ? e2 : e2+1)]
        (e1,e2) ∉ blacklist || # try another pair without adding to the # of attempts
            continue           # if (e1,e2) was already attempted
        nattempts += 1
        ## check that constraints are met
        p1 = getParent(edge1)
        p2 = getParent(edge2)
        constraintsmet = true
        for con in constraints
            if con.type == 2 # forbid going out of (edge1) or into (edge2) the species group
                if con.nodenum == p1.number || con.edgenum == p2.number
                    push!(blacklist, (e1,e2))
                    constraintsmet = false
                    break # of constraint loop
                end
            end
        end
        constraintsmet || continue # try another pair if constraints not met
        ## check for no 3-cycle, if no3cycle is requested
        if no3cycle && hybrid3cycle(edge1, edge2) # does not depend on which partner edge is chosen
            push!(blacklist, (e1,e2))
            continue
        end
        ## check for no hybrid ladder, if requested: edge2 cannot be a hybrid
        if nohybridladder && edge2.hybrid
            push!(blacklist, (e1,e2))
            continue
        end
        hybridpartnernew = (rand() > 0.2) # if true: partner hybrid = new edge above edge 2
        ## check that the new network will be a DAG: no directional conflict
        if net.numHybrids > 0 && directionalconflict(net, p1, edge2, hybridpartnernew)
            hybridpartnernew = !hybridpartnernew # try again with opposite
            if directionalconflict(net, p1, edge2, !hybridpartnernew)
                push!(blacklist, (e1,e2))
                continue
            end # else: switching hybridpartnernew worked
        end
        return addhybridedge!(net, edge1, edge2, hybridpartnernew)
    end
    # error("tried max number of attempts, none worked!")
    return nothing
end

"""
    addhybridedge!(net::HybridNetwork, edge1::Edge, edge2::Edge, hybridpartnernew::Bool)

Add hybridization to `net` coming from `edge1` going into `edge2`.
2 new nodes and 3 new edges are created: `edge1` are `edge2` are both cut into 2 edges,
and a new edge is created linking the 2 new "middle" nodes, pointing from `edge1` to `edge2`.
The new node in the middle of `edge1` is a tree node.
The new node in the middle of `edge2` is a hybrid node.
Its parent edges are the newly created hybrid edge (minor, with a random γ<0.5),
and either the newly edge "above" `edge2` if `hybridpartnernew=true`,
or the old `edge2` otherwise (which would reverse the direction of `edge2` and others).

Always called from the other function method, which performs a bunch of checks.
Updates `containRoot` attributes for edges below
new hybrid node if applicable.

`net` is modified and the new hybrid node (middle of the old `edge2`) is returned.

```jldoctest
julia> net = readTopology("(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));");

julia> hybnode = PhyloNetworks.addhybridedge!(net, net.edge[3], net.edge[9], true)
PhyloNetworks.Node:
number:12
hybrid node
attached to 3 edges, numbered: 9 22 23

julia> hybnode.edge
3-element Array{PhyloNetworks.Edge,1}:
PhyloNetworks.Edge:
number:9
length:-1.0
attached to 2 node(s) (parent first): 12 -7

PhyloNetworks.Edge:
number:22
length:-1.0
major hybrid edge with gamma=0.7768082985204838
attached to 2 node(s) (parent first): -6 12

PhyloNetworks.Edge:
number:23
length:0.0
minor hybrid edge with gamma=0.22319170147951617
attached to 2 node(s) (parent first): 11 12

julia> writeTopology(net)
"(((((#H1,(S6,S7)),(((S1,S4),(S5)#H1))#:::0.7768082985204838))#H2,((S8,S9),#:0.0::0.22319170147951617)),(#H2,S10));"

```
"""
function addhybridedge!(net::HybridNetwork, edge1::Edge, edge2::Edge, hybridpartnernew::Bool)
    newnode1_tree, edgeabovee1 = breakedge!(edge1, net) # new tree node
    newnode2_hybrid, edgeabovee2 = breakedge!(edge2, net) # new hybrid node
    newnode2_hybrid.hybrid = true
    pushHybrid!(net, newnode2_hybrid) # updates net.hybrid and net.numHybrids
    # new hybrid edge
    newgamma = rand()*0.5; # should be between 0 and 0.5 to create a minor hybrid edge
    hybrid_edge = Edge(maximum(e.number for e in net.edge) + 1, 0.0, true, newgamma, false) # number, length, hybrid, gamma, isMajor
    # partner edge: update hybrid status, γ and direction
    if hybridpartnernew
        edgeabovee2.hybrid = true
        edgeabovee2.gamma = 1.0 - newgamma
    else
        c2 = getChild(edge2) # child of edge2 before we switch its direction
        i2 = findfirst(isequal(c2), net.node)
        net.root = i2 # makes c2 the new root node
        edge2.hybrid = true
        edge2.gamma = 1.0 - newgamma
        edge2.isChild1 =  !edge2.isChild1 # reverse the direction of edge2
        edgeabovee2.isChild1 = !edgeabovee2.isChild1
    end
    # parse hyb names to find the next available. assignhybridnames! would do them all
    rx = r"^H(\d+)$"
    hnum = net.numHybrids # to name the new hybrid, potentially
    for n in net.node
        m = match(rx, n.name)
        if m !== nothing
            hi = parse(Int, m[1])
            hnum > hi || (hnum = hi+1)
        end
    end
    newnode2_hybrid.name = "H$hnum"
    setNode!(hybrid_edge, [newnode2_hybrid, newnode1_tree]) # [child node, parent node] to match isChild1=true
    setEdge!(newnode1_tree, hybrid_edge)
    setEdge!(newnode2_hybrid, hybrid_edge)
    pushEdge!(net, hybrid_edge)
    if hybridpartnernew
        # @debug "triple-check it's a DAG" directEdges!(net)
        norootbelow!(edge2)
    else
        directEdges!(net)
        norootbelow!(edgeabovee2)
    end
    return newnode2_hybrid
end

"""
    hybrid3cycle(edge1::Edge, edge2::Edge)

Check if proposed hybrid edge from `edge1` into `edge2` would create a 3 cycle,
that is, if `edge1` and `edge2` have a node in common.
(This move cannot create a 2-cycles because new nodes would be created in the
middle of edges 1 and 2.)
"""
function hybrid3cycle(edge1::Edge, edge2::Edge)
    !isempty(findall(in(edge1.node), edge2.node))
end

"""
    directionalconflict(net::HybridNetwork, parent::Node, edge::Edge,
                        hybridpartnernew::Bool)

Check if creating a hybrid edge down of `parent` node into the middle of `edge`
would create a directed cycle in `net`, i.e. not a DAG. The proposed hybrid
would go in the direction of `edge` down its child node if `hybridpartnernew`
is true. Otherwise, both halves of `edge` would have their direction reversed,
for the hybrid to go towards the original parent node of `edge`.
Does *not* modify the network.

Output: `true` if a conflict would arise (non-DAG), `false` if no conflict.
"""
function directionalconflict(net::HybridNetwork, parent::Node, edge2::Edge, hybridpartnernew::Bool)
    if hybridpartnernew # all edges would retain theirs directions: use isChild1 fields
        c2 = getChild(edge2)
        return parent === c2 || isdescendant(parent, c2)
    else # after hybrid addition, edge 2 would be reversed: "up" toward its own parent
        if !edge2.containRoot || edge2.hybrid # direction of edge2 cannot be reversed
            return true
        else # net would be a DAG with reversed directions, could even be rooted on edge2
            p2 = getParent(edge2)
            return parent === p2 || isdescendant_undirected(parent, p2, edge2)
        end
    end
end

"""
    isdirectionaldescendant(des::Node, anc::Node)

Check if, following only directional edges, `anc` node flows into `des` node.
"""
function isdirectionaldescendant(des::Node, anc::Node)
    visited = Int[]
    for e in anc.edge
        if isdirectionaldescendant!(visited, des, e)
            return true
        end
    end
    return false
end

"""
    isdirectionaldescendant!(visited::Vector{Int}, des::Node, e::Edge)

Return true if `des` node is directional descendant of `anc` node.
Does not use isChild1 attribute because isChild takes rooting as fixed
"""
function isdirectionaldescendant!(visited::Vector{Int}, des::Node, e::Edge)
    for n in e.node #check all nodes of e
        if n == des
            return true
        end
        if n.hybrid # only need to check this for hybrid nodes
            if n.number in visited
                return false  # n was already visited: exit.
            end
            push!(visited, n.number)
        end
        for ce in n.edge
            if ce == e # avoids stack overflow error caused by calling with same arguments
                return false # because e already checked
            elseif !ce.containRoot && n in ce.node # replaces n == getParent(ce) to avoid isChild1
                if isdirectionaldescendant!(visited, des, ce) return true; end
            end
        end
        return false
    end
end
