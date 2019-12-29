# functions to add hybridzation edges during optimization in all level networks

# add hybrid alternated with nni moves (see moves_semidirected.jl) as
# part of discrete traits full likelihood network optimization.

# uses parts of addHybrid.jl

"""
    addhybridedge!(net::HybridNetwork)

Randomly choose two edge indices. If they pass constraint and 3- cycle checks,
adds hybrid edge from edge 1 to edge 2. If needed, update containRoot attribute 
below new hybrid node.

If successful, return net, hybridedge, newhybridnode
If not, return nothing.
#? should we include an optional edgelength parameter here?
```jldoctest
julia> str_species_net = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));";

julia> species_net = readTopology(str_species_net);

julia> species_net, newhybridnode = PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[3], net_level1.edge[9]);

julia> newhybridnode #TODO

julia> newhybridnode.edge #TODO

julia> writeTopology(species_net) #TODO

```
"""
function addhybridedge!(net::HybridNetwork, constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    edgesfound = false
    while !edgesfound 
        e1, e2 = Random.randperm(length(net.edge))[1:2] # randomly chooses two edges without replacement
        edge1 = net.edge[e1]
        edge2 = net.edge[e2]
        edgesfound = true
        ## Check 1: Check that these edges are not the top of a clade or species group
        for con in constraints
            if con.edgenum == edge1.number || con.edgenum == edge2.number
                edgesfound = false
            end
        end
        ## Check 2: Check that the edges are interior
        if edge1.leaf || edge2.leaf
            edgesfound = false
        end
    end
    addhybridedge!(net, edge1, edge2)
end

function addhybridedge!(net::HybridNetwork, edge1::Edge, edge2::Edge)
    ## Check 3: Hybrid would not create directional conflict
    if net.numHybrids > 0 && !directionalconflict(net, edge1, edge2) 
        error("directional conflict: $edge1 is a directional descendant of $edge2.")
        #return nothing # directional conflict: e1 is directed and e1 is a descendant of e2
    ## Check 4: Hybrid would not create a 3-cycle
    elseif hybridwouldcreate3cycle(net, edge1, edge2)
        error("hybrid between $edge1 and $edge2 would create 3-cycle.")
        #return nothing
    else # add hybridization
        newnode1, edgeabovee1 = addnodeonedge!(edge1) #new nodes
        newnode2, edgeabovee2 = addnodeonedge!(edge2)
        # new hybrid edge
        newedge = Edge(maximum(e.number for e in net.edge) + 1) # isChild1 = true by default in edge creation        # connect new edge to node 1, node 2
        setNode!(newedge, [newnode1, newnode2])
        setEdge!(newnode1, newedge)
        setEdge!(newnode2, newedge)
        pushEdge!(net, newedge)

        #update node and edgehybrid statuses
        newnode2.hybrid = true
        newedge.hybrid = true
        edgeabovee2.hybrid = true

        # update edge gammas #? should we do something more sophisticated here?
        newedge.gamma = rand()*0.5;
        edgeabovee2.gamma = 1-newedge.gamma
        updateContainRoot!(net, newnode2);
        # updateInCycle!(net, newnode2); #? need this?
        # updateMajorHybrid!(net, newnode2); #? need this?
        return net, newnode2
    end
end

"""
    hybrid3cycle(net::HybridNetwork, edge1::Edge, edge2::Edge)

Check if proposed hybrid edge would create a 3 cycle. (Because we create new
nodes, this move cannot create a 2-cycle.)
"""
function hybridwouldcreate3cycle(net::HybridNetwork, edge1::Edge, edge2::Edge)
    # if node lists overlap, they share a node, so adding a hybrid between
    # them would create a 3-cycle.
    if !isempty(findall(in(edge1.node),edge2.node))
        return true
    else
        return false
    end
end

"""
    directionalconflict(net::HybridNetwork, edge1::Edge, edge2::Edge)

Check if proposed hybrid edge would create a non-DAG. Uses `isChild1` attribute.
"""
function directionalconflict(net::HybridNetwork, edge1::Edge, edge2::Edge)
    if !edge1.containRoot
        if isdirectionaldescendant(getChild(edge1), getChild(edge2))
            # edge 2 is a directional ancestor of edge 1 OR
            # edge 2 is a hybrid & this hybrid flows into E1 
            # (see "Checking for Directional Conflicts" case 2, case 6)
            return true
        elseif any([n.hybrid for n in PhyloNetworks.getChildren(PhyloNetworks.getChild(edge2))]) &&
            isdirectionaldescendant(getChild(edge1), getChild(edge2))
            # edge 2 is the parent of a hybrid & this hybrid flows into E1
            #? Do we want to allow this case? (see "Checking for Directional Conflicts" case 3, case 5)
            return true
        end
    else
        return false
    end
end

"""
    isdirectionaldescendant(des::Node, anc::Node)

Check if, following only directional edges, `anc` node flows into `des` node.
Uses `isChild1` attribute.
"""
function isdirectionaldescendant(des::Node, anc::Node)
    visited = Int[]
    for e in anc.edge
        if !e.containRoot && isdirectionaldescendant!(visited, des, e)
            return true
        end
    end
    return false
end

"""
    isdirectionaldescendant!(visited::Vector{Int}, des::Node, e::Edge)

Return true if `des` node is directional descendant of `anc` node. 
Uses `isChild1` attribute.
"""
function isdirectionaldescendant!(visited::Vector{Int}, des::Node, e::Edge)
    n = getChild(e)
    if n == des
        return true
    end
    if n.hybrid # only need to check this for hybrid nodes
        push!(visited, n.number)
        if n.number in visited
            return false  # n was already visited: exit. avoid infinite loop is isChild1 was bad.
        end
    end
    for ce in n.edge
        if !ce.containRoot && n == getParent(ce)
            if isdirectionaldescendant!(visited, des, ce) return true; end
        end
    end
    return false
end

# future todos: create deleteHybridedge.jl based on deleteHybrid.jl
