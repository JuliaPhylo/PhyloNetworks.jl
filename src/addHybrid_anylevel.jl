# functions to add hybridzation edges during optimization in all level networks

# add hybrid alternated with nni moves (see moves_semidirected.jl) as
# part of discrete traits full likelihood network optimization.

# uses parts of addHybrid.jl

"""
    addhybridedge!(net::HybridNetwork)

Randomly choose two edge indices. If they pass 2- and 3- cycle checks,
adds hybrid edge from edge 1 to edge 2. If needed, update containRoot attribute 
below new hybrid node.

If successful, return net, hybridedge, newhybridnode
If not, return nothing.

```jldoctest
julia> str_species_net = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));";

julia> species_net = readTopology(str_species_net);

julia> species_net, newhybridnode = PhyloNetworks.addhybridedge!(net_level1, net_level1.edge[3], net_level1.edge[9]);

julia> newhybridnode #TODO

julia> newhybridnode.edge #TODO

julia> writeTopology(species_net) #TODO

```
"""
function addhybridedge!(net::HybridNetwork)
    e1, e2 = Random.randperm(length(net.edge))[1:2] # randomly chooses two edges without replacement
    edge1 = net.edge[e1]
    edge2 = net.edge[e2]
    addhybridedge!(net, edge1, edge2)
end

function addhybridedge!(net::HybridNetwork, edge1::Edge, edge2::Edge)
    if (!edge1.containRoot && !edge2.containRoot) && (edge1 in descendants(edge2))
        return nothing # directional conflict: e1 is directed and e1 is a descendant of e2
    elseif hybridwouldcreate3cycle(net, edge1, edge2)
        return nothing
    else # add hybrid edge
        gamma = rand()*0.5; #? Cécile, should we use chooseEdgesGamma() function here? Or is this sufficient since we're not using a blacklist?
        edge3, edge4 = parameters4createHybrid!(edge1, edge2, net)
        newhybridnode = createHybrid!(edge1, edge2, edge3, edge4, net, gamma)
        updateInCycle!(net, newhybridnode);
        updateMajorHybrid!(net, newhybridnode);
        updateContainRoot!(net, newhybridnode);
        #updateAllNewHybrid!(newhybridnode, net, true, true, false) 
            # has parts we dont need/seem broken? ERROR: no method matching getOtherNode(::Nothing, ::PhyloNetworks.Node)
            # booleans: updatemajor (bool) to decide if we need to update major edge
            # allow = true allows extreme/very bad triangles, needed when reading
            # updatePart = true will update PArtition at this moment, it makes sense with a newly added hybrid
        #? Cécile, are there other update functions from addHybrid.jl we should use?
        return net, newhybridnode
    end
end

"""
    hybrid3cycle(net::HybridNetwork, e1::Edge, e2::Edge)

Check if proposed hybrid edge would create a 3 cycle. (Because we create new
nodes, this move cannot create a 2-cycle.)
"""
function hybridwouldcreate3cycle(net::HybridNetwork, e1::Edge, e2::Edge)
    # if node lists overlap, they share a node, so adding a hybrid between
    # them would create a 3-cycle.
    if !isempty(findall(in(e1.node),e2.node))
        return true
    else
        return false
    end
end

# future todos: create deleteHybrid_anylevel.jl based on deleteHybrid.jl
