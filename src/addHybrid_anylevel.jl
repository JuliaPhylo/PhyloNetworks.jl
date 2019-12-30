# functions to add hybridzation edges during optimization in all level networks

# add hybrid alternated with nni moves (see moves_semidirected.jl) as
# part of discrete traits full likelihood network optimization.

# uses parts of addHybrid.jl

"""
    addhybridedge!(net::HybridNetwork)

Randomly choose two edge indices. If they pass constraint and 3- cycle checks,
adds hybrid edge from edge 1 to edge 2 by calling next version of `adddhybridedge!`. 

If successful, return net, newhybridnode in tuple.
If not, return nothing. #? Should return nothing or an error?

semi-directed version modeled after approach in addHybrid.jl)
#? should we include an optional edgelength parameter here? 
#? Or perhaps its better to update branch lengths later after addhybrid step of optimization?
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

        if rand() > 0.5
            hybridpartnernew = true # the partner hybrid will be the new edge above edge 2 (returned by addnodeonedge)
        else
            hybridpartnernew = false #the partner hybrid will be existing edge "below" edge 2
        end
    addhybridedge!(net, edge1, edge2, hybridpartnernew)
end

"""
    addhybridedge!(net::HybridNetwork, edge1::Edge, edge2::Edge, hybridpartnernew::Bool)

Helper function for addhybridedge(net::HybridNework). Always called from above.
If needed, updates containRoot attribute below new hybrid node.

`hybridpartnernew` boolean indicates which half of `edge2` will be the partner hybrid edge 
to the newly-added hybrid edge. This partner hybrid edge will point toward the 
newly-created node on the middle of the original `edge2`.

If successful, return net, newhybridnode in tuple.

```jldoctest
julia> str_species_net = "(((S8,S9),((((S1,S4),(S5)#H1),(#H1,(S6,S7))))#H2),(#H2,S10));";

julia> species_net = readTopology(str_species_net);

julia> species_net, newhybridnode = PhyloNetworks.addhybridedge!(species_net, species_net.edge[3], species_net.edge[9], true);

julia> newhybridnode
PhyloNetworks.Node:
 number:12
 hybrid node
 attached to 3 edges, numbered: 9 22 23

 julia> newhybridnode.edge
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
 
 
 julia> writeTopology(species_net)
 "(((((#H1,(S6,S7)),(((S1,S4),(S5)#H1))#:::0.7768082985204838))#H2,((S8,S9),#:0.0::0.22319170147951617)),(#H2,S10));"

```
"""
function addhybridedge!(net::HybridNetwork, edge1::Edge, edge2::Edge, hybridpartnernew::Bool)
    ## Check 3: Hybrid would not create directional conflict
    if hybrid3cycle(net, edge1, edge2)
        error("hybrid between edge $(edge1.number) and edge $(edge2.number) would create 3-cycle.")
        #return nothing
    ## Check 4: Hybrid would not create a 3-cycle
    elseif net.numHybrids > 0 && directionalconflict(net, edge1, edge2, hybridpartnernew) 
        error("directional conflict: edge $(edge1.number) is a directional descendant of edge $(edge2.number).")
        #return nothing 
    else # add hybridization
        newnode1_tree, edgeabovee1 = addnodeonedge!(net, edge1, false) # new tree node
        newnode2_hybrid, edgeabovee2 = addnodeonedge!(net, edge2, true) # new hybrid node
        # new hybrid edge
        newgamma = rand()*0.5; #gamma should be between 0 and 0.5
        hybrid_edge = Edge(maximum(e.number for e in net.edge) + 1, 0.0, true, newgamma, newgamma>=0.5) # isChild1 = true by default in edge creation        
        #update node and edgehybrid statuses
        newnode2_hybrid.hybrid = true
        hybrid_edge.hybrid = true
        if hybridpartnernew
            edgeabovee2.hybrid = true
        else
            edge2.hybrid = true
        end
        # update edge gammas 
        edgeabovee2.gamma = 1 - newgamma
        updateContainRoot!(net, newnode2_hybrid);
        setNode!(hybrid_edge, [newnode1_tree, newnode2_hybrid])
        setEdge!(newnode1_tree, hybrid_edge)
        setEdge!(newnode2_hybrid, hybrid_edge)
        pushEdge!(net, hybrid_edge)
        
        # updateInCycle!(net, newnode2); #? need this?
        # updateMajorHybrid!(net, newnode2); #? need this?
        return net, newnode2_hybrid
    end
end

"""
    hybrid3cycle(net::HybridNetwork, edge1::Edge, edge2::Edge)

Check if proposed hybrid edge would create a 3 cycle. (Because we create new
nodes, this move cannot create a 2-cycle.)
"""
function hybrid3cycle(net::HybridNetwork, edge1::Edge, edge2::Edge)
    # if node lists overlap, they share a node, so adding a hybrid between
    # them would create a 3-cycle.
    if !isempty(findall(in(edge1.node),edge2.node))
        return true
    else
        return false
    end
end

"""
    directionalconflict(net::HybridNetwork, edge1::Edge, edge2::Edge, hybridpartnernew::Bool)

Check if proposed hybrid edge would create a non-DAG.

Answers question: Is edge 1 a descendant of the newly-created node on edge 2?

"""
function directionalconflict(net::HybridNetwork, edge1::Edge, edge2::Edge, hybridpartnernew::Bool)
    # (see "Checking for Directional Conflicts")
    if hybridpartnernew # after hybrid addition, edge 2 would flow toward isChild1(edge2)
        if !edge1.containRoot && isdirectionaldescendant(getParent(edge1), getChild(edge2))
            # we can trust getChild1 for edge1 because !edge1.containRoot #TODO confirm this
            return true
        else
            return false
        end
    else # after hybrid addition, edge 2 would flow toward getParent(edge2)
        if !edge1.containRoot && isdirectionaldescendant(getParent(edge1), getParent(edge2)) 
            # we can trust getChild1 for edge1 because !edge1.containRoot #TODO confirm this
            return true
        else
            return false
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
