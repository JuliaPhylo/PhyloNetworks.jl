# functions to add hybridizations during optimization of an any-level 
# semi-directed network

# addhybridedge! alternated with nni moves (see moves_semidirected.jl) and local
# optimization as part of discrete traits full likelihood network optimization.

"""
    addhybridedge!(net::HybridNetwork, treechild::Bool, constraints=TopologyConstraint[]::Vector{TopologyConstraint})

Randomly choose two edge indices. If they pass constraint and 3- cycle checks,
adds hybrid edge from edge 1 to edge 2 by calling next version of `adddhybridedge!`. 

Randomly decides whether the new the partner hybrid will be the new edge above `edge2`
or the partner hybrid will be existing edge below `edge 2`. This partner hybrid 
edge will point toward the newly-created node on the middle of the original `edge2`.

When testing the viability of an edge, it tries both `hybridpartnernew` equal to 
`true` and equal to `false`.

A version of addhybridization for semi-directed networks, modeled after 
approach in addHybrid.jl.

If successful, return net, newhybridnode in tuple.
If not, adds edge pairs to the blacklist and tries with a new set of edges.
"""
function addhybridedge!(net::HybridNetwork, treechild, constraints=TopologyConstraint[]::Vector{TopologyConstraint})
    edgesfound = false
    Array{Array{Int64,1},1}
    blacklist = Array{Array{Edge, 1}, 1}() # ordered Dict{Edge1,Edge2}
    while !edgesfound 
        if length(blacklist) == Int(factorial(length(net.edge))/factorial(length(net.edge)-2)) # all permutations
            break
        end
        e1, e2 = Random.randperm(length(net.edge))[1:2] # randomly chooses two edges without replacement
            # fixit: if we could restrict to only interior edges, this while loop would be much faster
        global edge1 = net.edge[e1]
        global edge2 = net.edge[e2]
        if !([edge1, edge2] in blacklist) # edges not in blacklist
            edgesfound = true
            ## Check 1: Check that these edges are not the top of a clade or species group
            for con in constraints
                if con.edgenum == edge1.number || con.edgenum == edge2.number
                    edgesfound = false
                end
            end
            global hybridpartnernew = (rand() > 0.5) # if true the partner hybrid will be new edge above edge 2
            ## Check 2: Check that the edges are interior
            if any([n.leaf for n in edge1.node]) || any([n.leaf for n in edge2.node])
                edgesfound = false
            ## Check 3: Hybrid would not create a 3-cycle 
            elseif hybrid3cycle(net, edge1, edge2)
                edgesfound = false
            ## Check 4: if treechild, edge2 cannot be a hybrid
            elseif treechild && edge2.hybrid
                edgesfound = false
            ## Check 5: Hybrid would not create directional conflict
            elseif net.numHybrids > 0 && directionalconflict(net, edge1, edge2, hybridpartnernew) 
                hybridpartnernew = !hybridpartnernew # try again with opposite
                if directionalconflict(net, edge1, edge2, hybridpartnernew)
                    edgesfound = false
                end
            end
            if !edgesfound # add edges to blacklist
                push!(blacklist, [edge1, edge2])
            end
        end
    end
    addhybridedge!(net, edge1, edge2, hybridpartnernew)
end

"""
    addhybridedge!(net::HybridNetwork, edge1::Edge, edge2::Edge, hybridpartnernew::Bool)

Adds hybridization to `net`. Helper function for `addhybridedge(net::HybridNework)`. 
Always called from above. Updates `containRoot` attributes for edges below 
new hybrid node if applicable.

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
    newnode2_hybrid.name = "H#$(net.numHybrids + 1)"
    # update edge gammas 
    edgeabovee2.gamma = 1 - newgamma
    updateContainRoot!(net, newnode2_hybrid);
    setNode!(hybrid_edge, [newnode1_tree, newnode2_hybrid])
    setEdge!(newnode1_tree, hybrid_edge)
    setEdge!(newnode2_hybrid, hybrid_edge)
    pushEdge!(net, hybrid_edge)
    #TODO need to push node to network?
        
    # updateInCycle!(net, newnode2); #? need this?
    # updateMajorHybrid!(net, newnode2); #? need this?
    return net, newnode2_hybrid
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
Asks if edge 1 a descendant of the newly-created node on edge 2.
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

"""
    removehybridedge!(net::HybridNetwork, node::Node, minor::Bool, blacklist::Bool)

Remove hybridization from network. Adapted from deleteHybrid in deleteHybrid.jl. 
Do not update branch lengths. 

If `minor = true`, deletes minor edge, else deletes major edge.

If `blacklist`, add an edge a blacklist.

Return net.

```jldoctest
#TODO
````
"""
function removehybridedge!(net::HybridNetwork, node::Node,minor::Bool, blacklist::Bool)
    node.hybrid || error("node $(node.number) has to be hybrid for deleteHybrid")
    if(minor)
        hybedge1,hybedge2,treeedge1 = hybridEdges(node);
        other1 = getOtherNode(hybedge1,node);
        other2 = getOtherNode(hybedge2,node);
        other3 =  getOtherNode(treeedge1,node);
        if(hybedge1.number > treeedge1.number)
            removeNode!(node,treeedge1);
            setNode!(treeedge1,other1);
            setEdge!(other1,treeedge1);
            removeEdge!(other1, hybedge1);
            deleteEdge!(net,hybedge1);
            #treeedge1.containRoot = (!treeedge1.containRoot || !hybedge1.containRoot) ? false : true #causes problems if hybrid.CR=false
            if(blacklist)
                println("put in blacklist edge $(treeedge1.number)")
                push!(net.blacklist, treeedge1.number)
            end
        else
            makeEdgeTree!(hybedge1,node)
            other1.hasHybEdge = false;
            removeNode!(node,hybedge1);
            setNode!(hybedge1,other3);
            setEdge!(other3,hybedge1);
            removeEdge!(other3,treeedge1);
            deleteEdge!(net,treeedge1);
            hybedge1.containRoot = (!treeedge1.containRoot || !hybedge1.containRoot) ? false : true
            if(blacklist)
                println("put in blacklist edge $(hybedge1.number)")
                push!(net.blacklist, hybedge1.number)
            end
        end
        hybindex = findfirst([e.hybrid for e in other2.edge]);
        hybindex != nothing || error("didn't find hybrid edge in other2")
        if(hybindex == 1)
            treeedge1 = other2.edge[2];
            treeedge2 = other2.edge[3];
        elseif(hybindex == 2)
            treeedge1 = other2.edge[1];
            treeedge2 = other2.edge[3];
        elseif(hybindex == 3)
            treeedge1 = other2.edge[1];
            treeedge2 = other2.edge[2];
        else
            error("strange node has more than three edges")
        end
        treenode1 = getOtherNode(treeedge1,other2);
        treenode2 = getOtherNode(treeedge2,other2);
        if(abs(treeedge1.number) > abs(treeedge2.number))
            removeNode!(other2,treeedge2);
            setNode!(treeedge2,treenode1);
            setEdge!(treenode1,treeedge2);
            removeEdge!(treenode1,treeedge1);
            deleteEdge!(net,treeedge1);
            treeedge2.containRoot = (!treeedge1.containRoot || !treeedge2.containRoot) ? false : true
            if(blacklist)
                println("put in blacklist edge $(treeedge2.number)")
                push!(net.blacklist, treeedge2.number)
            end
        else
            removeNode!(other2,treeedge1);
            setNode!(treeedge1,treenode2);
            setEdge!(treenode2,treeedge1);
            removeEdge!(treenode2,treeedge2);
            deleteEdge!(net,treeedge2);
            treeedge1.containRoot = (!treeedge1.containRoot || !treeedge2.containRoot) ? false : true
            if(blacklist)
                println("put in blacklist edge $(treeedge1.number)")
                push!(net.blacklist, treeedge1.number)
            end
        end
        #removeHybrid!(net,node);
        deleteNode!(net,node);
        deleteNode!(net,other2);
        deleteEdge!(net,hybedge2);
    else
        hybedge1,hybedge2,treeedge1 = hybridEdges(node);
        other1 = getOtherNode(hybedge1,node);
        other2 = getOtherNode(hybedge2,node);
        removeEdge!(other2,hybedge2)
        removeNode!(node,treeedge1)
        setEdge!(other2,treeedge1)
        setNode!(treeedge1,other2)
        #removeHybrid!(net,node)
        deleteNode!(net,node)
        deleteEdge!(net,hybedge1)
        deleteEdge!(net,hybedge2)
        removeEdge!(other1,hybedge1)
        size(other1.edge,1) == 2 || error("strange node $(other1.number) had 4 edges")
        if(abs(other1.edge[1].number) < abs(other1.edge[2].number))
            edge = other1.edge[1]
            otheredge = other1.edge[2]
        else
            edge = other1.edge[2]
            otheredge = other1.edge[1]
        end
        other3 =  getOtherNode(otheredge,other1);
        removeNode!(other1,edge)
        removeEdge!(other3,otheredge)
        setEdge!(other3,edge)
        setNode!(edge,other3)
        deleteNode!(net,other1)
        deleteEdge!(net,otheredge)
    end
end
#TODOs
# write a function to subset to only interior edges of a network
# test using deleteHybrid from deleteHybrid.jl to remove a hybridization
