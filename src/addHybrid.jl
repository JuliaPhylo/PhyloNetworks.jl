# functions to add hybridizations to a semi-directed network
# subject to topological constraints; no level restriction

"""
    addhybridedge!(
        [rng::AbstractRNG,]
        net::HybridNetwork,
        nohybridladder::Bool,
        no3cycle::Bool,
        constraints::Vector{TopologyConstraint}=TopologyConstraint[];
        maxattempts::Int=10,
        fixroot::Bool=false
    )

Randomly choose two edges in `net` then: add hybrid edge from edge 1 to edge 2
of length 0.01. The two halves of edge 1 (and of edge 2) have equal lengths.
The hybrid partner edge (top half of edge 2, if fixroot is true) will point
towards the newly-created node on the middle of the original edge 2,
and have an inheritance γ randomly uniformly chosen in (0,0.5).

If the resulting network is a DAG, satisfies the constraint(s),
does not contain any 3-cycle (if `no3cycle=true`), and does not have
a hybrid ladder (if `nohybridladder=true`) then the proposal is successful:
`net` is modified, and the function returns the newly created hybrid node and
newly created hybrid edge.

If the resulting network is not acceptable, then a new set of edges
is proposed (using a blacklist) until one is found acceptable, or until
a maximum number of attempts have been made (`maxattempts`).
If none of the attempted proposals are successful, `nothing` is returned
(without causing an error).

After a pair of edges is picked, the "top" half of edge2 is proposed
as the partner hybrid edge with probability 0.8 if `fixroot` is false,
(to avoid changing the direction of edge2 with more than 50% chance)
and with probability 1.0 if `fixroot` is true.
If this choice does not work and if `fixroot` is false,
the other half of edge2 is proposed as the partner hybrid edge.
Note that choosing the "bottom" half of edge2 as the partner edge
requires to flip the direction of edge 2, and to move the root accordingly
(to the original child of edge2).

# examples

```jldoctest
julia> net = readTopology("((S1,(((S2,(S3)#H1),(#H1,S4)))#H2),(#H2,S5));");

julia> using Random

julia> Random.seed!(170);

julia> PhyloNetworks.addhybridedge!(net, true, true)
(PhyloNetworks.Node:
 number:9
 name:H3
 hybrid node
 attached to 3 edges, numbered: 5 16 17
, PhyloNetworks.EdgeT{PhyloNetworks.Node}:
 number:17
 length:0.01
 minor hybrid edge with gamma=0.32771460911632916
 attached to 2 node(s) (parent first): 8 9
)

julia> writeTopology(net, round=true, digits=2)
"((S1,(((#H1,S4),((S2,(S3)#H1))#H3:::0.67))#H2),((#H2,S5),#H3:0.01::0.33));"

```
"""
function addhybridedge!(net::HybridNetwork, nolad::Bool, args...; kwargs...)
    addhybridedge!(Random.default_rng(), net, nolad, args...; kwargs...)
end
function addhybridedge!(
    rng::Random.AbstractRNG,
    net::HybridNetwork,
    nohybridladder::Bool,
    no3cycle::Bool,
    constraints::Vector{TopologyConstraint}=TopologyConstraint[];
    maxattempts::Int=10,
    fixroot::Bool=false
)
    all(con.type == 1 for con in constraints) || error("only type-1 constraints implemented so far")
    numedges = length(net.edge)
    blacklist = Set{Tuple{Int,Int}}()
    nmax_blacklist = numedges * (numedges-1) # all sets of edge1 -> edge2
    nattempts = 0
    while nattempts < maxattempts && length(blacklist) < nmax_blacklist
        e1 = Random.rand(rng, 1:numedges) # presumably faster than Random.randperm or Random.shuffle
        edge1 = net.edge[e1]
        e2 = Random.rand(rng, 1:(numedges-1)) # e2 must be different from e1: only numedges-1 options
        edge2 = net.edge[(e2<e1 ? e2 : e2+1)]
        (e1,e2) ∉ blacklist || # try another pair without adding to the # of attempts
            continue           # if (e1,e2) was already attempted
        nattempts += 1
        ## check that constraints are met
        p1 = getparent(edge1)
        p2 = getparent(edge2)
        constraintsmet = true
        for con in constraints
            if con.type == 1 # forbid going out of (edge1) or into (edge2) the species group
                if con.node === p1 || con.node === p2
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
        ## check for no hybrid ladder, if requested:
         # edge2 cannot be a hybrid edge or the child of a hybrid node
        if nohybridladder && (edge2.hybrid || getparent(edge2).hybrid)
            push!(blacklist, (e1,e2))
            continue
        end
        hybridpartnernew = (fixroot ? true : rand(rng) > 0.2) # if true: partner hybrid = new edge above edge 2
        ## check that the new network will be a DAG: no directional conflict
        if directionalconflict(p1, edge2, hybridpartnernew)
            if fixroot # don't try to change the direction of edge2
                push!(blacklist, (e1,e2))
                continue
            end # else: try harder: change direction of edge2 and move root
            hybridpartnernew = !hybridpartnernew # try again with opposite
            if directionalconflict(p1, edge2, hybridpartnernew)
                push!(blacklist, (e1,e2))
                continue
            end # else: switching hybridpartnernew worked
        end
        newgamma = rand(rng)/2; # in (0,.5) to create a minor hybrid edge
        return addhybridedge!(net, edge1, edge2, hybridpartnernew, 0.01, newgamma)
    end
    # if we get here: none of the max number of attempts worked - return nothing.
    # if an attempt had worked, we would have returned something.
    return nothing
end

"""
    addhybridedge!(
        net::HybridNetwork,
        edge1::Edge,
        edge2::Edge,
        hybridpartnernew::Bool,
        edgelength::Float64=-1.0,
        gamma::Float64=-1.0
    )

Add hybridization to `net` coming from `edge1` going into `edge2`.
2 new nodes and 3 new edges are created: `edge1` are `edge2` are both cut into 2 edges,
and a new edge is created linking the 2 new "middle" nodes, pointing from `edge1` to `edge2`.
The new node in the middle of `edge1` is a tree node.
The new node in the middle of `edge2` is a hybrid node.
Its parent edges are the newly created hybrid edge (with γ = gamma, missing by default),
and either the newly edge "above" `edge2` if `hybridpartnernew=true`,
or the old `edge2` otherwise (which would reverse the direction of `edge2` and others).

Should be called from the other method, which performs a bunch of checks.
Updates `containRoot` attributes for edges below the new hybrid node.

Output: new hybrid node (middle of the old `edge2`) and new hybrid edge.

# examples

```jldoctest
julia> net = readTopology("((S8,(((S1,(S5)#H1),(#H1,S6)))#H2),(#H2,S10));");

julia> hybnode, hybedge = PhyloNetworks.addhybridedge!(net, net.edge[13], net.edge[8], true, 0.0, 0.2)
(PhyloNetworks.Node:
 number:9
 name:H3
 hybrid node
 attached to 3 edges, numbered: 8 16 17
, PhyloNetworks.EdgeT{PhyloNetworks.Node}:
 number:17
 length:0.0
 minor hybrid edge with gamma=0.2
 attached to 2 node(s) (parent first): 8 9
)


julia> writeTopology(net)
"((S8,(((S1,(S5)#H1),((#H1,S6))#H3:::0.8))#H2),(#H2,(S10,#H3:0.0::0.2)));"

```
"""
function addhybridedge!(
    net::HybridNetwork,
    edge1::Edge,
    edge2::Edge,
    hybridpartnernew::Bool,
    edgelength::Float64=-1.0,
    gamma::Float64=-1.0
)
    gamma == -1.0 || (gamma <= 1.0 && gamma >= 0.0) || error("invalid γ to add a hybrid edge")
    gbar = (gamma == -1.0 ? -1.0 : 1.0 - gamma) # 1-gamma, with γ=-1 as missing
    newnode1_tree, edgeabovee1 = breakedge!(edge1, net) # new tree node
    newnode2_hybrid, edgeabovee2 = breakedge!(edge2, net) # new hybrid node
    newnode2_hybrid.hybrid = true
    pushHybrid!(net, newnode2_hybrid) # updates net.hybrid and net.numHybrids
    # new hybrid edge, minor if γ missing (-1)
    hybrid_edge = Edge(maximum(e.number for e in net.edge) + 1, edgelength, true, gamma, gamma>0.5) # number, length, hybrid, gamma, isMajor
    # partner edge: update hybrid status, γ and direction
    if hybridpartnernew
        edgeabovee2.hybrid = true
        edgeabovee2.gamma = gbar
        if gamma>0.5
            edgeabovee2.isMajor = false
        end
    else
        c2 = getchild(edge2) # child of edge2 before we switch its direction
        i2 = findfirst(isequal(c2), net.node)
        net.root = i2 # makes c2 the new root node
        edge2.hybrid = true
        edge2.gamma = gbar
        if gamma>0.5
            edge2.isMajor = false
        end
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
    return newnode2_hybrid, hybrid_edge
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
    directionalconflict(parent::Node, edge::Edge, hybridpartnernew::Bool)

Check if creating a hybrid edge down of `parent` node into the middle of `edge`
would create a directed cycle in `net`, i.e. not a DAG. The proposed hybrid
would go in the direction of `edge` down its child node if `hybridpartnernew`
is true. Otherwise, both halves of `edge` would have their direction reversed,
for the hybrid to go towards the original parent node of `edge`.
Does *not* modify the network.

Output: `true` if a conflict would arise (non-DAG), `false` if no conflict.
"""
function directionalconflict(parent::Node, edge2::Edge, hybridpartnernew::Bool)
    if hybridpartnernew # all edges would retain their directions: use isChild1 fields
        c2 = getchild(edge2)
        return parent === c2 || isdescendant(parent, c2)
    else # after hybrid addition, edge 2 would be reversed: "up" toward its own parent
        if !edge2.containRoot || edge2.hybrid
            return true # direction of edge2 cannot be reversed
        else # net would be a DAG with reversed directions, could even be rooted on edge2
            p2 = getparent(edge2)
            return parent === p2 || isdescendant_undirected(parent, p2, edge2)
        end
    end
end


# --- add alternative hybridizations found in bootstrap
"""
    addAlternativeHybridizations!(net::HybridNetwork, BSe::DataFrame;
                                  cutoff::Number=10, top::Int=3)

Modify the network `net` (the best estimated network) by adding some of
the hybridizations present in the bootstrap networks. By default, it will only
add hybrid edges with more than 10% bootstrap support (`cutoff`) and it will
only include the top 3 hybridizations (`top`) sorted by bootstrap support.

The dataframe `BSe` is also modified. In the original `BSe`,
supposedly obtained with `hybridclades_support`, hybrid edges that do not
appear in the best network have a missing number.
After hybrid edges from bootstrap networks are added,
`BSe` is modified to include the edge numbers of the newly added hybrid edges.
To distinguish hybrid edges present in the original network versus new edges,
an extra column of true/false values is also added to `BSe`, named "alternative",
with true for newly added edges absent from the original network.

The hybrid edges added to `net` are added as minor edges, to keep the underlying
major tree topology.

# example

```jldoctest
julia> bootnet = readMultiTopology(joinpath(dirname(pathof(PhyloNetworks)), "..","examples", "bootsnaq.out")); # vector of 10 networks

julia> bestnet = readTopology("((O,(E,#H7:::0.196):0.314):0.332,(((A)#H7:::0.804,B):10.0,(C,D):10.0):0.332);");

julia> BSn, BSe, BSc, BSgam, BSedgenum = hybridclades_support(bootnet, bestnet);

julia> BSe[1:6,[:edge,:hybrid_clade,:sister_clade,:BS_hybrid_edge]]
6×4 DataFrame
 Row │ edge     hybrid_clade  sister_clade  BS_hybrid_edge 
     │ Int64?   String        String        Float64        
─────┼─────────────────────────────────────────────────────
   1 │       7  H7            B                       33.0
   2 │       3  H7            E                       32.0
   3 │ missing  c_minus3      c_minus8                44.0
   4 │ missing  c_minus3      H7                      44.0
   5 │ missing  E             O                       12.0
   6 │ missing  c_minus6      c_minus8                 9.0

julia> PhyloNetworks.addAlternativeHybridizations!(bestnet, BSe)

julia> BSe[1:6,[:edge,:hybrid_clade,:sister_clade,:BS_hybrid_edge,:alternative]]
6×5 DataFrame
 Row │ edge     hybrid_clade  sister_clade  BS_hybrid_edge  alternative 
     │ Int64?   String        String        Float64         Bool        
─────┼──────────────────────────────────────────────────────────────────
   1 │       7  H7            B                       33.0        false
   2 │       3  H7            E                       32.0        false
   3 │      16  c_minus3      c_minus8                44.0         true
   4 │      19  c_minus3      H7                      44.0         true
   5 │      22  E             O                       12.0         true
   6 │ missing  c_minus6      c_minus8                 9.0        false

julia> # using PhyloPlots; plot(bestnet, edgelabel=BSe[:,[:edge,:BS_hybrid_edge]]);
```
"""
function addAlternativeHybridizations!(
    net::HybridNetwork,
    BSe::DataFrame;
    cutoff::Number=10,
    top::Int=3
)
    top > 0 || error("top must be greater than 0")
    BSe[!,:alternative] = falses(nrow(BSe))
    newBSe = subset(BSe,
        :BS_hybrid_edge => x -> x.> cutoff, :edge   => ByRow( ismissing),
        :hybrid => ByRow(!ismissing),       :sister => ByRow(!ismissing),
    )
    top = min(top,nrow(newBSe))
    if top==0
        @info "no alternative hybridizations with support > cutoff $cutoff%, so nothing added."
        return
    end
    for i in 1:top
        hybnum = newBSe[i,:hybrid]
        sisnum = newBSe[i,:sister]
        edgenum = addHybridBetweenClades!(net, hybnum, sisnum)
        if isnothing(edgenum)
          @warn "cannot add desired hybrid (BS=$(newBSe[i,:BS_hybrid_edge])): the network would have a directed cycle"
          continue
        end
        ind1 = findall(x->!ismissing(x) && x==hybnum, BSe[!,:hybrid])
        ind2 = findall(x->!ismissing(x) && x==sisnum, BSe[!,:sister])
        ind = intersect(ind1,ind2)
        BSe[ind,:edge] .= edgenum
        BSe[ind,:alternative] .= true
    end
end


"""
    addHybridBetweenClades!(net::HybridNetwork, hybnum::Number, sisnum::Number)

Modify `net` by adding a minor hybrid edge from "donor" to "recipient",
where "donor" is the major parent edge `e1` of node number `hybnum` and
"recipient" is the major parent edge `e2` of node number `sisnum`.
The new nodes are currently inserted at the middle of these parent edges.

If a hybrid edge from `e1` to `e2` would create a directed cycle in the network,
then this hybrid cannot be added.
In that case, the donor edge `e1` is moved up if its parent is a hybrid node,
to ensure that the sister clade to the new hybrid would be a desired (the
descendant taxa from `e1`) and a new attempt is made to create a hybrid edge.

Output: number of the new hybrid edge, or `nothing` if the desired hybridization
is not possible.

See also:
[`addhybridedge!`](@ref) (used by this method) and
[`directionalconflict`](@ref) to check that `net` would still be a DAG.
"""
function addHybridBetweenClades!(net::HybridNetwork, hybnum::Number, sisnum::Number)
    hybind = getIndexNode(hybnum,net)
    sisind = getIndexNode(sisnum,net)
    e1 = getparentedge(net.node[sisind]) # major parent edges
    e2 = getparentedge(net.node[hybind])
    p1 = getparent(e1)
    if directionalconflict(p1, e2, true) # then: first try to move the donor up
        # so long as the descendant taxa (= sister clade) remain the same
        while p1.hybrid
          e1 = getparentedge(p1) # major parent edge: same descendant taxa
          p1 = getparent(e1)
        end
        directionalconflict(p1, e2, true) && return nothing
    end
    hn, he = addhybridedge!(net, e1, e2, true) # he: missing length & gamma by default
    # ideally: add option "where" to breakedge!, used by addhybridedge!
    # so as to place the new nodes at the base of each clade.
    # currently: the new nodes are inserted at the middle of e1 and e2.
    return he.number
end