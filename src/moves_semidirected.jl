#=
move types:
 - NNI = nearest neighbor interchange, as defined in
   Gambette, van Iersel, Jones, Lafond, Pardi and Scornavacca (2017).
   Rearrangement moves on rooted phylogenetic networks.
   PLOS Computational Biology 13(8):e1005611.
 - change the direction of a hybrid edge
 - add or remove a hybrid edge
constrained to user-defined topological constraints (like an outgroup)

in `moves.jl`, functions are tailored to level-1 networks;
here, functions apply to semidirected networks of all levels.
=#

"""
Type for various topological constraints, such as:

1. a set of taxa forming a clade in the major tree
2. a set of individuals belonging to the same species
3. a set of taxa forming an outgroup clade or grade in the major tree
4. a set of taxa forming a clade in any one of the displayed trees
5. a set of taxa forming an outgroup clade or grade in any one of the displayed trees
"""
struct TopologyConstraint
    "type of constraint: 1, 2, or 3. Types 4 and 5 not implemented yet"
    type::UInt8
    "names of taxa in the constraint (clade / species / outgroup members)"
    taxonname::Set{String}
    """
    node numbers of taxa in the constraint (clade / species / outgroup members).
    warning: interpretation will be network-dependent.
    """
    taxonnum::Set{Int}
    "number of the stem edge of the constrained group"
    edgenum::Int
end
function TopologyConstraint(type::UInt8, taxonname::Set{String}, net::HybridNetwork)
    taxonnum = Set{Int}() # fixit: use the network instead
    edgenum = -1 # fixit: check that the network satisfies the constraint and return the associated edge
    TopologyConstraint(type, taxonname, taxonnum, edgenum)
end
# fixit: use these constraints in the functions below,
# to check that the proposed NNI move is acceptable

#= TODO
- new test file with example uses
- add option `no3cycle` to forbid moves that would create a 3-cycle
- add option to modify branch lengths during NNI:
  `nni!` that take detailed input, and output of same signature
- `unzip` (bool) to constrain branch lengths to 0 below a hybrid node:
  make this part of the optimization process. Or option for `nni!` functions?
- RR move, during optimization: avoid the re-optimization of the likelihood.
  accept the move, just change inheritance values to get same likelihood
  (still update direct / forward likelihoods and tree priors?)
- add options to forbid non-tree-child networks
- add moves that change the root position (and edge directions accordingly)
  without re-calculating the likelihood (accept the move): perhaps recalculate forward / direct likelihoods
  why: the root position affects the feasibility of the NNIs starting from a BR configuration
- to check if there is a path from n1 to n2: see if we can use `descendants` and `descendants!` in `compareNetworks.jl`
  check condition for getting a DAG in case the network is semi-rooted, BR case
=#

"""
    nni!(net::HybridNetwork, e::Edge, constraint::Vector{TopologyConstraint})

Attempt to perform a nearest neighbor interchange (NNI) around edge `e`,
randomly chosen among all possible NNIs (e.g 3, sometimes more depending on `e`)
satisfying the constraints, and such that the new network is a DAG.

Output: information indicating which NNI was performed, or `nothing` if all NNIs failed.
fixit: what information to return exactly? NNI number? info to undo?

# examples

```jldoctest
fixit: define a network, find suitable edge index i,
undoinfo = nni!(net, net.edge[i])
net
nni!(undoinfo...)
net
```
fixit: add example with constraint
"""
function nni!(net::HybridNetwork, e::Edge, constraint::Vector{TopologyConstraint})
    hybparent = getParent(e).hybrid
    nnirange = 0x01:nnimax(e) # 0x01 = 1 but UInt8 instead of Int
    nnis = Random.shuffle(nnirange)
    for nummove in nnis # iterate through all possible NNIs, but in random order
        moveinfo = nni!(net, e, nummove)
        moveinfo !== nothing || continue # to next possible NNI
        # fixit: if failure, make sure the network was not changed by nni!, else undo it
        if checknetwork(net, constraint)
            return moveinfo
        else nni!(moveinfo...); end # undo the previous NNI, try again
    end
    return nothing # if we get to this point, all NNIs failed
end

"""
    nnimax(e::Edge)

Return the number of NNI moves around edge `e`,
assuming that the network is semi-directed.
Return 0 if `e` is not internal, and more generally if either node
attached to `e` does not have 3 edges.

Output: UInt8 (e.g. `0x02`)
"""
function nnimax(e::Edge)
    length(e.node[1].edge) == 3 || return 0x00
    length(e.node[2].edge) == 3 || return 0x00
    # e.hybrid and hybrid parent: RR case, 2 possible NNIs only
    # e.hybrid and tree parent:   BR case, 3 or 6 NNIs if e is may contain the root
    # e not hybrid, hyb parent:   RB case, 4 NNIs
    # e not hybrid, tree parent:  BB case, 2 NNIs if directed, 8 if undirected
    hybparent = getParent(e).hybrid
    n = (e.hybrid ? (hybparent ? 0x02 : (e.containRoot ? 0x06 : 0x03)) : # RR & BR
                    (hybparent ? 0x04 : (e.containRoot ? 0x08 : 0x02)))  # RB & BB
    return n
end

"""
    edgerelation(e::Edge, node::Node, origin::Edge)

fixit: move this function into auxiliary.jl (with better name?)

Return a symbol:

- `:origin` if `e` is equal to `origin`, and otherwise:
- `:parent` if `e` is a parent of `node`,
- `:child` if `e` is a child of `node`

using the `isChild1` attribute of edges.
Useful when `e` iterates over all edges adjacent to `node` and when
`origin` is one of the edges adjacent to `node`,
to known the order in which these edges come.

example:
```julia
labs = [edgerelation(e, u, uv) for e in u.edge] # assuming u is a node of edge uv
parentindex = findfirst(isequal(:parent), labs) # could be 'nothing' if no parent
childindices = findall( isequal(:child), labs)  # vector. could be empty
```
"""
function edgerelation(ee::Edge, n::Node, origin::Edge)
    (ee===origin ? :origin : (n===getChild(ee) ? :parent : :child))
end

"""
    nni!(net::HybridNetwork, uv::Edge, nummove::UInt8)

Modify `net` with a nearest neighbor interchange (NNI) around edge `uv`.
Return the information necessary to undo the NNI, or `nothing` if the move
was not successful (such as if the resulting graph was not acyclic or if
the focus edge is adjacent to a polytomy).
`nummove` specifies which of the available NNIs is performed.
The network is modified even if the NNI failed. (fixit: check and update the comment accordingly!!)

rooted-NNI options according to Gambette et al. (2017), fig. 8:
* BB: 2 moves, both to BB, if directed edges. 8 moves if undirected.
* RR: 2 moves, both to RR.
* BR: 3 moves, 1 RB & 2 BRs, if directed. 6 moves if e is undirected.
* RB: 4 moves, all 4 BRs.
The extra options are due to assuming a semi-directed network, whereas
Gambette et al (2017) describe options for rooted networks.
On a semi-directed network, there might be a choice of how to direct
the edges that may contain the root, e.g. choice of e=uv versus vu, and
choice of labelling adjacent nodes as α/β (BB), or as α/γ (BR).

Warning: the edges' field `isChild1` is assumed to be correct,
according to the `net.root`.
"""
function nni!(net::HybridNetwork, uv::Edge, nummove::UInt8)
    nmovemax = nnimax(uv)
    if nmovemax == 0x00 # uv not internal, or polytomy e.g. species represented by polytomy with >2 indiv
        return nothing
    end
    nummove > 0x00 || error("nummove $(nummove) must be >0")
    nummove <= nmovemax || error("nummove $(nummove) must be <= $nmovemax for edge number $(uv.number)")
    u = getParent(uv)
    v = getChild(uv)
    ## TASK 1: grab the edges adjacent to uv: αu, βu, vγ, vδ and adjacent nodes
    # get edges αu & βu connected to u
    # α = u's major parent if it exists, β = u's last child or minor parent
    if u.hybrid
        αu = getMajorParentEdge(u)
        βu = getMinorParentEdge(u)
    else # u may not have any parent, e.g. if root node
        # pick αu = parent edge if possible, first edge of u (other than uv) otherwise
        labs = [edgerelation(e, u, uv) for e in u.edge]
        pi = findfirst(isequal(:parent), labs) # parent index: there should be 1 or 0
        ci = findall(  isequal(:child), labs)  # vector. there should be 1 or 2
        if pi === nothing
            length(ci) == 2 || error("node $(u.number) should have 2 children other than $(v.number)")
            pi = popfirst!(ci)
        end
        αu = u.edge[pi]
        βu = u.edge[ci[1]]
    end
    # get edges vδ & vγ connected to v
    # δ = v's last child, γ = other child or v's parent other than u
    labs = [edgerelation(e, v, uv) for e in v.edge]
    if v.hybrid # then v must have another parent edge, other than uv
        vγ = v.edge[findfirst(isequal(:parent), labs)] # γ = getParent(vδ)
        vδ = v.edge[findfirst(isequal(:child), labs)]
    else
        ci = findall(isequal(:child), labs)
        length(ci) == 2 || error("node $(v.number) should have 2 children")
        vγ = v.edge[ci[1]] # γ = getChild(vγ)
        vδ = v.edge[ci[2]]
    end
    ## TASK 2: semi-directed network swaps:
    ## swap u <-> v, α <-> β if undirected and according to move number
    if !u.hybrid && uv.containRoot # case BB or BR, with uv that may contain the root
        nmoves = 0x03 # number of moves in the rooted (directed) case
        if !v.hybrid # BB, uv and all adjacent edges are undirected: 8 moves
            if nummove > 0x04  # switch u & v with prob 1/2
                nummove -= 0x04  # now in 1,2,3,4
                (u,v) = (v,u)
                (αu,βu,vγ,vδ) = (vγ,vδ,αu,βu)
                # (α,αu,β,βu,γ,vγ,δ,vδ) = (γ,vγ,δ,vδ,α,αu,β,βu)
            end
            nmoves = 0x02
        end
        # now nmoves should be in 1, ..., 2*nmoves: 1:4 (BB) or 1:6 (BR)
        # next: switch α & β with probability 1/2
        if nummove > nmoves
            nummove -= nmoves # now nummoves in 1,2 (BB) or 1,2,3 (BR)
            (αu,βu) = (βu,αu)
        end
    end
    ## TASK 3: rooted network swaps, then do the NNI
    if u.hybrid # RR or RB cases
        # swap α & β: reduce NNIs from 2 to 1 (RR) or from 4 to 2 (RB)
        if nummove == 0x02 || nummove == 0x04
            (αu,βu) = (βu,αu)
        end
        if !v.hybrid && nummove > 0x02 # case RB, moves 3 or 4
            (vγ,vδ) = (vδ,vγ)
        end
        # detach β and graft onto vδ
        res = nni!(αu, u, uv, v, vδ)
    else # u = bifurcation: BB or BR cases
        if !v.hybrid # case BB: 2 options
            if nummove == 0x02
                (vγ,vδ) = (vδ,vγ)
            end
            # detach β and graft onto vδ
            res = nni!(αu, u, uv, v, vδ)
        else # case BR: 3 options
            if nummove == 0x01 # graft γ onto α
                # fixit: check that there's no directed path from child of u -> γ (if u has a parent),
                #        or no path from α -> γ (if u = root). If so: return nothing
                res = nni!(vδ, v, uv, u, αu)
            elseif nummove == 0x02 # graft γ onto β
                # fixit: check that there's no directed path from child of u -> γ (if u has a parent),
                #        or no path from β -> γ (if u = root). If so: return nothing
                res = nni!(vδ, v, uv, u, βu)
            else # nummove == 0x03
                if getChild(αu)===u # if α->u: graft δ onto α
                    # fixit: check no directed path from α -> γ
                    res = nni!(vγ, v, uv, u, αu)
                elseif getChild(βu)===u # if β->u: graft δ onto β
                    res = nni!(vγ, v, uv, u, βu)
                else
                    return nothing # fail: nummove 3 is impossible when u is the root: not DAG
                end
            end
        end
    end
    return res
end

"""
    nni!(αu, u, uv::Edge, v, vδ)
    nni!(αu,u,uv,v,vδ, flip::Bool, inner::Bool, indices)

Transform a network locally around the focus edge `uv` with
the following NNI, that detaches u-β and grafts it onto vδ:

    α - u -- v ------ δ
        |    |
        β    γ

    α ------ v -- u - δ
             |    |
             γ    β

**Warnings**:
- *No* check of assumed adjacencies
- Not implemented for cases that are not necessary thanks to symmetry,
  such as cases covered by `nni!(vδ, v, uv, u, αu)` or `nni!(βu, u, v, vγ)`.
  More specifically, these cases are not implemented (and not checked):
  * u not hybrid & v hybrid
  * u hybrid, v not hybrid, α -> u <- v -> δ
- Because of this, `nni(αu,u,uv,v,vδ, ...)` should not be used directly;
  use instead `nni!(net, uv, move_number)`.

Node numbers and edge numbers are not modified.
Edge `uv` keeps its direction unchanged *unless* the directions were
`α -> u -> v -> δ` or `α <- u <- v <- δ`,
in which case the direction of `uv` is flipped.

The second version is to undo the NNI more easily.
With the second version, the input has the same signature as the output.
If `output = nni!(input)`, then `nni!(output...)` is valid and undoes
the first operation.

fixit: describe how branch lengths are modified,
"""
function nni!(αu::Edge, u::Node, uv::Edge, v::Node, vδ::Edge)
    # find indices to detach αu from u, and to detach v from vδ
    u_in_αu = findfirst(n->n===u, αu.node)
    αu_in_u = findfirst(e->e===αu, u.edge)
    vδ_in_v = findfirst(e->e===vδ, v.edge)
    v_in_vδ = findfirst(n->n===v, vδ.node)
    # none of them should be 'nothing' --not checked
    αu_child = getChild(αu)
    uv_child = getChild(uv)
    vδ_child = getChild(vδ)
    # flip the direction of uv if α->u->v->δ or α<-u<-v<-δ
    flip = (αu_child === u && uv_child === v && vδ_child !== v) ||
             (αu_child !== u && uv_child === u && vδ_child === v)
    inner = (flip ? false : αu_child === u && vδ_child === v )
    # fixit: calculate new edge lengths here; give them as extra arguments to nni! below
    res = nni!(αu,u,uv,v,vδ, flip,inner, u_in_αu,αu_in_u,vδ_in_v,v_in_vδ)
    return res
end
function nni!(αu::Edge, u::Node, uv::Edge, v::Node, vδ::Edge,
              flip::Bool, inner::Bool,
              u_in_αu::Int, αu_in_u::Int, vδ_in_v::Int, v_in_vδ::Int)
              # extra arguments: new edge lengths
    # fixit: extract old edge lengths; update to new lengths; return old lengths as extra items
    ## TASK 1: do the NNI swap: to detach & re-attach at once
    (αu.node[u_in_αu], u.edge[αu_in_u], v.edge[vδ_in_v], vδ.node[v_in_vδ]) =
        (v, vδ, αu, u)
    ## TASK 2: update edges' directions
    if flip
        uv.isChild1 = !uv.isChild1
    end
    ## TASK 3: update hybrid status of nodes & edges, and edges containRoot field
    if u.hybrid && !v.hybrid
        if flip # RB, BR 1 or BR 2': flip hybrid status of uv and (αu or vδ)
            if uv.hybrid # BR 1 or BR 2'
                uv.hybrid = false
                vδ.hybrid = true
                norootbelow!(uv)
            else # RB
                uv.hybrid = true
                αu.hybrid = false
                if αu.containRoot
                    allowrootbelow!(v, αu)
                end
            end
        else # assumes α<-u or v<-δ --- but not checked
            # not implemented: α->u<-v->δ: do uv.hybrid=false and γv.hybrid=true
            if inner # BR 3 or 4': α->u<-v<-δ: flip hybrid status of αu and vδ
                αu.hybrid = false
                vδ.hybrid = true
                # containRoot: nothing to update
            else # BR 1' or 2: α<-u<-v->δ : hybrid status just fine
                norootbelow!(vδ)
                if uv.containRoot
                    allowrootbelow!(αu)
                end
            end
        end
    end
    # throw error if u not hybrid & v hybrid? (does not need to be implemented)
    return vδ,u,uv,v,αu, flip,inner, v_in_vδ,αu_in_u,vδ_in_v,u_in_αu
end


"""
    cladesviolated(network::HybridNetwork, cladeconstraints)

Check if network violates user-given clade constraints. Return false if passes,
return true if violates clades.
"""
function cladesviolated(net::HybridNetwork, cladeconstraints::Dict)
    return false
end
