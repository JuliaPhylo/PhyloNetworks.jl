"""
    CacheLengthLiNC

Type to store intermediate values during optimization of one branch length,
to limit time spent on garbage collection. Re-used across branches.
"""
struct CacheLengthLiNC
    """forward likelihood of descendants of focus edge: nstates x nsites x nrates x ntrees"""
    flike::Array{Float64,4}
    """likelihood of non-descendants of focus edge: direct * backwards.
       size: nstates x nsites x nrates x ntrees"""
    dblike::Array{Float64,4}
    """whether displayed trees have the focus edge (true) or not (false),
       which can happen for hybrid edges or on non-tree-child networks"""
    hase::Vector{Bool}
    "transition probability matrices: P(r*t) for each rate r"
    Prt::Vector{StaticArrays.MMatrix}
    "rQP(rt) for each rate r, for the gradient"
    rQP::Vector{StaticArrays.MMatrix}
    "site-specific gradient numerator for log-likelihood: nsites"
    glik::Vector{Float64}
    opt::NLopt.Opt
end
function CacheLengthLiNC(obj::SSM,
        ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
        maxeval::Int)
    k = nstates(obj.model)
    nt = size(obj._loglikcache, 1) # but: there might be fewer displayed trees
    # fixit: make other things with new nt?
    # move nt to last dimension? change the order of all loops?
    nr = length(obj.ratemodel.ratemultiplier)
    ns = obj.nsites
    Prt = [MMatrix{k,k, Float64}(undef) for i in 1:nr]
    rQP = [MMatrix{k,k, Float64}(undef) for i in 1:nr]
    opt = NLopt.Opt(:LD_SLSQP, 1)
    NLopt.ftol_rel!(opt,ftolRel) # relative criterion
    NLopt.ftol_abs!(opt,ftolAbs) # absolute criterion
    NLopt.xtol_rel!(opt,xtolRel)
    NLopt.xtol_abs!(opt,xtolAbs)
    NLopt.maxeval!(opt, maxeval) # max number of iterations
    opt.lower_bounds = 0.0 # fixit: change to BLmin
    # fixit: add upper bound
    CacheLengthLiNC(
        Array{Float64}(undef, (k, ns, nr, nt)),
        Array{Float64}(undef, (k, ns, nr, nt)),
        Vector{Bool}(undef, nt),
        Prt,
        rQP,
        Vector{Float64}(undef, ns),
        opt,
    )
end

"""
    updatecache_edge!(lcache::CacheLengthLiNC, obj::SSM, focusedge)

Update fields in `lcache`, to correspond to the forward likelihood
at the child node of the focus edge,
backwards (at parent node) x directional (at sister edges) likelihoods,
and keep track of which displayed trees do have the focus edge.
These don't change if the length of the focus edge is modified.
These quantities are then rescaled on the log scale (to get a max of 0)
and exponentiated to get back to the probability scale.

Assumptions: the following fields of `obj` are up-to-date:
- log transition probabilities in `.logtrans`
- displayed trees, with nodes preordered
- prior log tree weights in `.priorltw`

Output: constant used to rescale each site, on the log scale.
"""
function updatecache_edge!(lcache::CacheLengthLiNC, obj::SSM, focusedge)
    flike  = lcache.flike
    dblike = lcache.dblike
    hase   = lcache.hase
    tree   = obj.displayedtree
    k, ns, nr, nt = size(flike)
    nt = length(tree) # could be less than dimensions in lcache
    # ! a sister edge in network may be absent in a displayed tree
    # but: edge & node numbers are the same in net and in trees
    v = getParent(focusedge)
    vnum = v.number
    unum = getChild(focusedge).number
    enum = focusedge.number
    snum = Int[]
    for e in v.edge
        if e !== focusedge && v == getParent(e) # then e sister to focus edge
            push!(snum, e.number)
        end
    end
    nsis = length(snum) # normally 1, 0 below hybrid node, 2 at root or at polytomies
    hassis = Array{Bool}(undef, (nsis,nt)) # is sister s in tree t?
    @inbounds for it in 1:nt
        hase[it] = any(x -> x.number == enum, tree[it].edge)
        for j in 1:nsis
            hassis[j,it] = any(x -> x.number == snum[j], tree[it].edge)
        end
    end
    ftmp = obj.forwardlik
    btmp = obj.backwardlik
    dtmp = obj.directlik
    clik = obj._loglikcache # conditional likelihood: conditional on tree & rate
    fill!(flike, -Inf); fill!(dblike, -Inf)
    for it in 1:nt
        if hase[it] # we do need flike and dblike
          for ir in 1:nr for is in 1:ns
            clik[it,ir,is] = discrete_corelikelihood_trait!(obj, it,is,ir, ftmp,dtmp)
            discrete_backwardlikelihood_trait!(obj, it,ir, btmp,dtmp)
            flike[:,is,ir,it]  .= ftmp[:,unum]
            dblike[:,is,ir,it] .= btmp[:,vnum]
            for isis in 1:nsis
              @inbounds if hassis[isis,it]
                dblike[:,is,ir,it] .+= dtmp[:,snum[isis]] # sum on log scale
              end
            end
          end; end
        else # tree does't have the focus edge: use clik only
          for ir in 1:nr for is in 1:ns
            clik[it,ir,is] = discrete_corelikelihood_trait!(obj, it,is,ir, ftmp,dtmp)
          end; end
        end
    end
    cf = maximum(flike) # unused values from trees without focus edge have -Inf
    cg = maximum(dblike)
    flike  .= exp.(flike .- cf)
    dblike .= exp.(dblike .- cg)
    cfg = cf + cg # lcache.cadjust = cf + cg
    clik .= exp.(clik .- cfg)
    return cfg
end

"""
optimizelocalBL_LiNC! to replace current version

fixit: add lcache as argument, create it within phyLiNC.

fixit: edit docstrings.
assume: displayed trees up-to-date with nodes preordered,
and priorltw up-to-date.
"""
function optimizelocalBL_LiNC!(obj::SSM, focusedge::Edge, # lcache::CacheLengthLiNC,
        ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64,
        maxeval=1000::Int)
    lcache = CacheLengthLiNC(obj, ftolRel,ftolAbs,xtolRel,xtolAbs, maxeval)
    neighboredges = adjacentedges(focusedge)
    # remove constrained edges: to keep network unzipped
    if !isempty(obj.net.hybrid)
        for i in length(neighboredges):-1:1
            @inbounds e = neighboredges[i]
            if getParent(e).hybrid
                deleteat!(neighboredges, i) # shallow copy!!
            end
        end
    end
    # assume: correct preordered displayed trees & tree weights
    # startlik = discrete_corelikelihood!(obj)
    update_logtrans(obj) # fixit: waste of resources?? avoid in some cases?
    qmat  = Q(obj.model)
    ltw = exp.(obj.priorltw)
    for e in neighboredges
        optimizelength_LiNC!(obj, e, lcache, qmat, ltw)
    end
    return neighboredges # excludes constrained edges
end

"""
    optimizelength_LiNC!(obj::SSM, edge::Edge, cache::CacheLengthLiNC,
                         Qmatrix, treeweights)

Optimize the length of a single `edge` using a gradient-based optimization.
Output: fixit

Warnings:
- displayed trees are assumed up-to-date, with nodes preordered
- no check that `edge` is not below a reticulation (the unzipped canonical
  version should have a length of 0 below reticulations).
"""
function optimizelength_LiNC!(obj::SSM, focusedge::Edge,
                              lcache::CacheLengthLiNC, qmat, ltw)
    startlik = obj.loglik
    cfg = updatecache_edge!(lcache, obj, focusedge)
    flike  = lcache.flike
    dblike = lcache.dblike
    hase   = lcache.hase
    Prt    = lcache.Prt
    rQP    = lcache.rQP
    glik   = lcache.glik   # glik = d(ulik)/dt
    tree   = obj.displayedtree
    k, ns, nr, nt = size(flike)
    nt = length(tree) # could be less than dimension in lcache
    rates = obj.ratemodel.ratemultiplier
    ulik  = obj._sitecache   # unconditional likelihood
    clik  = obj._loglikcache # conditional likelihood
    adjustment = obj.totalsiteweight * (cfg - log(nr))

    # use site weights, if any
    wsum = (isnothing(obj.siteweight) ? sum : x -> sum(x .* obj.siteweight))

    # objective to maximize: overall log-likelihood
    function objective(t::Vector, grad::Vector)
        len = t[1] # candidate branch length
        @inbounds for ir in 1:nr
            P!(Prt[ir], obj.model, len * rates[ir])
            lmul!(rates[ir], mul!(rQP[ir], qmat, Prt[ir])) # in-place multiplication
        end
        # integrate over trees & rates
        fill!(ulik, 0.0); fill!(glik, 0.0)
        for it in 1:nt
          if hase[it] # tree has edge: contributes to gradient
            for is in 1:ns
                u=0.0; g=0.0
                for ir in 1:nr
                    # fixit: dot(x,A,y) requires at least Julia 1.4
                    u += dot(dblike[:,is,ir,it], Prt[ir], flike[:,is,ir,it])
                    g += dot(dblike[:,is,ir,it], rQP[ir], flike[:,is,ir,it])
                end
                ulik[is] += u * ltw[it]
                glik[is] += g * ltw[it]
            end
          else # tree doesn't have focus edge: use clik only, gradient=0
            for is in 1:ns
                u=0.0
                for ir in 1:nr
                    u += clik[it,ir,is] # already exp-ed inside updatecache_edge!
                end
                ulik[is] += u * ltw[it]
            end
          end
        end
        # product over sites (or sum on log scale)
        loglik = wsum(log.(ulik))
        if length(grad) > 0 # finish gradient calculation
            glik ./= ulik
            grad[1] = wsum(glik)
            #@show grad[1]
        end
        #@show loglik + adjustment
        return loglik
    end
    # ll = objective([focusedge.length], [0.0]) + adjustment
    # ll ≈ startlik || @warn "ll = $ll is different from starting likelihood $startlik"
    optBL = lcache.opt
    NLopt.max_objective!(optBL, objective)
    fmax, xmax, ret = NLopt.optimize(optBL, [focusedge.length])
    @debug "BL: got $(round(fmax; digits=5)) at BL = $(round.(xmax, digits=5)) after $(optBL.numevals) iterations (return code $(ret))"
    newlik = fmax + adjustment
    if startlik > newlik
        @debug "starting likelihood better than after optimization. Skipping branch length update."
        return nothing
    end
    focusedge.length = xmax[1]
    obj.loglik = newlik
    update_logtrans(obj, focusedge) # obj ready for optimizing another edge length
    return nothing
end

#=
function optimizelength_LiNC!(obj::SSM, focusedge::Edge, lcache::CacheLengthLiNC)
    fun = objective(t) # returns vector [loglik,gradient,hessian]
    f = x::Vector -> -fun(x[1])[1] # minimize -loglik
    function g!(g,x::Vector)
        g = [-fun(x[1])[2]]
    end
    function h!(h,x::Vector)
        h = [-fun(x[1])[3]]
    end
    x0 = [focusedge.length]
    df = TwiceDifferentiable(f, g!, h!, x0)
    lx = [0.]; ux = [Inf];
    dfc = TwiceDifferentiableConstraints(lx, ux)
    res = Optim.optimize(df, dfc, x0, IPNewton())
end
=#

# fixit: docstrings
# edge should be an edge in the network
function update_logtrans(obj::SSM, edge::Edge)
    rates = obj.ratemodel.ratemultiplier
    enum = edge.number
    len = edge.length
    for i in 1:length(rates)
        pmat = view(obj.logtrans, :,:,enum,i)
        @inbounds pmat .= log.(P!(pmat, obj.model, len * rates[i]))
    end
end
