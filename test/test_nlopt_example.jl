# small test with one tree quartet for the NLopt
# Claudia January 2015

treal = 0.1
obsCF=[1-2/3*exp(-treal),1/3*exp(-treal),1/3*exp(-treal)]
using NLopt

type ToyQuartet
    t::Float64
end

function objective(x::Float64,obsCF::Vector{Float64})
    length(obsCF) == 3 || error("obsCF must be size 3")
    val = obsCF[1]*log(1-2/3*exp(-x)) + obsCF[2]*log(1/3*exp(-x)) + obsCF[3]*log(1/3*exp(-x))
    return -val
end

function optimizar(obsCF::Vector{Float64}, q::ToyQuartet)
    t = q.t
    k = length(t)
    opt = NLopt.Opt(:LN_BOBYQA,k) # :LD_MMA if use gradient
    # criterion based on prof Bates code
    NLopt.ftol_rel!(opt,1e-12) # relative criterion
    NLopt.ftol_abs!(opt,1e-8) # absolute critetion
    NLopt.xtol_abs!(opt,1e-10) # criterion on parameter value changes
    NLopt.lower_bounds!(opt, zeros(k))
    NLopt.upper_bounds!(opt, Inf)
    count = 0
    function obj(x::Vector{Float64},g::Vector{Float64}) # added g::Vector{Float64} for gradient, ow error
        count += 1
        println("t is $(t) initially, x: $(x)")
        t = deepcopy(x)
        println("t is now $(t), x: $(x)")
        val = objective(x[1],obsCF)
        println("f_$count: $(round(val, digits=5)), x: $(x)")
        return val
    end
    NLopt.min_objective!(opt,obj)
    fmin, xmin, ret = NLopt.optimize(opt,[t])
    return fmin,xmin
end

q = ToyQuartet(1.0)
fmin,xmin = optimizar(obsCF,q)
