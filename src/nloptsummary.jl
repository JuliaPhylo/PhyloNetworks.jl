"""
    OptSummary{T<:AbstractFloat}

Summary of an `NLopt` optimization. Idea and code taken from
[`MixedModels`](https://github.com/JuliaStats/MixedModels.jl).
`T` is the type of the function argument(s) and of the function value.
"""
mutable struct OptSummary{T<:AbstractFloat}
    "copy of initial param values in the optimization"
    initial::Vector{T}
    "lower bounds on the parameter values"
    lowerbd::Vector{T}
    "relative tolerance on the function value f(x), as in NLopt"
    ftol_rel::T
    "absolute tolerance on the function value f(x), as in NLopt"
    ftol_abs::T
    "relative tolerance on the argument value x, as in NLopt"
    xtol_rel::T
    "absolute tolerance on the argument value x"
    xtol_abs::Vector{T}
    "initial step sizes, as in NLopt"
    initial_step::Vector{T}
    "maximum number of function evaluations, as in NLopt"
    maxfeval::Int
    "copy of the final parameter values from the optimization"
    final::Vector{T}
    "final value of the objective function"
    fmin::T
    "number of function evaluations"
    feval::Int
    "name of the optimizer used, as a `Symbol`"
    algorithm::Symbol
    "return value, as a `Symbol`"
    returnvalue::Symbol
end
function OptSummary(initial::Vector{T}, lowerbd::Vector{T}, algorithm::Symbol;
        ftol_rel::T = zero(T),
        ftol_abs::T = zero(T),
        xtol_rel::T = zero(T),
        xtol_abs::Vector{T} = fill(zero, length(initial)),
        initial_step::Vector{T} = T[]) where {T<:AbstractFloat}
    OptSummary(initial, lowerbd,
        ftol_rel, ftol_abs, xtol_rel, xtol_abs,
        initial_step,
        -1, # maxeval: 0 or negative for no limit
        copy(initial), # final x
        T(Inf), # final value (will be finite after minimization)
        -1, # number of evals
        algorithm,
        :FAILURE)
end

function NLopt.Opt(optsum::OptSummary)
    lb = optsum.lowerbd
    np = length(lb) # number of parameters to optimize
    opt = NLopt.Opt(optsum.algorithm, np)
    NLopt.ftol_rel!(opt, optsum.ftol_rel) # relative criterion on objective
    NLopt.ftol_abs!(opt, optsum.ftol_abs) # absolute criterion on objective
    NLopt.xtol_rel!(opt, optsum.xtol_rel) # relative criterion on parameters
    NLopt.xtol_abs!(opt, optsum.xtol_abs) # absolute criterion on parameters
    NLopt.lower_bounds!(opt, lb)
    NLopt.maxeval!(opt, optsum.maxfeval)
    if isempty(optsum.initial_step)
        optsum.initial_step = NLopt.initial_step(opt, similar(lb))
    else
        NLopt.initial_step!(opt, optsum.initial_step)
    end
    return opt
end
