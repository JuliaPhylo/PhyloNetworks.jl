using Documenter, PhyloNetworks

makedocs()

# Versions of `mkdocs` and `mkdocs-material` specified manually to avoid conflicts. 
# See here: https://discourse.julialang.org/t/mkdocs-material-in-documenter/13764/3
# To be kept in mind: those versions might evolve in the future.
deploydocs(
    deps   = Deps.pip("pygments", "mkdocs==0.17.5", "mkdocs-material==2.9.4", "python-markdown-math"),
    repo = "github.com/crsl4/PhyloNetworks.jl.git",
    julia  = "0.6",
    osname = "linux"
)
