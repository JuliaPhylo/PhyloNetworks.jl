using Documenter, PhyloNetworks

makedocs()

deploydocs(
    deps   = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math"),
    repo = "github.com/crsl4/PhyloNetworks.jl.git",
    julia  = "0.5",
    osname = "osx"
)
