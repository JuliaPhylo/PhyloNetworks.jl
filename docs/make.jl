using Documenter, PhyloNetworks

makedocs()

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/crsl4/PhyloNetworks.jl.git",
    julia  = "0.4",
    osname = "osx"
)
