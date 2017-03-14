using Documenter, PhyloNetworks

include(Pkg.dir("PhyloNetworks","docs","src", "man", "src", "make_weave.jl"))

makedocs()

deploydocs(
    deps   = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math"),
    repo = "github.com/pbastide/PhyloNetworks.jl.git",
    julia  = "0.5",
    osname = "osx"
)
