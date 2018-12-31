using Documenter, DocumenterMarkdown

using Pkg
Pkg.add(PackageSpec(name="PhyloPlots", rev="master"))

using PhyloNetworks

makedocs(sitename = "PhyloNetworks.jl",
         modules = [PhyloNetworks], # to list methods from PhyloNetworks only, not from Base etc.
         format = Markdown(),
         Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true") # to make it easier when building locally, if using HTML (not markdown) format
        )

# Versions of `mkdocs` and `mkdocs-material` specified manually to avoid conflicts. 
# See here: https://discourse.julialang.org/t/mkdocs-material-in-documenter/13764/3
# To be kept in mind: those versions might evolve in the future.
deploydocs(
    repo = "github.com/cecileane/PhyloNetworks.jl.git",
    #deps= Deps.pip("pygments", "mkdocs==0.17.5", "mkdocs-material==2.9.4", "python-markdown-math"),
    deps = Deps.pip("pygments", "mkdocs==1.0.4",  "mkdocs-material==3.2.0", "python-markdown-math"),
    # deps = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math"),
    make = () -> run(`mkdocs build`),
    target = "site" # which files get copied to gh-pages
)
