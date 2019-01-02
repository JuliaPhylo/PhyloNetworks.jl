using Documenter, DocumenterMarkdown

using Pkg
Pkg.add(PackageSpec(name="PhyloPlots", rev="master"))

using PhyloNetworks

makedocs(
    sitename = "PhyloNetworks.jl",
    authors = "Claudia Solís-Lemus, Cécile Ané, Paul Bastide and contributors.",
    modules = [PhyloNetworks], # to list methods from PhyloNetworks only, not from Base etc.
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"), # easier local build
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Installation" => "man/installation.md",
            "Input Data for SNaQ" => "man/inputdata.md",
            "TICR pipeline" => "man/ticr_howtogetQuartetCFs.md",
            "Network estimation and display" => "man/snaq_plot.md",
            "Network comparison and manipulation" => "man/dist_reroot.md",
            "Candidate Networks" => "man/fixednetworkoptim.md",
            "Extract Expected CFs" => "man/expectedCFs.md",
            "Bootstrap" => "man/bootstrap.md",
            "Multiple Alleles" => "man/multiplealleles.md",
            "Continuous Trait Evolution" => "man/trait_tree.md",
            "Parsimony on networks" => "man/parsimony.md",
        ],
        "Library" => [
            "Public" => "lib/public.md",
            "Internals" => "lib/internals.md",
        ]
    ],
)

# Versions of `mkdocs` and `mkdocs-material` specified manually to avoid conflicts. 
# See here: https://discourse.julialang.org/t/mkdocs-material-in-documenter/13764/3
# To be kept in mind: those versions might evolve in the future.
deploydocs(
    repo = "github.com/cecileane/PhyloNetworks.jl.git"
    #,deps= Deps.pip("pygments", "mkdocs==0.17.5", "mkdocs-material==2.9.4", "python-markdown-math"),
    #make = () -> run(`mkdocs build`),
    #target = "site" # which files get copied to gh-pages
)
