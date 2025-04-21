using Documenter

# this installs the dev version of PhyloPlots for compatibility.
using Pkg
Pkg.add(PackageSpec(name="PhyloPlots", rev="master"))


using PhyloNetworks

# interlink with other packages, for @ref calls to become "external ref"
using DocumenterInterLinks
links = InterLinks(
    "PhyloPlots" => "https://juliaphylo.github.io/PhyloPlots.jl/stable/objects.inv",
    "SNaQ" => "https://juliaphylo.github.io/SNaQ.jl/stable/objects.inv",
    "QGoF" => "https://juliaphylo.github.io/QuartetNetworkGoodnessFit.jl/stable/objects.inv",
    # "PhyloTraits"=> "https://juliaphylo.github.io/PhyloTraits.jl/stable/objects.inv",
)
# default loading of interlinked packages in all docstring examples
DocMeta.setdocmeta!(PhyloNetworks, :DocTestSetup,
    :(using PhyloNetworks, PhyloPlots); # , PhyloTraits, SNaQ
    recursive=true)
using PhyloPlots # to trigger any precompilation warning outside jldoctests

makedocs(
    sitename = "PhyloNetworks.jl",
    authors = "Claudia Solís-Lemus, Cécile Ané, Paul Bastide and contributors.",
    modules = [PhyloNetworks], # to list methods from PhyloNetworks only, not from Base etc.
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true", # easier local build
        size_threshold = 600 * 2^10,
        size_threshold_warn = 500 * 2^10, # 600 KiB
        canonical="https://juliaphylo.github.io/PhyloNetworks.jl/stable/",
        edit_link="master",
    ),
    # exception, so warning-only for :cross_references. List all others:
    warnonly = Documenter.except(:autodocs_block, :missing_docs, # :cross_references,
        :docs_block,
        :doctest, :eval_block, :example_block, :footnote, :linkcheck_remotes,
        :linkcheck, :meta_block, :parse_error, :setup_block),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Installation" => "man/installation.md",
            "Introduction" => "man/introduction.md",
            "Network manipulation" => "man/netmanipulation.md",
            "Network visualization" => "man/net_plot.md",
            "Network comparison and manipulation" => "man/dist_reroot.md",
            "Network support" => "man/network_support.md",
            "Parsimony on networks" => "man/parsimony.md",
            "Neighbour Joining" => "man/nj.md",
            "Pairwise and quartet data" => "man/expecteddata.md",
        ],
        "Library" => [
            "Public" => "lib/public.md",
            "Internals" => "lib/internals.md",
        ]
    ],
)

deploydocs(
    repo = "github.com/JuliaPhylo/PhyloNetworks.jl.git",
    push_preview = true,
    devbranch = "master",
)
