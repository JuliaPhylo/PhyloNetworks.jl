using Weave

Pkg.dir("PhyloNetworks","docs","src", "man", "src", "make_weave.jl")

weave(Pkg.dir("PhyloNetworks","docs","src", "man", "src","trait_tree.jmd"),
      informat="markdown",
      out_path = Pkg.dir("PhyloNetworks","docs","src", "man"),
      fig_ext = ".png",
      doctype = "github")

# weave("docs/src/man/src/trait_tree.jmd", informat="markdown",
    #   out_path = "docs/src/man/", doctype = "github")
