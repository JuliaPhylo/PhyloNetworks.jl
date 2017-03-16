using Weave, Fontconfig, Cairo

set_chunk_defaults(Dict{Symbol, Any}(:results => "hidden", :eval => false))

weave(Pkg.dir("PhyloNetworks","docs","src", "man", "src","trait_tree.jmd"),
      informat="markdown",
      out_path = Pkg.dir("PhyloNetworks","docs","src", "man"),
      fig_path = "../assets/figures",
      fig_ext = ".png",
      doctype = "github")

