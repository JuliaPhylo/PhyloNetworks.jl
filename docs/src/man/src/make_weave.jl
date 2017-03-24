using Weave, Fontconfig, Cairo

set_chunk_defaults(Dict{Symbol, Any}(:results => "hidden", :eval => false))

weave(Pkg.dir("PhyloNetworks","docs","src", "man", "src","trait_tree.jmd"),
      informat="markdown",
      out_path = Pkg.dir("PhyloNetworks","docs","src", "man"),
      fig_path = "../assets/figures",
      fig_ext = ".png",
      doctype = "github")

## Delete figures of not built
draw_fig = get(ENV, "DRAW_FIG", "")
draw_fig = !(draw_fig == "false")
travis_tag = get(ENV, "TRAVIS_TAG", "")
# If tagged version (to stable), always draw figures
if travis_tag != ""
    draw_fig = true
end
if !draw_fig
    # rm the figures if not drawn.
    rm(Pkg.dir("PhyloNetworks","docs","src", "assets", "figures");
       force=true, recursive=true)
end
