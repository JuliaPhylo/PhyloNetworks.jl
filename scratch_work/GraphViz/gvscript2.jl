using GraphViz

println("Opening .dot file and storing as variable")
dot = open("graph.dot","r") do io Graph(io) end
println(".dot file stored as variable")
println("Calling layout")
GraphViz.layout!(dot,engine="dot")
println("Starting file write.")
open("graph.svg","w") do f
	GraphViz.writemime(f, MIME"image/svg+xml"(),dot)
end #do