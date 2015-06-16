using GraphViz
g = open("graph.dot","r") do io Graph(io) end
GraphViz.layout!(g)
GraphViz.render_x11(GraphViz.Context(),g)