#Scratch file (building graphs)
#Creating graphs using CFnetworks functions and connecting with external packages
#John Spaw 6-8-15

#include("types.jl")
#include("functions.jl")

e1 = Edge(1,5.0);
e4 = Edge(4,2.0);

n1 = Node(1, false, false, [e1, e4]);
n2 = Node(2, false, false, [e1]);
n5 = Node(5, true, false, [e4]);


setNode!(e1,[n1,n2]);
setNode!(e4,[n1,n5]);


test_net=HybridNetwork([n1,n2,n5],[e1,e4]);

printNodes(test_net)

for i=1:2
    println(" ")
end

printEdges(test_net)

for i=1:2
    println(" ")
end

using LightGraphs
g = LightGraphs.Graph(5)
add_edge!(g,1,2)
add_edge!(g,1,3)

am = full(adjacency_matrix(g))

println(am)

using GraphLayout
loc_x, loc_y = GraphLayout.layout_spring_adj(am)
draw_layout_adj(am, loc_x, loc_y, filename = "Testgraph.svg")


