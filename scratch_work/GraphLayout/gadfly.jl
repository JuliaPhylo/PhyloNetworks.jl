# Plot the dependency graph for Gadfly.jl
using GraphLayout

# Load adjacency matrix from file
fp = open("gadfly.txt","r")
adj_matrix = readdlm(fp,',')
close(fp)

# Load node labels
fp = open("gadfly_names.txt","r")
labels = map(chomp,readlines(fp))
close(fp)

# Calculate a layout
srand(4)
println("Layout...")
loc_x, loc_y = layout_spring_adj(adj_matrix)

# Draw it
println("Drawing...")
draw_layout_adj(adj_matrix, loc_x, loc_y, labels=labels, filename="gadfly.svg")


adj_list = Vector{Int}[]
for i in 1:size(adj_matrix,1)
    new_list = Int[]
    for j in 1:size(adj_matrix,2)
        if adj_matrix[i,j] != zero(eltype(adj_matrix))
            push!(new_list,j)
        end
    end
    push!(adj_list, new_list)
end

GraphLayout.layout_tree(adj_list, labels, cycles=false, filename="gadfly_tree.svg")