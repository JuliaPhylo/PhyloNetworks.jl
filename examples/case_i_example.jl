# Case I example: Bad diamond II
# Claudia January 2015


# types in "types.jl"
if !isdefined(:individualtest) individualtest = false; end

if(individualtest)
    include("../src/types.jl")
    include("../src/functions.jl")
end

tree = "((((8,10))#H1,7),6,(4,#H1));" # Case I Bad diamond II
#f = open("prueba_tree.txt","w")
#write(f,tree)
#close(f)
net = readTopologyUpdate(tree)

#printEdges(net)
#printNodes(net)
