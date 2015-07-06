#John Spaw
#This script includes all of the important base files for starting a Julia session

include("types.jl")
include("functions.jl")
include("scratch_work/Misc/test_graphs.jl")
include("visualization/traverseEdges.jl")
include("visualization/finalDraw.jl")

net = create_g3();

