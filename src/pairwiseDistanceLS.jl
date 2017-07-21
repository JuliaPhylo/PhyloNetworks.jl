"""
    setBLfromDivergenceTimes!(net::HybridNetwork, divergenceTimes::Dict{Int64,Float64})

Modify branch lengths in `net` to match divergence times given in `divergenceTimes`.

# Examples

"""


function setBLfromDivergenceTimes!(net::HybridNetwork, divergenceTimes::Dict{Int64,Float64})
    for e in net.edge
        parentTime = divergenceTimes[e.node[e.isChild1 ? 2 : 1].number]
        childTime = divergenceTimes[e.node[e.isChild1 ? 1 : 2].number]
        e.length = parentTime - childTime
    end
end            

"""

    parwiseTaxonDistanceMatrix(net::HybridNetwork, [checkPreorder=true::Bool])

Create an n-by-n matrix (where n=# taxa in `net`) of pairwise distances between taxa in `net`.

"""

function pairwiseTaxonDistanceMatrix(net::HybridNetwork;
                                     checkPreorder=true::Bool)
    recursionPreOrderMatrixOnly(net,
                      checkPreorder,
                      initsharedPathMatrix,
                      updateRootSharedPathMatrix!,
                      updateTreePairwiseTaxonDistanceMatrix!,
                      updateHybridPairwiseTaxonDistanceMatrix!,
                      "b")
end

function recursionPreOrderMatrixOnly(net::HybridNetwork,
                           checkPreorder=true::Bool,
                           init=identity::Function,
                           updateRoot=identity::Function,
                           updateTree=identity::Function,
                           updateHybrid=identity::Function,
                           indexation="b"::AbstractString,
                           params...)
    net.isRooted || error("net needs to be rooted for preorder recursion")
    if(checkPreorder)
        preorder!(net)
    end
    M = recursionPreOrder(net.nodes_changed, init, updateRoot, updateTree, updateHybrid, params)
end

function updateTreePairwiseTaxonDistanceMatrix!(V::Matrix,
                                     i::Int,
                                     parentIndex::Int,
                                     edge::PhyloNetworks.Edge,
                                     params)
    for j in 1:(i-1)
        V[i,j] = V[parentIndex,j]+edge.length
        V[j,i] = V[i,j]
    end
    V[i,i] = 0.0
end

function updateHybridPairwiseTaxonDistanceMatrix!(V::Matrix,
                                       i::Int,
                                       parentIndex1::Int,
                                       parentIndex2::Int,
                                       edge1::PhyloNetworks.Edge,
                                       edge2::PhyloNetworks.Edge,
                                       params)
    for j in 1:(i-1)
        V[i,j] = edge1.gamma*(edge1.length+V[parentIndex1,j])+edge2.gamma*(edge2.length+V[parentIndex2,j])
        V[j,i] = V[i,j]
    end
end

"""
    pairwiseDistanceLSscore(net::HybridNetwork, divergenceTimes::Dict{Int64,Float64}, dnaDistances::Array{Float64,2})

Produce mismatch between the network distances and the observed distances from DNA (`dnaDistances`).
"""

function pairwiseDistanceLSscore(net::HybridNetwork, divergenceTimes::Dict{Int64,Float64}, dnaDistances::Array{Float64,2})
    setBLfromDivergenceTimes!(net, divergenceTimes)
    networkDistances = pairwiseTaxonDistanceMatrix(net;
                                checkPreorder=true)
    distanceMismatch = sqrt(sum(sum((dnaDistances-networkDistances).^2)))
    return distanceMismatch
end
