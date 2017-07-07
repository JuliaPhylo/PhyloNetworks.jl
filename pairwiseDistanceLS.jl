"""
Sets edges of HybridNetwork based on values in divergenceTimes dictionary;
keys must be node numbers;
values must be integers representing divergence time (mya) at a given node  
"""

function setBLfromDivergenceTimes(net::HybridNetwork, divergenceTimes::Dict{Int64,Int64})
    for e in net.edge
        parent = 0
        child = 0
        for time in keys(divergenceTime)
            if e.node[e.isChild1 ? 2 : 1].number == time
                parent = divergenceTime[time]
            end
            if e.node[e.isChild1 ? 1 : 2].number == time
                child = divergenceTime[time]
            end
        end
        e.length = parent - child
    end
end            

"""
Creates a MatrixTopologicalOrder of distances from each node to all other nodes
"""

function pairwiseTaxonDistanceMatrix(net::HybridNetwork;
                                     checkPreorder=true::Bool)
    recursionPreOrder(net,
                      checkPreorder,
                      initPairwiseTaxonDistanceMatrix,
                      updateRootPairwiseTaxonDistanceMatrix!,
                      updateTreePairwiseTaxonDistanceMatrix!,
                      updateHybridPairwiseTaxonDistanceMatrix!,
                      "b")
end

"""
initializes taxon distance matrix
"""
function initPairwiseTaxonDistanceMatrix(nodes::Vector{PhyloNetworks.Node}, params)
    n = length(nodes)
    return(zeros(Float64,n,n))
end

function updateRootPairwiseTaxonDistanceMatrix!(V::Matrix, i::Int, params)
    return
end

"""
calculates distance from tree nodes
"""

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

"""
calculates distance from hybrid nodes
"""

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
calcutates mismatch between pairwise distances observed between DNA sequences and pairwised distances observed from a network
"""

function pairwiseDistanceMismatch(networkDistances::PhyloNetworks.MatrixTopologicalOrder, dnaDistances::PhyloNetworks.MatrixTopologicalOrder)
    distanceMismatch = sqrt(sum(sum((dnaDistances.V-networkDistances.V).^2)))
    return distanceMismatch
end

"""
calcutates mismatch between pairwise distances observed between DNA sequences and pairwised distances observed from a network given a HybridNetwork, a dictionary of divergence times at each node, and a pairwise DNA sequences distance matrix
"""

function pairwiseDistanceLSscore(net::HybridNetwork, divergenceTimes::Dict{Int64,Int64}, dnaDistances::PhyloNetworks.MatrixTopologicalOrder)
    setBLfromDivergenceTimes(net::HybridNetwork, divergenceTime::Dict{Int64,Int64})
    pairwiseTaxonDistanceMatrix(net::HybridNetwork;
                                checkPreorder=true::Bool)
    distanceMismatch = pairwiseDistanceMismatch(networkDistances::PhyloNetworks.MatrixTopologicalOrder, dnaDistances::PhyloNetworks.MatrixTopologicalOrder)
    return distanceMismatch
end

