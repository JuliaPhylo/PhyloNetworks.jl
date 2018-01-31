# parsimony tools: by William Sparks, 2017

"""
`parsimonyBottomUp!(node, states, score)`

Bottom-up phase (from tips to root) of the Fitch algorithm:
assign sets of character states to internal nodes based on
character states at tips. Polytomies okay.
Assumes a *tree* (no reticulation) and correct isChild1 attribute.

output: dictionary with state sets and most parsimonious score
"""

function parsimonyBottomUp!{T}(node::Node, possibleStates::Dict{Int64,Set{T}}, parsimonyscore::Array{Int64,1})
    node.leaf && return # change nothing if leaf
    childrenStates = Set{T}[] # state sets for the 2 (or more) children
    for e in node.edge
        if e.node[e.isChild1 ? 1 : 2] == node continue; end
        # excluded parent edges only: assuming tree here
        child = getOtherNode(e, node)
        parsimonyBottomUp!(child, possibleStates, parsimonyscore)
        if haskey(possibleStates, child.number) # false if missing data
            push!(childrenStates, possibleStates[child.number])
        end
    end
    if length(childrenStates)==0 return; end # change nothing if no data below

    votes = Dict{T,Int}() # number of children that vote for a given state
    for set in childrenStates
      for i in set
        if haskey(votes,i)
          votes[i] += 1
        else
          votes[i] = 1
        end
      end
    end
    mv = maximum(values(votes))
    filter!((k,v) -> v==mv, votes) # keep states with max votes
    possibleStates[node.number] = Set(keys(votes))
    parsimonyscore[1] += length(childrenStates) - mv # extra cost
end

"""
`parsimonyTopDown!(node, states)`

Top-down phase (root to tips) of the Fitch algorithm:
constrains character states at internal nodes based on
the state of the root. Assumes a *tree*: no reticulation.

output: dictionary with state sets
"""

function parsimonyTopDown!{T}(node::Node, possibleStates::Dict{Int64,Set{T}})
    for e in node.edge
        child = e.node[e.isChild1 ? 1 : 2]
        if child == node continue; end # exclude parent edges
        if child.leaf continue; end    # no changing the state of tips
        commonState = intersect(possibleStates[node.number], possibleStates[child.number])
        if !isempty(commonState)
            possibleStates[child.number] = Set(commonState) # restrain states to those with cost 0
        end
        parsimonyTopDown!(child, possibleStates)
    end
    return possibleStates # updated dictionary of sets of character states
end

"""
`parsimonySummary(tree, nodestates)`

summarize character states at nodes, assuming a *tree*
"""

function parsimonySummary{T}(tree::HybridNetwork, nodestates::Dict{Int64,Set{T}})
    println("node number => character states on tree ",
            writeTopology(tree,di=true,round=true,digits=1))
    for n in tree.node
        if !haskey(nodestates, n.number) continue; end
        print(n.number)
        if n.name != "" print(" (",n.name,")"); end
        if n == tree.node[tree.root] print(" (root)"); end
        println(": ", sort(collect(nodestates[n.number])))
    end
end

"""
`parsimonyDiscrete(net, tipdata)`

Calculate the most parsimonious (MP) score of a network given
a discrete character at the tips.
The softwired parsimony concept is used: where the number of state
transitions is minimized over all trees displayed in the network.
Tip data can be given in a data frame, in which case the taxon names
are to appear in column 1 or in a column named "taxon" or "species", and
trait values are to appear in column 2 or in a column named "trait".
Alternatively, tip data can be given as a dictionary taxon => trait.

also return the union of all optimized character states
at each internal node as obtained by Fitch algorithm,
where the union is taken over displayed trees with the MP score.
"""

function parsimonyDiscrete{T}(net::HybridNetwork, tips::Dict{String,T})
    # T = type of characters. Typically Int if data are binary 0-1
    # initialize dictionary: node number -> admissible character states
    possibleStates = Dict{Int64,Set{T}}()
    for l in net.leaf
        if haskey(tips, l.name)
            possibleStates[l.number] = Set(tips[l.name])
        end
    end
    charset = union(possibleStates) # fixit
    # assign this set to all tips with no data

    directEdges!(net) # parsimonyBottomUp! uses isChild1 attributes
    trees = displayedTrees(net, 0.0) # all displayed trees
    mpscore = Int[] # one score for each tree
    statesets = Dict{Int64,Set{T}}[] # one state set dict per tree
    for tree in trees
        statedict = deepcopy(possibleStates)
        parsimonyscore = [0] # initialization, mutable
        parsimonyBottomUp!(tree.node[tree.root], statedict, parsimonyscore)
        push!(mpscore, parsimonyscore[1])
        push!(statesets, statedict)
    end
    mps = findmin(mpscore)[1] # MP score
    mpt = find(x -> x==mps, mpscore) # indices of all trees with MP score
    statedictUnion = statesets[mpt[1]] # later: union over all MP trees
    println("parsimony score: ", mps)
    for i in mpt # top down calculation for best trees only
        parsimonyTopDown!(trees[i].node[trees[i].root], statesets[i])
        parsimonySummary(trees[i], statesets[i])
        if i == mpt[1] continue; end
        for n in keys(statesets[i])
          if haskey(statedictUnion, n) # degree-2 nodes absent from trees
            union!(statedictUnion[n], statesets[i][n])
          else
            statedictUnion[n] = statesets[i][n]
          end
        end
        # fixit: for each hybrid edge, count the number of MP trees that have it,
        #        to return information on which hybrid edge is most parsimonious
    end
    return mps, statedictUnion
end

function parsimonyDiscrete(net::HybridNetwork, dat::DataFrame)
    i = findfirst(DataFrames.names(dat), :taxon)
    if i==0 i = findfirst(DataFrames.names(dat), :species); end
    if i==0 i=1; end # first column if not column named "taxon" or "species"
    j = findfirst(DataFrames.names(dat), :trait)
    if j==0 j=2; end
    if i==j
        error("""expecting taxon names in column 'taxon', or 'species' or column 1,
              and trait values in column 'trait' or column 2.""")
    end
    tips = Dict{String,eltypes(dat)[j]}()
    for r in 1:nrow(dat)
        if DataFrames.isna(dat[r,j]) continue; end
        tips[dat[r,i]] = dat[r,j]
    end
    parsimonyDiscrete(net,tips)
end


function generateSwitchingVectors(n, switchingVectors::Array{Bool,2} ,prefix=""::String) #fixit avoid to create all the vectors using myabe products and iterators
    if (!n)
        i=1
        switching = Vector{Bool}(0)
        while i <= length(prefix)
            if prefix.ch[i]=='0'
                push!(switching, false)
            else
                push!(switching, true)
            end
            i+=1
        end
       push!(switchingVectors, switching)
       return nothing
    end
    generateSwitchingVectors(n-1, switchingVectors, prefix + '0')
    generateSwitchingVectors(n-1, switchingVectors, prefix + '1')
    return nothing
end


"""
`parsimonyBottomUpWeight!(node, states, w, scores)`

Computing the MP scores (one per each assigemenet assignment of the root state)
of a swicthing as described in Algorithm 1 in the following paper:
Fischer, M., van Iersel, L., Kelk, S., Scornavacca, C. On computing the Maximum
Parsimony score of a phylogenetic network. SIDMA, 29 (1), pp 559 - 585.
Polytomies okay.

Assumes a *switching* (ie correct fromBadDiamnodI field) and correct isChild1 field.

"""

function parsimonyBottomUpWeight!{T}(node::Node, tips::Dict{String,Set{T}},
    charset::Set{T}, w::Array{Int64,2}, parsimonyscore::Array{Int64,2})

    for s in charset
        if node.leaf
            if s== tips[node.name]
                parsimonyscore[node.number][s] = 0
            else
                parsimonyscore[node.number][s] = Inf
            end
        else
            parsimonyscore[node.number][s] = w[node.number][s]
            if ! node.isExtBadTriangle # isExtBadTriangle=dummy leaf in the algorithm
                for e in node.edge
                    bestMin =Inf
                    if e.fromBadDiamnodI || e.node[e.isChild1 ? 1 : 2] == node continue; end #fromBadDiamnodI= off edge
                    son=getOtherNode(e, node)
                    parsimonyBottomUpWeight!(son, possibleStates, charset, w, parsimonyscore)
                    for sf in charset # best assignement for the son
                        min= parsimonyscore[son.number][sf]
                        if s !=sf
                            min +=1 #fixit (celine) code it with cost of change delta(sf,s)
                        end
                        bestMin = min(min, bestMin)
                    end
                    parsimonyscore[node.number][s] += bestMin  # add best assignement for the son to the PS of the parent
                end
            end
        end
    end
    return nothing
end


"""
`parsimonyDiscreteDP(net, tipdata)`

Calculate the most parsimonious (MP) score of a network given
a discrete character at the tips using the dynamic programming algorithm given
in the paper
Fischer, M., van Iersel, L., Kelk, S., Scornavacca, C. On computing the Maximum
Parsimony score of a phylogenetic network. SIDMA, 29 (1), pp 559 - 585.
The softwired parsimony concept is used: where the number of state
transitions is minimized over all trees displayed in the network.
Tip data can be given in a data frame, in which case the taxon names
are to appear in column 1 or in a column named "taxon" or "species", and
trait values are to appear in column 2 or in a column named "trait".
Alternatively, tip data can be given as a dictionary taxon => trait.

"""

function parsimonyDiscreteDP{T}(net::HybridNetwork, tips::Dict{String,T})
    # T = type of characters. Typically Int if data are binary 0-1
    # initialize dictionary: node number -> admissible character states

    charset = union(tips) # fixit
    # assign this set to all tips with no data

    w = zeros(Int,(length(net.node), length(charset)))
    parsimonyscore = zeros(Int,(length(net.node), length(charset)))

    blobroots, majorEdges, minorEdges = blobInfo!(net) # calls directEdges! by default: sets isChild1

    bcnumber=1;
    for r in blobroots
        mpscoreSwitchings = Array{Int64,2}
        switchingVectors = Vector{Bool}(0)(0)
        generateSwitchingVectors(length(majorEdges[bcnumber]), switchingVectors ) # length(majorEdges) is the level of our biconnected component

        for switching in switchingVectors
            switchingEdge=1;
            while switchingEdge <= length(majorEdges[bcnumber]) #assigning on off to the hybid edges of the blob w.r.t. the switching[switchingEdge]
                if switching[switchingEdge]
                    majorEdges[bcnumber][switchingEdge].fromBadDiamnodI=true
                    minorEdges[bcnumber][switchingEdge].fromBadDiamnodI=false
                else
                    majorEdges[bcnumber][switchingEdge].fromBadDiamnodI=false
                    minorEdges[bcnumber][switchingEdge].fromBadDiamnodI=true
                end
                switchingEdge+=1
            end
            push!(mpscoreSwitchings, parsimonyBottomUpWeight!{T}(r, charset, tips, w, parsimonyscore))
        end

        for s in charset # best assignement for the son
            switchingNumber=1;
            bestMin=Inf
            while switchingNumber <= length(mpscoreSwitchings)
                min= mpscoreSwitchings[switchingNumber][s]
                bestMin = min(min, bestMin)
                switchingNumber+=1
            end
            w[r.number][s] += bestMin  # add best assignement for the son to the PS of the parent
        end
        bcnumber+=1
    end

    bestMin=Inf
    for s in charset
        min= w[net.root.number][s]
        bestMin = min(min, bestMin)
    end
    return  bestMin
end



function parsimonyDiscreteDP(net::HybridNetwork, dat::DataFrame)
    i = findfirst(DataFrames.names(dat), :taxon)
    if i==0 i = findfirst(DataFrames.names(dat), :species); end
    if i==0 i=1; end # first column if not column named "taxon" or "species"
    j = findfirst(DataFrames.names(dat), :trait)
    if j==0 j=2; end
    if i==j
        error("""expecting taxon names in column 'taxon', or 'species' or column 1,
              and trait values in column 'trait' or column 2.""")
    end
    tips = Dict{String,eltypes(dat)[j]}()
    for r in 1:nrow(dat)
        if DataFrames.isna(dat[r,j]) continue; end
        tips[dat[r,i]] = dat[r,j]
    end
    parsimonyDiscreteDP(net,tips)
end
