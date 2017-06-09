using PhyloNetworks
import PhyloNetworks.getOtherNode
import PhyloNetworks.getIndex

net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1)

"""
pushes leaf numbers and character states to nodeStates dictionary
"""

function initializeStates(tips::Dict{String,Int64}, possibleStates::Dict{Int64,Any})
    for l in net.leaf
        for k in keys(tips)
            l.name == k ? possibleStates[l.number] = Set(tips[k]) : continue # if the name of a leaf node matches a key in the dictionary of tip states, assign key and corresponding value to possibleStates dictionary
        end
    end
end

"""
performs the bottom-up phase of the fitch algorithm: assigns character states to internal nodes based on character states at tips; calculates parsimony score
"""

function parsimonyBottomUp!(node::PhyloNetworks.Node, possibleStates::Dict{Int64,Any}, parsimonyscore::Array{Int64,1})
    children = [] #initialize array of child nodes
    childrenStates = [] #initialize array of values for child keys

    for e in node.edge #iterate across all edges for a given node
        if e.node[e.isChild1 ? 1 : 2] != node #exclude parent edges
            append!(children, [getOtherNode(e, node)]) #push child nodes to children array
        end
    end

    for child in children # iterate over children array
        if haskey(possibleStates, child.number) # if all child nodes already have a character state assigned:
            append!(childrenStates, possibleStates[child.number]) # append the states of the child nodes to children values array
        else # if one or more children do not have states assigned:
            parsimonyBottomUp!(child, possibleStates, parsimonyscore) #recursively call parsimony function on children
            append!(childrenStates, possibleStates[child.number])
        end
    end

    if node.leaf == false
        votes = Dict{Int,Int}()
        for set in childrenStates
            for i in set
                if haskey(votes,i)
                    votes[i] += 1
                else
                    votes[i] = 1
                end
            end
        end

        mv = 0 # max number of votes
        for state in keys(votes)
            if votes[state]>mv
                mv = votes[state]
            end
        end

        filter!((k,v) -> v==mv, votes) # filters d: keep states with mv votes
        possibleStates[node.number] = Set(keys(votes)) # {1,2}: set of states for the parent node
        parsimonyscore[1] += length(childrenStates) - mv # 2: parsimony cost
    end

    if haskey(possibleStates, net.node[net.root].number) # do nothing if the current node is also the root
        return possibleStates, parsimonyscore[1]
    end
    return possibleStates, parsimonyscore[1]
end

"""
performs top-down phase of Fitch Algorithm; assigns character states to internal nodes based on arbitrary state of the root; calculates parsimony
"""

function parsimonyTopDown!(node::PhyloNetworks.Node, possibleStates::Dict{Int64,Any}, treecost::Array{Int64,1})
    for e in node.edge
        if e.node[e.isChild1 ? 1 : 2] != node
            # continue if current edge is not the parent of current node
            child = getOtherNode(e, node) # define child nodes of current node
            childState = possibleStates[child.number]
            parentState = possibleStates[node.number]

            commonState = intersect(parentState, childState) # define commonstate as the intersection of the values of all children; see function listintersect
            if isempty(commonState) # if parent and child nodes have no common character state
                treecost[1] += 1 # increase parsimony score of tree
            else # if parent and child nodes have one or more states in common
                if !child.leaf
                    possibleStates[child.number] = Set(commonState) # change character state of child to that of parent; do not change character states for leaves
                end
            end

            if child.leaf
                continue
            else
                parsimonyTopDown!(child, possibleStates, treecost) # recursively call function on child of current node
            end
        end
    end
    return possibleStates, treecost[1] # return updated character states dictionary corresponding parsimony score
end

"""
summarize character states at nodes and parsimony scores
"""

function parsimonySummary(tree::HybridNetwork, rootSpecificStates::Dict{Int64,Any}, parsimonyscore::Array{Int64,1}, treecost::Array{Int64,1})
    println("\nparsimonyscore: ", parsimonyscore[1]) #print parsimony score for a given tree
    println("Node Number => Character States \n")
    for (k,v) in rootSpecificStates #iterate over keys and values
        if k == net.node[net.root].number
            println("root:") #print key, value pair for root node
            println(k => v)
        end
    end

    println("\ntips:") #print k,v for tips
    for (k,v) in rootSpecificStates
        for l in net.leaf
            if l.number == k
                println(k => v)
            end
        end
    end

    println("\ninternal nodes:") #print k,v for internal nodes
    for (k,v) in rootSpecificStates
        if k != net.node[net.root].number && issubset(k, [l.number for l in net.leaf]) == false
            println(k => v)
        end
    end
    println("\ntreecost: ", treecost[1]) #print tree cost for a given root state
end

"""
calculates parsimony score and optimized character states for a binary character on a PhyloNetworks.HybridNetwork given the character states at the tips
"""

function parsimony(net::HybridNetwork, tips::Dict{String,Int64})
    trees = displayedTrees(net, 0.0) #create list of all possible trees in a given network

    for tree in trees #iterate over all trees in a network
        possibleStates = Dict{Int64,Any}() # initialize character states dictionary
        parsimonyscore = [0] # initialize parsimony score

        directEdges!(net) #orders hybrid edges; necesary before preorder!, cladewiseorder!
        preorder!(net)  #creates nodes_changed vector: nodes arranged for preorder traversal

        initializeStates(tips, possibleStates) #create dictionary for states of all nodes

        parsimonyBottomUp!(tree.node[tree.root], possibleStates, parsimonyscore) #calculate parsimony score and possible states at a given node

        for state in possibleStates[tree.node[tree.root].number] #iterate over possible root states
            rootSpecificStates = possibleStates #create a new dictionary that will contain optimized states for a given state at the root
            rootSpecificStates[tree.node[tree.root].number] = Set(state) # set root state
            treecost = [0] #initialize treecost

            parsimonyTopDown!(tree.node[tree.root], rootSpecificStates, treecost) #calculate optimized states and relative cost for a given root state on a given tree

            parsimonySummary(tree::HybridNetwork, rootSpecificStates::Dict, parsimonyscore::Array{Int64,1}, treecost::Array{Int64,1}) #summarize results
        end
    end
end

parsimony(net, tips)
