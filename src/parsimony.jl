
"""
`parsimonyBottomUpFitch!(node, states, score)`

Bottom-up phase (from tips to root) of the Fitch algorithm:
assign sets of character states to internal nodes based on
character states at tips. Polytomies okay.
Assumes a *tree* (no reticulation) and correct isChild1 attribute.

output: dictionary with state sets and most parsimonious score
"""

function parsimonyBottomUpFitch!(node::Node, possibleStates::Dict{Int64,Set{T}}, parsimonyscore::Array{Int64,1}) where {T}
    node.leaf && return # change nothing if leaf
    childrenStates = Set{T}[] # state sets for the 2 (or more) children
    for e in node.edge
        if e.node[e.isChild1 ? 1 : 2] == node continue; end
        # excluded parent edges only: assuming tree here
        child = getOtherNode(e, node)
        parsimonyBottomUpFitch!(child, possibleStates, parsimonyscore)
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
`parsimonyTopDownFitch!(node, states)`

Top-down phase (root to tips) of the Fitch algorithm:
constrains character states at internal nodes based on
the state of the root. Assumes a *tree*: no reticulation.

output: dictionary with state sets
"""

function parsimonyTopDownFitch!(node::Node, possibleStates::Dict{Int64,Set{T}}) where {T}
    for e in node.edge
        child = e.node[e.isChild1 ? 1 : 2]
        if child == node continue; end # exclude parent edges
        if child.leaf continue; end    # no changing the state of tips
        commonState = intersect(possibleStates[node.number], possibleStates[child.number])
        if !isempty(commonState)
            possibleStates[child.number] = Set(commonState) # restrain states to those with cost 0
        end
        parsimonyTopDownFitch!(child, possibleStates)
    end
    return possibleStates # updated dictionary of sets of character states
end

"""
`parsimonySummaryFitch(tree, nodestates)`

summarize character states at nodes, assuming a *tree*
"""

function parsimonySummaryFitch(tree::HybridNetwork, nodestates::Dict{Int64,Set{T}}) where {T}
    println("node number => character states on tree ",
            writeTopology(tree,di=true,round=true,digits=1))
    for n in tree.node
        haskey(nodestates, n.number) || continue
        print(n.number)
        if n.name != "" print(" (",n.name,")"); end
        if n == tree.node[tree.root] print(" (root)"); end
        println(": ", sort(collect(nodestates[n.number])))
    end
end

"""
    parsimonyDiscreteFitch(net, tipdata)

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

function parsimonyDiscreteFitch(net::HybridNetwork, tips::Dict{String,T}) where {T}
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

    directEdges!(net) # parsimonyBottomUpFitch! uses isChild1 attributes
    trees = displayedTrees(net, 0.0) # all displayed trees
    mpscore = Int[] # one score for each tree
    statesets = Dict{Int64,Set{T}}[] # one state set dict per tree
    for tree in trees
        statedict = deepcopy(possibleStates)
        parsimonyscore = [0] # initialization, mutable
        parsimonyBottomUpFitch!(tree.node[tree.root], statedict, parsimonyscore)
        push!(mpscore, parsimonyscore[1])
        push!(statesets, statedict)
    end
    mps = findmin(mpscore)[1] # MP score
    mpt = find(x -> x==mps, mpscore) # indices of all trees with MP score
    statedictUnion = statesets[mpt[1]] # later: union over all MP trees
    # println("parsimony score: ", mps)
    for i in mpt # top down calculation for best trees only
        parsimonyTopDownFitch!(trees[i].node[trees[i].root], statesets[i])
        parsimonySummaryFitch(trees[i], statesets[i])
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

function parsimonyDiscreteFitch(net::HybridNetwork, dat::DataFrame)
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
    parsimonyDiscreteFitch(net,tips)
end

"""
`parsimonyBottomUpSoftwired!(node, blobroot, states, w, scores)`

Computing the MP scores (one for each assignment of the root state)
of a swicthing as described in Algorithm 1 in the following paper:

Fischer, M., van Iersel, L., Kelk, S., Scornavacca, C. (2015).
On computing the Maximum Parsimony score of a phylogenetic network.
SIAM J. Discrete Math., 29(1):559-585.

Assumes a *switching* (ie correct `fromBadDiamondI` field) and correct isChild1 field.
The field `isExtBadTriangle` is used to know which nodes are at the root of a blob.
"""

function parsimonyBottomUpSoftwired!(node::Node, blobroot::Node, nchar::Integer,
    w::AbstractArray, parsimonyscore::AbstractArray)

    #println("entering with node $(node.number)")
    parsimonyscore[node.number,:] = w[node.number,:] # at all nodes to re-initialize between switchings
    if node.leaf || (node.isExtBadTriangle && node != blobroot)
        return nothing # isExtBadTriangle=dummy leaf: root of another blob
    end
    for e in node.edge
        if (e.hybrid && e.fromBadDiamondI) || getChild(e) == node continue; end # fromBadDiamondI= edge switched off
        son = getChild(e)
        parsimonyBottomUpSoftwired!(son, blobroot, nchar, w, parsimonyscore)
    end
    for s in 1:nchar
        for e in node.edge
            if (e.hybrid && e.fromBadDiamondI) || getChild(e) == node continue; end
            son = getChild(e)
            bestMin = Inf # best score from to this one child starting from s at the node
            for sf in 1:nchar # best assignement for the son
                minpars = parsimonyscore[son.number,sf]
                if s != sf
                    minpars +=1 # could be variable cost of change delta(sf,s)
                end
                bestMin = min(minpars, bestMin)
                end
            parsimonyscore[node.number,s] += bestMin  # add best assignement for the son to the PS of the parent
        end
    end
    #println("end of node $(node.number)")
    #@show parsimonyscore
    return nothing
end


"""
    parsimonySoftwired(net, tipdata)
    parsimonySoftwired(net, species, sequences)

Calculate the most parsimonious (MP) score of a network given
a discrete character at the tips using the dynamic programming algorithm given
in Fischer et al. (2015).
The softwired parsimony concept is used: where the number of state
transitions is minimized over all trees displayed in the network.

Data can given in one of the following:
- `tipdata`: data frame for a single trait, in which case the taxon names
are to appear in column 1 or in a column named "taxon" or "species", and
trait values are to appear in column 2 or in a column named "trait".
- `tipdata`: dictionary taxon => state, for a single trait.
- `species`: array of strings, and `sequences`: array of sequences,
   in the order corresponding to the order of species names.

# References

1. Fischer, M., van Iersel, L., Kelk, S., Scornavacca, C. (2015).
   On computing the Maximum Parsimony score of a phylogenetic network.
   SIAM J. Discrete Math., 29(1):559-585.
"""

function parsimonySoftwired(net::HybridNetwork, tips::Dict{String,T}) where {T}
    # T = type of characters. Typically Int if data are binary 0-1
    species = Array{String}(0)
    dat = Vector{Vector{T}}(0)
    for (k,v) in tips
        push!(species, k)
        push!(dat, [v])
    end
    parsimonySoftwired(net, species, dat)
end

function parsimonySoftwired(net::HybridNetwork, dat::DataFrame)
    i = findfirst(DataFrames.names(dat), :taxon)
    if i==0 i = findfirst(DataFrames.names(dat), :species); end
    if i==0 i=1; end # first column if not column named "taxon" or "species"
    j = findfirst(DataFrames.names(dat), :trait)
    if j==0 j=2; end
    if i==j
        error("""expecting taxon names in column 'taxon', or 'species' or column 1,
              and trait values in column 'trait' or column 2.""")
    end
    tips = Dict{String,eltypes(dat)[j]}() # fixit: give the species and trait arrays, not dictionary
    for r in 1:nrow(dat)
        if DataFrames.isna(dat[r,j]) continue; end
        tips[dat[r,i]] = dat[r,j]
    end
    parsimonySoftwired(net,tips)
end

function parsimonySoftwired(net::HybridNetwork, species=Array{String},
    sequenceData=AbstractArray)

    resetNodeNumbers!(net)
    nsites = length(sequenceData[1])
    nspecies = length(species) # no check to see if == length(sequenceData)
    nchar = 0 # this is to allocate the needed memory for w and
    # parsimonyscore matrices below. Only the first few columns will be
    # used for sites that have fewer than the maximum number of states.
    for i in 1:nsites
        nchar = max(nchar, length(unique([s[i] for s in sequenceData])))
    end
    w = zeros(Float64,(length(net.node), nchar))
    parsimonyscore = zeros(Float64,(length(net.node), nchar))
    tips = Dict{String, typeof(sequenceData[1][1])}() # to re-use memory later (?)
    checkGap = eltype(sequenceData) == BioSequences.BioSequence
    sequenceType = eltype(sequenceData[1])
    blobroots, majorEdges, minorEdges = blobInfo(net) # calls directEdges!: sets isChild1

    score = 0.0
    #allscores = Vector{Float64}(0)
    for isite in 1:nsites
      for sp in 1:nspecies
        nu = sequenceData[sp][isite] # nucleotide
        if checkGap &&
            (isgap(nu) ||
             (sequenceType==BioSymbols.DNA && nu == BioSymbols.DNA_N) ||
             (sequenceType==BioSymbols.RNA && nu == BioSymbols.RNA_N) ||
             (sequenceType==BioSymbols.AminoAcid && nu == BioSymbols.AA_N) )
            delete!(tips, species[sp])
        else
            tips[species[sp]] = nu # fixit: allow for ambiguous nucleotides?
        end
      end # tips now contains data for 1 site
      charset = union(values(tips))
      nchari = length(charset) # possibly < nchar
      nchari > 0 || continue
      initializeWeightsFromLeavesSoftwired!(w, net, tips, charset) # data are now in w
      fill!(parsimonyscore, 0.0) # reinitialize parsimonyscore too
      # @show tips; @show w

      for bcnumber in 1:length(blobroots)
        r = blobroots[bcnumber]
        nhyb = length(majorEdges[bcnumber])
        # println("site $isite, r.number = $(r.number), nhyb=$(nhyb)")
        mpscoreSwitchings = Array{Float64}(2^nhyb, nchari) # grabs memory
        # fixit: move that outside to avoid grabbing memory over and over again
        iswitch = 0
        for switching in IterTools.product([[true, false] for i=1:nhyb]...)
            # next: modify the `fromBadDiamondI` of hybrid edges in the blob:
            # switching[h] = pick the major parent of hybrid h if true, pick minor if false
            iswitch += 1
            #@show switching
            for h in 1:nhyb
                majorEdges[bcnumber][h].fromBadDiamondI =  switching[h]
                minorEdges[bcnumber][h].fromBadDiamondI = !switching[h]
            end
            parsimonyBottomUpSoftwired!(r, r, nchari, w, parsimonyscore) # updates parsimonyscore
            for s in 1:nchari
              mpscoreSwitchings[iswitch,s] = parsimonyscore[r.number,s]
            end
        end

        for i in 1:nchari
            # add best assignement for the son to the PS of the parent
            w[r.number,i] += minimum(mpscoreSwitchings[:,i])
        end
        bcnumber+=1
      end

      bestMin = minimum(w[net.node[net.root].number, 1:nchari])
      #if bestMin>0 println("site $isite, score $bestMin"); end
      #push!(allscores, bestMin)
      score += bestMin
    end
    return score #, allscores
end

"""
    initializeWeightsFromLeavesSoftwired!(w, net, tips, charset)

Modify weight in w: to Inf for w[n, s] if the "tips" data
has a state different from s at node number n.
Assumes that w was initialized to 0 for the leaves.
"""
function initializeWeightsFromLeavesSoftwired!(w::AbstractArray, net::HybridNetwork, tips, charset)
    fill!(w, 0.0)
    for node in net.node
        node.leaf || continue
        haskey(tips, node.name) || continue
        for i in 1:length(charset)
            s = charset[i]
            if s != tips[node.name]
                w[node.number,i] = Inf
            end
        end
        #@show w[node.number,:]
    end
end

function readFastaToSequenceDict(filename::String)
    reader = BioSequences.FASTA.Reader(open(filename))
    #dat = Dict{String, }()
    sequences = Array{BioSequences.BioSequence}(0)
    species = Array{String}(0)
    for record in reader
        #push!(dat, FASTA.identifier(record) => sequence(record))
        push!(sequences, sequence(record))
        push!(species, FASTA.identifier(record))
    end
    seqlengths = [length(s) for s in sequences] # values(dat)
    nsites, tmp = extrema(seqlengths)
    nsites == tmp || error("sequences not of same lengths: from $nsites to $tmp sites")
    """
    fixit:
    - deleteat!(seq::BioSequence, i::Integer) to delete all but 1 site with identical patterns
    - calculate and output the weights
    - also delete invariant sites
    - clever: change AACC to 0011 because score AACC=AATT=CCTT etc. Sort, keep uniques.
    """
    return species, sequences
end

"""
parsimonyGF(net::HybridNetwork, species, sequenceData)

GF = general framework
"""

function parsimonyGF(net::HybridNetwork, tips::Dict{String,T}) where {T}
    # T = type of characters. Typically Int if data are binary 0-1
    species = Array{String}(0)
    dat = Vector{Vector{T}}(0)
    for (k,v) in tips
        push!(species, k)
        push!(dat, [v])
    end
    parsimonyGF(net, species, dat)
end

function parsimonyGF(net::HybridNetwork, species=Array{String},
    sequenceData=AbstractArray, criterion=:softwired::Symbol)

    resetNodeNumbers!(net)
    rootnumber = net.node[net.root].number
    nsites = length(sequenceData[1])
    nspecies = length(species) # no check to see if == length(sequenceData)
    nstates = 0 # this is to allocate the needed memory for w and
    # parsimonyscore matrices below. Only the first few columns will be
    # used for sites that have fewer than the maximum number of states.
    for i in 1:nsites
        nstates = max(nstates, length(unique([s[i] for s in sequenceData])))
    end
    if criterion == :softwired # lineage states: ∅, 1, 2, 3, ..., nstates
        lineagestate = [Set{Int}()] # first element: empty set ∅
        for i in 1:nstates
            push!(lineagestate, Set(i))
        end
        nchar = length(lineagestate)
        allowedAtRoot = 2:nchar # ∅ is not allowed at root of network
        costfunction = function(finalset, parentsets) # parentsets = array of sets: 0, 1, 2
            if length(finalset) > sum([length(s) for s in parentsets])
                return Inf
            else # calculate final set minus (union of parent sets)
                missingstates = finalset
                for s in parentsets
                    missingstates = setdiff(missingstates, s)
                end
                return length(missingstates)
            end
        end
        costmatrix1 = Array{Float64}(nchar,nchar) # i = 1 parent set, j = final set
        for i in 1:nchar
            for j in 1:nchar
                costmatrix1[i,j] = costfunction(lineagestate[j], [lineagestate[i]])
            end
        end
        costmatrix2 = Vector{Array{Float64}}(nchar)
        # costmatrix2[k][i,j]: from parents k,i to final set j
        for k in 1:nchar
            costmatrix2[k] = Array{Float64}(nchar,nchar)
            for i in 1:nchar
                for j in 1:nchar
                    costmatrix2[k][i,j] = costfunction(lineagestate[j], [lineagestate[i], lineagestate[k]])
                end
            end
        end
    elseif criterion == :parental
        error("parental parsimony not implemented yet")
    else
        error("criterion $criterion unknown")
    end

    w = zeros(Float64,(length(net.node), nchar))
    parsimonyscore = zeros(Float64,(length(net.node), nchar))
    tips = Dict{String, typeof(sequenceData[1][1])}() # to re-use memory later (?)
    checkGap = eltype(sequenceData) == BioSequences.BioSequence
    sequenceType = eltype(sequenceData[1])
    blobroots, majorEdges, minorEdges = blobInfo(net) # calls directEdges!: sets isChild1
    # fixit: use trivial biconnected components, and compare running time
        # pick 1 parent node (the minor parent arbitrarily) for each hybrid, then
    # "cut" both children edges of that parent: mark its `fromBadDiamondI` = false
    for e in net.edge
        e.fromBadDiamondI = false # don't cut by default
    end
    guessedparent = Vector{Vector{Node}}(0)
    for melist in minorEdges # minor edge list for one single blob + loop over blobs
        guessedparentBlob = Vector{Node}(0)
        for e in melist
            p = getParent(e)
            if p ∉ guessedparentBlob
                push!(guessedparentBlob, p)
            end
            for e2 in p.edge
                p == getParent(e2) || continue
                e2.fromBadDiamondI = true # cut the edge: not followed in recursive call
            end
        end
        push!(guessedparent, guessedparentBlob)
    end

    maxguessedParents = maximum([length(me) for me in guessedparent])
    guessStates = Array{Float64}(nchar^maxguessedParents, nchar) # grabs memory

    #@show guessedparent

    score = 0.0
    #allscores = Vector{Float64}(0)
    for isite in 1:nsites
      for sp in 1:nspecies
        nu = sequenceData[sp][isite] # nucleotide
        if checkGap &&
            (isgap(nu) ||
             (sequenceType==BioSymbols.DNA && nu == BioSymbols.DNA_N) ||
             (sequenceType==BioSymbols.RNA && nu == BioSymbols.RNA_N) ||
             (sequenceType==BioSymbols.AminoAcid && nu == BioSymbols.AA_N) )
            delete!(tips, species[sp])
        else
            tips[species[sp]] = nu # fixit: allow for ambiguous nucleotides?
        end
      end # tips now contains data for 1 site
      stateset = union(values(tips))
      nstatesi = length(stateset) # possibly < nstates
      nstatesi > 0 || continue
      if criterion == :softwired
        nchari = nstatesi + 1
      end
      initializeWeightsFromLeaves!(w, net, tips, stateset, criterion) # data are now in w
      # @show tips; @show w

      for bcnumber in 1:length(blobroots)
        r = blobroots[bcnumber]
        nhyb = length(majorEdges[bcnumber])
        # println("site $isite, r.number = $(r.number), nhyb=$(nhyb)")
        iguess = 0
        for guesses in IterTools.product([1:nchari for i=1:length(guessedparent[bcnumber])]...)
            iguess += 1
            fill!(parsimonyscore, 0.0) # reinitialize parsimonyscore too
            #@show guesses
            for pind in 1:nhyb
                p = guessedparent[bcnumber][pind] # detached parent of hybrid with index pind in blob number pcnumber
                for s in 1:nchari
                    w[p.number, s] = Inf
                end
                w[p.number, guesses[pind]] = 0.0
            end
            parsimonyBottomUpGF!(r, r, nchari, w, parsimonyscore, costmatrix1, costmatrix2)
            # recursion above: updates parsimonyscore
            #@show nchari
            for s in 1:nchari
                guessStates[iguess,s] = parsimonyscore[r.number,s]
            end
            #@show guessStates

        end
        # now iguess = total number of guess assignments for that blob
        for i in 1:nchari
            # add best assignement for the son to the PS of the parent
            w[r.number,i] += minimum(guessStates[1:iguess,i])
        end
        bcnumber+=1
      end

      bestMin = Inf
      for s in allowedAtRoot
        s <= nchari || break
        bestMin = min(bestMin, w[rootnumber, s])
      end
      #if bestMin>0 println("site $isite, score $bestMin"); end
      #push!(allscores, bestMin)
      score += bestMin
    end
    return score #, allscores
end

"""
    initializeWeightsFromLeaves!(w, net, tips, stateset, criterion)

Modify weight in w: to Inf for w[n, i] if the "tips" data
has a state different from the lineage state of index i at node number n.
Assumes that w was initialized to 0 for the leaves.

criterion: should be one of `:softwired`, `:parental` or `:hardwired`.
- softwired parsimony: lineage states are in this order: ∅,{1},{2},{3},...,{nstates}
"""
function initializeWeightsFromLeaves!(w::AbstractArray, net::HybridNetwork, tips, stateset,
        criterion::Symbol)
    fill!(w, 0.0)
    if criterion == :softwired
        nchar = length(stateset) + 1
    elseif criterion == :parental
        error("parental parsimony not implemented yet")
    else
        error("unknown criterion: $criterion")
    end
    for node in net.node
        node.leaf || continue
        w[node.number,1] = Inf # ∅ comes first
        haskey(tips, node.name) || continue
        for i in 2:nchar
            s = stateset[i-1] # state s correspond to lineage at index s+1
            if s != tips[node.name]
                w[node.number,i] = Inf
            end
        end
        #@show w[node.number,:]
    end
end

"""
    `parsimonyBottomUpGF!(node, blobroot, states, w, scores,
        costmatrix1, costmatrix2)`

Computing the MP scores (one for each assignment of the blob root state)
given the descendants of a blob, conditional on the states at predefined parents
of hybrids in the blobs (one parent per hybrid) as described in

Leo Van Iersel, Mark Jones, Celine Scornavacca (2017).
Improved Maximum Parsimony Models for Phylogenetic Networks,
Systematic Biology,
(https://doi.org/10.1093/sysbio/syx094).

Assumes a set of state *guesses*, ie correct initialization of `w` for
predefined hybrid parents, and correct `fromBadDiamondI` field for the children
edges of these predefined parents. For edges that are cut `fromBadDiamondI` is true.

The field `isExtBadTriangle` is used to know which nodes are at the root of a blob.
The field isChild1 is used (and assumed correct).

fixit: explain last arguments
"""

function parsimonyBottomUpGF!(node::Node, blobroot::Node, nchar::Integer,
    w::AbstractArray, parsimonyscore::AbstractArray,
    costmatrix1::AbstractArray, costmatrix2::AbstractArray)

    println("entering with node $(node.number)")
    parsimonyscore[node.number,1:nchar] = w[node.number,1:nchar] # at all nodes to re-initialize between guesses
    if !node.leaf && (!node.isExtBadTriangle || node == blobroot)
    # isExtBadTriangle=dummy leaf: root of another blob
    for e in node.edge # post-order traversal according to major tree: detached edges were minor.
        if !e.isMajor || getChild(e) == node continue; end # Even if we didn't visit one parent (yet),
        son = getChild(e) # that parent is a minor parent with an assigned guessed state.
        parsimonyBottomUpGF!(son, blobroot, nchar, w, parsimonyscore, costmatrix1, costmatrix2)
    end
    # check to see if "node" has guessed value: by checking to see if all its children edges were cut
    cutparent = false
    for e in node.edge
        if getChild(e) == node continue; end
        cutparent = cutparent || e.fromBadDiamondI # if fromBadDiamondI=true: edge is cut
    end
    if !cutparent # look at best assignment of children, to score each assignment at node
        for e in node.edge # avoid edges that were cut: those for which fromBadDiamondI is true
            if getChild(e) == node continue; end
            son = getChild(e)
            bestpars = [Inf for s in 1:nchar] # best score, so far, for state s at node.
            # calculate cost to go from each s (and from parents' guesses, stored in w) to son state:
            if son.hybrid
                # find potential other parent of son, detached with guessed states
                p2 = getMinorParent(son)
                k = findfirst(w[p2.number, 1:nchar], 0.0) # guess made for parent p2
                for sfinal in 1:nchar
                    pars = parsimonyscore[son.number, sfinal] + costmatrix2[k][1:nchar,sfinal]
                    for s in 1:nchar
                        bestpars[s] = min(bestpars[s], pars[s])
                    end
                end
            else # son has no guessed parent: has only 1 parent: "node"
                for sfinal in 1:nchar
                    pars = parsimonyscore[son.number, sfinal] + costmatrix1[1:nchar,sfinal]
                    for s in 1:nchar
                        bestpars[s] = min(bestpars[s], pars[s])
                    end
                end
            end
            parsimonyscore[node.number,1:nchar] += bestpars # add score from best assignement for this son
        end
    end
    end # of if: not leaf, not root of another blob
    println("almost the end of node $(node.number)")
    #@show parsimonyscore
    #@show w
    node != blobroot ||  return nothing
    # if blob root has detached parents only, these are paid for in the blob
    # in which this blob root is a leaf.
    # pay now for re-attaching the guessed nodes to their children
    detachedParents = Vector{Node}(0)
    for e in node.edge
        if getParent(e) == node continue; end
        e.fromBadDiamondI || return nothing # break if one parent edge was not cut
        push!(detachedParents, getParent(e))
    end
    # if we get here, it means that "node" is the root of the blob
    # or the root of a tree after detaching the guessed parents from all their children.
    #@show detachedParents

    if length(detachedParents) == 0 return nothing # nothing to add

    elseif length(detachedParents) == 1
        k = findfirst(w[detachedParents[1].number, 1:nchar], 0.0) # guess made for parent p1
        pars = parsimonyscore[node.number, 1:nchar] + costmatrix1[k,1:nchar]
    elseif length(detachedParents) == 2
        k1 = findfirst(w[detachedParents[1].number, 1:nchar], 0.0)
        k2 = findfirst(w[detachedParents[2].number, 1:nchar], 0.0)
        pars = parsimonyscore[node.number, 1:nchar] + costmatrix2[k1][k2,1:nchar]
    end
    cost = minimum(pars)
    for s in 1:nchar
        parsimonyscore[blobroot.number, s] += cost
    end

    println("the end of node $(node.number)")
    #@show parsimonyscore
    #@show w

    println("end of node $(node.number)")
    return nothing
end
