
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
    filter!(p -> p.second==mv, votes) # keep states with max votes
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
    mpt = findall(x -> x==mps, mpscore) # indices of all trees with MP score
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
    i = findfirst(isequal(:taxon), DataFrames.names(dat))
    if i===nothing i = findfirst(isequal(:species), DataFrames.names(dat)); end
    if i===nothing i=1; end # first column if no column named "taxon" or "species"
    j = findfirst(isequal(:trait), DataFrames.names(dat))
    if j===nothing j=2; end
    if i==j
        error("""expecting taxon names in column 'taxon', or 'species' or column 1,
              and trait values in column 'trait' or column 2.""")
    end
    tips = Dict{String,eltypes(dat)[j]}()
    for r in 1:nrow(dat)
        if ismissing(dat[r,j]) continue; end
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
a discrete character at the tips.
The softwired parsimony concept is used: where the number of state
transitions is minimized over all trees displayed in the network.

Data can given in one of the following:
- `tipdata`: data frame for a single trait, in which case the taxon names
   are to appear in column 1 or in a column named "taxon" or "species", and
   trait values are to appear in column 2 or in a column named "trait".
- `tipdata`: dictionary taxon => state, for a single trait.
- `species`: array of strings, and `sequences`: array of sequences,
   in the order corresponding to the order of species names.

## algorithm

The dynamic programming algorithm by Fischer et al. (2015) is used.
The function loops over all the displayed subtrees within a blob
(biconnected component), so its complexity is of the order of
`n * m * c^2 * 2^level` where `n` is the number of tips,
`m` the number of traits, `c` the number of states, and `level`
is the level of the network: the maximum number of hybridizations
within a blob.

See [`parsimonyGF`](@ref) for a different algorithm, slower but
extendable to other parsimony criteria.

## references

1. Fischer, M., van Iersel, L., Kelk, S., Scornavacca, C. (2015).
   On computing the Maximum Parsimony score of a phylogenetic network.
   SIAM J. Discrete Math., 29(1):559-585.
"""
function parsimonySoftwired(net::HybridNetwork, tips::Dict{String,T}) where {T}
    # T = type of characters. Typically Int if data are binary 0-1
    species = String[]
    dat = Vector{T}[]
    for (k,v) in tips
        push!(species, k)
        push!(dat, [v])
    end
    parsimonySoftwired(net, species, dat)
end

function parsimonySoftwired(net::HybridNetwork, dat::DataFrame)
    i = findfirst(isequal(:taxon), DataFrames.names(dat))
    if i===nothing i = findfirst(isequal(:species), DataFrames.names(dat)); end
    if i===nothing i=1; end # first column if no column named "taxon" or "species"
    j = findfirst(isequal(:trait), DataFrames.names(dat))
    if j===nothing j=2; end
    if i==j
        error("""expecting taxon names in column 'taxon', or 'species' or column 1,
              and trait values in column 'trait' or column 2.""")
    end
    innet = findall(in(tipLabels(net)), dat[i]) # species in the network
    species = dat[innet, i]
    tips = dat[j][innet]
    indna = findall(ismissing, tips) # species with missing data
    deleteat!(species, indna)
    deleteat!(tips,    indna)
    parsimonySoftwired(net,tips)
end

function parsimonySoftwired(net::HybridNetwork, species::Array{String},
    sequenceData::AbstractArray)

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
    #allscores = Float64[]
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
        mpscoreSwitchings = Array{Float64}(undef, 2^nhyb, nchari) # grabs memory
        # fixit: move that outside to avoid grabbing memory over and over again
        iswitch = 0
        perms = nhyb == 0 ? [()] : Iterators.product([[true, false] for i=1:nhyb]...)
        for switching in perms
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

function readFastaToArray(filename::String)
    reader = BioSequences.FASTA.Reader(open(filename))
    #dat = Dict{String, }()
    sequences = Array{BioSequences.BioSequence}(undef, 0)
    species = String[]
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
    readfastatodna(filename::String)

Read a fasta file to a dataframe containing a column for each site.
Calculate weights and remove matching site patterns to reduce matrix dimension.

Return a tuple containing:
1. data frame of BioSequence DNA sequences, with taxon names in column 1
   followed by a column for each site pattern, in columns 2-npatterns;
2. array of weights, one weight for each of the site columns.
   The length of the weight vector is equal to npatterns.
"""
function readfastatodna(fastafile::String, countPatterns=false::Bool)
    reader = BioSequences.FASTA.Reader(open(fastafile))
    siteList = Vector{Vector}(undef, 0) #array of arrays, one array for each site (8 in example)
    species = String[]
    firstspecies = Bool(true)
    nsites = 0
    for record in reader #by species (row)
        if firstspecies
            nsites = length(sequence(record))
            for site in 1:nsites # initialize an array for each site
                push!(siteList, Vector{BioSequences.DNA}(undef, 0))
            end
            firstspecies = false
        end
        push!(species, FASTA.identifier(record)) #adds species name to end of species vector
        length(sequence(record)) == nsites || error("sequences of different length: current sequences is ", 
            length(sequence(record)), " long while first sequence is ", nsites, " long")
        for site in 1:nsites
            push!(siteList[site], sequence(record)[site])
        end
    end

    # calculate weights and remove matching site patterns to reduce dimension of siteList
    weights = ones(Float64, nsites)
    if countPatterns
        for c in nsites:-1:1
            target = siteList[c]
            for comparison in 1:(c-1)
                if target == siteList[comparison] 
                    weights[comparison] += weights[c] #add c's weight to comparison's weight
                    deleteat!(siteList, c) # delete c in siteList
                    deleteat!(weights, c) #delete c in weights
                    break
                end
            end
        end
    end

    #create dat here
    dat = DataFrame(siteList)
    insertcols!(dat, 1, taxon = species)
    return (dat, weights)
end

"""
    readCSVtoArray(dat::DataFrame)
    readCSVtoArray(filename::String)

Read a CSV table containing both species names and data,
create two separate arrays: one for the species names,
a second for the data, in a format that [`parsimonyGF`](@ref) needs.

Warning:
- it will try to find a column 'taxon' or 'species' for the taxon names.
  If none found, it will assume the taxon names are in column 1.
- will use all other columns as characters
"""
function readCSVtoArray(dat::DataFrame)
    i = findfirst(isequal(:taxon), DataFrames.names(dat))
    if i===nothing i = findfirst(isequal(:species), DataFrames.names(dat)); end
    if i===nothing
        @warn "expecting taxon names in column 'taxon', or 'species', so will assume column 1"
        i = 1
    end

    species = String[]
    for d in dat[!,i]
        push!(species,string(d))
    end

    ind = deleteat!(collect(1:size(dat,2)),i) ##character columns
    seq = Vector{Vector{Any}}(undef, 0)
    for j=1:size(dat,1) ##for every taxon:
        v = Vector{Any}(undef, 0)
        for ii in ind
            if ismissing(dat[j,ii])
                push!(v,missing)
            else
                push!(v,dat[j,ii])
            end
        end
        push!(seq,v)
    end
    return species,seq
end

function readCSVtoArray(filename::String)
    dat = CSV.read(filename)
    readCSVtoArray(dat)
end

"""
    parsimonyGF(net, tip_dictionary, criterion=:softwired)
    parsimonyGF(net, species, sequenceData, criterion=:softwired)

Calculate the most parsimonious score of a network given
discrete characters at the tips using a general framework
(Van Iersel et al. 2018) allowing for various parsimony criteria:
softwired (default), hardwired, parental etc.
Only softwired is implemented at the moment.

Data can given in one of the following:
- `tipdata`: data frame for a single trait, in which case the taxon names
   are to appear in column 1 or in a column named "taxon" or "species", and
   trait values are to appear in column 2 or in a column named "trait".
- `tipdata`: dictionary taxon => state, for a single trait.
- `species`: array of strings, and `sequences`: array of sequences,
    in the order corresponding to the order of species names.

## algorithm

The complexity of the algorithm is exponential in the level of the
network, that is, the maximum number of hybridizations in a single
blob, or biconnected component (Fischer et al. 2015).
The function loops over all the state assignments of the minor parent
of each hybrid node within a blob, so its complexity is of the order of
`n * m * c^2 * c^level` where `n` is the number of tips,
`m` the number of traits and `c` the number of states.

See [`parsimonySoftwired`](@ref) for a faster algorithm, but
solving the softwired criterion only.

## references

1. Leo Van Iersel, Mark Jones, Celine Scornavacca (2017).
   Improved Maximum Parsimony Models for Phylogenetic Networks,
   Systematic Biology,
   (https://doi.org/10.1093/sysbio/syx094).

2. Fischer, M., van Iersel, L., Kelk, S., Scornavacca, C. (2015).
   On computing the Maximum Parsimony score of a phylogenetic network.
   SIAM J. Discrete Math., 29(1):559-585.

Use the recursive helper function [`parsimonyBottomUpGF!`](@ref).
Use the fields `isChild1`,
`isExtBadTriangle` to know which nodes are at the root of a blob, and
`fromBadDiamondI` to know which edges are cut (below the minor parent of each hybrid).
"""
function parsimonyGF(net::HybridNetwork, tips::Dict{String,T},
                     criterion=:softwired::Symbol) where {T}
    # T = type of characters. Typically Int if data are binary 0-1
    species = String[]
    dat = Vector{T}[]
    for (k,v) in tips
        push!(species, k)
        push!(dat, [v])
    end
    parsimonyGF(net, species, dat, criterion)
end


function parsimonyGF(net::HybridNetwork, species::Array{String},
    sequenceData::AbstractArray, criterion=:softwired::Symbol)

    resetNodeNumbers!(net) # direct edges and checks pre-order by default
    rootnumber = net.node[net.root].number
    nsites = length(sequenceData[1])
    nspecies = length(species) # no check to see if == length(sequenceData)
    nstates = 0 # this is to allocate the needed memory for w and
    # parsimonyscore matrices below. Only the first few columns will be
    # used for sites that have fewer than the maximum number of states.
    for i in 1:nsites
        nstatesi = length(unique([s[i] for s in sequenceData]))
        if nstates < nstatesi
           nstates = nstatesi
        end
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
        costmatrix1 = Array{Float64}(undef, nchar,nchar) # i = 1 parent set, j = final set
        for i in 1:nchar
            for j in 1:nchar
                costmatrix1[i,j] = costfunction(lineagestate[j], [lineagestate[i]])
            end
        end
        costmatrix2 = Vector{Array{Float64}}(undef, nchar)
        # costmatrix2[k][i,j]: from parents k,i to final set j
        for k in 1:nchar
            costmatrix2[k] = Array{Float64}(undef, nchar,nchar)
            for i in 1:nchar
                for j in 1:nchar
                    costmatrix2[k][i,j] = costfunction(lineagestate[j], [lineagestate[i], lineagestate[k]])
                end
            end
        end
    elseif criterion == :parental
        error("parental parsimony not implemented yet")
        # costfunctionP = function(finalset, parentsets) # parentsets = array of sets: 0, 1, 2
        #     # Input: node v in network N with at most one parent u not in P, set S in Y, set Sprime 
        #     #   in Y, assignment fprime.
        #     # Return: cost on eduges entering v for any lineage function f that extends fprime and 
        #     #   assigns f(u) = S (if u exists) and f(v) = Sprime
        #     if node v is root of Network N
        #         return 0
        #     else
        #         if length(finalset) > sum([length(s) for s in parentsets])
        #         return Inf
        #         else # calculate final set minus (union of parent sets)
        #             missingstates = finalset
        #         for s in parentsets
        #             missingstates = setdiff(missingstates, s)
        #         end
        #         return length(missingstates)
        #     end
        # end
        # for fprime in set of fprimes
        #     for i in 1:r
        #         for each vertex u (in reverse topology ordering) and s in Y
        #             if u is a leaf in X
        #                 H[u,S] = 0*(S==alpha(u)) + inf(S!=alpha(u))
        #             end
        #             if u in P
        #                 H[u,S] = 0*(S==fprime(u)) + inf(S!=fprime(u))
        #             end
        #             if u has one child v in T
        #                 H[u,S] = min(H(v, Sprime)) + costfunctionP(v, s, sprime, fprime)
        #             end
        #             if u has two children v1,v2 in Ti
        #                 H[u,S] = min(H(v1, Sprime)) + costfunctionP(v1, s, sprime, fprime) +
        #                         min(H(v2, Sprime)) + costfunctionP(v2, s, sprime, fprime)
        #             end
        #         end
        #     end
        # end
        # optfprime = min(H[rho1, (j)]) + min(costfunctionP(rho1, null, S, fprime) + H(rho[i], S))
        # #note: its not clear where we optimize fprime ^
        # return optfprime
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
    guessedparent = Vector{Node}[]
    for melist in minorEdges # minor edge list for one single blob + loop over blobs
        guessedparentBlob = Node[]
        for e in melist
            p = getParent(e)
            if p ∉ guessedparentBlob
              push!(guessedparentBlob, p)
              for e2 in p.edge
                p == getParent(e2) || continue
                e2.fromBadDiamondI = true # cut the edge: not followed in recursive call
              end
            end
        end
        push!(guessedparent, guessedparentBlob)
    end
    maxguessedParents = maximum([length(me) for me in guessedparent])
    guessStates = Vector{Float64}(undef, nchar) # grabs memory
    # will contain the best score for a blob root starting at some state s,
    # from best guess so far (different guesses are okay for different s)

    # determine which nodes are the root of trees after we cut off edges
    # below guessed parents. For this, at each node, set its
    # inCycle = 0 if the node has (at least) 1 non-cut edge, and
    # inCycle = # of detached parents if the node has detached parents only.
    # Do this now to avoid re-doing it at each site.
    # `inCycle` was used earlier to find blobs.
    for n in net.node n.inCycle=0; end
    # first: calculate inCycle = # of detached parents
    for e in net.edge
        e.fromBadDiamondI || continue # to next edge if current edge not cut
        n = getChild(e)
        n.inCycle += 1
    end
    # second: make inCycle = 0 if node has a non-detached parent
    for n in net.node
        n.inCycle > 0 || continue # to next node if inCycle = 0 already
        for e in n.edge
            n == getChild(e) || continue
            !e.fromBadDiamondI || continue
            n.inCycle = 0 # if e parent of n and e not cut: make inCycle 0
            break
        end
    end

    score = 0.0
    #allscores = Float64[]
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
      #fill!(parsimonyscore, 0.0) # reinitialize parsimonyscore too
      # @show tips; @show w

      for bcnumber in 1:length(blobroots)
        r = blobroots[bcnumber]
        nhyb = length(majorEdges[bcnumber])
        #println("site $isite, r.number = $(r.number), nhyb=$(nhyb)")
        firstguess = true
        gplen = length(guessedparent[bcnumber])
        perms = gplen == 0 ? [()] : Iterators.product([1:nchari for i=1:gplen]...)
        for guesses in perms
            #@show guesses
            for pind in 1:nhyb
                p = guessedparent[bcnumber][pind] # detached parent of hybrid with index pind in blob number pcnumber
                # fixit: memory greedy below. check what was changed only.
                for s in 1:nchari
                    w[p.number, s] = Inf
                end
                w[p.number, guesses[pind]] = 0.0
            end
            parsimonyBottomUpGF!(r, r, nchari, w, parsimonyscore, costmatrix1, costmatrix2)
            # recursion above: updates parsimonyscore
            #@show nchari
            if firstguess
              for s in 1:nchari
                guessStates[s] = parsimonyscore[r.number,s]
              end
              firstguess = false
            else
              for s in 1:nchari
                if guessStates[s] > parsimonyscore[r.number,s]
                   guessStates[s] = parsimonyscore[r.number,s]
                end
              end
            end
            #@show guessStates
        end
        for i in 1:nchari
            # add best assignement for the son to the PS of the parent
            w[r.number,i] += guessStates[i]
        end
        bcnumber+=1
      end

      bestMin = Inf
      for s in allowedAtRoot
        s <= nchari || break
        if bestMin > w[rootnumber, s]
           bestMin = w[rootnumber, s]
        end
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
    `parsimonyBottomUpGF!(node, blobroot, nchar, w, scores,
        costmatrix1, costmatrix2)`

Compute the MP scores (one for each assignment of the blob root state)
given the descendants of a blob, conditional on the states at predefined parents
of hybrids in the blobs (one parent per hybrid) as described in

Leo Van Iersel, Mark Jones, Celine Scornavacca (2017).
Improved Maximum Parsimony Models for Phylogenetic Networks,
Systematic Biology,
(https://doi.org/10.1093/sysbio/syx094).

Assumes a set of state *guesses*, ie correct initialization of `w` for
predefined hybrid parents, and correct `fromBadDiamondI` field for the children
edges of these predefined parents. `fromBadDiamondI` is true for edges that are cut.

The field `isExtBadTriangle` is used to know which nodes are at the root of a blob.
The field `isChild1` is used (and assumed correct).
Field `inCycle` is assumed to store the # of detached parents (with guessed states)

- `nchar`: number of characters considered at internal lineages.
  For softwired parsimony, this is # states + 1, because characters
  at internal nodes are ∅, {1}, {2}, etc.
  For parental parsimony, this is 2^#states -1, because characters
  are all sets on {1,2,...} except for the empty set ∅.
- `costmatrix1`[i,j] and `costmatrix2`[k][i,j]: 2d array and vector of 2d arrays
  containing the cost of going to character j starting from character i
  when the end node has a single parent, or the cost of a child node
  having character j when its parents have characters k and i.
  These cost matrices are pre-computed depending on the parsimony criterion
  (softwired, hardwired, parental etc.)

used by [`parsimonyGF`](@ref).
"""
function parsimonyBottomUpGF!(node::Node, blobroot::Node, nchar::Integer,
    w::AbstractArray, parsimonyscore::AbstractArray,
    costmatrix1::AbstractArray, costmatrix2::AbstractArray)

    #println("entering with node $(node.number)")
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
        if e.fromBadDiamondI # if true: one child edge is cut, so all are cut
            cutparent = true
            break
        end
    end
    if !cutparent # look at best assignment of children, to score each assignment at node
        for e in node.edge # avoid edges that were cut: those for which fromBadDiamondI is true
            son = getChild(e)
            if son == node continue; end
            bestpars = [Inf for s in 1:nchar] # best score, so far, for state s at node.
            # calculate cost to go from each s (and from parents' guesses, stored in w) to son state:
            if son.hybrid
                # find potential other parent of son, detached with guessed states
                p2 = getMinorParent(son)
                k = findfirst(isequal(0.0), w[p2.number, 1:nchar]) # guess made for parent p2
                for sfinal in 1:nchar
                    pars = parsimonyscore[son.number, sfinal] .+ costmatrix2[k][1:nchar,sfinal]
                    for s in 1:nchar
                        if bestpars[s] > pars[s]
                           bestpars[s] = pars[s]
                        end
                    end
                end
            else # son has no guessed parent: has only 1 parent: "node"
                for sfinal in 1:nchar
                    pars = parsimonyscore[son.number, sfinal] .+ costmatrix1[1:nchar,sfinal]
                    for s in 1:nchar
                        if bestpars[s] > pars[s]
                           bestpars[s] = pars[s]
                        end
                    end
                end
            end
            parsimonyscore[node.number,1:nchar] .+= bestpars # add score from best assignement for this son
        end
    end
    end # of if: not leaf, not root of another blob
    #println("almost the end of node $(node.number)")
    #@show parsimonyscore
    #@show w
    node != blobroot ||  return nothing
    # if blob root has detached parents only, these are paid for in the blob
    # in which this blob root is a leaf.
    node.inCycle > 0 || return nothing # inCycle = # of detached parents
    # if we get here, it means that "node" is not root of current blob, but
    # is root of a tree after detaching all guessed parents from their children.
    # pay now for re-attaching the guessed parents to node
    cost = 0.0 # variable external to 'for' loops below
    if node.inCycle == 1 # 1 detached parent, no non-detached parents
        for e in node.edge
            par = getParent(e)
            if par == node continue; end
            # now 'par' is the single detached parent
            k = findfirst(isequal(0.0), w[par.number, 1:nchar]) # guess at parent
            cost = minimum(parsimonyscore[node.number, 1:nchar] +
                            costmatrix1[k,1:nchar])
            #println("node $(node.number), parent $(par.number), guess k=$k")
            break # out of for loop
        end
    else # node.inCycle should be 2: 2 detached parents
        k1 = 0 # guess made on first detached parent
        for e in node.edge
            par = getParent(e)
            if par == node continue; end
            # now 'par' is one of the 2 guessed parents
            if k1 == 0
                k1 = findfirst(isequal(0.0), w[par.number, 1:nchar]) # guess at 1st parent
            else
                k2 = findfirst(isequal(0.0), w[par.number, 1:nchar]) # guess at 2nd parent
                cost = minimum(parsimonyscore[node.number, 1:nchar] +
                               costmatrix2[k1][k2,1:nchar])
                break
            end
        end
    end # now 'cost' is up-to-date
    for s in 1:nchar
        parsimonyscore[blobroot.number, s] += cost
    end
    #@show parsimonyscore
    #@show w
    #println("end of node $(node.number)")
    return nothing
end


## --------------------------------------------------
## Search for most parsimonious network

## start the search from (or near) topology currT,
## .loglik will save now parsimony score
## fixit: we are using the functions for level1 network, we need to include level-k networks
## this will be done by changing the current proposedTop! function
## the current search/moves functions are the same for snaq, so they need a "good" semi-directed
## network, so currently we propose new topology from this set of networks, and simply root to
## compute the parsimony score (we do not keep rooted networks as the search objects),
## that is, newT and currT are unrooted through the whole algorithm (because computing the parsimony destroys inCycle)
## criterion=softwired by default
## at this stage, currT is compatible with the outgroup, but still unrooted (for moves functions to work)
## tolAbs: could be set to 0.1 and removed from list of arguments. up to 0.5 should work,
##         because parsimony scores are integers, but Float64 to extend cost function later perhaps
"""
Road map for various functions behind maxParsimonyNet

    maxParsimonyNet
    maxParsimonyNetRun1
    maxParsimonyNetRun1!

All return their optimized network. Only maxParsimonyNet returns a rooted network
(though all functions guarantee that the returned networks agree with the outgroup).

- maxParsimonyNet calls maxParsimonyNetRun1 per run, after a read(write(.)) of the starting network
  (to ensure level-1 and semi-directedness).
- maxParsimonyNetRun1 will make a copy of the topology, and will call findStartingTopology!
  to modify the topology according to random NNI/move origin/move target moves. It then calls maxParsimonyNetRun1!
  on the modified network
- maxParsimonyNetRun1! proposes new network with various moves (same moves as snaq), and stops when it finds the
  most parsimonious network, using [`parsimonyGF`](@ref).

None of these functions allow for multiple alleles yet.

Note that the search algorithm keeps two HybridNetworks at a time: currT (current topology) and newT (proposed topology).
Both are kept unrooted (semi-directed), otherwise the moves in proposedTop! function fail.
We only root the topologies to calculate the parsimony, so we create a rooted copy (currTr, newTr) to compute parsimony
score in this copied topology. We do not root and calculate parsimony score in the original HybridNetworks objects (currT,newT)
because the computation of the parsimony score overwrites the inCycle attribute of the Nodes, which messes with
the search moves.

Extensions:
- other criteria: hardwired, parental (only softwired implemented now)
- remove level-1 restriction: this will involve changing the proposedTop! function to use rSPR or rNNI moves
  (instead of the level-1 moves coded for snaq!).
  We need:
  - functions for rSPR and rNNI moves
  - create new proposedTop! function (proposedRootedTop?) to choose from the rSPR/rNNI/other moves
  - have the search in maxParsimonuNetRun1! to have rooted currT and rooted newT, instead of keeping semi-directed objects
    (currT, newT), only to root for the parsimony score (currTr, newTr)
- outgroup is currently String or Node number, but it would be good if it allowed Edge number as an option too.
  Not sure best way to distinguish between Node number and Edge number, which is why left as Node number for now.
"""
function maxParsimonyNetRun1!(currT::HybridNetwork, tolAbs::Float64, Nfail::Integer, df::DataFrame, hmax::Integer,
                              logfile::IO, writelog::Bool, outgroup::Union{AbstractString,Integer},
                              criterion=:softwired::Symbol)
    tolAbs >= 0 || error("tolAbs must be greater than zero: $(tolAbs)")
    Nfail > 0 || error("Nfail must be greater than zero: $(Nfail)")
    @debug begin printEverything(currT); "printed everything" end
    CHECKNET && checkNet(currT)
    count = 0
    movescount = zeros(Int,18) #1:6 number of times moved proposed, 7:12 number of times success move (no intersecting cycles, etc.), 13:18 accepted by parsimony score
    movesfail = zeros(Int,6) #count of failed moves for current topology
    failures = 0
    stillmoves = true
    Nmov = zeros(Int,6)
    species, traits = readCSVtoArray(df)
    currTr = deepcopy(currT)
    rootatnode!(currTr,outgroup)
    currT.loglik = parsimonyGF(currTr,species,traits,criterion)
    absDiff = tolAbs + 1
    newT = deepcopy(currT)
    writelog && write(logfile, "\nBegins heuristic search of most parsimonious network------\n")
    while absDiff > tolAbs && failures < Nfail && stillmoves
        count += 1
        calculateNmov!(newT,Nmov)
        move = whichMove(newT,hmax,movesfail,Nmov)
        if move != :none
            newT0 = deepcopy(newT) ## to go back if proposed topology conflicts with the outgroup
            flag = proposedTop!(move,newT,true, count,10, movescount,movesfail,false) #N=10 because with 1 it never finds an edge for nni
            newTr = deepcopy(newT) ##rooted version only to compute parsimony
            try
                rootatnode!(newTr,outgroup; verbose=false)
            catch
                flag = false
            end

            if flag #no need else in general because newT always undone if failed, but needed for conflicts with root
                @debug "successful move and correct root placement"
                accepted = false
                newT.loglik = parsimonyGF(newTr,species,traits,criterion)
                accepted = (newT.loglik < currT.loglik && abs(newT.loglik-currT.loglik) > tolAbs) ? true : false
                #newT better parsimony score: need to check for error or keeps jumping back and forth
                if accepted
                    absDiff = abs(newT.loglik - currT.loglik)
                    currT = deepcopy(newT)
                    failures = 0
                    movescount[move2int[move]+12] += 1
                    movesfail = zeros(Int,6) #count of failed moves for current topology
                else
                    failures += 1
                    movesfail[move2int[move]] += 1
                    newT = deepcopy(currT)
                end
            else
                @debug "unsuccessful move or incorrect root placement"
                newT = newT0 ## not counting errors in outgroup as failures, maybe we should
            end
        else
            stillmoves = false
        end
    end
    assignhybridnames!(newT)
    if absDiff <= tolAbs
        writelog && write(logfile,"\nSTOPPED by absolute difference criteria")
    elseif !stillmoves
        writelog && write(logfile,"\nSTOPPED for not having more moves to propose: movesfail $(movesfail), Nmov $(Nmov)")
    else
        writelog && write(logfile,"\nSTOPPED by number of failures criteria")
    end
    writelog && write(logfile,"\nEND: found minimizer topology at step $(count) (failures: $(failures)) with parsimony score=$(round(newT.loglik, digits=5))")
    writelog && printCounts(movescount,zeros(Int,13),logfile) ## zeroes in lieu of movesgamma, not used in parsimony
    setBLGammaParsimony!(newT)
    return newT
end


## find the maximum parsimony network;
## transform the starting topology first
## does not allow multiple alleles
@doc (@doc maxParsimonyNetRun1!) maxParsimonyNetRun1
function maxParsimonyNetRun1(currT0::HybridNetwork, df::DataFrame, Nfail::Integer, tolAbs::Float64,
                                      hmax::Integer,seed::Integer,logfile::IO, writelog::Bool, probST::Float64,
                                      outgroup::Union{AbstractString,Integer}, criterion=:softwired::Symbol)
    Random.seed!(seed)
    currT = findStartingTopology!(currT0, probST, false,writelog, logfile, outgroup=outgroup)
    net = maxParsimonyNetRun1!(currT, tolAbs, Nfail, df, hmax,logfile,writelog, outgroup, criterion)
    return net
end


## find the most parsimonious network over multiple runs
## no multiple alleles for now;
## if rootname not defined, it does not save output files
## fixit: now it only works if currT0 is tree, or level-1 network
## also, throws an error if outgroup not compatible with starting network
## (instead of choosing another root)
"""
    maxParsimonyNet(T::HybridNetwork, df::DataFrame)

Search for the most parsimonious network (or tree).
A level-1 network is assumed.
`df` should be a data frame containing the species names in column 1,
or in a column named `species` or `taxon`. Trait data are assumed to be
in all other columns. The search starts from topology `T`,
which can be a tree or a network with no more than `hmax` hybrid nodes
(see optional arguments below for `hmax`).

Output:

- estimated network in file `.out` (also in `.log`): best network overall and list of
  networks from each individual run.
- if any error occurred, file `.err` provides information (seed) to reproduce the error.

Optional arguments include

- hmax: maximum number of hybridizations allowed (default 1)
- runs: number of starting points for the search (default 10);
  each starting point is `T` with probability `probST`=0.3 or a
  modification of `T` otherwise (using a NNI move, or a hybrid edge
  direction change)
- Nfail: number of failures (proposed networks with equal or worse score)
  before the search is aborted. 75 by default: this is quite small,
  which is okay for a first trial. Larger values are recommended.
- outgroup: outgroup taxon.
            It can be a taxon name (String) or Node number (Integer).
            If none provided, or if the outgroup conflicts the starting
            topology, the function returns an error
- filename: root name for the output files. Default is "mp". If empty (""),
  files are *not* created, progress log goes to the screen only (standard out).
- seed: seed to replicate a given search
- criterion: parsimony score could be hardwired, softwired (default) or parental. Currently,
             only softwired is implemented

# References

1. Leo Van Iersel, Mark Jones, Celine Scornavacca (2017).
   Improved Maximum Parsimony Models for Phylogenetic Networks,
   Systematic Biology,
   (https://doi.org/10.1093/sysbio/syx094).

2. Fischer, M., van Iersel, L., Kelk, S., Scornavacca, C. (2015).
   On computing the Maximum Parsimony score of a phylogenetic network.
   SIAM J. Discrete Math., 29(1):559-585.

For a roadmap of the functions inside maxParsimonyNet, see [`maxParsimonyNetRun1!`](@ref).
"""
function maxParsimonyNet(currT::HybridNetwork, df::DataFrame;
    tolAbs=fAbs::Float64, Nfail=numFails::Integer,
    hmax=1::Integer, runs=10::Integer, outgroup="none"::Union{AbstractString,Integer},
    rootname="mp"::AbstractString, seed=0::Integer, probST=0.3::Float64,
    criterion=:softwired::Symbol)

    currT0 = readTopologyUpdate(writeTopologyLevel1(currT)) # update all level-1 things
    flag = checkNet(currT0,true) # light checking only
    flag && error("starting topology suspected not level-1")

    writelog = true
    writelog_1proc = false
    if (rootname != "")
        julialog = string(rootname,".log")
        logfile = open(julialog,"w")
        juliaout = string(rootname,".out")
        if Distributed.nprocs() == 1
            writelog_1proc = true
            juliaerr = string(rootname,".err")
            errfile = open(juliaerr,"w")
        end
    else
      writelog = false
      logfile = stdout # used in call to optTopRun1!
    end
    str = """optimization of topology using:
              hmax = $(hmax),
              max number of failed proposals = $(Nfail).
             """

    ## need to root the network in a good place before parsimony
    ## throw an error if outgroup conflicts with starting topology,
    ## instead of choosing another outgroup
    outgroup == "none" && error("provide a sensible outgroup for parsimony score computation")
    rootatnode!(currT,outgroup) ##currT (not currT0) bc we just want to check that starting topology
                                ##is not in conflict with outgroup,
                                ##but we need a semi-directed level-1 "good" network (currT0) for search

    str *= (writelog ? "rootname for files: $(rootname)\n" : "no output files\n")
    str *= "BEGIN: $(runs) runs on starting tree $(writeTopology(currT0))\n"
    if Distributed.nprocs()>1
        str *= "       using $(Distributed.nprocs()) processors\n"
    end
    str *= "       with parsimony criterion: $(string(criterion))\n"
    if (writelog)
        write(logfile,str)
        flush(logfile)
    end
    print(stdout,str)
    print(stdout, Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s") * "\n")
    # if 1 proc: time printed to logfile at start of every run, not here.

    if seed == 0
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) #better seed based on clock
    end
    if  writelog
        write(logfile,"\nmain seed $(seed)\n")
        flush(logfile)
    else print(stdout,"\nmain seed $(seed)\n"); end
    Random.seed!(seed)
    seeds = [seed;round.(Integer,floor.(rand(runs-1)*100000))]
    if writelog && !writelog_1proc
        for i in 1:runs # workers won't write to logfile
            write(logfile, "seed: $(seeds[i]) for run $(i)\n")
        end
        flush(logfile)
    end

    tstart = time_ns()
    bestnet = Distributed.pmap(1:runs) do i # for i in 1:runs
        logstr = "seed: $(seeds[i]) for run $(i), $(Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s"))\n"
        print(stdout, logstr)
        msg = "\nBEGIN Max Parsimony search for run $(i), seed $(seeds[i]) and hmax $(hmax)"
        if writelog_1proc # workers can't write on streams opened by master
            write(logfile, logstr * msg)
            flush(logfile)
        end
        GC.gc();
        try
            best = maxParsimonyNetRun1(currT0, df, Nfail, tolAbs, hmax, seeds[i],logfile,writelog_1proc,
                                                probST,outgroup,criterion);
            logstr *= "\nFINISHED Max $(string(criterion)) parsimony for run $(i), parsimony of best: $(best.loglik)\n"
            if writelog_1proc
                logstr = writeTopology(best)
                logstr *= "\n---------------------\n"
                write(logfile, logstr)
                flush(logfile)
            end
            return best
        catch(err)
            msg = "\nERROR found on Max Parsimony for run $(i) seed $(seeds[i]): $(err)\n"
            logstr = msg * "\n---------------------\n"
            if writelog_1proc
                write(logfile, logstr)
                flush(logfile)
                write(errfile, msg)
                flush(errfile)
            end
            @warn msg # returns: nothing
        end
    end
    tend = time_ns() # in nanoseconds
    telapsed = round(convert(Int64, tend-tstart) * 1e-9, digits=2) # in seconds
    writelog_1proc && close(errfile)
    msg = "\n" * Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s")
    if writelog
        write(logfile, msg)
    end
    filter!(n -> n !== nothing, bestnet) # remove "nothing", failed runs
    if length(bestnet)>0
        ind = sortperm([n.loglik for n in bestnet])
        bestnet = bestnet[ind]
        maxNet = bestnet[1]::HybridNetwork # tell type to compiler
    else
        error("all runs failed")
    end

    rootatnode!(maxNet,outgroup)
    writelog &&
    write(logfile,"\nMaxNet is $(writeTopology(maxNet)) \nwith $(string(criterion)) parsimony score $(maxNet.loglik)\n")
    print(stdout,"\nMaxNet is $(writeTopology(maxNet)) \nwith $(string(criterion)) parsimony score $(maxNet.loglik)\n")

    s = writelog ? open(juliaout,"w") : stdout
    str = writeTopology(maxNet) * """
     $(string(criterion)) parsimony score = $(maxNet.loglik)
     Dendroscope: $(writeTopology(maxNet,di=true))
     Elapsed time: $(telapsed) seconds, $(runs) attempted runs
    -------
    List of estimated networks for all runs (sorted by $(string(criterion)) parsimony score; the smaller, the better):
    """
    for n in bestnet
        str *= " "
        str *= writeTopology(rootatnode!(n,outgroup))
        str *= ", with $(string(criterion)) parsimony $(n.loglik)\n"
    end
    str *= "-------\n"
    write(s,str);
    writelog && close(s) # to close juliaout file (but not stdout!)
    writelog && close(logfile)

    return maxNet
end


## unused function:
## function to check if currT agrees with the outgroup
function correctRootPlace(currT::HybridNetwork, outgroup::AbstractString)
    try
        checkRootPlace!(currT,outgroup=outgroup)
    catch err
        if isa(err, RootMismatch)
            println("RootMismatch: ", err.msg,
                    """\nThe starting topology has hybrid edges that are incompatible with the desired outgroup.
                    """)
        else
            println("error trying to reroot: ", err.msg);
        end
        return false
    end
    return true
end

