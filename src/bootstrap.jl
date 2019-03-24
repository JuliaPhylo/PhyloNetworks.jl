# julia functions for bootstrap
# Claudia October 2015
# Cecile April 2016

"""
    readBootstrapTrees(listfile; relative2listfile=true)

Read the list of file names in `listfile`, then read all the trees in each of
these files. Output: vector of vectors of trees (networks with h>0 allowed).

`listfile` should be the name of a file containing the path/name to multiple
bootstrap files, one on each line (no header). Each named bootstrap file should
contain multiple trees, one per line (such as bootstrap trees from a single gene).

The path/name to each bootstrap file should be relative to `listfile`.
Otherwise, use option `relative2listfile=false`, in which case the file names
are interpreted as usual: relative to the user's current directory
if not given as absolute paths.
"""
function readBootstrapTrees(filelist::AbstractString; relative2listfile=true::Bool)
    filelistdir = dirname(filelist)
    bootfiles = CSV.read(filelist, header=false, types=Dict(1=>String))
    size(bootfiles)[2] > 0 ||
        error("there should be a column in file $filelist: with a single bootstrap file name on each row (no header)")
    ngenes = size(bootfiles)[1]
    bf = (relative2listfile ? joinpath.(filelistdir, bootfiles[1]) : bootfiles[1])
    treelists = Array{Vector{HybridNetwork}}(undef, ngenes)
    for igene in 1:ngenes
        treelists[igene] = readMultiTopology(bf[igene])
        print("read $igene/$ngenes bootstrap tree files\r") # using \r for better progress display
    end
    return treelists
end

"""
    sampleBootstrapTrees(vector of tree lists; seed=0::Integer, generesampling=false, row=0)
    sampleBootstrapTrees!(tree list, vector of tree lists; seed=0::Integer, generesampling=false, row=0)

Sample bootstrap gene trees, 1 tree per gene.
Set the seed with keyword argument `seed`, which is 0 by default.
When `seed=0`, the actual seed is set using the clock.
Assumes a vector of vectors of networks (see `readBootstrapTrees`),
each one of length 1 or more (error if one vector is empty, tested in `bootsnaq`).

- site resampling: always, from sampling one bootstrap tree from each given list.
  This tree is sampled at **random** unless `row>0` (see below).
- gene resampling: if `generesampling=true` (default is false),
  genes (i.e. lists) are sampled with replacement.
- `row=i`: samples the ith bootstrap tree for each gene.
  `row` is turned back to 0 if gene resampling is true.

output: one vector of trees. the modifying function (!) modifies the input tree list and returns it.
"""
function sampleBootstrapTrees(trees::Vector{Vector{HybridNetwork}};
                              seed=0::Integer, generesampling=false::Bool, row=0::Integer)
    bootTrees = Array{HybridNetwork}(undef, length(trees))
    sampleBootstrapTrees!(bootTrees, trees, seed=seed, generesampling=generesampling, row=row)
end

function sampleBootstrapTrees!(bootTrees::Vector{HybridNetwork}, trees::Vector{Vector{HybridNetwork}};
                              seed=0::Integer, generesampling=false::Bool, row=0::Integer)
    numgen = length(trees) ## number of genes
    numgen>0 || error("needs at least 1 array of trees")
    numgen <= length(bootTrees) || error("the input tree list needs to be of length $numgen at least")
    if (generesampling) row=0; end
    if row==0
      if seed == 0
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) #better seed based on clock
        println("using seed $(seed) for bootstrap trees")
      end
      Random.seed!(seed)
      if generesampling
        indxg = sample(1:numgen, numgen) # default is with replacement. good!
      end
    end
    for g in 1:numgen
        ig = (generesampling ? indxg[g] : g )
        if row==0
          indxt = sample(1:length(trees[ig]),1)[1]
        else
          indxt = row
          length(trees[ig]) >= row || error("gene $g has fewer than $row bootstrap trees.")
        end
        bootTrees[g] = trees[ig][indxt]
    end
    return bootTrees
end


"""
    sampleCFfromCI(data frame, seed=0)
    sampleCFfromCI!(data frame, seed=0)

Read a data frame containing CFs and their credibility intervals, and
sample new obsCF uniformly within the CIs.
These CFs are then rescaled to sum up to 1 for each 4-taxon sets.
Return a data frame with taxon labels in first 4 columns, sampled obs CFs in columns 5-7
and credibility intervals in columns 8-13.

- The non-modifying function creates a new data frame (with re-ordered columns) and returns it.
  If `seed=-1`, the new df is a deep copy of the input df, with no call to the random
  number generator. Otherwise, `seed` is passed to the modifying function.
- The modifying function overwrites the input data frame with the sampled CFs and returns it.
  If `seed=0`, the random generator is seeded from the clock. Otherwise the random generator
  is seeded using `seed`.

Warning: the modifying version does *not* check the data frame: assumes correct columns.

optional argument: `delim=','` by default: how columns are delimited.
"""
function sampleCFfromCI(df::DataFrame, seed=0::Integer)
    @debug "order of columns should be: t1,t2,t3,t4,cf1234,cf1324,cf1423,cf1234LO,cf1234HI,..."
    size(df,2) == 13 || size(df,2) == 14 || @warn "sampleCFfromCI function assumes table from TICR: CF, CFlo, CFhi"
    obsCFcol = [findfirst(isequal(:CF12_34), DataFrames.names(df)),
                findfirst(isequal(:CF13_24), DataFrames.names(df)),
                findfirst(isequal(:CF14_23), DataFrames.names(df))]
    nothing ∉ obsCFcol || error("""CF columns were not found: should be named like 'CF12_34'""")
    obsCFcol == [5,8,11] ||
        @warn """CF columns were found, but not in the expected columns.
                Lower/upper bounds of credibility intervals assumed in columns 6,7, 9,10 and 12,13."""
    colsTa = [1,2,3,4]        # column numbers for taxon names
    colsCI = [6,7,9,10,12,13] # for lower/upper CI bounds
    length(findall(in(obsCFcol), colsTa)) ==0 ||
        error("CFs found in columns 1-4 where taxon labels are expected")
    length(findall(in(obsCFcol), colsCI)) ==0 ||
        error("CFs found in columns where credibility intervals are expected")
    newdf = deepcopy(df[ [colsTa; obsCFcol; colsCI] ])
    if seed==-1
      return newdf
    else
      return sampleCFfromCI!(newdf::DataFrame, seed)
    end
end

function sampleCFfromCI!(df::DataFrame, seed=0::Integer)
    if seed == 0
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) #better seed based on clock
        println("using seed $(seed) for bootstrap table")
    end
    Random.seed!(seed)
    for i in 1:size(df,1)
        c1 = (df[i, 9]-df[i, 8])*rand()+df[i, 8]
        c2 = (df[i,11]-df[i,10])*rand()+df[i,10]
        c3 = (df[i,13]-df[i,12])*rand()+df[i,12]
        suma = c1+c2+c3
        df[5][i] = c1/suma
        df[6][i] = c2/suma
        df[7][i] = c3/suma
    end
    return df
end

sampleCFfromCI(file::AbstractString; delim=','::Char,seed=0::Integer) = sampleCFfromCI(CSV.read(file, delim=delim),seed)

# function that will do bootstrap of snaq estimation in series
# it repeats optTopRuns nrep times
# it has the same arguments as optTopRuns except for:
# - need data table of CF with conf intervals (instead of d DataCF),
#   or vector of vector of HybridNetworks
# - new argument nrep: number of bootstrap replicates (default 10)
# - new argument runs2: percentage of bootstrap replicates to start in the best network, by default 0.25
# - new argument bestNet: to start the optimization. if prcnet>0.0 and bestNet is not input as argument from a previous run, it will estimate it inside
# - quartetfile if it was used in original data ("none" if all quartets used)
# recall: optTopRuns! does *not* modify its input starting network
function optTopRunsBoot(currT0::HybridNetwork, data::Union{DataFrame,Vector{Vector{HybridNetwork}}},
                        hmax::Integer, liktolAbs::Float64, Nfail::Integer, ftolRel::Float64,ftolAbs::Float64,xtolRel::Float64,xtolAbs::Float64,
                        verbose::Bool, closeN::Bool, Nmov0::Vector{Int},
                        runs1::Integer, outgroup::AbstractString, filename::AbstractString, seed::Integer, probST::Float64,
                        nrep::Integer, runs2::Integer, bestNet::HybridNetwork, quartetfile::AbstractString)
    println("BOOTSTRAP OF SNAQ ESTIMATION")
    writelog = true
    if filename != ""
        logfile = open(string(filename,".log"),"w")
        write(logfile, "BOOTSTRAP OF SNAQ ESTIMATION \n")
    else
        writelog = false
        logfile = stdout
    end

    inputastrees = isa(data, Vector{Vector{HybridNetwork}})
    inputastrees || isa(data, DataFrame) ||
        error("Input data not recognized: $(typeof(data))")

    if runs1>0 && runs2>0
        str = """Will use this network as starting topology for $runs1 run(s) for each bootstrap replicate:
                 $(writeTopologyLevel1(currT0))
                 and this other network for $runs2 run(s):
                 $(writeTopologyLevel1(bestNet))
                 """
        writelog && write(logfile, str)
        print(str)
    end

    if inputastrees # allocate memory, to be re-used later
        newtrees = sampleBootstrapTrees(data, row=1)
        newd = readTrees2CF(newtrees, quartetfile=quartetfile, writeTab=false, writeSummary=false)
        taxa = unionTaxa(newtrees)
    else
        newdf = sampleCFfromCI(data, -1) # column names check, newdf has obsCF in columns 5-7
        # seed=-1: deep copy only, no rand()
        newd = readTableCF!(newdf, collect(1:7)) # allocate memory for DataCF object
    end
    if runs1>0 && isTree(currT0) # get rough first estimate of branch lengths in startnet
        updateBL!(currT0, newd)
    end

    if seed == 0
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) #better seed based on clock
    end
    println("main seed $(seed)")
    writelog && write(logfile,"\nmain seed $(seed)\n")
    writelog && flush(logfile)
    Random.seed!(seed)
    seedsData = round.(Int,floor.(rand(nrep)*100000)) # seeds to sample bootstrap data
    if runs1>0
        seeds = round.(Int,floor.(rand(nrep)*100000)) # seeds for all optimizations from currT0
    end
    if runs2>0
      seedsOtherNet = round.(Int,floor.(rand(nrep)*100000)) # for runs starting from other net
    end

    bootNet = HybridNetwork[]

    writelog && write(logfile,"\nBEGIN: $(nrep) replicates\n$(Libc.strftime(time()))\n")
    writelog && flush(logfile)
    for i in 1:nrep
        str = "\nbegin replicate $(i)\nbootstrap data simulation: seed $(seedsData[i])\n"
        writelog && write(logfile, str)
        print(str)
        if !inputastrees
            sampleCFfromCI!(newdf, seedsData[i])
            readTableCF!(newd, newdf, [5,6,7])
        else
            sampleBootstrapTrees!(newtrees, data, seed=seedsData[i])
            calculateObsCFAll!(newd,taxa) # do not use readTrees2CF: to save memory and gc time
        end
        if runs1>0
            str = "estimation, $runs1 run" * (runs1>1 ? "s" : "") * ": seed $(seeds[i])\n"
            writelog && write(logfile, str)
            print(str)
            rootname = ""
            @debug begin rootname = string(filename,"_",i);
                         "rootname set to $rootname"; end
            net1 = optTopRuns!(currT0, liktolAbs, Nfail, newd, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0, runs1, outgroup,
                               rootname,seeds[i],probST)
            if runs2==0
                net = net1
            end
        end
        if runs2>0
            str = "estimation, $runs2 run" * (runs2>1 ? "s" : "") * " starting from other net: seed $(seedsOtherNet[i])\n"
            writelog && write(logfile, str)
            print(str)
            rootname = ""
            @debug begin rootname = string(filename,"_",i,"_startNet2");
                         "rootname set to $rootname"; end
            net2 = optTopRuns!(bestNet, liktolAbs, Nfail, newd, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0, runs2, outgroup,
                               rootname,seedsOtherNet[i],probST)
            if runs1==0
                net = net2
            end
        end
        if runs1>0 && runs2>0
            net = (net1.loglik < net2.loglik ? net1 : net2)
        end
        writelog && flush(logfile)

        push!(bootNet, net)
        str = (outgroup=="none" ? writeTopologyLevel1(net) : writeTopologyLevel1(net,outgroup))
        if writelog
            write(logfile, str)
            write(logfile,"\n")
            flush(logfile)
        end
        println(str) # net also printed by each optTopRuns! but separately from the 2 starting points
    end # of the nrep bootstrap replicates
    writelog && close(logfile)

    if writelog
      s = open(string(filename,".out"),"w")
      for n in bootNet
        if outgroup == "none"
            write(s,"$(writeTopologyLevel1(n))\n")
        else
            write(s,"$(writeTopologyLevel1(n,outgroup))\n")
        end
        # "with -loglik $(n.loglik)" not printed: not comparable across bootstrap networks
      end
      close(s)
    end
    return bootNet
end

# like snaq, only calls optTopRunsBoot
# undocumented arguments: closeN, Nmov0
"""
    bootsnaq(T::HybridNetwork, df::DataFrame)
    bootsnaq(T::HybridNetwork, vector of tree lists)

Bootstrap analysis for SNaQ.
Bootstrap data can be quartet concordance factors (CF),
drawn from sampling uniformly in their credibility intervals,
as given in the data frame `df`.
Alternatively, bootstrap data can be gene trees sampled from
a vector of tree lists: one list of bootstrap trees per locus
(see `readBootstrapTrees` to generate this,
from a file containing a list of bootstrap files: one per locus).

From each bootstrap replicate, a network
is estimated with snaq!, with a search starting from topology `T`.
Optional arguments include the following, with default values in parentheses:

- `hmax` (1): max number of reticulations in the estimated networks
- `nrep` (10): number of bootstrap replicates.
- `runs` (10): number of independent optimization runs for each replicate
- `filename` ("bootsnaq"): root name for output files. No output files if "".
- `seed` (0 to get a random seed from the clock): seed for random number generator
- `otherNet` (empty): another starting topology so that each replicate will start prcnet% runs on otherNet and (1-prcnet)% runs on `T`
- `prcnet` (0): percentage of runs starting on `otherNet`; error if different than 0.0, and otherNet not specified.
- `ftolRel`, `ftolAbs`, `xtolRel`, `xtolAbs`, `liktolAbs`, `Nfail`,
  `probST`, `verbose`, `outgroup`: see `snaq!`, same defaults.

If `T` is a tree, its branch lengths are first optimized roughly with [`updateBL!`](@ref)
(by using the average CF of all quartets defining each branch and calculating the coalescent units
corresponding to this quartet CF).
If `T` has one or more reticulations, its branch lengths are taken as is to start the search.
The branch lengths of `otherNet` are always taken as is to start the search.
"""
function bootsnaq(startnet::HybridNetwork, data::Union{DataFrame,Vector{Vector{HybridNetwork}}};
                  hmax=1::Integer, liktolAbs=likAbs::Float64, Nfail=numFails::Integer,
                  ftolRel=fRel::Float64, ftolAbs=fAbs::Float64, xtolRel=xRel::Float64, xtolAbs=xAbs::Float64,
                  verbose=false::Bool, closeN=true::Bool, Nmov0=numMoves::Vector{Int},
                  runs=10::Integer, outgroup="none"::AbstractString, filename="bootsnaq"::AbstractString,
                  seed=0::Integer, probST=0.3::Float64, nrep=10::Integer, prcnet=0.0::Float64,
                  otherNet=HybridNetwork()::HybridNetwork, quartetfile="none"::AbstractString)

    inputastrees = isa(data, Vector{Vector{HybridNetwork}})
    inputastrees || isa(data, DataFrame) ||
        error("Input data not recognized: $(typeof(data))")

    if !inputastrees
    (DataFrames.names(data)[[6,7,9,10,12,13]] == [:CF12_34_lo,:CF12_34_hi,:CF13_24_lo,:CF13_24_hi,:CF14_23_lo,:CF14_23_hi]) ||
      @warn """assume table with CI from TICR: CFlo, CFhi in columns 6,7; 9,10; and 12,13.
              Found different column names: $(DataFrames.names(data)[[6,7,9,10,12,13]])"""
    else # check 1+ genes, each with 1+ trees, all with h=0.
        ngenes = length(data)
        ngenes > 0 || error("empty list of bootstrap trees (0 genes)")
        for igene in 1:ngenes
            btr = data[igene]
            length(btr) > 0 || error("no bootstrap trees for $(igene)th gene")
            for itree in 1:length(btr)
                btr[itree].numHybrids == 0 || error("network $itree is not a tree for $(igene)th gene")
            end
        end
    end
    prcnet >= 0 || error("percentage of times to use the best network as starting topology should be positive: $(prcnet)")
    prcnet = (prcnet <= 1.0) ? prcnet : prcnet/100
    runs2 = round(Int, runs*prcnet)            # runs starting from otherNet
    runs1 = runs - runs2                       # runs starting from startnet

    if runs1>0
        startnet=readTopologyLevel1(writeTopologyLevel1(startnet)) # does not modify startnet outside
        flag = checkNet(startnet,true) # light checking only
        flag && error("starting topology suspected not level-1")
        try
            checkNet(startnet)
        catch err
            println("starting topology is not of level 1:")
            rethrow(err)
        end
    end
    runs2 == 0 || otherNet.numTaxa > 0 ||
        error("""otherNet not given and prcnet>0. Please set prcnet to 0 to start optimizations
                from the same network always, or else provide an other network "otherNet"
                to start the optimization from this other network in pcrnet % of runs.""")
    if runs2 > 0
        otherNet=readTopologyLevel1(writeTopologyLevel1(otherNet))
        flag = checkNet(otherNet,true) # light checking only
        flag && error("starting topology 'otherNet' suspected not level-1")
        try
            checkNet(otherNet)
        catch err
            println("starting topology 'otherNet' not a level 1 network:")
            rethrow(err)
        end
    end

    # for multiple alleles: expand into two leaves quartets like sp1 sp1 sp2 sp3.
    if (@isdefined originald) && !isempty(originald.repSpecies) ## not defined if treefile empty, but not needed
        expandLeaves!(originald.repSpecies,startnet)
    end

    optTopRunsBoot(startnet,data,hmax, liktolAbs, Nfail,ftolRel, ftolAbs, xtolRel, xtolAbs,
                   verbose, closeN, Nmov0, runs1, outgroup, filename,
                   seed, probST, nrep, runs2, otherNet, quartetfile)
end

"""
`treeEdgesBootstrap(boot_net::Vector{HybridNetwork}, ref_net::HybridNetwork)`

Read a list of bootstrap networks (`boot_net`) and a reference network (`ref_net`),
and calculate the bootstrap support of the tree edges in the reference network.
All minor hybrid edges (γ<0.5) are removed to extract the major tree from
each network. All remaining edges are tree edges, each associated with a bipartition.

output:

- a data frame with one row per tree edge and two columns: edge number, bootstrap support
  (as a percentage)
- the major tree from the reference network, where minor hybrid edges (with γ<0.5)
  have been removed.
"""
function treeEdgesBootstrap(net::Vector{HybridNetwork}, net0::HybridNetwork)
    # estimated network, major tree and matrix
    S = tipLabels(net0)
    tree0 =majorTree(net0)
    M0 = tree2Matrix(tree0,S, rooted=false)

    M = Matrix[]
    tree = HybridNetwork[]
    for n in net
        t = majorTree(n)
        push!(tree,t)
        mm = tree2Matrix(t,S, rooted=false)
        push!(M,mm)
    end

    df = DataFrame(edgeNumber=Int[], proportion=Float64[])

    for i in 1:size(M0,1) #rows in M0: internal edges
        cnt = 0 #count
        for j in 1:length(M) #for every M
            for k in 1:size(M[j],1) #check every row in M
                if M0[i,2:end] == M[j][k,2:end] || M0[i,2:end] == map(x->(x+1)%2,M[j][k,2:end]) #same row
                    #println("found same row: $(M0[i,2:end]) and $(M[j][k,2:end])")
                    cnt += 1
                    break
                end
            end
        end
        push!(df,[M0[i,1] cnt*100/length(M)])
        # fixit later, when edges have an attribute for bootstrap support:
        # set BS to cnt/length(M) for the edge numbered M0[i,1].
    end
    @info """edge numbers in the data frame correspond to the current edge numbers in the network.
       If the network is modified, the edge numbers in the (modified) network might not correspond
       to those in the bootstrap table. Plot the bootstrap values onto the current network with
       plot(network_name, edgeLabel=bootstrap_table_name)"""
    return df, tree0
end


"""
    hybridDetection(net::Vector{HybridNetwork}, net1::HybridNetwork, outgroup::AbstractString)

function can only compare hybrid nodes in networks that have the same underlying major tree
also, need to root all networks in the same place, and the root has to be compatible with the
direction of the hybrid edges

it computes the rooted hardwired distance between networks, the root matters.
input: vector of bootstrap networks (net), estimated network (net1), outgroup

returns

- a matrix with one row per bootstrap network, and 2*number of hybrids in net1,
column i corresponds to whether hybrid i (net1.hybrid[i]) is found in the bootstrap network,
column 2i+1 corresponds to the estimated gamma on the bootstrap network
(0.0 if hybrid not found).
To know the order of hybrids, print net1.hybrid[i] i=1,...,num of hybrids

- list of discrepant trees (trees not matching the main tree in net1)
"""
function hybridDetection(net::Vector{HybridNetwork}, net1::HybridNetwork, outgroup::AbstractString)
    tree1 = majorTree(net1)
    rootnet1 = deepcopy(net1)
    rootatnode!(rootnet1,outgroup)

    # HF for "hybrid found?"
    HFmat = zeros(length(net),net1.numHybrids*2)

    # discrepant trees
    discTrees = HybridNetwork[]
    # major trees
    majorTrees = HybridNetwork[]

    i = 1
    for n in net
        tree = majorTree(n)
        push!(majorTrees,tree)
        RFmajor = hardwiredClusterDistance(tree, tree1, false)
        if RFmajor != 0
            push!(discTrees,tree)
            i+=1
            continue # skip replicate if major tree is *not* correct
        end

        found = zeros(Bool,   net1.numHybrids) # repeats false
        gamma = zeros(Float64,net1.numHybrids) # repeats 0.0
        # re-root estimated network if not rooted correctly
        reroot = true
        if length(n.node[n.root].edge) == 2 # check if root connects to correct outgroup
            for e in n.node[n.root].edge
                for node in e.node
                    if node.name == outgroup
                        reroot = false
                        break
                    end
                end
                if (!reroot) break end
            end
        end
        !reroot || println("Will need to reroot the estimated network...")
        for trueh = 1:net1.numHybrids
            netT = deepcopy(rootnet1)
            displayedNetworkAt!(netT, netT.hybrid[trueh]) # bug: need correct attributes to re-root later...
            for esth = 1:n.numHybrids
                netE = deepcopy(n)
                displayedNetworkAt!(netE, netE.hybrid[esth])
                if reroot
                    rootatnode!(netE, outgroup) # if re-rooting is not possible,
                end                         # then the hybridization doesn't match.
                if hardwiredClusterDistance(netT, netE, true) == 0 # true: rooted
                    found[trueh] = true
                    node = netE.hybrid[1]
                    edges = hybridEdges(node)
                    edges[2]. hybrid || error("edge should be hybrid")
                    !edges[2]. isMajor || error("edge should be minor hybrid")
                    edges[2].gamma <= 0.5 || error("gamma should be less than 0.5")
                    gamma[trueh] = edges[2].gamma
                    break # to exit loop over esth
                end
            end
        end
        HFmat[i,1:net1.numHybrids] = found
        HFmat[i,(net1.numHybrids+1):end] = gamma
        i+=1
    end
    treeMatch=size(HFmat[sum(HFmat,2).>0,:],1) #number of bootstrap trees that match tree1
    println("$(treeMatch) out of $(length(net)) bootstrap major trees match with the major tree in the estimated network")
    println("order of hybrids:")
    for h in net1.hybrid
        println("$(h.name)")
    end
    return HFmat,discTrees
end


"""
`hybridBootstrapSupport(boot_net::Vector{HybridNetwork}, ref_net::HybridNetwork; rooted=false)`

Match hybrid nodes in a reference network with those in an array of networks,
like bootstrap networks.
All networks must be fully resolved, and on the same taxon set.
If `rooted=true`, all networks are assumed to have been properly rooted beforehand.
Otherwise, the origin of each hybrid edge is considered as an unrooted bipartition (default).

Two hybrid edges in two networks are said to match if they share the same "hybrid" clade
(or recipient) and the same "donor clade", which is a sister to the hybrid clade in the network.
Since a hybrid clade has 2 parent edges, it is sister to two clades simultaneously: one is its
major sister (following the major hybrid edge with γ>0.5) and one is its minor sister (following
the major hybrid edge with γ<0.5).

To calculate these hybrid and sister clades at a given hybrid node, all other hybrid edges are
first removed from the network. Then, the hybrid clade is the hardwired cluster (descendants) of
either hybrid edge and major/minor clade is the hardwired cluster of the sibling edge of the
major/minor hybrid parent.
If `rooted=false`, sister clades are considered as bipartitions.

Output:

1. a "node" data frame (see below)
2. an "edge" data frame (see below)
3. a "clade" data frame to describe the make up of all clades found as hybrids or sisters,
   starting with a column `taxa` that lists all taxa. All other columns correspond to a given
   clade and contain true/false values. `true` means that a given taxon belongs in a given clade.
   For a clade named `H1`, for instance, and if the data frame was named `cla`, the
   list of taxa in this clade can be obtained with `cla[:taxa][cla[:H1]]`.
4. an array of gamma values, with one row for each bootstrap network and two columns (major/minor)
   for each hybrid edge in the reference network. If this hybrid edge was found in the bootstrap network
   (i.e. same hybrid and sister clades, after removal of all other hybrid nodes),
   its bootstrap gamma value is recorded here. Otherwise, the gamma entry is 0.0.
5. a vector with the number of each hybrid edge in the reference network,
   in the same order as for the columns in the array of gamma values above.

The "node" data frame has one row per clade and 9 columns giving:

   - **clade**: the clade's name, like the taxon name (if a hybrid is a single taxon) or
     the hybrid tag (like 'H1') in the reference network
   - **node**: the node number in the reference network. missing if the clade is not in this network.
   - **hybridnode**: typically the same node number as above, except for hybrid clades in the
     reference network. For those, the hybrid node number is listed here.
   - **edge**: number of the parent edge, parent to the node in column 2,
     if found in the ref network. missing otherwise.
   - **BS_hybrid**: percentage of bootstrap networks in which the clade is found to be a hybrid clade.
   - **BS_sister**: percentage of bootstrap networks in which the clade is found to be sister to
     some hybrid clade (sum of the next 2 columns)
   - **BS_major_sister**: percentage of bootstrap networks in which the clade is found to be the
     major sister to some hybrid clade
   - **BS_minor_sister**: same as previous, but minor
   - **BS_hybrid_samesisters**: percentage of bootstrap networks in which the clade is found to be
     a hybrid and with the same set of sister clades as in the reference network.
     Applies to hybrid clades found in the reference network only, missing for all other clades.

The "edge" data frame has one row for each pair of clades, and 8 columns:

  - **edge**: hybrid edge number, if the edge appears in the reference network. missing otherwise.
  - **hybrid_clade**: name of the clade found to be a hybrid, descendent of 'edge'
  - **hybrid**: node number of that clade, if it appears in the reference network. missing otherwise.
  - **sister_clade**: name of the clade that is sister to 'edge', i.e. be sister to a hybrid
  - **sister**: node number of that clade, if in the ref network.
  - **BS_hybrid_edge**: percentage of bootstrap networks in which 'edge' is found to be a hybrid
     edge, i.e. when the clade in the 'hybrid' column is found to be a hybrid and the clade in
     the 'sister' column is one of its sisters.
  - **BS_major**: percentage of bootstrap networks in which 'edge' is found to be a major hybrid
     edge, i.e. when 'hybrid' is found to be a hybrid clade and 'sister' is found to be its
     major sister.
  - **BS_minor**: same as previous, but minor
"""
function hybridBootstrapSupport(nets::Vector{HybridNetwork}, refnet::HybridNetwork;
         rooted=false::Bool)
    numNets = length(nets)
    numNets>0 || error("there aren't any test (bootstrap) networks")
    numHybs = refnet.numHybrids
    numHybs>0 || error("there aren't any hybrid in reference network")
    try directEdges!(refnet)
    catch err
        if isa(err, RootMismatch)
            err.msg *= "\nPlease change the root in reference network (see rootatnode! or rootonedge!)"
        end
        rethrow(err)
    end
    taxa = tipLabels(refnet)
    ntax = length(taxa)

    # extract hardwired clusters of each tree edge in major tree of reference net,
    # and of each hybrid edge after all other hybrids are removed following the major edge.
    clade = Vector{Bool}[] # list all clades in reference and in bootstrap networks
    treenode = Int[]     # node number for each clade: of tree node if found in reference net
    treeedge = Int[]     # number of tree edge corresponding to the clade (parent of treenode)
    leafname = AbstractString[] # "" if internal node, leaf name if leaf
    # clade, treenode, treeedge: same size. leafname and hybparent: same size initially only
    hybind = Int[]       # indices in 'clade', that appear as hybrid clades in reference
    hybnode = Int[]      # for those hybrid clades: number of hybrid node
    majsisedge = Int[]   # for those hybrid clades: number of major sister edge (in ref net)
    minsisedge = Int[]   #                                    minor
    majsisind = Int[]    # for those hybrid clades: index of major sister clade in 'clade'
    minsisind = Int[]    #                                   minor
    # hybind, hybnode, majsis*, minsis*: same size = number of hybrid nodes in reference network

    if !rooted && length(refnet.node[refnet.root].edge)==2 && any(e -> e.hybrid, refnet.node[refnet.root].edge)
        refnet = deepcopy(refnet) # new binding inside function
        fuseedgesat!(refnet.root, refnet)
    end # issues otherwise: correct tree edge for root bipartitions, find sister clades, ...

    reftre = majorTree(refnet)
    skipone = (!rooted && length(reftre.node[reftre.root].edge)<3) # not count same bipartition twice
    for pe in reftre.edge
        hwc = hardwiredCluster(pe,taxa) # not very efficient, but human readable
        if skipone && refnet.node[refnet.root] ≡ getParent(pe) && sum(hwc)>1
            skipone = false             # wrong algo for trivial 2-taxon rooted tree (A,B);
            println("skip edge $(pe.number)")
        else
            push!(clade, hwc)
            push!(treeedge, pe.number)
            cn = getChild(pe) # child node of pe
            push!(treenode, cn.number)
            push!(leafname, (cn.leaf ? cn.name : ""))
        end
    end
    hybparent = zeros(Int,length(clade)) # 0 if has no hybrid parent node. Index in hybnode otherwise.
    for trueh = 1:numHybs
        net0 = deepcopy(refnet)
        displayedNetworkAt!(net0, net0.hybrid[trueh]) # removes all minor hybrid edges but one
        hn = net0.hybrid[1]
        hemaj, hemin, ce = hybridEdges(hn) # assumes no polytomy at hybrid nodes and correct node.hybrid
        (hemin.hybrid && !hemin.isMajor) || error("edge should be hybrid and minor")
        (hemaj.hybrid &&  hemaj.isMajor) || error("edge should be hybrid and major")
        ic = findfirst(isequal(ce.number), treeedge)
        ic !== nothing || error("hybrid node $(hn.number): child edge not found in major tree")
        hybparent[ic] = trueh
        push!(hybind,ic)
        push!(hybnode, hn.number)
        push!(majsisedge, hemaj.number)
        push!(minsisedge, hemin.number)
        for sis in ["min","maj"]
          he = (sis=="min" ? hemin : hemaj)
          pn = getParent(he) # parent node of sister (origin of gene flow if minor)
          atroot = (!rooted && pn ≡ net0.node[net0.root]) # polytomy at root pn of degree 3: will exclude one child edge
          hwc = zeros(Bool,ntax) # new binding each time. pushed to clade below.
          for ce in pn.edge    # important if polytomy
            if ce ≢ he && pn ≡ getParent(ce)
                hw = hardwiredCluster(ce,taxa)
                if atroot && any(hw & clade[ic]) # sister clade intersects child clade
                    (hw & clade[ic]) == clade[ic] ||
                        @warn "weird clusters at the root in reference, hybrid node $(hn.number)"
                else
                    hwc .|= hw
                end
            end
          end
          i = findfirst(isequal(hwc), clade)
          if (!rooted && i===nothing) i = findfirst(isequal(.!hwc), clade) end
          i !== nothing || error(string("hyb node $(hn.number): ",sis,"or clade not found in main tree"))
          if (sis=="min") push!(minsisind, i)
          else            push!(majsisind, i)
          end
          if sis=="min"  # need to get clade not in main tree: hybrid + minor sister
            pe = nothing # looking for the (or one) parent edge of pn
            for ce in pn.edge
              if pn ≡ getChild(ce)
                  pe=ce
                  break
              end
            end
            # pe == nothing if minor hybrid is at root and rooted=true. Will just miss edge and node number
            # for that clade, but in that case, the minor sister might have been assigned that edge anyway...
            if pe != nothing
              hwc = hardwiredCluster(pe,taxa)
              i = findfirst(isequal(hwc), clade) # i>0: (hybrid + minor sister) can be in main tree if
              # hybrid origin is ancestral, i.e. hybrid clade is nested within minor sister.
              if i===nothing
                push!(clade, hwc)
                push!(treenode, pn.number)
                push!(treeedge, pe.number)
                push!(hybparent, 0)
                push!(leafname, (pn.leaf ? pn.name : ""))
              end
            end
          end
        end
    end
    # for cl in clade @show taxa[cl]; end; @show treenode; @show treeedge; @show hybparent; @show leafname
    # @show hybind; @show hybnode; @show majsisedge; @show minsisedge; @show majsisind; @show minsisind

    # node-associated summaries:
    nclades = length(clade)
    nh = length(hybnode)
    nedges = 2*nh
    BShyb        = zeros(Float64, nclades) # clade = hybrid
    BSmajsis     = zeros(Float64, nclades) # clade = major sister of some hybrid
    BSminsis     = zeros(Float64, nclades) # clade = minor sister of some hybrid
    # edge-associated summaries, i.e. associated to a pair of clades (hyb,sis):
    hybcladei   = repeat(hybind, inner=[2]) # indices in 'clade'
    siscladei   = Array{Int}(undef, nedges) # edge order: (major then minor) for all hybrids
    edgenum     = Array{Int}(undef, nedges)
    for i=1:nh
        siscladei[2*i-1] = majsisind[i]
        siscladei[2*i]   = minsisind[i]
        edgenum[2*i-1] = majsisedge[i]
        edgenum[2*i]   = minsisedge[i]
    end # less desirable order: [majsisind; minsisind], corresponds to outer=[2] for hybcladei
    BShybmajsis = zeros(Float64, nedges) # one clade = hybrid, one clade = major sister
    BShybminsis = zeros(Float64, nedges) # one clade = hybrid, one clade = minor sister
    # network*edge-associated values: only keep those for edges in ref net
    gamma = zeros(Float64,numNets,nedges)
    # node-associate 3-way partition summary:
    BShybsamesis = zeros(Float64, nh) # same hybrid & sister set as in ref net

    nextnum = max(maximum([n.number for n in refnet.edge]),
                  maximum([n.number for n in refnet.node]) ) + 1

    for i = 1:numNets
        net = nets[i]
        length(tipLabels(net))==ntax || error("networks have non-matching taxon sets")
        try directEdges!(net) # make sure the root is admissible
        catch err
          if isa(err, RootMismatch)
            err.msg *= "\nPlease change the root in test network (see rootatnode! or rootatedge!)"
          end
          rethrow(err)
        end
        for esth = 1:net.numHybrids   # try to match estimated hybrid edge
            hwcPar = zeros(Bool,ntax)   # minor sister clade
            hwcChi = zeros(Bool,ntax)   #       child  clade
            hwcSib = zeros(Bool,ntax)   # major sister clade
            net1 = deepcopy(net)
            displayedNetworkAt!(net1, net1.hybrid[esth])
            hn = net1.hybrid[1]
            if !rooted && length(net1.node[net1.root].edge)==2 && any(e -> e.hybrid, net1.node[net1.root].edge)
                fuseedgesat!(net1.root, net1)
            end
            hemaj, hemin, ce = hybridEdges(hn) # assumes no polytomy at hybrid node
            (hemin.hybrid && !hemin.isMajor) || error("edge should be hybrid and minor")
            (hemaj.hybrid &&  hemaj.isMajor) || error("edge should be hybrid and major")
            hardwiredCluster!(hwcChi,hemin,taxa)
            for sis in ["min","maj"]
                he = (sis=="min" ? hemin : hemaj)
                pn = getParent(he) # parent of hybrid edge
                atroot = (!rooted && pn ≡ net1.node[net1.root])
                # if at root: exclude the child edge in the same cycle as he.
                # its cluster includes hwcChi. all other child edges do not interest hwcChi.
                # if (atroot) @show i; @warn "$(sis)or edge is at the root!"; end
                for ce in pn.edge
                  if ce ≢ he && pn ≡ getParent(ce)
                    hwc = hardwiredCluster(ce,taxa)
                    if !atroot || sum(hwc .& hwcChi) == 0 # empty intersection
                      if (sis=="maj") hwcSib .|= hwc;
                      else            hwcPar .|= hwc; end
                    elseif (hwc .& hwcChi) != hwcChi
                        @warn "weird clusters at the root. bootstrap net i=$i, hybrid $(net.hybrid[esth].name)"
                    end
                  end
                end # will use complement too: test network may be rooted differently
            end
            # @show taxa[hwcChi]; @show taxa[hwcPar]
            if all(hwcPar) || all(hwcSib) || all(.!hwcPar) || all(.!hwcSib)
                @warn "parent or sibling cluster is full or empty. bootstrap net i=$i, hybrid $(net.hybrid[esth].name)"
            end

            ihyb = findfirst(isequal(hwcChi), clade)
            newhyb = (ihyb===nothing) # hwcChi not found in clade list
            if newhyb
              push!(clade, hwcChi)
              push!(treenode, nextnum)
              push!(treeedge, nextnum)
              push!(BShyb,    1.0)
              push!(BSmajsis, 0.0)
              push!(BSminsis, 0.0)
              ihyb = length(clade) # was nothing; converted to nex available integer
              nextnum += 1
            else
              BShyb[ihyb] += 1.0
            end

            iSmaj = findfirst(isequal(hwcSib), clade)
            iSmin = findfirst(isequal(hwcPar), clade)
            if (!rooted && iSmaj===nothing) iSmaj = findfirst(isequal(.!hwcSib), clade) end
            if (!rooted && iSmin===nothing) iSmin = findfirst(isequal(.!hwcPar), clade) end
            newmaj = (iSmaj===nothing)
            if newmaj
              push!(clade, hwcSib)
              push!(treenode, nextnum)
              push!(treeedge, nextnum)
              push!(BShyb,    0.0)
              push!(BSmajsis, 1.0)
              push!(BSminsis, 0.0)
              iSmaj = length(clade) # was nothing; now integer
              nextnum += 1
            else
              BSmajsis[iSmaj] += 1.0
            end
            newmin = (iSmin===nothing)
            if newmin
              push!(clade, hwcPar)
              push!(treenode, nextnum)
              push!(treeedge, nextnum)
              push!(BShyb,    0.0)
              push!(BSmajsis, 0.0)
              push!(BSminsis, 1.0)
              iSmin = length(clade) # was nothing; now integer
              nextnum += 1
            else
              BSminsis[iSmin] += 1.0
            end
            samehyb = (newhyb ? Int[] : findall(isequal(ihyb), hybcladei))
            for sis in ["maj","min"]
                newpair = newhyb || (sis=="min" ? newmin : newmaj)
                iSsis = (sis=="min" ? iSmin : iSmaj)
                if !newpair # hyb and sis clades were already found, but not sure if together
                    iish = findfirst(isequal(iSsis), siscladei[samehyb])
                    if iish !== nothing # pair was indeed found
                        ipair = samehyb[iish]
                        if (sis=="min") BShybminsis[ipair] += 1.0
                        else            BShybmajsis[ipair] += 1.0 end
                        if ipair <= nedges # pair = edge in reference network
                            gamma[i,ipair] = (sis=="min" ? hemin.gamma : hemaj.gamma)
                        end
                    else newpair = true; end
                end
                if newpair
                    push!(hybcladei, ihyb)
                    push!(siscladei, iSsis)
                    push!(BShybminsis, (sis=="min" ? 1.0 : 0.0))
                    push!(BShybmajsis, (sis=="min" ? 0.0 : 1.0))
                end
            end
            if ihyb<=nclades && hybparent[ihyb]>0 # hyb clade match in ref net
                th = hybparent[ihyb]
                if ((iSmaj==majsisind[th] && iSmin==minsisind[th]) ||
                    (iSmin==majsisind[th] && iSmaj==minsisind[th]) )
                    BShybsamesis[th] += 1
                end
            end
        end # loop over hybrid nodes
    end # loop over bootstrap nets

    fac = 100/numNets
    BShyb    *= fac # size: length(clade)
    BSmajsis *= fac
    BSminsis *= fac
    BShybmajsis  *= fac # size: length of hybcladei or siscladei
    BShybminsis  *= fac
    BShybsamesis *= fac # size: length(hybnode) = nedges/2

    # combine results into 2 dataframes: node & edge summaries
    # and 1 dataframe to describe clades
    # first: detect nodes with all BS=0
    keepc = ones(Bool,length(clade)) # keep clade h in output tables?
    for h=nclades:-1:1
        if BShyb[h]==0.0 && BSmajsis[h]==0.0 && BSminsis[h]==0.0
            keepc[h] = false
            deleteat!(BShyb,h); deleteat!(BSmajsis,h); deleteat!(BSminsis,h)
        end
    end
    nkeepc = sum(keepc)
    # clade descriptions
    resCluster = DataFrame(taxa=taxa)
    cladestr = Array{String}(undef, length(clade))
    rowh = 1
    for h=1:length(clade)
        nn = treenode[h] # node number in ref net
        cladestr[h] = string("c_", (nn<0 ? "minus" : ""), abs(nn))
        # column name for a clade at node 5: "c_5". At node -5: "c_minus5"
        # because symbol :c_-5 causes an error.
        if h <= nclades &&  leafname[h] != "" # replace "c5" by leaf name, e.g. "taxon1"
            cladestr[h] = leafname[h]
        end
        if h <= nclades &&  hybparent[h]>0 # replace "c5" or leaf name by hybrid name, e.g. "H1"
            na = refnet.hybrid[hybparent[h]].name
            cladestr[h] = (na=="" ? string("H", refnet.hybrid[hybparent[h]].number) : replace(na, r"^#" => ""))
        end
        if keepc[h]
            rowh += 1
            insertcols!(resCluster, rowh, Symbol(cladestr[h]) => clade[h])
        end
    end
    # node summaries
    resNode = DataFrame(clade=cladestr[keepc],
                        node=allowmissing(treenode[keepc]),
                        hybridnode=allowmissing(treenode[keepc]),
                        edge=allowmissing(treeedge[keepc]),
                        BS_hybrid=BShyb, BS_sister = BSmajsis + BSminsis,
                        BS_major_sister=BSmajsis, BS_minor_sister=BSminsis,
                        BS_hybrid_samesisters=Vector{Union{Missing,Float64}}(undef, nkeepc))
    rowh = 1
    for h=1:length(clade)
        if h <= nclades && keepc[h] && hybparent[h]>0
            resNode[:hybridnode][rowh]            =      hybnode[hybparent[h]]
            resNode[:BS_hybrid_samesisters][rowh] = BShybsamesis[hybparent[h]]
        elseif keepc[h]
            resNode[:BS_hybrid_samesisters][rowh] = missing
        end
        if h>nclades # clade *not* in the reference network
            resNode[:node][rowh] = missing
            resNode[:hybridnode][rowh] = missing
            resNode[:edge][rowh] = missing
        end
        if keepc[h]  rowh += 1; end
    end
    insertcols!(resNode, 10, :BS_all => resNode[:BS_hybrid]+resNode[:BS_sister])
    sort!(resNode, [:BS_all,:BS_hybrid]; rev=true)
    deletecols!(resNode, :BS_all)
    # edge summaries
    resEdge = DataFrame(edge = Vector{Union{Int, Missing}}(undef, length(hybcladei)),
                        hybrid_clade=cladestr[hybcladei],
                        hybrid=Vector{Union{Int, Missing}}(treenode[hybcladei]),
                        sister_clade=cladestr[siscladei],
                        sister=Vector{Union{Int, Missing}}(treenode[siscladei]),
                        BS_hybrid_edge = BShybmajsis+BShybminsis,
                        BS_major=BShybmajsis, BS_minor=BShybminsis)
    for i=1:length(hybcladei)
        h = hybcladei[i]
        if h <= nclades && hybparent[h]>0
            resEdge[:hybrid][i] = hybnode[hybparent[h]]
        end
        if h>nclades            resEdge[:hybrid][i]=missing; end
        if siscladei[i]>nclades resEdge[:sister][i]=missing; end
        if i <= nedges
             resEdge[:edge][i] = edgenum[i]
        else resEdge[:edge][i] = missing
        end
    end
    o = [1:nedges; sortperm(resEdge[:BS_hybrid_edge][nedges+1:length(hybcladei)],rev=true) .+ nedges]
    return resNode, resEdge[o,:], resCluster, gamma, edgenum
end

"""
    summarizeHFdf(HFmat::Matrix)

Summarize data frame output from [`hybridDetection`](@ref).
Output: dataframe with one row per hybrid, and 5 columns:

- hybrid index (order from estimated network, see [`hybridDetection`](@ref),
- number of bootstrap trees that match the
  underlying tree of estimated network
- number of bootstrap networks that have the hybrid
- mean estimated gamma in the bootstrap networks that have the hybrid
- sd estimated gamma in the bootstrap networks that have the hybrid also

last row has index -1, and the third column has the number of networks
that have all hybrids (hybrid index, mean gamma, sd gamma are
meaningless in this last row)
"""
function summarizeHFdf(HFmat::Matrix)
    HFmat2 = HFmat[sum(HFmat,2) .>0,:]
    gt = size(HFmat2,1)
    total = size(HFmat,1)
    numH = round(Int,size(HFmat,2)/2)
    df = DataFrame(hybrid=Int[],goodTrees=Float64[],netWithHybrid=Float64[],meanGamma=Float64[], sdGamma=Float64[])
    for i in 1:numH
        mat = HFmat2[HFmat2[:,i] .> 0, :]
        n = size(mat,1)
        g = mean(mat[:,round(Int,numH+i)])
        s = std(mat[:,round(Int,numH+i)])
        push!(df,[i gt n g s])
    end
    which = Bool[]
    for i in 1:size(HFmat2,1)
        push!(which,sum(HFmat2[i,1:numH]) == numH)
    end
    push!(df, [-1 gt sum(which) -1.0 -1.0])
    return df
end
