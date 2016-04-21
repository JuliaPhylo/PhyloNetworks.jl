# julia functions for bootstrap
# Claudia October 2015
# Cecile April 2016

# function to read a CF table with CI
# and sample new obsCF
# the function returns the new dataframe
# warning: for some reason UTF8String is needed instead of AbstractString
# seed: to choose the rand() for the table
function bootstrapCFtable(df::DataFrame;seed=0::Int)
    warn("bootstrapCFtable function assumes table from TICR: CF, CFlo, CFhi")
    DEBUG && warn("order of columns should be: t1,t2,t3,t4,cf1234,cf1324,cf1423,cf1234LO,cf1234HI,...")
    size(df,2) == 13 || size(df,2) == 14 || warn("bootstrapCFtable function assumes table from TICR: CF, CFlo, CFhi")
    newdf = DataFrames.DataFrame(t1=UTF8String[],t2=UTF8String[],t3=UTF8String[],t4=UTF8String[],CF12_34=0.,CF13_24=0.,CF14_23=0.)
    if(seed == 0)
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) #better seed based on clock
    end
    println("using seed $(seed) for bootstrap table")
    srand(seed)
    for(i in 1:size(df,1))
        c1 = (df[i,7]-df[i,6])*rand()+df[i,6] #fixit: check this is uniform
        c2 = (df[i,10]-df[i,9])*rand()+df[i,9]
        c3 = (df[i,13]-df[i,12])*rand()+df[i,12]
        suma = c1+c2+c3
        c1 = c1/suma
        c2 = c2/suma
        c3 = c3/suma
        append!(newdf,DataFrame(t1=convert(UTF8String,string(df[i,1])), t2=convert(UTF8String,string(df[i,2])), t3=convert(UTF8String,string(df[i,3])), t4=convert(UTF8String,string(df[i,4])), CF12_34=c1,CF13_24=c2, CF14_23=c3))
    end
    return newdf
end

bootstrapCFtable(file::AbstractString;sep=','::Char,seed=0::Int) = bootstrapCFtable(readtable(file,separator=sep),seed=seed)


# function that will do bootstrap of snaq estimation in series
# it repeats optTopRuns nrep times
# it has the same arguments as optTopRuns except for:
# - need df table of CF with conf intervals (instead of d DataCF)
# - new argument nrep: number of bootstrap replicates (default 10)
# - new argument prcnet: percentage of bootstrap replicates to start in the best network, by default 0.25
# - new argument bestNet: to start the optimization. if prcnet>0.0 and bestNet is not input as argument from a previous run, it will estimate it inside
# returns vector of HybridNetworks, bestNet is the first one, and the other nrep networks after
# I believe no ! needed because we need a clean copy of currT in each replicate, so deepcopied
function optTopRunsBoot(currT0::HybridNetwork, df::DataFrame, hmax::Int64, M::Number, Nfail::Int64,ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64, verbose::Bool, closeN::Bool, Nmov0::Vector{Int64}, runs::Int64, outgroup::AbstractString, filename::AbstractString, returnNet::Bool, seed::Int64, probST::Float64, nrep::Int64, prcnet::Float64, bestNet::HybridNetwork)
    prcnet >= 0 || error("percentage of times to use the best network as starting topology should be positive: $(prcnet)")
    prcnet = (prcnet <= 1.0) ? prcnet : prcnet/100
    println("BOOTSTRAP OF SNAQ ESTIMATION")
    julialog = string(filename,".log")
    logfile = open(julialog,"w")
    write(logfile, "BOOTSTRAP OF SNAQ ESTIMATION \n")

    if(seed == 0)
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) #better seed based on clock
    end
    write(logfile,"\nmain seed $(seed)\n")
    flush(logfile)
    srand(seed)
    seeds = [seed;round(Int64,floor(rand(nrep)*100000))] #seeds for all runs
    bootNet = HybridNetwork[]

    if(prcnet > 0.0)
        write(logfile, "Starting topology: will use the best network $(prcnet*100) percent of times \n")
        if(bestNet.numTaxa == 0)
            write(logfile, "bestNet not given as input, so we need to estimate it before doing bootstrap\n")
            println("bestNet not input, so we need to estimate it before doing bootstrap in order to use it as starting topology in $(prcnet*100) percent of times")
            d = readTableCF(df)
            startnet=deepcopy(currT0)
            bestNet = optTopRuns!(startnet, M, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0, runs, outgroup, "bestNet", true,seeds[1],probST)
            push!(bootNet, deepcopy(bestNet))
        else
            write(logfile,"bestNet input: $(writeTopology(bestNet))\n to use in $(prcnet *100) percent of times as starting topology")
            println("bestNet input: $(writeTopology(bestNet))\n to use in $(prcnet *100) percent of times as starting topology")
            push!(bootNet, deepcopy(bestNet))
        end
    end

    write(logfile,"\nBEGIN: $(nrep) replicates")
    write(logfile,"\n$(Libc.strftime(time()))")
    flush(logfile)

    for(i in 1:nrep)
        write(logfile,"\n begin replicate $(i) with seed $(seeds[i+1])---------\n")
        println("\nbegin replicate $(i) with seed $(seeds[i+1])\n")
        newdf = bootstrapCFtable(df, seed=seeds[i+1])
#        writetable(string("CFtable",i,".csv"),newdf)
        newd = readTableCF(newdf)
        if(i/nrep <= prcnet)
            write(logfile,"\nStarting topology: best network")
            startnet = deepcopy(bestNet)
        else
            write(logfile,"\nStarting topology: same starting tree as original optimization\n")
            startnet=deepcopy(currT0)
        end
        flush(logfile)
        net = optTopRuns!(startnet, M, Nfail, newd, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0, runs, outgroup, string(filename,"_",i), true,seeds[i+1],probST)
        push!(bootNet, deepcopy(net))
        if(outgroup == "none")
            write(logfile,writeTopology(net)) #no outgroup
        else
            write(logfile,writeTopology(net,outgroup)) #outgroup
        end
        flush(logfile)
    end #end nrep
    close(logfile)

    s = open(string(filename,".out"),"w")
    for(n in bootNet)
        if(outgroup == "none")
            write(s,"\n $(writeTopology(n)), with -loglik $(n.loglik)")
        else
            write(s,"\n $(writeTopology(n,outgroup)), with -loglik $(n.loglik)")
        end
    end
    close(s)
    if(returnNet)
        return bootNet
    end
end

# like snaq, only calls optTopRunsBoot
# will later decide which to call depending on nproc()
function bootsnaq(currT0::HybridNetwork, df::DataFrame; hmax=1::Int64, M=multiplier::Number, Nfail=numFails::Int64,ftolRel=fRel::Float64, ftolAbs=fAbs::Float64, xtolRel=xRel::Float64, xtolAbs=xAbs::Float64, verbose=false::Bool, closeN=true::Bool, Nmov0=numMoves::Vector{Int64}, runs=10::Int64, outgroup="none"::AbstractString, filename="bootsnaq"::AbstractString, returnNet=true::Bool, seed=0::Int64, probST=0.3::Float64, nrep=10::Int64, prcnet=0.25::Float64, bestNet=HybridNetwork()::HybridNetwork)
    startnet=deepcopy(currT0)
    if(nprocs() > 1) #more than 1 processor, still not working
        error("bootsnaq not implemented for parallelization yet")
        optTopRunsBootParallel(startnet, df, hmax, M, Nfail,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0, runs, outgroup, filename, returnNet, seed, probST, nrep, prcnet, bestNet)
    else
        optTopRunsBoot(startnet, df, hmax, M, Nfail,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0, runs, outgroup, filename, returnNet, seed, probST, nrep, prcnet, bestNet)
    end
end


# same as optTopRunsBoot but for many processors in parallel
# warning: still not debugged
function optTopRunsBootParallel(currT0::HybridNetwork, df::DataFrame, hmax::Int64, M::Number, Nfail::Int64,ftolRel::Float64, ftolAbs::Float64, xtolRel::Float64, xtolAbs::Float64, verbose::Bool, closeN::Bool, Nmov0::Vector{Int64}, runs::Int64, outgroup::AbstractString, filename::AbstractString, returnNet::Bool, seed::Int64, probST::Float64, nrep::Int64, prcnet::Float64, bestNet::HybridNetwork)
    warn("bootsnaq function not debugged yet")
    prcnet > 0 || error("percentage of times to use the best network as starting topology should be positive: $(prcnet)")
    prcnet = (prcnet <= 1.0) ? prcnet : prcnet/100
    println("BOOTSTRAP OF SNAQ ESTIMATION")
    # shared arrays, variables
    dfS = convert2SharedArray(df) #fixit: need to code
    intS = convert2SharedArray(hmax,M,Nfail,runs,Nmov0) #fixit: need to code
    floatS = convert2SharedArray(ftolRel,ftolAbs,xtolRel,xtolAbs,probST,prcnet) #fixit: need to code
    @everywhere verbose,closeN,outgroup,filename

    # split of replicates
    nrep_proc = floor(nrep/nworkers())
    nrep_missing = nrep - nrep_proc * nworkers()
    bootNet = HybridNetwork[]

    # seeds
    if(seed == 0)
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) #better seed based on clock
    end
    srand(seed)
    seeds = [parse(Int,floor(rand(nworkers())*100000))] #seeds for all workers


    @sync begin
        i = 1
        for p in workers()
            addrep = nrep_missing > 0 ? 1 : 0
            net = @async remotecall_fetch(p,loc_bootsnaq,dfS, intS, floatS, currT, bestNet, seeds[i], nrep_proc+addrep)
            push!(bootNet,net)
            nrep_missing -= addrep
            i += 1
        end
    end
end

# function to the each local run of bootsnaq
# in optTopRunsBootParallel
# dfS = shared array with CF table
# intS = shared array with integer arguments: hmax, M, Nfail, runs, Nmov0
# floatS= shared array with float arguments: ftolRel, ftolAbs, xtolRel, xtolAbs, probST, prcnet
# currT, bestNet are starting topology and best network
# seed= to start the procedure, if seed=0, then clock used
# nrep= number of replicates in this particular processor
# other parameters of optTopRunsBoot are global by @everywhere in optTopRunsBootParallel
function loc_bootsnaq(dfS::SharedArray, intS::SharedArray, floatS::SharedArray, currT::HybridNetwork, bestNet::HybridNetwork, seed::Int, nrep::Int)
    df = bootstrapCFtable(dfS) #fixit: need to code this
    optTopRunsBoot(currT,df,intS[1],intS[2],floatS[1],floatS[2],floatS[3],floatS[4],verbose,closeN,intS[5:end],intS[4],outgroup,string(filename,seed),true,seed,floatS[5],nrep,floatS[6],bestNet)
end


# function to take a DataFrame and convert to SharedArray
# fixit: does not work, S is filled with zeros, but also how to pass the strings of taxon names??
function convert2SharedArray(df::DataFrame)
    error("convert2SharedArray not working, should not be called")
    S = SharedArray(Float64,size(df))
    for(i in size(df,1))
        for(j in size(df,2))
            S[i,j] = df[i,j]
        end
    end
    return S
end

# function that reads the list of bootstrap networks (net), and the estimated network (net0)
# and calculates the bootstrap support of the tree edges in the estimated network
# it returns a data frame with one row per tree edge, and two columns: edge number, bootstrap support
"""
`treeEdgesBootstrap(net::Vector{HybridNetwork}, net0::HybridNetwork)`

read an array of bootstrap networks (net) and a reference network (net0),
and calculates the bootstrap support of the tree edges in the reference network.

return a data frame with one row per tree edge and two columns: edge number, bootstrap support
"""
function treeEdgesBootstrap(net::Vector{HybridNetwork}, net0::HybridNetwork)
    # estimated network, major tree and matrix
    S = tipLabels(net0)
    tree0 =majorTree(net0)
    M0 = tree2Matrix(tree0,S)

    M = Matrix[]
    tree = HybridNetwork[]
    for(n in net)
        t = majorTree(n)
        push!(tree,t)
        mm = tree2Matrix(t,S)
        push!(M,mm)
    end

    df = DataFrame(edge=Int64[], bs=Float64[])

    for(i in 1:size(M0,1)) #rows in M0: internal edges
        cnt = 0 #count
        for(j in 1:length(M)) #for every M
            for(k in 1:size(M[j],1)) #check every row in M
                if(M0[i,2:end] == M[j][k,2:end] || M0[i,2:end] == map(x->(x+1)%2,M[j][k,2:end])) #same row
                    #println("found same row: $(M0[i,2:end]) and $(M[j][k,2:end])")
                    cnt += 1
                    break
                end
            end
        end
        push!(df,[M0[i,1] cnt/length(M)])
    end
    return df, tree0
end


# function can only compare hybrid nodes in networks that have the same underlying major tree
# also, need to root all networks in the same place, and the root has to be compatible with the
# direction of the hybrid edges
# it computes the rooted hardwired distance between networks, the root matters.
# input: vector of bootstrap networks (net), estimated network (net1), outgroup
# returns 1)a matrix with one row per bootstrap network, and 2*number of hybrids in net1,
# 2) list of discrepant trees (trees not matching the main tree in net1)
# column i corresponds to whether hybrid i (net1.hybrid[i]) is found in the bootstrap network,
# column 2i+1 corresponds to the estimated gamma on the bootstrap network (0.0 if hybrid not found)
# to know the order of hybrids, print net1.hybrid[i] i=1,...,num of hybrids
"""
`hybridDetection(net::Vector{HybridNetwork}, net1::HybridNetwork, outgroup::AbstractString)`

function can only compare hybrid nodes in networks that have the same underlying major tree
also, need to root all networks in the same place, and the root has to be compatible with the
direction of the hybrid edges

it computes the rooted hardwired distance between networks, the root matters.
input: vector of bootstrap networks (net), estimated network (net1), outgroup

returns

- a matrix with one row per bootstrap network, and 2*number of hybrids in net1,
column i corresponds to whether hybrid i (net1.hybrid[i]) is found in the bootstrap network,
column 2i+1 corresponds to the estimated gamma on the bootstrap network (0.0 if hybrid not found)

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
    for(n in net)
        tree = majorTree(n)
        push!(majorTrees,tree)
        RFmajor = hardwiredClusterDistance(tree, tree1, false)
        if(RFmajor != 0)
            push!(discTrees,tree)
            i+=1
            continue # skip replicate if major tree is *not* correct
        end

        found = rep(false,net1.numHybrids)
        gamma = rep(0.0,net1.numHybrids)
        # re-root estimated network if not rooted correctly
        reroot = true
        if (length(n.node[n.root].edge) == 2) # check if root connects to correct outgroup
            for (e in n.node[n.root].edge)
                for (node in e.node)
                    if (node.name == outgroup)
                        reroot = false
                        break
                    end
                end
                if (!reroot) break end
            end
        end
        !reroot || println("Will need to reroot the estimated network...")
        for (trueh = 1:net1.numHybrids)
            netT = deepcopy(rootnet1)
            displayedNetworkAt!(netT, netT.hybrid[trueh]) # bug: need correct attributes to re-root later...
            for (esth = 1:n.numHybrids)
                netE = deepcopy(n)
                displayedNetworkAt!(netE, netE.hybrid[esth])
                if (reroot)
                    rootatnode!(netE, outgroup) # if re-rooting is not possible,
                end                         # then the hybridization doesn't match.
                if (hardwiredClusterDistance(netT, netE, true) == 0) # true: rooted
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
    for(h in net1.hybrid)
        println("$(h.name)")
    end
    return HFmat,discTrees
end


"""
`hybridBootstrapFrequency(boot_net::Vector{HybridNetwork}, ref_net::HybridNetwork;rooted=false::Bool)`

Match hybrid nodes in a reference network with those in an array of networks, 
like bootstrap networks.
All networks must be fully resolved, and on the same taxon set.
If `rooted=true`, all networks are assumed to have been properly rooted beforehand.
Otherwise, the origin of each hybrid edge is considered as an unrooted bipartition (default).

Two hybrid edges in two networks are said to match if they share the same "recipient clade" and
the same "donor clade". To calculate these clades, all other hybrid edges are first removed from
both networks. Then, the recipient clade is the hardwired cluster (descendants) of the hybrid
edge and the donor clade (origin) is the hardwired cluster of its sibling edge.
If `rooted=false`, this cluster is considered a bipartition.

- For a **introgression** or **gene flow** event, it makes most sense
  to compare the *minor* hybrid parental edges only, across the two networks.
- For a full **hybridization** event that resulted in a new hybrid
  species, then it makes most sense to compare both hybrid parental edges across the two
  networks, and to allows matching of the minor edge in one network to the major edge in the
  other network (and vice versa).

Output:

- a data frame with one row per hybrid node in the reference network and 3 columns giving:

  1. the hybrid tag
  2. the proportion of bootstrap networks in which the minor hybrid edge matches that in the
  reference network (for matching introgression events)
  3. the proportion of bootstrap networks in which the set of hybrid edges {minor,major}
  matches that in the reference network (for matching hybridization events)

- a data frame with one row for each "recipient" clade. The first rows correspond to the
  hybrids in the reference network. Other rows are from other hybrids found in the bootstrap networks.
  The second columns gives the proportion of bootstrap networks in which a minor hybrid edge had
  a matching recipient cluster.
- a similar data frame with one row for each "donor" clade (or minor origin for hybridization
  events), and the proportion of bootstrap networks in which a minor hybrid edge had a matching
  donor clade.
- a similar data frame with one row for each "sibling" clade (or major origin for hybridization
  events), to summarize the donor clades of the major hybrid edges in the bootstrap networks.
- a data frame to describe all these clusters: for the "recipient", "donor" and "sibling" clades.
  The first column that lists all taxa. Each clade is described by a column of true/false values.
  `true` indicates that a taxon is a descendant of the recipient/donor/sibling lineage.
- an array of gamma values, with one row for each bootstrap network and one column for each hybrid
  in the reference network. The entry is the gamma value of the minor hybrid edge if it was
  found in the network (matching with either criterion: introgression or hybridization),
  or 0.0 if it was not found.
- an array with the edge number of each minor hybrid edge in the reference network.

The arrays with gammas and edge numbers do not have column names, but their columns correspond
to hybrids in the same order in which they are listed in the first data frames,
which do have column names.
"""
function hybridBootstrapFrequency(nets::Vector{HybridNetwork}, refnet::HybridNetwork;
         rooted=false::Bool)
    numNets = length(nets)
    numNets>0 || error("there aren't any test (bootstrap) networks")
    numHybs = refnet.numHybrids
    numHybs>0 || error("there aren't any hybrid in reference network")
    try directEdges!(refnet)
    catch err
        if isa(err, RootMismatch)
            err.msg *= "\nPlease change the root in reference network (see rootatnode! or rootatedge!)"
        end
        rethrow(err)
    end
    taxa = tipLabels(refnet)
    ntax = length(taxa)

    # extract hardwired clusters (origin and recipient) of each hybrid in reference net
    edgeRef = Int64[]
    hwcRefChi = Vector{Bool}[] # recipient (gene flow), hybrid clade (hybridization)
    hwcRefPar = Vector{Bool}[] # donor (gene flow), or minor parent (hybridization)
    hwcRefSib = Vector{Bool}[] # sibling of recipient (gene flow), or major parent (hybridization)
    hwc = zeros(Bool,ntax) # temporary but constant memory storage

    for (trueh = 1:numHybs)
        net0 = deepcopy(refnet)
        displayedNetworkAt!(net0, net0.hybrid[trueh]) # removes all minor hybrid edges but one
        hn = net0.hybrid[1]
        he = hybridEdges(hn)[2] # assumes no polytomy at hybrid nodes, and correct node.hybrid
        (he.hybrid && !he.isMajor) || error("edge should be hybrid and minor")
        push!(edgeRef, he.number)
        push!(hwcRefChi, hardwiredCluster(he,taxa))
        pn = he.node[he.isChild1?2:1] # parent node: origin of hybrid
        push!(hwcRefPar, zeros(Bool,ntax))
        for (ce in pn.edge)   # important if polytomy
            if (ce!=he && pn == ce.node[ce.isChild1?2:1])
                hardwiredCluster!(hwcRefPar[trueh],ce,taxa)
            end
        end
        he = hybridEdges(hn)[1]
        (he.hybrid && he.isMajor) || error("edge should be hybrid and major")
        pn = he.node[he.isChild1?2:1] # major parent node
        push!(hwcRefSib, zeros(Bool,ntax))
        for (ce in pn.edge)
            if (ce!=he && pn == ce.node[ce.isChild1?2:1])
                hardwiredCluster!(hwcRefSib[trueh],ce,taxa)
            end
        end
    end
    # for (h=1:numHybs) @show taxa[hwcRefChi[h]]; @show taxa[hwcRefPar[h]];  @show taxa[hwcRefSib[h]]; end
    # @show edgeRef

    HFintrog = zeros(Bool,numNets,numHybs) # hybrid found? introgression: match donor & recipient, minor edge
    HFhybrid = zeros(Bool,numNets,numHybs) # for full hybridization: match hybrid & both parents
    # one row per minor hybrid edge in reference net
    gamma = zeros(Float64,numNets,numHybs)
    HPfreq = zeros(Float64,numHybs) # all hybrid parents (origins), frequencies only
    HCfreq = zeros(Float64,numHybs) # all hybrid children (recipients)
    HSfreq = zeros(Float64,numHybs) # all        siblings (or major donor)
    # 2 columns: parent found, child found

    for (i = 1:numNets)
        net = nets[i]
        length(tipLabels(net))==ntax || error("networks have non-matching taxon sets")
        try directEdges!(net) # make sure the root is admissible
        catch err
          if isa(err, RootMismatch)
            err.msg *= "\nPlease change the root in test network (see rootatnode! or rootatedge!)"
          end
          rethrow(err)
        end
        for (esth = 1:net.numHybrids)   # try to match estimated hybrid edge
            hwcPar = zeros(Bool,ntax)
            hwcChi = zeros(Bool,ntax)
            hwcSib = zeros(Bool,ntax)
            net1 = deepcopy(net)
            displayedNetworkAt!(net1, net1.hybrid[esth])
            hn = net1.hybrid[1]
            he = hybridEdges(hn)[2] # assumes no polytomy at hybrid node
            (he.hybrid && !he.isMajor) || error("edge should be hybrid and minor")
            hardwiredCluster!(hwcChi,he,taxa)
            pn = he.node[he.isChild1?2:1] # minor parent node: origin of gene flow
            for (ce in pn.edge)
                if (ce!=he && pn == ce.node[ce.isChild1?2:1])
                    hardwiredCluster!(hwcPar,ce,taxa)
                end
            end # will use complement too: test network may be rooted differently
            hem= hybridEdges(hn)[1] # major parent edge
            (hem.hybrid && hem.isMajor) || error("edge should be hybrid and major")
            pn = hem.node[hem.isChild1?2:1] # major parent
            for (ce in pn.edge)
                if (ce!=hem && pn == ce.node[ce.isChild1?2:1])
                    hardwiredCluster!(hwcSib,ce,taxa)
                end
            end
            # @show taxa[hwcChi]; @show taxa[hwcPar]

            hc = findfirst(hwcRefChi, hwcChi)
            hp = findfirst(hwcRefPar, hwcPar) # 0 if not found
            hs = findfirst(hwcRefSib, hwcSib)
            if (!rooted && hp==0) hp = findfirst(hwcRefPar, !hwcPar) end
            if (!rooted && hs==0) hs = findfirst(hwcRefSib, !hwcSib) end
            if hc==0
                push!(hwcRefChi, hwcChi)
                push!(HCfreq, 1.0)
            else HCfreq[hc] += 1.0; end
            if hp==0
                push!(hwcRefPar, hwcPar)
                push!(HPfreq, 1.0)
            else HPfreq[hp] += 1.0 end
            if hs==0
                push!(hwcRefSib, hwcSib)
                push!(HSfreq, 1.0)
            else HSfreq[hs] += 1.0 end
            if (hc>0 && hc<=numHybs) # same hybrid clade
                if hp==hc # same minor parent
                    HFintrog[i,hc] = true
                    gamma[i,hc] = he.gamma # he is minor edge at this point
                    if hs==hc # same major parent
                        HFhybrid[i,hc] = true
                    end
                else # check if match if minor & major are flipped
                    hp = findfirst(hwcRefSib, hwcPar)
                    hs = findfirst(hwcRefPar, hwcSib)
                    if (!rooted && hp==0) hp = findfirst(hwcRefSib, !hwcPar) end
                    if (!rooted && hs==0) hp = findfirst(hwcRefPar, !hwcSib) end
                    if hp==hc && hs==hc
                        HFhybrid[i,hc] = true
                        gamma[i,hc] = hem.gamma # major hem matches minor edge in ref net
                    end
                end
            end
        end
    end
    Hintrogfreq = Array(Float64,numHybs) # overal % detection
    Hhybridfreq = Array(Float64,numHybs)
    for trueh=1:refnet.numHybrids
        Hintrogfreq[trueh] = sum(HFintrog[:,trueh])/numNets
        Hhybridfreq[trueh] = sum(HFhybrid[:,trueh])/numNets
    end
    HPfreq /= numNets
    HCfreq /= numNets
    HSfreq /= numNets

    # change matrices to dataframes, initialize with NA, add column names
    # push! to add rows, insert! to add column

    resCluster = DataFrame(taxa=taxa)
    resfreq  = DataFrame(hybrid=AbstractString[],prop_introgression=Float64[],prop_hybridization=Float64[])
    resCfreq = DataFrame(recipient=AbstractString[],proportion=Float64[])
    resPfreq = DataFrame(donor=AbstractString[],proportion=Float64[])
    resSfreq = DataFrame(sibling=AbstractString[],proportion=Float64[])
    for h=1:length(hwcRefChi)
        str = (h<=numHybs ? string("recipient",refnet.hybrid[h].name) : string("recipient",h-numHybs))
        insert!(resCluster, h+1, hwcRefChi[h], symbol(str))
        push!(resCfreq,[str, HCfreq[h]])
        if (h <= refnet.numHybrids)
            push!(resfreq,[refnet.hybrid[h].name, Hintrogfreq[h], Hhybridfreq[h]])
        end
    end
    h1 = length(hwcRefChi)+1; h2 = length(hwcRefChi) - 2*numHybs
    for h=1:length(hwcRefPar)
        str = (h<=numHybs ? string("donor",refnet.hybrid[h].name) : string("donor",h+h2))
        insert!(resCluster, h+h1, hwcRefPar[h], symbol(str))
        push!(resPfreq,[str, HPfreq[h]])
    end
    h1 += length(hwcRefPar); h2 += length(hwcRefPar) - numHybs
    for h=1:length(hwcRefSib)
        str = (h<=numHybs ? string("sibling",refnet.hybrid[h].name) : string("sibling",h+h2))
        insert!(resCluster, h+h1, hwcRefSib[h], symbol(str))
        push!(resSfreq,[str, HSfreq[h]])
    end
    return resfreq, resCfreq, resPfreq, resSfreq, resCluster, gamma, edgeRef
end



# function to summarize df output from hybridDetection input: HFdf
# (see hybridDetection) returns dataframe with one row per hybrid, and
# 5 columns: - hybrid index (order from estimated network, see
# hybridDetection), - number of bootstrap trees that match the
# underlying tree of estimated network, - number of bootstrap networks
# that have the hybrid - mean estimated gamma in the bootstrap
# networks that have the hybrid - sd estimated gamma in the bootstrap
# networks that have the hybrid also, last row has index -1, and the
# third column has the number of networks that have all hybrids
# (hybrid index, mean gamma, sd gamma are meaningless in this last
# row)
"""
`summarizeHFdf(HFmat::Matrix)`

function to summarize df output from hybridDetection input: HFdf
(see hybridDetection) returns dataframe with one row per hybrid, and
5 columns:

- hybrid index (order from estimated network, see
hybridDetection),

- number of bootstrap trees that match the
underlying tree of estimated network, - number of bootstrap networks
that have the hybrid

- mean estimated gamma in the bootstrap
networks that have the hybrid

- sd estimated gamma in the bootstrap
networks that have the hybrid also

last row has index -1, and the third column has the number of networks
that have all hybrids (hybrid index, mean gamma, sd gamma are
meaningless in this last row)
"""
function summarizeHFdf(HFmat::Matrix)
    HFmat2 = HFmat[sum(HFmat,2) .>0,:]
    gt = size(HFmat2,1)
    total = size(HFmat,1)
    numH = round(Int,size(HFmat,2)/2)
    df = DataFrame(hybrid=Int64[],goodTrees=Float64[],netWithHybrid=Float64[],meanGamma=Float64[], sdGamma=Float64[])
    for(i in 1:numH)
        mat = HFmat2[HFmat2[:,i] .> 0, :]
        n = size(mat,1)
        g = mean(mat[:,round(Int,numH+i)])
        s = std(mat[:,round(Int,numH+i)])
        push!(df,[i gt n g s])
    end
    which = Bool[]
    for(i in 1:size(HFmat2,1))
        push!(which,sum(HFmat2[i,1:numH]) == numH)
    end
    push!(df, [-1 gt sum(which) -1.0 -1.0])
    return df
end
