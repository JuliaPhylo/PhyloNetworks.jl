# julia functions for bootstrap
# Claudia October 2015

# function to read a CF table with CI
# and sample new obsCF
# the function returns the new dataframe
# warning: for some reason UTF8String is needed instead of AbstractString
# seed: to choose the rand() for the table
function bootstrapCFtable(df::DataFrame;seed=0::Int)
    warn("bootstrapCFtable function assumes table from TICR: CF, CFlo, CFhi")
    DEBUG && warn("order of columns should be: t1,t2,t3,t4,cf1234,cf1324,cf1423,cf1234LO,cf1234HI,...")
    size(df,2) == 13 || warn("Dataframe should have 7 columns: 4taxa, 3CF*3")
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
        c1 = c1/(c1+c2+c3)
        c2 = c2/(c1+c2+c3)
        c3 = c3/(c1+c2+c3)
        append!(newdf,DataFrame(t1=convert(UTF8String,string(df[i,1])), t2=convert(UTF8String,string(df[i,2])), t3=convert(UTF8String,string(df[i,3])), t4=convert(UTF8String,string(df[i,4])), CF12_34=c1,CF13_24=c2, CF14_23=c3))
    end
    return newdf
end

bootstrapCFtable(file::AbstractString;sep=','::Char) = bootstrapCFtable(readtable(file,separator=sep))


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
    warn("bootsnaq function not debugged yet")
    prcnet >= 0 || error("percentage of times to use the best network as starting topology should be positive: $(prcnet)")
    prcnet = (prcnet <= 1.0) ? prcnet : prcnet/100
    println("BOOTSTRAP OF SNAQ ESTIMATION")
    julialog = string(filename,".log")
    logfile = open(julialog,"w")
    write(logfile, "BOOTSTRAP OF SNAQ ESTIMATION \n")

    if(seed == 0)
        t = time()/1e9
        a = split(string(t),".")
        seed = int(a[2][end-4:end]) #better seed based on clock
    end
    write(logfile,"\nmain seed $(seed)\n")
    flush(logfile)
    srand(seed)
    seeds = [seed;int(floor(rand(nrep)*100000))] #seeds for all runs
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

    bootNet = HybridNetwork[];

    write(logfile,"\nBEGIN: $(nrep) replicates")
    write(logfile,"\n$(Libc.strftime(time()))")
    flush(logfile)

    for(i in 1:nrep)
        write(logfile,"\n begin replicate $(i) with seed $(seeds[i+1])---------\n")
        println("\nbegin replicate $(i) with seed $(seeds[i+1])\n")
        newdf = bootstrapCFtable(df)
        writetable(string("CFtable",i,".csv"),newdf)
        newd = readTableCF(newdf)
        if(i/nrep <= prcnet)
            write(logfile,"\nStarting topology: best network")
            startnet = deepcopy(bestNet)
        else
            write(logfile,"\nStarting topology: same starting tree as original optimization")
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
    warn("bootsnaq function not debugged yet")
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
        seed = int(a[2][end-4:end]) #better seed based on clock
    end
    srand(seed)
    seeds = [int(floor(rand(nworkers())*100000))] #seeds for all workers


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
    S = SharedArray(Float64,size(df))
    for(i in size(df,1))
        for(j in size(df,2))
            S[i,j] = df[i,j]
        end
    end
    return S
end
