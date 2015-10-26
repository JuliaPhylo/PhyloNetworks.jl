# julia functions for bootstrap
# Claudia October 2015

# function to read a CF table with CI
# and sample new obsCF
# the function returns the new dataframe
# for some reason UTF8String is needed instead of AbstractString
function bootstrapCFtable(df::DataFrame)
    DEBUG && warn("order of columns should be: t1,t2,t3,t4,cf1234,cf1324,cf1423,cf1234LO,cf1234HI,...")
    size(df,2) == 13 || warn("Dataframe should have 7 columns: 4taxa, 3CF*3")
    newdf = DataFrames.DataFrame(t1=UTF8String[],t2=UTF8String[],t3=UTF8String[],t4=UTF8String[],CF1234=0.,CF1324=0.,CF1423=0.)
    for(i in 1:size(df,1))
        c1 = (df[i,9]-df[i,8])*0.5*randn()+df[i,5]
        c2 = (df[i,11]-df[i,10])*0.5*randn()+df[i,6]
        c3 = (df[i,13]-df[i,12])*0.5*randn()+df[i,7]
        c1 = max(0.0,c1)
        c2 = max(0.0,c2)
        c3 = max(0.0,c3)
        c1 = min(1.0,c1)
        c2 = min(1.0,c2)
        c3 = min(1.0,c3)
        append!(newdf,DataFrame(t1=convert(UTF8String,string(df[i,1])), t2=convert(UTF8String,string(df[i,2])), t3=convert(UTF8String,string(df[i,3])), t4=convert(UTF8String,string(df[i,4])), CF1234=c1,CF1324=c2, CF1423=c3))
    end
    return newdf
end

bootstrapCFtable(file::AbstractString;sep=','::Char) = bootstrapCFtable(readtable(file,separator=sep))


# function that will do bootstrap of snaq estimation
# it has the same arguments as snaq except for:
# - need df table of CF with conf intervals (instead of d DataCF)
# - new argument nrep: number of bootstrap replicates (default 10)
# - new argument prcnet: percentage of bootstrap replicates to start in the best network, by default 0.25
# - new argument bestNet: to start the optimization. if prcnet>0.0 and bestNet is not input as argument from a previous run, it will estimate it inside
function bootsnaq(currT0::HybridNetwork, df::DataFrame; hmax=1::Int64, M=multiplier::Number, Nfail=numFails::Int64,ftolRel=fRel::Float64, ftolAbs=fAbs::Float64, xtolRel=xRel::Float64, xtolAbs=xAbs::Float64, verbose=false::Bool, closeN=true::Bool, Nmov0=numMoves::Vector{Int64}, runs=10::Int64, outgroup="none"::AbstractString, filename="bootsnaq_main"::AbstractString, returnNet=true::Bool, seed=0::Int64, probST=0.3::Float64, nrep=10::Int64, prcnet=0.25::Float64, bestNet=HybridNetwork()::HybridNetwork)
    prcnet > 0 || error("percentage of times to use the best network as starting topology should be positive: $(prcnet)")
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
    seeds = [seed,int(floor(rand(nrep)*100000))]

    if(prcnet > 0.0)
        write(logfile, "Starting topology: will use the best network $(prcnet*100) percent of times \n")
        if(bestNet.numTaxa == 0)
            write(logfile, "bestNet not input, so we need to estimate it before doing bootstrap\n")
            println("bestNet not input, so we need to estimate it before doing bootstrap in order to use it as starting topology in $(prcnet*100) percent of times")
            d = readTableCF(df)
            startnet=deepcopy(currT0)
            bestNet = optTopRuns!(startnet, M, Nfail, d, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0, runs, outgroup, "bestNet", true,seeds[1],probST)
        else
            write(logfile,"bestNet input: $(writeTopology(bestNet))\n to use in $(prcnet *100) percent of times as starting topology")
            println("bestNet input: $(writeTopology(bestNet))\n to use in $(prcnet *100) percent of times as starting topology")
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
        net = optTopRuns!(startnet, M, Nfail, newd, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0, runs, outgroup, string("bootsnaq",i), true,seeds[i+1],probST)
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








