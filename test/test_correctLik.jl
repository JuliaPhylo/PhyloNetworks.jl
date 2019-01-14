# test to see if the likelihood is correctly calculated
# and if the networks are correctly estimated
# Claudia August 2015

# -------------------5taxon tree------------------

PhyloNetworks.CHECKNET || error("need CHECKNET==true in PhyloNetworks to test snaq in test_correctLik.jl")

#df = readtable("Tree_output.txt")
df=DataFrame(t1=["6","6","10","6","6"],
             t2=["7","7","7","10","7"],
             t3=["4","10","4","4","4"],
             t4=["8","8","8","8","10"],
             obsCF12=[0.2729102510259939, 0.3967750546426937, 0.30161247267865315, 0.24693940689390592, 0.2729102510259939],
             obsCF13=[0.45417949794801216, 0.30161247267865315, 0.30161247267865315, 0.5061211862121882, 0.45417949794801216],
             obsCF14=[0.2729102510259939, 0.30161247267865315, 0.3967750546426937, 0.24693940689390592, 0.2729102510259939])
d = readTableCF(df)
@test_throws ErrorException PhyloNetworks.writeExpCF(d)
@test writeTableCF(d) == rename(df, [:obsCF12 => :CF12_34, :obsCF13 => :CF13_24, :obsCF14 => :CF14_23])
@test tipLabels(d) ==  ["4","6","7","8","10"]
@test_logs PhyloNetworks.descData(d, devnull)

df[:ngenes] = [10,10,10,10,20]
allowmissing!(df, :ngenes)
d = readTableCF(df)
df[:ngenes][1] = missing; d.quartet[1].ngenes = -1.0
newdf = writeTableCF(d)
@test newdf[1:7] == rename(df, [:obsCF12 => :CF12_34, :obsCF13 => :CF13_24, :obsCF14 => :CF14_23])[1:7]
@test ismissing(newdf[:ngenes][1])
@test newdf[:ngenes][2:end] == df[:ngenes][2:end]

# starting tree:
tree = "((6,4),(7,8),10);"
currT = readTopologyLevel1(tree);
#printEdges(currT)

@testset "correct pseudo likelihood and snaq" begin
@testset "lik on tree" begin
extractQuartet!(currT,d)
calculateExpCFAll!(d)
tmp = (@test_logs PhyloNetworks.writeExpCF(d))
for i in [5,7] for j in 2:5 @test tmp[i][j] ≈ 0.12262648039048077; end end
for j in 2:5 @test tmp[6][j] ≈ 0.7547470392190385; end
lik = logPseudoLik(d)
@test lik ≈ 193.7812623319291
#estTree = optTopRun1!(currT,d,0,5454) # issue with printCounts, TravisCI?
#@test estTree.loglik ≈ 0.0 atol=1e-8
#println("passed optTopRun1! on tree")
end

# ------------------5taxon network 1 hybridization: Case H-----------------
# starting topology: Case G
global tree = "((((6:0.1,4:1.5)1:0.2,(7)11#H1)5:0.1,(11#H1,8)),10:0.1);" # Case G
global currT = readTopologyLevel1(tree);
# real network: Case H
global df=DataFrame(t1=["6","6","10","6","6"],t2=["7","7","7","10","7"],t3=["4","10","4","4","4"],t4=["8","8","8","8","10"],CF1234=[0.13002257237728915, 0.36936019721747243, 0.34692592933269173, 0.12051951084152591, 0.11095702789935982], CF1324=[0.7399548552454217, 0.28371387344983595, 0.28371387344983595, 0.7589609783169482, 0.7780859442012804],CF1423=[0.13002257237728915, 0.34692592933269173, 0.36936019721747243, 0.12051951084152591, 0.11095702789935982])
global d = readTableCF(df)

@testset "lik of network" begin
extractQuartet!(currT,d)
calculateExpCFAll!(d)
lik = logPseudoLik(d)
@test lik ≈ 50.17161079450669
end

@testset "network estimation h=1" begin
estNet = optTopRun1!(currT, 0.01,75, d,1, 1e-5,1e-6,1e-3,1e-4,
                     false,true,Int[], 5454, stdout,false,0.3)
# topology, likAbs,Nfail, data,hmax, fRel,fAbs,xRel,xAbs,
# verbose,closeN,numMoves, seed, logfile,writelog,probST,sout)
@test estNet.loglik ≈ 0.002165 atol=5.0e-6
# 0.00216 < estNet.loglik < 0.00217 || Base.error("not correct estimated network")
end

@testset "snaq! in serial and in parallel" begin
  global tree = readTopology("((((6:0.1,4:1.5),9)1:0.1,8),10:0.1);")
  @test_throws ErrorException snaq!(tree, d) # some taxa are in quartets, not in tree
  originalstdout = stdout
  redirect_stdout(open("/dev/null", "w")) # not portable to Windows
  global net = readTopology("((((6:0.1,4:1.5)1:0.2,((7,60))11#H1)5:0.1,(11#H1,8)),10:0.1);")
  @test_logs (:warn, r"^these taxa will be deleted") snaq!(net, d, # taxon "60" in net: not in quartets
    hmax=1, runs=1, Nfail=1, seed=1234, ftolRel=1e-2,ftolAbs=1e-2,xtolAbs=1e-2,xtolRel=1e-2)
  global n1 = snaq!(currT, d, hmax=1, runs=1, Nfail=1, seed=1234,
             ftolRel=1e-2,ftolAbs=1e-2,xtolAbs=1e-2,xtolRel=1e-2,
             verbose=true)
  addprocs(1)
  @everywhere using PhyloNetworks
  global n2 = snaq!(currT, d, hmax=1, runs=2, Nfail=1, seed=1234,
             ftolRel=1e-2,ftolAbs=1e-2,xtolAbs=1e-2,xtolRel=1e-2)
  redirect_stdout(originalstdout)
  rmprocs(workers())
  @test writeTopology(n1, round=true)==writeTopology(n2, round=true)
  @test n1.loglik == n2.loglik
  rm("snaq.out")
  rm("snaq.networks")
  rm("snaq.log") # .log and .err should be git-ignored, but still
  rm("snaq.err")
end
end
