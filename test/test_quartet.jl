PN = PhyloNetworks

@testset "QuartetT basic tests" begin
# basic tests for QuartetT type
nCk = PN.nchoose1234(5)
oneQ = PN.QuartetT(1,3,4,6, [.92,.04,.04, 100], nCk)
@test string(oneQ) == "4-taxon set number 8; taxon numbers: 1,3,4,6\ndata: [0.92, 0.04, 0.04, 100.0]"
end

@testset "convert QuartetT vector to table" begin
SVector = PN.StaticArrays.SVector
MVector = PN.StaticArrays.MVector
MMatrix = PN.StaticArrays.MMatrix
# without taxon names, length-3 vector of boolean data
qlist = [
  PN.QuartetT(1, SVector{4}(1,2,3,4), MVector{3,Bool}(0,1,0)),
  PN.QuartetT(2, SVector{4}(1,2,3,5), MVector{3,Bool}(0,0,1)),
]
@test tablequartetCF(qlist) == (qind=[1,2],t1=[1,1],t2=[2,2],t3=[3,3],
    t4=[4,5],CF12_34=[false,false],CF13_24=[true,false],CF14_23=[false,true])
# with taxon names
@test tablequartetCF(qlist, ["a","b","c","d","e"]) == (
  qind=[1,2],t1=["a","a"],t2=["b","b"],t3=["c","c"],t4=["d","e"],
  CF12_34=[false,false],CF13_24=[true,false],CF14_23=[false,true])
# 3x2 data, Int8-valued
qlist = [PN.QuartetT(1, SVector{4}(1,2,3,4), MMatrix{3,2,Int8}(1:6))]
@test tablequartetCF(qlist) == (qind=[1],t1=[1],t2=[2],t3=[3],t4=[4],
    CF12_34=Int8[1], CF13_24=Int8[2], CF14_23=Int8[3],
    V2_12_34=Int8[4],V2_13_24=Int8[5],V2_14_23=Int8[6])
# custom column names
@test_logs (:error, r"^'colnames'") tablequartetCF(qlist, colnames=["v1","v2"])
@test tablequartetCF(qlist, colnames=["v1","v2","v3","w1","w2","w3"]) ==
  (qind=[1],t1=[1],t2=[2],t3=[3],t4=[4], v1=Int8[1],v2=Int8[2],v3=Int8[3],
    w1=Int8[4],w2=Int8[5],w3=Int8[6])
# quartet data: length-4 vectors, second 4-taxon set with 0 genes
qlist = [
  PN.QuartetT(1, SVector{4}(1,2,3,4), MVector{4,Int8}(0,1,0,10)),
  PN.QuartetT(2, SVector{4}(1,2,3,5), MVector{4,Int8}(0,0,1, 0)),
  PN.QuartetT(3, SVector{4}(1,2,4,5), MVector{4,Int8}(1,0,0, 5)),
]
@test tablequartetCF(qlist; keepQwithoutgenes=false) ==
  (qind=[1,3],t1=[1,1],t2=[2,2],t3=[3,4], t4=[4,5],
  CF12_34=[0,1],CF13_24=[1,0],CF14_23=[0,0], ngenes=[10,5])
end

if false # was used to time countquartetsintrees vs readTrees2CF under PN ≤ v0.16
# under PN v1: countquartetsintrees returns a NamedTuple, not a DataFrame,
#              and readTrees2CF was moved to SNaQ.
dir = "/Users/ane/Documents/private/concordance/quartetNetwork/multiind/data"
treefile = joinpath(dir, "raxml_1387_sample_5species4alleles.tre")
tree = readmultinewick(treefile); # 1387 trees
# extrema([t.numtaxa for t in tree]) # 4-16 taxa in each
@time df1 = tablequartetCF(countquartetsintrees(tree)...)
# 0.139761 seconds (900.12 k allocations: 52.000 MiB, 11.52% gc time). 3876×8 DataFrames.DataFrame
@time df2 = tablequartetCF(readTrees2CF(tree, writeTab=false, writeSummary=false))
# 13.154085 seconds (84.86 M allocations: 10.010 GiB, 7.21% gc time).  3876×8 DataFrames.DataFrame
df12 = innerjoin(df1, df2, on=[:t1,:t2,:t3,:t4], makeunique=true)
@test all([df12[:,4+i] == df12[:,8+i] for i in 1:4])
# using BenchmarkTools
# @benchmark countquartetsintrees(tree)
#=
memory estimate:  49.08 MiB
allocs estimate:  788504
--------------
minimum time:     110.907 ms (10.89% GC)
median time:      121.782 ms (11.20% GC)
mean time:        125.725 ms (14.00% GC)
maximum time:     227.248 ms (54.16% GC)
--------------
samples:          40
evals/sample:     1
=#
# @benchmark readTrees2CF(tree, writeTab=false, writeSummary=false)
#=
BenchmarkTools.Trial:
  memory estimate:  10.01 GiB
  allocs estimate:  84744840
  --------------
  minimum time:     13.547 s (8.63% GC)
  median time:      13.547 s (8.63% GC)
  mean time:        13.547 s (8.63% GC)
  maximum time:     13.547 s (8.63% GC)
  --------------
  samples:          1
  evals/sample:     1
=#
mappingfile = joinpath(dir, "strain2bin_map.csv")

taxonmap = DataFrame(CSV.File(mappingfile); copycols=false) # 110×3 DataFrames.DataFrame
taxonmap = Dict(taxonmap[i,:allele] => taxonmap[i,:species] for i in 1:110)
@time df1 = tablequartetCF(countquartetsintrees(tree, taxonmap; weight_byallele=true)...)
# 0.119289 seconds (698.57 k allocations: 43.305 MiB, 17.40% gc time). 5×8 DataFrames.DataFrame
## larger examples: 98 to 110 taxa, 1387 trees
tree = readmultinewick(joinpath(dir, "raxml_1387.tre")) # 1387 trees, 98-110 taxa in each
@time df1 = tablequartetCF(countquartetsintrees(tree)...)
# 1219.94 seconds = 20.3 min (298.45 M allocations: 8.509 GiB, 0.79% gc time). 5773185×8 DataFrames.DataFrame
## mid-size example: to be able to run the slower algorithm and compare times
tree = readmultinewick(joinpath(dir, "raxml_1387_sample_13species4alleles.tre")); # 1387 trees, 19-40 taxa in each
@time df1 = tablequartetCF(countquartetsintrees(tree)...)
# 5.639443 seconds (16.84 M allocations: 568.496 MiB, 7.50% gc time). 292825×8 DataFrame
# ~ 600 times faster
@time df2 = tablequartetCF(readTrees2CF(tree, writeTab=false, writeSummary=false))
# 3365.783672 seconds = 50.1 min (13.43 G allocations: 2.665 TiB, 21.33% gc time). 292825×8 DataFrame
df12 = innerjoin(df1, df2, on=[:t1,:t2,:t3,:t4], makeunique=true)
@test df12[!,8] ≈ df12[!,12] # number of genes
hasdata = map(iszero, df12[!,8]) # sum: 34 four-taxon sets have data for 0 genes
df12[hasdata,5:12] # countquartetsintrees gives 0s, readTrees2CF gives NaN
@test all([df12[.!hasdata,4+i] ≈ df12[.!hasdata,8+i] for i in 1:3]) # true. yeah!
end

@testset "countquartetsintrees" begin
sixtreestr = ["(E,((A,B),(C,D)),O);","(((A,B),(C,D)),(E,O));","(A,B,((C,D),(E,O)));",
              "(B,((C,D),(E,O)));","((C,D),(A,(B,E)),O);","((C,D),(A,B,E),O);"]
sixtrees = readnewick.(sixtreestr)
nt = tablequartetCF(countquartetsintrees(sixtrees)...)
@test nt == (qind = 1:15,
  t1 = ["A","A","A","A","B","A","A","A","B","A","A","B","A","B","C"],
  t2 = ["B","B","B","C","C","B","B","C","C","B","C","C","D","D","D"],
  t3 = ["C","C","D","D","D","C","D","D","D","E","E","E","E","E","E"],
  t4 = ["D","E","E","E","E","O","O","O","O","O","O","O","O","O","O"],
  CF12_34 = [1,.75,.75, 0, 0, 1, 1, 0, 0,.75,.6,2/3,.6, 2/3, 1],
  CF13_24 = [0,.25,.25, 0, 0, 0, 0, 0, 0, 0, .4,1/3,.4, 1/3, 0],
  CF14_23 = [0,  0,  0, 1, 1, 0, 0, 1, 1,.25, 0,  0, 0, 0, 0],
  ngenes = [5, 4, 4, 5, 6, 5, 5, 5, 6, 4, 5, 6, 5, 6, 6])
o = [1,2,4,7,11,3,5,8,12,6,9,13,10,14,15]
q,t = countquartetsintrees(sixtrees, Dict("A"=>"AB", "B"=>"AB"); showprogressbar=false);
nt = tablequartetCF(q,t)
@test nt[:CF12_34] ≈ [0,0,2/3,2/3,1]
@test nt[:CF13_24] ≈ [0,0,1/3,1/3,0]
@test nt[:CF14_23] ≈ [1.,1,0,0,0]
@test nt[:ngenes]  ≈ [6.,6,6,6,6]
# again, but weight each allele
q,t = countquartetsintrees(sixtrees, Dict("A"=>"AB", "B"=>"AB"); weight_byallele=true, showprogressbar=false);
nt = tablequartetCF(q,t)
@test nt[:CF12_34] ≈ [0,0,7/11,7/11,1]
@test nt[:CF13_24] ≈ [0,0,4/11,4/11,0]
@test nt[:CF14_23] ≈ [1.,1,0,0,0]
@test nt[:ngenes]  ≈ [11.,11,11,11,6]
end
