## Claudia: commented these tests out because they were not part of runtests.jl
## when I added the test for readNexusTrees. So, I was not sure if they were
## not included for a reason (slow test?)

## @testset "test: read CF data" begin
## # d=readInputData("1.ms"); #tableCF0.txt
## d = (@test_logs PhyloNetworks.readInputData("1.ms",:all,0,["1","2","3","4","5","6"],false,"none",false,false);) #tableCF0.txt
## d = PhyloNetworks.readInputData("1.ms",:all,0,false); # writes summaryTreesQuartets.txt
## d=readInputData("1.ms",:rand,10); #tableCF3.txt
## d=readInputData("1.ms",[1,2,3,4,5]); #tableCF4.txt

## d=readInputData("1.ms","allQuartets.txt"); #tableCF1.txt
## d=readInputData("1.ms","allQuartets.txt",:rand,10); #tableCF2.txt
## d=readInputData("1.ms","allQuartets.txt",true,"try4.txt"); #try4.txt

## descData(d)


## d = readTrees2CF("1.ms");
## d = readTrees2CF("1.ms",filename="try.txt");
## d = readTrees2CF("1.ms","allQuartets.txt",filename="hola.txt");
## d = readTrees2CF("1.ms","allQuartets.txt",whichQ=:rand,numQ=10,filename="hola.txt");

## d= readTableCF("tableCF.txt");
## d2= readTableCF("tableCFbad.txt")
## end


@testset "test: reading nexus file" begin
nexusfile = joinpath(@__DIR__, "..", "examples", "test.nex")
# nexusfile = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","test.nex")
vnet = readNexusTrees(nexusfile);
@test length(vnet) == 10
@test length(vnet[10].edge) == 10
@test vnet[10].edge[7].length ≈ 0.00035
vnet = readNexusTrees(nexusfile, PhyloNetworks.readTopologyUpdate, false, false);
@test length(vnet) == 10
@test length(vnet[10].edge) == 9
@test vnet[10].edge[7].length ≈ 0.00035
end

@testset "test: calculate quartet CF from input gene trees" begin
sixtreestr = ["(E,((A,B),(C,D)),O);","(((A,B),(C,D)),(E,O));","(A,B,((C,D),(E,O)));",
              "(B,((C,D),(E,O)));","((C,D),(A,(B,E)),O);","((C,D),(A,B,E),O);"]
sixtrees = readTopology.(sixtreestr)
df1 = writeTableCF(countquartetsintrees(sixtrees)...)
df2 = writeTableCF(readTrees2CF(sixtrees, writeTab=false, writeSummary=false))
o = [1,2,4,7,11,3,5,8,12,6,9,13,10,14,15]
@test df1 == df2[o,:]
q,t = countquartetsintrees(sixtrees, Dict("A"=>"AB", "B"=>"AB"); showprogressbar=false);
df1 = writeTableCF(q,t)
@test df1[!,:CF12_34] ≈ [0,0,2/3,2/3,1]
@test df1[!,:CF13_24] ≈ [0,0,1/3,1/3,0]
@test df1[!,:CF14_23] ≈ [1.,1,0,0,0]
@test df1[!,:ngenes]  ≈ [6.,6,6,6,6]
# again, but weight each allele
q,t = countquartetsintrees(sixtrees, Dict("A"=>"AB", "B"=>"AB"); weight_byallele=true, showprogressbar=false);
df1 = writeTableCF(q,t)
@test df1[!,:CF12_34] ≈ [0,0,7/11,7/11,1]
@test df1[!,:CF13_24] ≈ [0,0,4/11,4/11,0]
@test df1[!,:CF14_23] ≈ [1.,1,0,0,0]
@test df1[!,:ngenes]  ≈ [11.,11,11,11,6]
# different averaging: first across each set of 4 alleles, then
# across sets of 4 alleles that map to the same set of 4 species
CSV.write("tmp_qCF.csv", df2)
CSV.write("tmp_map.csv", DataFrame(allele = ["A","B"], species = ["AB","AB"]))
df2_byallele = (@test_logs (:warn, r"not all alleles were mapped") mapAllelesCFtable("tmp_map.csv", "tmp_qCF.csv"))
rm.(["tmp_qCF.csv","tmp_map.csv"]);
q = readTableCF!(df2_byallele);
df2 = writeTableCF(q) # 45×8 DataFrames.DataFrame
# df2[7:11,:]
# df12 = join(df1, df2, on=[:t1,:t2,:t3,:t4], makeunique=true)
# all([df12[:,4+i] ≈ df12[:,8+i] for i in 1:4]) # false: because averaging done differently by the 2 functions
@test df2 == DataFrame(
  t1=["AB","AB","AB","AB","AB","AB","AB","AB","AB","AB","C"],
  t2=["AB__2","AB__2","AB__2","AB__2","AB__2","AB__2","C","C","C","D","D"],
  t3=["C","C","C","D","D","E","D","D","E","E","E"],
  t4=["D","E","O","E","O","O","E","O","O","O","O"],
  CF12_34=[1.0,0.75,1.0,0.75,1.0,0.75,0.0,0.0,(3/5+4/6)/2,(3/5+4/6)/2,1.0],
  CF13_24=[0,0.25,0,0.25,0,0,0,0,(2/5+2/6)/2,(2/5+2/6)/2,0.0],
  CF14_23=[0,0,0,0,0,0.25,1,1,0,0,0],
  ngenes=Union{Missing,Float64}[5,4,5,4,5,4,5.5,5.5,5.5,5.5,6]
)
end

if false # was used to time `countquartetsintrees` vs `readTrees2CF`
dir = "/Users/ane/Documents/private/concordance/quartetNetwork/multiind/data"
treefile = joinpath(dir, "raxml_1387_sample_5species4alleles.tre")
tree = readMultiTopology(treefile); # 1387 trees
# extrema([t.numTaxa for t in tree]) # 4-16 taxa in each
@time df1 = writeTableCF(countquartetsintrees(tree)...)
# 0.139761 seconds (900.12 k allocations: 52.000 MiB, 11.52% gc time). 3876×8 DataFrames.DataFrame
@time df2 = writeTableCF(readTrees2CF(tree, writeTab=false, writeSummary=false))
# 13.154085 seconds (84.86 M allocations: 10.010 GiB, 7.21% gc time).  3876×8 DataFrames.DataFrame
df12 = join(df1, df2, on=[:t1,:t2,:t3,:t4], makeunique=true)
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
using CSV
taxonmap = CSV.read(mappingfile) # 110×3 DataFrames.DataFrame
taxonmap = Dict(taxonmap[i,:allele] => taxonmap[i,:species] for i in 1:110)
@time df1 = writeTableCF(countquartetsintrees(tree, taxonmap; weight_byallele=true)...)
# 0.119289 seconds (698.57 k allocations: 43.305 MiB, 17.40% gc time). 5×8 DataFrames.DataFrame
## larger examples: 98 to 110 taxa, 1387 trees
tree = readMultiTopology(joinpath(dir, "raxml_1387.tre")) # 1387 trees, 98-110 taxa in each
@time df1 = writeTableCF(countquartetsintrees(tree)...)
# 1219.94 seconds = 20.3 min (298.45 M allocations: 8.509 GiB, 0.79% gc time). 5773185×8 DataFrames.DataFrame
## mid-size example: to be able to run the slower algorithm and compare times
tree = readMultiTopology(joinpath(dir, "raxml_1387_sample_13species4alleles.tre")); # 1387 trees, 19-40 taxa in each
@time df1 = writeTableCF(countquartetsintrees(tree)...)
# 5.639443 seconds (16.84 M allocations: 568.496 MiB, 7.50% gc time). 292825×8 DataFrame
# ~ 600 times faster
@time df2 = writeTableCF(readTrees2CF(tree, writeTab=false, writeSummary=false))
# 3365.783672 seconds = 50.1 min (13.43 G allocations: 2.665 TiB, 21.33% gc time). 292825×8 DataFrame
df12 = join(df1, df2, on=[:t1,:t2,:t3,:t4], makeunique=true)
@test df12[!,8] ≈ df12[!,12] # number of genes
hasdata = map(iszero, df12[!,8]) # sum: 34 four-taxon sets have data for 0 genes
df12[hasdata,5:12] # countquartetsintrees gives 0s, readTrees2CF gives NaN
@test all([df12[.!hasdata,4+i] ≈ df12[.!hasdata,8+i] for i in 1:3]) # true. yeah!
end
