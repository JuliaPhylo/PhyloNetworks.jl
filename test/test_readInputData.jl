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
df1 = writeTableCF(observedquartetCF(sixtrees)...)
df2 = writeTableCF(readTrees2CF(sixtrees, writeTab=false, writeSummary=false))
o = [1,2,4,7,11,3,5,8,12,6,9,13,10,14,15]
@test df1 == df2[o,:]
end

if false # code used to time the different implementations
dir = "/Users/ane/Documents/private/concordance/quartetNetwork/multiind/data"
treefile1 = joinpath(dir, "raxml_1387_sample_5species4alleles.tre")
tree1 = readMultiTopology(treefile1); # 1387 trees
# extrema([t.numTaxa for t in tree1]) # 4-16 taxa in each
@time df1 = writeTableCF(observedquartetCF(tree1)...)
# 0.139761 seconds (900.12 k allocations: 52.000 MiB, 11.52% gc time). 3876×8 DataFrames.DataFrame
@time df2 = writeTableCF(readTrees2CF(tree1, writeTab=false, writeSummary=false))
# 13.154085 seconds (84.86 M allocations: 10.010 GiB, 7.21% gc time).  3876×8 DataFrames.DataFrame
df12 = join(df1, df2, on=[:t1,:t2,:t3,:t4], makeunique=true)
@test all([df12[:,4+i] == df12[:,8+i] for i in 1:4])
# using BenchmarkTools
# @benchmark observedquartetCF(tree1)
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
# @benchmark readTrees2CF(tree1, writeTab=false, writeSummary=false)
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
@time df11 = writeTableCF(observedquartetCF(tree1, :all, taxonmap)...)
# 0.119289 seconds (698.57 k allocations: 43.305 MiB, 17.40% gc time). 5×8 DataFrames.DataFrame
CSV.write("qCF_tmp.csv", df2)
@time df22_byallele = mapAllelesCFtable(mappingfile, "qCF_tmp.csv")
# 0.042906 seconds (63.41 k allocations: 10.254 MiB, 15.92% gc time). 3876×8 DataFrames.DataFrame
@time q = readTableCF!(df22_byallele);
# 0.036276 seconds (202.10 k allocations: 30.620 MiB, 17.11% gc time)
df22 = writeTableCF(q) # 45×8 DataFrames.DataFrame
df12 = join(df11, df22, on=[:t1,:t2,:t3,:t4], makeunique=true)
@test all([df12[:,4+i] == df12[:,8+i] for i in 1:4]) # false: because of gene weight = 1 versus allele weight = 1
end
