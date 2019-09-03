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
