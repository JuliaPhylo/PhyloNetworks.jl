@testset "testing TICR" begin
global df

@testset "ticr! on data frame, on tree" begin
truenet1 = readTopology("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
setGamma!(truenet1.edge[9],0.0);
truenet2 = deepcopy(truenet1);
df = CSV.read(joinpath(@__DIR__,"..","examples","buckyCF.csv"));
#df = CSV.read(joinpath(dirname(pathof(PhyloNetworks)), "..","examples","buckyCF.csv"));
#without optimizing branch lengthes
result1 = ticr!(truenet1,df,false);
@test result1[2] ≈ 25.962962962962965463 # chi-squared statistic obtained from R
@test result1[1] ≈ 9.7092282251534852702e-06 # p-value obtained from R
@test result1[3] ≈ 48.152697007372566418 # pseudo log-lik obtained from R
@test result1[4] ≈ 10.576940922426542713 atol=1e-5 # alpha obtained from R
#with optimizing branch lengthes
result2 = ticr!(truenet2,df,true);
setGamma!(result2[6].edge[9],0.0);
# fixit: hard-code the network result2[6] to get predictable branch lengths in the tests below
result3 = ticr!(result2[6],df,false);
@test result3[2] ≈ 25.962962962962965463 # chi-squared statistic obtained from R
@test result3[1] ≈ 9.7092282251534852702e-06 # p-value obtained from R
@test result3[3] ≈ 54.449883693197676848 atol=3e-1 # pseudo log-lik obtained from R
@test result3[4] ≈ 20.694991969052374259 atol=6e-1 # alpha obtained from R
end

@testset "ticr! on data frame, on network" begin
truenet3 = readTopology("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
#without optimizing branch lengths
netresult1 = ticr!(truenet3,df,false);
#setGamma!(netresult1[6].edge[9],0.0);
#netresult2 = ticr!(netresult1[6],netCF1,false);
@test netresult1[2] ≈ 2.851851851851852 # chi-squared statistic
@test netresult1[1] ≈ 0.41503515532593677 # p-value
@test netresult1[3] ≈ 68.03708830981597 # pseudo log-lik
@test netresult1[4] ≈ 29.34808731515701 atol=1e-5 # alpha
end

end
