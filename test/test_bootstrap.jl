# test the functions in src/bootstrap.jl

exdir = joinpath(@__DIR__,"..","examples")
# exdir = joinpath(dirname(pathof(PhyloNetworks)), "..","examples")

@testset "testing hybridclades_support" begin
bestnet = readTopology(joinpath(exdir,"fish2hyb.net"));
bootnet = readmultitopology(joinpath(exdir,"fish3hyb_20boostrap.net"));
# issues with bootstrap networks 12, 21, 42, 96
# plot(bootnet[20], showedgenumber=true)
# include(string(home, "bootstrap.jl"))
resn, rese, resc, gam, edgenum = hybridclades_support(bootnet,bestnet);
#@show resn; @show rese; showall(gam); @show edgenum; resc
# plot(bestnet, shownodenumber=true);

@test resn[1:2,:clade] == ["H26","H25"]
@test resn[1:2,:BS_hybrid_samesisters] == [25.0,100.0]
@test resn[!,:BS_hybrid] == [100.0,100,0,0,0,75,0,0,0,0,0,0,5,5,5,5,5,0,0,0]
@test resn[!,:BS_minor_sister] == [0.0,0,100,0,0,5,10,70,75,25,5,5,0,0,0,0,0,0,0,5]
@test resn[!,:BS_major_sister] == [0.0,0,0,100,100,0,70,10,0,0,5,5,0,0,0,0,0,5,5,0]
@test rese[2,:BS_minor] == 25.0  # BS of introgression for H26
@test rese[4,:BS_minor] == 100.0 # BS of introgression for H25
@test resc[!,:taxa]==["Xgordoni","Xmeyeri","Xcouchianus","Xvariatus","Xevelynae","Xxiphidium",
 "Xmilleri","Xandersi","Xmaculatus","Xhellerii","Xalvarezi","Xmayae","Xsignum","Xclemenciae_F2",
 "Xmonticolus","Xmontezumae","Xnezahuacoyotl","Xbirchmanni_GARC","Xmalinche_CHIC2","Xcortezi",
 "Xcontinens","Xpygmaeus","Xnigrensis","Xmultilineatus"]
@test resc[!,:taxa][resc[!,:H26]] == ["Xnezahuacoyotl"]
@test resc[!,:taxa][resc[!,:H25]] == ["Xmontezumae","Xnezahuacoyotl","Xbirchmanni_GARC","Xmalinche_CHIC2","Xcortezi","Xcontinens","Xpygmaeus","Xnigrensis","Xmultilineatus"]
@test resc[!,:taxa][resc[!,:c_minus27]] == ["Xnigrensis","Xmultilineatus"] # minor sis of H26
@test resc[!,:taxa][resc[!,:Xxiphidium]] == ["Xxiphidium"] # minor sis of H25
@test resc[!,:taxa][resc[!,:Xsignum]] == ["Xsignum"] # donor8 previously
@test resc[!,:taxa][resc[!,:c_minus24]] == ["Xcontinens","Xpygmaeus","Xnigrensis","Xmultilineatus"] # donor7
@test resc[!,:taxa][resc[!,:Xmontezumae]] == ["Xmontezumae"] # major sis of H26. Below: major sis of H25
@test resc[!,:taxa][resc[!,:c_minus12]] == ["Xhellerii","Xalvarezi","Xmayae","Xsignum","Xclemenciae_F2","Xmonticolus"]
@test gam[:,2] == [.0,.0,.192,.0,.0,.0,.0,.0,.193,.0,.184,.193,.0,.0,.0,.0,.0,.193,.0,.0]
@test gam[:,4] == [.165,.166,.165,.166,.165,.165,.166,.165,.165,.166,.164,.166,.166,.165,.165,.165,.166,.165,.166,.166]
@test edgenum ==[25,39,43,7]

PhyloNetworks.addAlternativeHybridizations!(bestnet, rese, cutoff=4)
@test rese[[5,10,11],:edge] == [54,57,60]
@test ismissing(rese[6,:edge])
@test length(bestnet.hybrid) == 5
end # of testset, hybridclades_support

