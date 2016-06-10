# test the functions in src/bootstrap.jl
# assumes types.jl and functions.jl have been included already, or using PhyloNetworks

if !isdefined(:doalltests) doalltests = false; end
scriptfile = @__FILE__
exdir = (scriptfile==nothing ? joinpath("..","examples") :
        joinpath(dirname(dirname(@__FILE__)), "examples"))

info("Testing hybridBootstrapSupport")
bestnet = readTopology(joinpath(exdir,"fish2hyb.net"));
bootnet = readMultiTopology(joinpath(exdir,"fish3hyb_20boostrap.net"));
# issues with bootstrap networks 12, 21, 42, 96
# plot(bootnet[21], showEdgeNumber=true, showEdgeLength=false)
# include(string(home, "bootstrap.jl"))
resn, rese, resc, gam, edgenum = hybridBootstrapSupport(bootnet,bestnet);
#@show resn; @show rese; showall(gam); @show edgenum

@test resn[:clade][1:2] == ["H26","H25"]
@test resn[:BS_hybrid_samesisters][1:2] == [25.0,100.0]
@test resn[:BS_hybrid] == [100.0,100,0,0,0,75,0,0,0,0,0,0,5,5,5,5,5,0,0,0]
@test resn[:BS_minor_sister] == [0.0,0,100,0,0,5,10,70,75,25,5,5,0,0,0,0,0,0,0,5]
@test resn[:BS_major_sister] == [0.0,0,0,100,100,0,70,10,0,0,5,5,0,0,0,0,0,5,5,0]
@test rese[:BS_minor][2] == 25.0  # BS of introgression for H26
@test rese[:BS_minor][4] == 100.0 # BS of introgression for H25
@test resc[:taxa]==["Xgordoni","Xmeyeri","Xcouchianus","Xvariatus","Xevelynae","Xxiphidium",
 "Xmilleri","Xandersi","Xmaculatus","Xhellerii","Xalvarezi","Xmayae","Xsignum","Xclemenciae_F2",
 "Xmonticolus","Xmontezumae","Xnezahuacoyotl","Xbirchmanni_GARC","Xmalinche_CHIC2","Xcortezi",
 "Xcontinens","Xpygmaeus","Xnigrensis","Xmultilineatus"]
@test resc[:taxa][resc[:H26]] == ["Xnezahuacoyotl"]
@test resc[:taxa][resc[:H25]] == ["Xmontezumae","Xnezahuacoyotl","Xbirchmanni_GARC","Xmalinche_CHIC2","Xcortezi","Xcontinens","Xpygmaeus","Xnigrensis","Xmultilineatus"]
@test resc[:taxa][resc[:c_minus26]] == ["Xnigrensis","Xmultilineatus"] # minor sis of H26
@test resc[:taxa][resc[:Xxiphidium]] == ["Xxiphidium"] # minor sis of H25
@test resc[:taxa][resc[:Xsignum]] == ["Xsignum"] # donor8 previously
@test resc[:taxa][resc[:c_minus23]] == ["Xcontinens","Xpygmaeus","Xnigrensis","Xmultilineatus"] # donor7
@test resc[:taxa][resc[:Xmontezumae]] == ["Xmontezumae"] # major sis of H26. Below: major sis of H25
@test resc[:taxa][resc[:c_minus11]] == ["Xhellerii","Xalvarezi","Xmayae","Xsignum","Xclemenciae_F2","Xmonticolus"]
@test gam[:,2] == [.0,.0,.192,.0,.0,.0,.0,.0,.193,.0,.184,.193,.0,.0,.0,.0,.0,.193,.0,.0]
@test gam[:,4] == [.165,.166,.165,.166,.165,.165,.166,.165,.165,.166,.164,.166,.166,.165,.165,.165,.166,.165,.166,.166]
@test edgenum ==[25,39,43,7]
