# test the functions in src/bootstrap.jl
# assumes types.jl and functions.jl have been included already, or using PhyloNetworks

if !isdefined(:doalltests) doalltests = false; end
scriptfile = @__FILE__
exdir = (scriptfile==nothing ? joinpath("..","examples") :
        joinpath(dirname(dirname(@__FILE__)), "examples"))

info("Testing hybridBootstrapFrequency")
bestnet = readTopology(joinpath(exdir,"fish2hyb.net"));
bootnet = readMultiTopology(joinpath(exdir,"fish3hyb_20boostrap.net"));
# issues with bootstrap networks 12, 21, 42, 96
# plot(bootnet[21], showEdgeNumber=true, showEdgeLength=false)
# include(string(home, "bootstrap.jl"))
f, fr, fd, fs, clade, gam, edgenum = hybridBootstrapFrequency(bootnet,bestnet);
@test f[:hybrid] == ["H26","H25"]
@test f[:prop_introgression] == [0.25,1.0]
@test f[:prop_hybridization] == [0.25,1.0]
@test fr[:proportion] == [1.0,1.0,0.75,0.05,0.05,0.05,0.05,0.05]
@test fd[:proportion] == [.25,1.0,0.75,0.7,0.05,0.1,0.05,0.05,0.05]
@test fs[:proportion] == [1.0,1.0,.7,.05,.1,.05,.05,.05]
# show(f); show(fr); show(fd); show(fs)
@test clade[:taxa]==["Xgordoni","Xmeyeri","Xcouchianus","Xvariatus","Xevelynae","Xxiphidium",
 "Xmilleri","Xandersi","Xmaculatus","Xhellerii","Xalvarezi","Xmayae","Xsignum","Xclemenciae_F2",
 "Xmonticolus","Xmontezumae","Xnezahuacoyotl","Xbirchmanni_GARC","Xmalinche_CHIC2","Xcortezi",
 "Xcontinens","Xpygmaeus","Xnigrensis","Xmultilineatus"]
@test clade[:recipientH26]==[rep(false,16);true;rep(false,7)]
@test clade[:recipientH25]==[rep(false,15);rep(true,9)]
@test clade[:donorH26]==[rep(false,22);rep(true,2)]
@test clade[:donorH25]==[rep(false,5);true;rep(false,18)]
@test clade[:donor7]  ==[rep(false,20);rep(true,4)]
@test clade[:donor8]  ==[rep(false,12);true;rep(false,11)]
@test clade[:siblingH26]==[rep(false,15);true;rep(false,8)]
@test clade[:siblingH25]==[rep(false,9);rep(true,6);rep(false,9)]
@test gam[:,1] == [.0,.0,.192,.0,.0,.0,.0,.0,.193,.0,.184,.193,.0,.0,.0,.0,.0,.193,.0,.0]
@test gam[:,2] == [.165,.166,.165,.166,.165,.165,.166,.165,.165,.166,.164,.166,.166,.165,.165,.165,.166,.165,.166,.166]
@test edgenum ==[39,7]