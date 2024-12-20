@testset "calibrate with distances: 3-cycle, non-identifiable example" begin
net2 = readnewick("((((D:0.1,C:0.2):1.5,(B:0.1)#H1:0.9::0.7):0.1,(#H1:0.01::0.3,A:0.3):0.8):0.1);")
dAll = pairwisetaxondistancematrix(net2, keepInternal=true)
@test dAll ≈ [0.0 .8 1.1 .1 .943 1.043 1.6 1.8 1.7
 0.8   0.0   0.3   0.9   1.263 1.363 2.4   2.6   2.5
 1.1   0.3   0.0   1.2   1.563 1.663 2.7   2.9   2.8
 0.1   0.9   1.2   0.0   0.903 1.003 1.5   1.7   1.6
 0.943 1.263 1.563 0.903 0.0   0.1   2.403 2.603 2.503
 1.043 1.363 1.663 1.003 0.1   0.0   2.503 2.703 2.603
 1.6   2.4   2.7   1.5   2.403 2.503 0.0   0.2   0.1
 1.8   2.6   2.9   1.7   2.603 2.703 0.2   0.0   0.3
 1.7   2.5   2.8   1.6   2.503 2.603 0.1   0.3   0.0]
taxa = [l.name for l in net2.leaf]
net2distances = pairwisetaxondistancematrix(net2)
@test net2distances ≈ [.0 .3 2.603 2.8
 0.3   0.0   2.703 2.9
 2.603 2.703 0.0   1.663
 2.8   2.9   1.663 0.0]
g = PhyloNetworks.pairwisetaxondistance_gradient(net2);
@test size(g) == (9,9,9)
na = getnodeages(net2);
@test na ≈ [1.7,0.11,0.0,1.6,0.1,0.0,0.1,0.0,0.0]
g = PhyloNetworks.pairwisetaxondistance_gradient(net2, nodeAges=na);
@test size(g) == (9,9,9)
@test g[:,:,1] ≈ [0 1 1 1 1 1 1 1 1;1 0 0 2 1.4 1.4 2 2 2; 1 0 0 2 1.4 1.4 2 2 2
 1 2 2 0 .6 .6 0 0 0; 1 1.4 1.4 .6 0 0 .6 .6 .6; 1 1.4 1.4 .6 0 0 .6 .6 .6
 1 2 2 0 .6 .6 0 0 0; 1 2 2 0 .6 .6 0 0 0; 1 2 2 0 .6 .6 0 0 0]

global net
net  = readnewick("((#H1:0.06::0.3,A:0.6):1.3,(B:0.1)#H1:0.7::0.7,(C,D):1.4);");
# same topology, unrooted, different BL, leaves ordered differently
# here: branch lengths not identifiable, even if minor fixed to 0
calibratefrompairwisedistances!(net, net2distances, taxa, #verbose=true,
  forceMinorLength0=true, ultrametric=false)
est = pairwisetaxondistancematrix(net, checkpreorder=false)
@test est ≈ [0 1.663 2.9 2.8; 1.663 0 2.703 2.603; 2.9 2.703 0 .3;
             2.8 2.603 .3 0] atol=.00002
# deep    LN_BOBYQA: 146 evaluations
# got 0.0 at [0.00536, 1.32091, 0.14155, 0.84494, 0.2, 0.1, 1.37373] after 152 iterations (returned FTOL_REACHED)
# shallow LN_BOBYQA: 79 evaluations
# deep    LD_MMA: 47 evaluations
# got 0.0 at [0.20465, 1.03622, 0.19367, 0.77048, 0.2, 0.1, 1.45914] after 47 iterations (returned FTOL_REACHED)
# second call, starting from solution:
# got 0.0 at [0.00536, 1.32092, 0.14155, 0.84493, 0.2, 0.1, 1.37373] after 29 iterations (returned FTOL_REACHED)

net = readnewick("((#H1:0.06::0.3,A:0.6):1.3,(B:0.1)#H1:0.7::0.7,(C,D):1.4);");
calibratefrompairwisedistances!(net, net2distances, taxa,
  verbose=false, forceMinorLength0=false, ultrametric=false)
est = pairwisetaxondistancematrix(net)
@test est ≈ [0 1.663 2.9 2.8; 1.663 0 2.703 2.603; 2.9 2.703 0 .3;
             2.8 2.603 .3 0] atol=.00002
end

@testset "calibrate with distances: 5-cycle example (when unrooted)" begin
# most lengths identifiable when we force
# minor hybrid to length 0 and ultrametric network
# not identifiable otherwise
global net
net = readnewick("((Ag:3.0,(#H1:0.0::0.2,Ak:2.5):0.5):0.5,(((((Az:0.2,Ag2:0.2):1.3,As:1.5):1.0)#H1:0.3::0.8,Ap:2.8):0.2,Ar:3.0):0.5);");
taxa = [l.name for l in net.leaf];
netdist = pairwisetaxondistancematrix(net);
@test getnodeages(net) ≈ [3.5,3,0,2.8,0,3,2.5,0,2.5,1.5,0,.2,0,0,0]
net = readnewick("((((((Ag2,Az),As))#H1:::0.8,Ap),Ar),(Ag,(#H1:::0.2,Ak)));");
# same topology, no BL, leaves ordered differently
o = [4,3,5,6,7,1,2]; # to get leaves in same order as in taxa. inverse: [6,7,2,1,3,4,5]
calibratefrompairwisedistances!(net, netdist, taxa,
  forceMinorLength0=true, ultrametric=false)
# got 0.0 at [0.19998, 0.19998, 1.30003, 1.5, 0.6802, 0.69964, 2.79995, 0.20002, 2.99995, 0.50006, 2.9999, 2.49952, 0.50049, 0.50006] after 1000 iterations (returned MAXEVAL_REACHED)
est = pairwisetaxondistancematrix(net, checkpreorder=false)
@test est ≈ netdist[o,o] atol=.0005
# most edge lengths are correct, except for 5 edges:
# - the 2 edges at the root: expected (the network should be unrooted)
#   when the optimization uses LN_BOBYQA
# - the 3 edges connecting to the hybrid node: lack of identifiability?

net = readnewick("((((((Ag2,Az),As))#H1:::0.8,Ap),Ar),(Ag,(#H1:::0.2,Ak)));");
# same topology, no BL, leaf ordered differently
for e in net.edge e.length=1.0; end
preorder!(net)
na = getnodeages(net);
@test na ≈ [6.,1,4,0,0,5,0,4,0,3,2,0,1,0,0]
pairwisetaxondistancematrix(net, nodeAges=na); # update edge lengths in net
@test [e.length for e in net.edge] ≈ [1.,1,1,2,1,1,4,1,5,1,1,1,4,-3,5]
g = PhyloNetworks.pairwisetaxondistance_gradient(net, nodeAges=na);
net.edge[11].length = 4.; # to respect constraints
net.edge[14].length = 0.;
net.edge[15].length = 2.;
net1 = deepcopy(net);
fmin, na = calibratefrompairwisedistances!(net1, netdist, taxa, #verbose=true,
  forceMinorLength0=false, ultrametric=true)
# got 0.0 at [3.5, 3.0, 2.50001, 3.0, 2.8, 1.95797, 1.5, 0.2] after 58 iterations (returned FTOL_REACHED)
est = pairwisetaxondistancematrix(net1)
@test est ≈ netdist[o,o] atol=.0005
# still not identifiable: all good but the 3 edges adjacent to hybrid node
# only the age of the hybrid node is wrong: zipped down nearer tips

calibratefrompairwisedistances!(net, netdist, taxa, #verbose=true,
  forceMinorLength0=true, ultrametric=true, NLoptMethod=:LN_COBYLA)
est = pairwisetaxondistancematrix(net)
@test est ≈ netdist[o,o] atol=.06
# hard problem? reached max # iterations
calibratefrompairwisedistances!(net, netdist, taxa, #verbose=true,
  forceMinorLength0=true, ultrametric=true)
est = pairwisetaxondistancematrix(net)
@test est ≈ netdist[o,o] atol=.06
est = [e.length for e in net.edge]
el = [.2,.2,1.3,1.5,1,.3,2.8,.2,3,.5,3,0,2.5,.5,.5]
@test est ≈ el atol=0.06
end

@testset "starting BL by fitting data-estimated distances" begin
trait_weights = [3,2,1,1,1]
trait_dat = Vector{Vector{Union{Missings.Missing,Int}}}([
  [4,3,4,1,3], # taxon 1
  [4,3,1,1,3],
  [4,3,1,1,2],
  [4,3,1,missing,missing],
  [4,3,1,4,2], # taxon 5
])
net = readnewick("((t1,t2),(t3,t4),t5);");
# printnodes(net) # leaf named "ti" does have node.number i
PhyloNetworks.startingBL!(net, trait_dat, trait_weights)
# ((t1:0.159,t2:0.0):0.099,(t3:0.018,t4:0.0):0.0,t5:0.146);
@test all([e.length for e in net.edge] .≈
  [0.1586697132288586, 0.0001, 0.09935344568620236, 0.018499865670972223, 0.0001, 0.0001, 0.14620582339229848])
end
