runall = false;
@testset "Testing traitLikDiscrete" begin
global net, n1, n2, d
@testset "Testing Substitution Models, P and Q matrices" begin

m1 = BinaryTraitSubstitutionModel(1.0, 2.0);
@test_logs show(devnull, m1)
m1 = BinaryTraitSubstitutionModel(1.0,2.0, ["carnivory", "non-carnivory"]);
@test nStates(m1)==2
@test PhyloNetworks.nparams(m1)==2
@test_logs show(devnull, m1)
@test_throws ErrorException PhyloNetworks.BinaryTraitSubstitutionModel(-1.0,2.0)
m2 = EqualRatesSubstitutionModel(4, 3.0);
@test nStates(m2)==4
@test PhyloNetworks.nparams(m2)==1
m2 = EqualRatesSubstitutionModel(4, 3.0, ["S1","S2","S3","S4"]);
@test_logs show(devnull, m2)
m3 = TwoBinaryTraitSubstitutionModel([2.0,1.2,1.1,2.2,1.0,3.1,2.0,1.1],
["carnivory", "noncarnivory", "wet", "dry"]);
@test_logs show(devnull, m3)
@test nStates(m3)==4
@test PhyloNetworks.nparams(m3)==8

@test Q(m1) == SMatrix{2,2}(-1.0, 2.0, 1.0, -2.0)
@test Q(m2) == SMatrix{4,4}(-9.0, 3, 3, 3, 3, -9, 3, 3, 3, 3, -9, 3, 3, 3, 3, -9)
@test Q(m3) ≈ SMatrix{4,4}(-3.0, 3.1, 1.2, 0.0, 1.0, -4.2, 0.0, 2.2, 2.0, 0.0, -3.2, 1.1, 0.0, 1.1, 2.0, -3.3) atol=1e-4
@test P(m1, 0.5) ≈ SMatrix{2,2}(0.7410433867161432,0.5179132265677134,0.2589566132838567,0.4820867734322865) atol=1e-15
@test P(m2, 0.1) ≈ SMatrix{4,4}(0.4758956589341516,0.17470144702194956,0.17470144702194962,0.17470144702194945,0.17470144702194956,0.47589565893415153,0.17470144702194967,0.17470144702194945,0.17470144702194962,0.17470144702194967,0.4758956589341516,0.17470144702194945,0.17470144702194945,0.17470144702194945,0.17470144702194945,0.4758956589341518) atol=1e-15
@test P(m3, 0.5) ≈ SMatrix{4,4}(0.39839916380463375,0.36847565707393248,0.23055614536582461,0.22576141081414305,0.1545971371424259,0.25768553619230444,0.14816051303688715,0.24300762748855972,0.29194735222005136,0.20198250750421617,0.35349416558860958,0.20267178083856716,0.15505634683288913,0.17185629922954704,0.26778917600867863,0.32855918085873009) atol=1e-15
@test P(m1, [0.02,0.01]) ≈ StaticArrays.SArray{Tuple{2,2},Float64,2,4}[[0.980588 0.0194118; 0.0388236 0.961176], [0.990149 0.00985149; 0.019703 0.980297]] atol=1e-6
@test P(m2, [0.02,0.01]) ≈ Array{Float64,2}[[0.839971 0.053343 0.053343 0.053343; 0.053343 0.839971 0.053343 0.053343; 0.053343 0.053343 0.839971 0.053343; 0.053343 0.053343 0.053343 0.839971], [0.91519 0.0282699 0.0282699 0.0282699; 0.0282699 0.91519 0.0282699 0.0282699; 0.0282699 0.0282699 0.91519 0.0282699; 0.0282699 0.0282699 0.0282699 0.91519]] atol=1e-6

end

@testset "Testing random discrete trait simulation" begin

m1 = BinaryTraitSubstitutionModel(1.0,2.0, ["carnivory", "non-carnivory"]);
m2 = EqualRatesSubstitutionModel(4, 3.0, ["S1","S2","S3","S4"]);
# on a single branch
Random.seed!(12345);
@test randomTrait(m1, 0.2, [1,2,1,2,2]) == [1,2,1,1,2]
Random.seed!(12345);
@test randomTrait(m2, 0.05, [1,3,4,2,1]) == [1,3,4,2,1]
# on a network
net = readTopology("(A:1.0,(B:1.0,(C:1.0,D:1.0):1.0):1.0);")
Random.seed!(12345);
a,b = randomTrait(m1, net)
@test a == [1 2 1 1 1 1 2]
@test b == ["-2", "-3", "-4", "D", "C", "B", "A"]
if runall
    for e in net.edge e.length = 10.0; end
    @time a,b = randomTrait(m1, net; ntraits=100000) # ~ 0.014 seconds
    mean(a[:,1]) # expect 1.5 at root
    mean(a[:,2]) # expect 1.333 at other nodes
    @time a,b = randomTrait(m2, net; ntraits=100000) # ~ 0.02 seconds
    length([x for x in a[:,1] if x==4])/length(a[:,1]) # expect 0.25
    length([x for x in a[:,2] if x==4])/length(a[:,2])
    length([x for x in a[:,3] if x==4])/length(a[:,3])
    length([x for x in a[:,4] if x==4])/length(a[:,4])
    length([x for x in a[:,5] if x==4])/length(a[:,5])
    length([x for x in a[:,6] if x==4])/length(a[:,6])
    length([x for x in a[:,7] if x==4])/length(a[:,7]) # expect 0.25
end

net2 = readTopology("(((A:4.0,(B:1.0)#H1:1.1::0.9):0.5,(C:0.6,#H1:1.0):1.0):3.0,D:5.0);")
Random.seed!(12345);
a,b = randomTrait(m1, net2; keepInternal=false)
@test a == [1  1  1  2]
@test b == ["D", "C", "B", "A"]
Random.seed!(12345);
a,b = randomTrait(m1, net2; keepInternal=true)
@test a == [1  2  1  1  1  1  1  1  1]
@test b == ["-2", "D", "-3", "-6", "C", "-4", "#H1", "B", "A"]
if runall
    for e in net2.edge
        if e.hybrid 
            e.length = 0.0
        end
    end
    a,b = randomTrait(m1, net2; ntraits=100000)
    # plot(net2, showNodeNumber=true) shows: H1 listed 7th, parents listed 4th and 6th
    c = map( != , a[:, 4],a[:, 6] ); # traits when parents have different traits
    n1 = sum(map( ==, a[c,7],a[c,6] )) # 39644 traits: hybrid ≠ major parent
    n2 = sum(map( ==, a[c,7],a[c,4] )) #  4401 traits: hybrid ≠ minor parent
    n1/sum(c) # expected 0.9
    n2/sum(c) # expected 0.1
    for e in net2.edge
        e.length = 0.0
    end
    net2.edge[4].length = 10.0
    a,b = randomTrait(m1, net2; ntraits=100000);
    a[:, 1] == a[:, 2]  # true: root = leaf D, as expected
    a[:, 1] == a[:, 5]  # true: root = leaf C
    mean(a[:, 6]) # expected 1.3333
    a[:, 6] == a[:, 9] # true: major hybrid parent node = leaf A
end

end

@testset "Test discrete likelihood, fixed topology" begin

# test on a tree
#=
likelihood calculated in R using a fixed Q matrix, first with ace() then
with fitDiscrete(), then with fitMK(). problem: they give different results,
see http://blog.phytools.org/2015/09/the-difference-between-different.html
- ace: misses log(#states) in its log-likelihood
- fitDiscrete in geiger: uses empirical prior at root, not stationary dist,
  but "lik" object is very flexible
- fitMk is correct. also great for 2 correlated binary traits

library(ape)
mytree = read.tree(text = "(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);")
states = c(1,1,2,2)
names(states)  = mytree$tip.label
fitER = ace(states, mytree, model="ER", type="discrete")
print(fitER$loglik, digits=17) # log-likelihood = -1.9706530878326345
print(fitER$loglik - log(2), digits=17) #         -2.6638002683925799
print(fitER$rates, digits=17)  # rates = 0.3743971742794559
print(fitER$lik.anc, digits=17)# posterior probs of states at nodes: 3x2 matrix (3 internal nodes, 2 states)

library(geiger)
fitER = fitDiscrete(mytree, states, model="ER")
print(fitER$opt$q12, digits=17) # rates = 0.36836216513047726
print(fitER$opt$lnL, digits=17) # log-likelihood = -2.6626566310743804
lik = fitER$lik
lik(0.3743971742794559, root="given",root.p=c(.5,.5)) # -2.6638002630818232: same as ace + log(2)

library(phytools)
Q2 = matrix(c(-1,1,1,-1),2,2)*fitER$opt$q12
fit2 = fitMk(mytree, states, model="ER", fixedQ=Q2)
print(fit2$logLik, digits=17) # log-likelihood = -2.6638637960257574

fitER = fitDiscrete(mytree, states, model="ARD")
lik = fitER$lik
Q = c(0.29885191850718751, 0.38944304456937912) # q12, q21
lik(Q, root="given", root.p=Q[2:1]/sum(Q)) # -2.6457428692377234
lik(Q, root="flat") # -2.6447321523303113
Q = c(0.2, 0.3) # q12, q21
lik(Q, root="flat") # -2.6754091090953693 .1,.7: -3.3291679800706073
optim(Q, lik, lower=1e-8, control=list(fnscale=-1), root="flat")
# rates = 0.29993140042699212 0.38882902905265493 loglik=-2.6447247349802496

states=c(1,2,1); names(states)=c("A","B","D")
fitER = fitDiscrete(mytree, states, model="ARD"); lik = fitER$lik
lik(Q, root="flat") # -2.1207856874033491
=#

net = readTopology("(A:3.0,(B:2.0,(C:1.0,D:1.0):1.0):1.0);");
tips = Dict("A" => "lo", "B" => "lo", "C" => "hi", "D" => "hi");
m1 = EqualRatesSubstitutionModel(2,0.36836216513047726, ["lo", "hi"]);
fit1 = (@test_logs fitDiscrete(net, m1, tips; fixedparam=true));
@test_logs show(devnull, fit1)
@test loglikelihood(fit1) ≈ -2.6638637960257574
species = ["G","C","A","B","D"]
dat1 = DataFrame(trait = ["hi","hi","lo","lo","hi"], species = species)
m2 = BinaryTraitSubstitutionModel(0.2, 0.3, ["lo", "hi"])
fit2 = (@test_logs fitDiscrete(net, m2, dat1; fixedparam=true))
@test fit2.trait == [[1],[1],[2],[2]]
@test loglikelihood(fit2) ≈ -2.6754091090953693
originalstdout = stdout
redirect_stdout(open("/dev/null", "w"))
fit2 = @test_logs fitDiscrete(net, m2, dat1; verbose=true) # 65 iterations
redirect_stdout(originalstdout)
@test fit2.model.rate ≈ [0.29993140042699212, 0.38882902905265493] atol=2e-4
@test loglikelihood(fit2) ≈ -2.6447247349802496 atol=2e-4
m2.rate = [0.2, 0.3];
dat2 = DataFrame(trait1= ["hi","hi","lo","lo","hi"], trait2=["hi",missing,"lo","hi","lo"]);
fit3 = (@test_logs fitDiscrete(net, m2, species, dat2; fixedparam=true))
@test fit3.loglik ≈ (-2.6754091090953693 - 2.1207856874033491)
PhyloNetworks.fit!(fit3; fixedparam=false)
@test fit3.model.rate ≈ [0.3245640354187991, 0.5079501745877728]
fit3.net = readTopology("(A,(B,(C,D):1.0):1.0);"); # no branch lengths
@test_throws ErrorException PhyloNetworks.fit!(fit3; fixedparam=true)
# if optimized, NLOpt catches the error (due to negative branch lengths)
# and stops with no error, return code FORCED_STOP

# test on a network, 1 hybridization
net = readTopology("(((A:4.0,(B:1.0)#H1:1.1::0.9):0.5,(C:0.6,#H1:1.0::0.1):1.0):3.0,D:5.0);")
# function below used to check that simulation proportions == likelihood
m1 = BinaryTraitSubstitutionModel([1.0, 2.0], [1,2]) # model.label = model.index
function traitprobabilities(model, net, ntraits=10)
    res, lab = randomTrait(model, net; ntraits=ntraits)
    tips = findall(in(tipLabels(net)), lab) # indices of tips: columns in res
    dat = DataFrame(species = lab[tips])
    tmp = StatsBase.countmap([res[i,tips] for i in 1:ntraits])
    i = 0
    prop = Float64[]
    for (k,v) in tmp
        i += 1
        dat[Symbol("x",i)] = k
        push!(prop, v/ntraits)
    end
    npatterns = i
    lik = Float64[]
    for i in 1:npatterns
        fit = fitDiscrete(net, model, dat[[:species, Symbol("x",i)]]; fixedparam=true)
        push!(lik, fit.loglik)
    end
    return dat, prop, lik
end
#=
using PhyloNetworks, StatsBase, DataFrames
d, p, ll = traitprobabilities(m1, net, 100000000);
all(isapprox.(log.(p), ll, atol=1e-3)) # true
hcat(log.(p), ll)
 -1.62173  -1.62184
 -3.00805  -3.00807
 -4.39506  -4.39436
 -3.00747  -3.0082 
 -3.70119  -3.70121
 -3.00759  -3.0082 
 -2.31516  -2.31505
 -2.31554  -2.31499
 -3.0083   -3.0082 
 -3.008    -3.0082 
 -2.31475  -2.31505
 -3.702    -3.70135
 -3.00836  -3.00813
 -3.70033  -3.70121
 -2.31546  -2.31499
 -3.70124  -3.70135
=#
d = DataFrame(species=["D","C","B","A"], x1=[1,1,1,1], x2=[1,2,2,1], x3=[2,2,2,2], x4=[1,1,2,2],
    x5=[2,2,2,1], x6=[2,2,1,1], x7=[1,1,2,1], x8=[2,1,1,1], x9=[2,1,2,1], x10=[1,2,1,2],
    x11=[1,2,1,1], x12=[2,2,1,2], x13=[2,1,1,2], x14=[1,2,2,2], x15=[1,1,1,2], x16=[2,1,2,2])
lik = Float64[]
for i in 1:16
    fit = fitDiscrete(net, m1, d[[:species, Symbol("x",i)]]; fixedparam=true)
    push!(lik, fit.loglik)
end
@test lik ≈ [-1.6218387598967712, -3.008066347196894, -4.3943604143403245, -3.008199100743402,
    -3.70121329832901, -3.0081981601869483, -2.315051933868397, -2.314985711030534,
    -3.0081988850020873, -3.0081983709272504, -2.3150512090547584, -3.70134532205944,
    -3.008132923628349, -3.7012134632082083, -2.3149859724945876, -3.7013460518770915]
fit1 = fitDiscrete(net, m1, d[[:species, :x6]])

# with parameter estimation
net = readTopology("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);")
m1 = BinaryTraitSubstitutionModel([1.0, 1.0], ["lo", "hi"])
dat = DataFrame(species=["C","A","B","D"], trait=["hi","lo","lo","hi"])
fit1 = fitDiscrete(net, m1, dat; fixedparam=true)
@test fit1.loglik ≈ -2.77132013004859
PhyloNetworks.fit!(fit1; fixedparam=false)
@test fit1.model.rate ≈ [0.2722263130324768, 0.34981109618902395] atol=1e-4
@test fit1.loglik ≈ -2.727701700695741
# for information only: function used locally to check for correct parameter estimation
function simulateManyTraits_estimate(ntraits)
    m1 = BinaryTraitSubstitutionModel([1.0, 0.5], [1,2])
    res, lab = randomTrait(m1, net; ntraits=ntraits)
    tips = findall(in(tipLabels(net)), lab) # indices of tips: columns in res
    dat = DataFrame(transpose(res[:,tips])); species = lab[tips]
    return fitDiscrete(net, m1, species, dat)
end
# simulateManyTraits_estimate(100000)
# α=1.1124637623451075, β=0.5604529225895175, loglik=-25587.1  with ntraits=10000
# α=0.9801472136310236, β=0.4891696992781437, loglik=-255755.6 with ntraits=100000
# time with ntraits=100000: 907.2s = 15min 7s (one single processor, no binning of traits with same pattern)

# ancestral state reconstruction - fixit!!
fit1.model.rate[1] = 0.2722263130324768;
fit1.model.rate[2] = 0.34981109618902395;
@test_throws ErrorException ancestralStateReconstruction(fit1, 4) # 1 trait, not 4: error
asr = ancestralStateReconstruction(fit1)
@test names(asr) == [:nodenumber, :nodelabel, :lo, :hi]
@test asr[:nodenumber] == collect(1:9)
@test asr[:nodelabel] == ["A","B","C","D","5","6","7","8","#H1"]
@test asr[:lo] ≈ [1.,1.,0.,0., 0.28602239466671175, 0.31945742289603263,
    0.16855042517785512, 0.7673588716207436, 0.7827758475866091]
@test asr[:hi] ≈ [0.,0.,1.,1.,0.713977605333288, 0.6805425771039674,
    0.8314495748221447, 0.23264112837925616, 0.21722415241339132]
@test fit1.postltw ≈ [-0.08356534477069566, -2.5236181051014333]
end # end of testset, fixed topology

end # of nested testsets
