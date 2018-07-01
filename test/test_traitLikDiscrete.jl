runall = false;

@testset "Testing Substitution Models, P and Q matrices" begin

m1 = BinaryTraitSubstitutionModel(1.0, 2.0);
@test_nowarn show(DevNull, m1)
m1 = BinaryTraitSubstitutionModel(1.0,2.0, ["carnivory", "non-carnivory"]);
@test nStates(m1)==2
@test PhyloNetworks.nparams(m1)==2
@test_nowarn show(DevNull, m1)
@test_throws ErrorException PhyloNetworks.BinaryTraitSubstitutionModel(-1.0,2.0)
m2 = EqualRatesSubstitutionModel(4, 3.0);
@test nStates(m2)==4
@test PhyloNetworks.nparams(m2)==1
m2 = EqualRatesSubstitutionModel(4, 3.0, ["S1","S2","S3","S4"]);
@test_nowarn show(DevNull, m2)
m3 = TwoBinaryTraitSubstitutionModel([2.0,1.2,1.1,2.2,1.0,3.1,2.0,1.1],
["carnivory", "noncarnivory", "wet", "dry"]);
@test_nowarn show(DevNull, m3)
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

info("\ton a single branch")
srand(12345);
@test randomTrait(m1, 0.2, [1,2,1,2,2]) == [1,2,1,1,2]
srand(12345);
@test randomTrait(m2, 0.05, [1,3,4,2,1]) == [1,3,4,2,1]

info("\ton a network")
net = readTopology("(A:1.0,(B:1.0,(C:1.0,D:1.0):1.0):1.0);")
srand(12345);
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
srand(12345);
a,b = randomTrait(m1, net2; keepInternal=false)
@test a == [1  1  1  2]
@test b == ["D", "C", "B", "A"]
srand(12345);
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
