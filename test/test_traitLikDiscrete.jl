using PhyloNetworks
using Base.Test
using StaticArrays

runall = false;

m1 = PhyloNetworks.BinaryTraitSubstitutionModel(1.0, 2.0) 
@test_throws ErrorException PhyloNetworks.BinaryTraitSubstitutionModel(-1.0,2.0)m1
@show m1
m2 = EqualRatesSubstitutionModel(4, 3.0)
@show m2
@test PhyloNetworks.nStates(m1) == 2
@test PhyloNetworks.nStates(m2) == 4
@test Q(m1) == SMatrix{2,2}(-1.0, 2.0, 1.0, -2.0)
@test Q(m2) == SMatrix{4,4}(-9.0, 3, 3, 3, 3, -9, 3, 3, 3, 3, -9, 3, 3, 3, 3, -9)
@test P(m1, 0.5) ≈ SMatrix{2,2}(0.7410433867161432,0.5179132265677134,0.2589566132838567,0.4820867734322865) atol=1e-15
@test P(m2, 0.1) ≈ SMatrix{4,4}(0.4758956589341516,0.17470144702194956,0.17470144702194962,0.17470144702194945,0.17470144702194956,0.47589565893415153,0.17470144702194967,0.17470144702194945,0.17470144702194962,0.17470144702194967,0.4758956589341516,0.17470144702194945,0.17470144702194945,0.17470144702194945,0.17470144702194945,0.4758956589341518) atol=1e-15
@test P(m1, [0.02,0.01]) ≈ StaticArrays.SArray{Tuple{2,2},Float64,2,4}[[0.980588 0.0194118; 0.0388236 0.961176], [0.990149 0.00985149; 0.019703 0.980297]] atol=1e-6
@test P(m2, [0.02,0.01]) ≈ Array{Float64,2}[[0.839971 0.053343 0.053343 0.053343; 0.053343 0.839971 0.053343 0.053343; 0.053343 0.053343 0.839971 0.053343; 0.053343 0.053343 0.053343 0.839971], [0.91519 0.0282699 0.0282699 0.0282699; 0.0282699 0.91519 0.0282699 0.0282699; 0.0282699 0.0282699 0.91519 0.0282699; 0.0282699 0.0282699 0.0282699 0.91519]] atol=1e-6

srand(12345);
@test randomTrait(m1, 0.2, [1,2,1,2,2]) == [1,2,1,1,2]
@test randomTrait(m2, 0.05, [1,3,4,2,1]) == [1,3,3,4,1]
a = randomTrait(m1, 100.0, fill(1, (10000,)))

net = readTopology("(A:1.0,(B:1.0,(C:1.0,D:1.0):1.0):1.0);")

srand(12345);
a,b = randomTrait(m1, net)
@test a == [1 2 1 1 1 1 2]
@test b == ["-2", "-3", "-4", "D", "C", "B", "A"]
if runall
 plot(net, :RCall, showNodeNumber=true);
 for e in net.edge e.length = 10.0; end
 @time a,b = randomTrait(m1, net; ntraits=100000)
 mean(a[:,1]) # expect 1.5 at root
 mean(a[:,2]) # expect 1.333 at other nodes
end

net2 = readTopology("(((A:4.0,(B:1.0)#H1:1.1::0.9):0.5,(C:0.6,#H1:1.0):1.0):3.0,D:5.0);")

a,b = randomTrait(m1, net2; keepInternal=false)
@test a == [1  1  1  1]
@test b == [ "-2","D","-3","-6","C","-4","#H1","B","A"]
if runall
    for e in net2.edge
        if e.hybrid 
            e.length = 0.0
        end
    end
    a,b = randomTrait(m1, net2; ntraits=100000)
    plot(net2, :RCall, showNodeNumber=true) # H1 listed 7th, Parents listed 4th and 6th
    c = map( != , a[:, 4],a[:, 6] )
    n1 = sum(map( ==, a[c,7],a[c,6] ))
    n2 = sum(map( ==, a[c,7],a[c,4] ))
    n1/sum(c) # expected 0.9
    n2/sum(c) # expected 0.1
    for e in net2.edge
        e.length = 0.0
    end
    net2.edge[4].length = 10.0
    a,b = randomTrait(m1, net2; ntraits=100000)
    a[:, 1] == a[:, 2]  # check if root == leaf D
    a[:, 1] == a[:, 5]  # check if root == leaf C
    mean(a[:, 6]) # expected 1.3333
    a[:, 6] == a[:, 9] # check if major hybrid parent node == leaf A
end