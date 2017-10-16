using PhyloNetworks
using Base.Test

m1 = BinaryTraitSubstitutionModel(1.0, 2.0)

@show m1

m2 = EqualRatesSubstitutionModel(4, 3.0)

@show m2

@test nStates(m1) == 2

@test nStates(m2) == 4

Q(m1)

Q(m2)

P(m1, 4.0)

P(m2, 4.0)

randomTrait(m1, 4.0, [1,2,1,2,2])

randomTrait(m2, 4.0, [1,3,4,2,1])

#Fixit test on example tree, network
