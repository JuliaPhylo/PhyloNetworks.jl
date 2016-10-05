 # test for trait evolution
 # Claudia November 2015
## modified by Paul Bastide
## still not an automatic test function, needs work


using PhyloNetworks
using GLM

#include("../src/types.jl")
#include("../src/functions.jl")

tree= "(A,((B,#H1),(C,(D)#H1)));"

net=readTopologyLevel1(tree)
printEdges(net)

## Re-root the tree so that it matches my example
rootatnode!(net, "A")
printEdges(net)

## Preorder
directEdges!(net)
preorder!(net)

## V matrix
V1 = sharedPathMatrix(net)

## By hand V matrix
l = ones(1, 9)
l[1] = 0.5; l[9] = 0.5
gamma = ones(1, 9); gamma[3] = 0.1; gamma[7] = 0.9
V2 = [0 0 0 0 0 0 0 0 0;
    0 l[1] 0 0 0 0 0 0 0;
    0 0 l[9] l[9] l[9] l[9] l[9] l[9] l[9];
    0 0 0 l[9]+l[4] l[9] l[9]+l[4] l[9]+gamma[3]*l[4] l[9] l[9]+gamma[3]*l[4];
    0 0 0 0 l[9]+l[8] l[9] l[9]+gamma[7]*l[8] l[9]+l[8] l[9]+gamma[7]*l[8];
    0 0 0 0 0 l[9]+l[4]+l[2] l[9]+gamma[3]l[4] l[9] l[9]+gamma[3]l[4];
    0 0 0 0 0 0 l[9]+gamma[3]^2*(l[3]+l[4])+gamma[7]^2*(l[7]+l[8]) l[9]+gamma[7]*l[8] l[9]+gamma[3]^2*(l[3]+l[4])+gamma[7]^2*(l[7]+l[8]);
    0 0 0 0 0 0 0 l[9]+l[8]+l[5] l[9]+gamma[7]*l[8];
    0 0 0 0 0 0 0 0 l[9]+gamma[3]^2*(l[3]+l[4])+gamma[7]^2*(l[7]+l[8])+l[6]]

for i in 1:9
    for j in 1:9
        if i > j
            V2[i,j] = V2[j,i]
        end
    end
end

## Seem identical
V1[:All]
V2

## Return false (floating errors ?)
V1[:All] == V2

## A hanful of values are not strictly equal
V1[:All] .== V2

## But accounting for floating errors, the to matrices are equal
isapprox(V1[:All], V2) # Frobenius norm of the difference is < tol

eq = zeros(9, 9)
for i in 1:9
    for j in 1:9
        eq[i,j] = isapprox(V1.V[i,j], V2[i,j]) # every couple of terms are approx equal
    end
end
eq
