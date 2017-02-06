# test for trait evolution
# Claudia November 2015
## modified by Paul Bastide

function test_show(x)
    io = IOBuffer()
    show(io, x)
end 

tree_str= "(A:0.5,((B:1,#H1:1::0.1):1,(C:1,(D:1)#H1:1::0.9):1):0.5);"
net = readTopology(tree_str)
preorder!(net)

## V matrix
V1 = sharedPathMatrix(net)
test_show(V1)
@show V1

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
# Make it symetric
for i in 1:9
    for j in 1:9
        if i > j
            V2[i,j] = V2[j,i]
        end
    end
end

nodesV2 = [-2, 1, -3, -4, -5, 2, 3, 4, 5] # root was number 6 before: with readTopologyLevel1 + rootatnode
ind = indexin(V1.nodeNumbersTopOrder, nodesV2)
V2 = V2[ind, ind]
test_show(V2)

@test_approx_eq V1[:All] V2

########################
## Incidence Matrix Test
########################
tree_str= "(A:0.5,((B:1,#H1:1::0.4):1,(C:1,(D:1)#H1:1::0.6):1):0.5);"
net = readTopology(tree_str)
preorder!(net)

T = incidenceMatrix(net)

T2 =  [1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
       1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
       1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
       1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0
       1.0  1.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0
       1.0  1.0  0.6  0.0  0.4  1.0  0.0  0.0  0.0
       1.0  1.0  0.6  0.0  0.4  1.0  1.0  0.0  0.0
       1.0  1.0  0.0  0.0  1.0  0.0  0.0  1.0  0.0
       1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0]

@test_approx_eq T[:All] T2
