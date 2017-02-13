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

##################
## New formula
#################

net = readTopology("(((Ag:5,(#H1:1::0.056,((Ak:2,(E:1,#H2:1::0.004):1):1,(M:2)#H2:1::0.996):1):1):1,(((((Az:1,Ag2:1):1,As:2):1)#H1:1::0.944,Ap:4):1,Ar:5):1):1,(P:4,20:4):3,165:7);")
preorder!(net)
V1 = sharedPathMatrix(net)

hyb = net.hybrid[2]

## Find child edge
hyb_edges = [e.hybrid for e in hyb.edge]
child_edge = hyb.edge[!hyb_edges][1]
child = child_edge.isChild1 ? child_edge.node[1] : child_edge[17].node[2]
## Find number of child edge
nodeNumbersTopOrder = [n.number for n in net.nodes_changed]
p = indexin([child.number], nodeNumbersTopOrder)
## Find parents edges and node numbers
par_edge_1 = hyb.edge[hyb_edges][1]
a_n = PhyloNetworks.getOtherNode(par_edge_1, child)
a = indexin([a_n.number], nodeNumbersTopOrder)
par_edge_2 = hyb.edge[hyb_edges][2]
b_n = PhyloNetworks.getOtherNode(par_edge_2, child)
b = indexin([b_n.number], nodeNumbersTopOrder)

## Gamma value
gam = par_edge_1.gamma

## Two extreme networks
par_edge_1.gamma = 1
par_edge_2.gamma = 0
V_t_1 = sharedPathMatrix(net)

par_edge_1.gamma = 0
par_edge_2.gamma = 1
V_t_2 = sharedPathMatrix(net)

## Descendant indicatrice matrix
des = PhyloNetworks.descendants(par_edge_1)
mask = indexin(des, V_t_1.nodeNumbersTopOrder)
D = zeros(net.numNodes, net.numNodes)
D[mask, mask] = 1.0

## Formula
V2 = gam*V_t_1[:All] + (1-gam)*V_t_2[:All] - gam*(1-gam) * (V_t_1.V[p, p] - V_t_1.V[a, b] + V_t_2.V[p, p] - V_t_2.V[a, b]) .* D

@test_approx_eq V1[:All] V2

