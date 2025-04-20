@testset "expected f-stat" begin

# 7 tips, 2 blobs: a 3-cycle (h=1) and a 4-blob (h=2) not galled
# rooted, 4 non-external cut edges (3 when semidirected)
# has a degree-2 node
nwk = "((((a1:0.24)#H2:0.58::0.7,(((#H2:0.62):0.49,(a21:0.14,a22:0.86):0.51):0.29)#H1:0.19::0.6):0.88,(#H1:0.82,a3:0.89):0.33):0.17,((#H3:0.67,((b1:0.87)#H3:0.86::0.8,b2:0.42):0.69):0.18,b3:0.66):0.43);"
net = readnewick(nwk)
f2 = PhyloNetworks.expectedf2matrix!(net)
#= manual calculations:
m = PhyloNetworks.descendenceweight(net)[:tips]
Ω = m * Diagonal([e.length for e in net.edge]) * transpose(m) # covariance
@test Matrix(vcv(net)) ≈ Ω
f2m = -2Ω
for i in 1:7, j in 1:7
    f2m[i,j] += Ω[i,i] + Ω[j,j]
end
=#
f2m = [
0.0      1.608868 2.328868 2.495188 4.023188 3.244388 2.614388;
1.608868 0.0      1.0      2.4652   4.178    3.3992   2.7692;
2.328868 1.0      0.0      3.1852   4.898    4.1192   3.4892;
2.495188 2.4652   3.1852   0.0      3.8888   3.11     2.48;
4.023188 4.178    4.898    3.8888   0.0      1.8948   2.7288;
3.244388 3.3992   4.1192   3.11     1.8948   0.0      1.95;
2.614388 2.7692   3.4892   2.48     2.7288   1.95     0.0
]
@test f2 ≈ f2m

end
