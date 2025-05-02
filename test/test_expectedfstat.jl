@testset "expected f2, f3, f4" begin

# 7 tips, 2 blobs: a 3-cycle (h=1) and a 4-blob (h=2) not galled
# rooted, 4 non-external cut edges (3 when semidirected)
# has a degree-2 node
nwk = "((((a1:0.24)#H2:0.58::0.7,(((#H2:0.62):0.49,(a21:0.14,a22:0.86):0.51):0.29)#H1:0.19::0.6):0.88,(#H1:0.82,a3:0.89):0.33):0.17,((#H3:0.67,((b1:0.87)#H3:0.86::0.8,b2:0.42):0.69):0.18,b3:0.66):0.43);"
net = readnewick(nwk)
f2 = expectedf2matrix(net)
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

@test_throws "reference taxon bx" PhyloNetworks.expectedf3matrix(net, "bx"; preorder=false)
f3 = PhyloNetworks.expectedf3matrix(net, "b3"; preorder=false)
@test f3 ≈ [
0       1.88736 1.88736 1.2996 0.66  0.66  0
1.88736 0       2.6292  1.392  0.66  0.66  0
1.88736 2.6292  0       1.392  0.66  0.66  0
1.2996  1.392   1.392   0      0.66  0.66  0
0.66    0.66    0.66    0.66   0     1.392 0
0.66    0.66    0.66    0.66   1.392 0     0
0       0       0       0      0     0     0
]

originalstdout = stdout
redirect_stdout(devnull) # to hide progress bar
f4,t = expectedf4table(net, preorder=false)
redirect_stdout(originalstdout)
nt = tablequartetf4(f4, t)
# DataFrame(nt) # splits with a 0: aa|bb, ab3|b1b2, a21a22|xy
@test keys(nt) == (:qind, :t1, :t2, :t3, :t4, :f4_12_34, :f4_13_42, :f4_14_23)
@test nt[:f4_12_34] ≈ [-.64944,-.74184,-.0924,-.0924,0,-.74184,-0.0924,-.0924,
  0,0,0,0,0,0,0,-.74184,-.0924,-.0924,
  0,0,0,0,0,0,0,0,0,0,
  0,0,0,           -.732,-.732,-.732,-.732]
@test nt[:f4_13_42] ≈ [.64944,.74184,-.49536,-.49536,-1.2372,.74184,-.49536,-.49536,
 -1.2372,-1.95936,-1.95936,-2.7012,-1.3716,-1.464,-1.464,.74184,-.49536,-.49536,
 -1.2372,-1.22736,-1.22736,-1.9692,-.6396,-.732,-.732,-1.22736,-1.22736,-1.9692,
 -.6396,-.732,-.732,.732, .732, .732, .732]
@test nt[:f4_14_23] ≈ [0,0,.58776,.58776,1.2372,0,.58776,.58776,
  1.2372, 1.95936, 1.95936, 2.7012, 1.3716, 1.464, 1.464,0,.58776,.58776,
  1.2372, 1.22736, 1.22736, 1.9692, .6396, .732, .732, 1.22736, 1.22736, 1.9692,
  .6396, .732, .732,0,0,0,0]
end
