@testset "readTopology Hybrid Node isMajor Tests" begin
    global n1, n2, n3
	# n1: see issue #44, both hybrid parent edges were minor
	n1 = readTopology("(A,((B,#H1),(C,(D)#H1:::0.2)));")
	@test writeTopology(n1) == "(A,((B,(D)#H1:::0.8),(C,#H1:::0.2)));"
	hybridParents = [x for x in n1.hybrid[1].edge if x.hybrid]
	# Exclusive or to check exactly one parent edge is major
	@test hybridParents[1].isMajor ⊻ hybridParents[2].isMajor
	# Gammas must sum to 1 and be non-negative for hybrids
	@test hybridParents[1].gamma + hybridParents[2].gamma == 1.0
	@test hybridParents[1].gamma >= 0.0
	@test hybridParents[2].gamma >= 0.0
	# n2 previously failed: both hybrid parent edges were major
	n2 = readTopology("(A,((B,#H1:::0.5),(C,(D)#H1:::0.5)));")
	@test writeTopology(n2) == "(A,((B,#H1:::0.5),(C,(D)#H1:::0.5)));"
	hybridParents = [x for x in n2.hybrid[1].edge if x.hybrid]
	# Exclusive or to check exactly one parent edge is major
	@test hybridParents[1].isMajor ⊻ hybridParents[2].isMajor
	# Gammas must sum to 1 and be non-negative for hybrids
	@test hybridParents[1].gamma + hybridParents[2].gamma == 1.0
	@test hybridParents[1].gamma >= 0.0
	@test hybridParents[2].gamma >= 0.0
	n3 = readTopology("(A,((C,(D)#H1:::0.5),(B,#H1:::0.5)));")
	@test writeTopology(n3) =="(A,((C,(D)#H1:::0.5),(B,#H1:::0.5)));"
end

@testset "parsing extended newick" begin
# second reticulation has minor before major, issue #70:
net1 = (@test_logs readTopology("(((lo,#H3),#H4),((sp)#H3,(mu)#H4));"));
net3 = (@test_logs readTopology("(((sp)#H3,(mu)#H4),((lo,#H3),#H4));"));
@test hardwiredClusterDistance(net1, net3, true) == 0
@test_logs readTopology("((((((((((((((Ae_ca1,Ae_ca2),Ae_ca3))#H1,#H2),(((Ae_um1,Ae_um2),Ae_um3),#H1)),((Ae_co1,Ae_co),(((Ae_un1,Ae_un2),Ae_un3),Ae_un4))),(((Ae_ta1,Ae_ta2),(Ae_ta3,Ae_ta4)),((((((((Ae_lo1,Ae_lo2),Ae_lo3),(Ae_sh1,Ae_sh2)),((Ae_bi1,Ae_bi2),Ae_bi3)),((Ae_se1,Ae_se2),Ae_se3)))#H2,#H3),#H4))),(((T_bo1,(T_bo2,T_bo3)),T_bo4),((T_ur1,T_ur2),(T_ur3,T_ur4)))),(((((Ae_sp1,Ae_sp2),Ae_sp3),Ae_sp4))#H3,((((Ae_mu1,Ae_mu2),Ae_mu3),Ae_mu4))#H4))),Ta_ca),S_va),Er_bo),H_vu);");
@test_throws ErrorException readTopology("(((lo,#H3),#H4);") # "Tree ended while reading in subtree ..."
@test_throws ErrorException readTopology("((lo,#H3),);") # "Expected beginning of subtree but read ..."
end
