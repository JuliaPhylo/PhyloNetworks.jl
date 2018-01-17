@testset "readTopology Hybrid Node isMajor Tests" begin
	# n1 previously failed: both hybrid parent edges were minor
	n1 = readTopology("(A,((B,#H1),(C,(D)#H1:::0.2)));")
	hybridParents = [x for x in n1.hybrid[1].edge if x.hybrid]
	# Exclusive or to check exactly one parent edge is major
	@test hybridParents[1].isMajor ^ !hybridParents[2].isMajor
	# Gammas must sum to 1 and be non-negative for hybrids
	@test hybridParents[1].gamma + hybridParents[2].gamma == 1
	@test hybridParents[1].gamma >= 0
	@test hybridParents[2].gamma >= 0
	# n2 previously failed: both hybrid parent edges were major
	n2 = readTopology("(A,((B,#H1:::0.5),(C,(D)#H1:::0.5)));")
	hybridParents = [x for x in n2.hybrid[1].edge if x.hybrid]
	# Exclusive or to check exactly one parent edge is major
	@test hybridParents[1].isMajor ^ !hybridParents[2].isMajor
	# Gammas must sum to 1 and be non-negative for hybrids
	@test hybridParents[1].gamma + hybridParents[2].gamma == 1
	@test hybridParents[1].gamma >= 0
	@test hybridParents[2].gamma >= 0
end