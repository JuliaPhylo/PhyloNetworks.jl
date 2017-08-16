@testset "Testing branch length distance LS score" begin

# on a tree:
net = readTopology("(A,(B,(C,D)));")
net2 = readTopology("(((A,B),C),D);")
divergenceTimes = Dict(1 => 0.0, 2 => 0.0, 3 => 0.0, 4 => 0.0, -4 => 10.0, -3 => 15.0, -2 => 25.0)
setBLfromDivergenceTimes!(net2, divergenceTimes)
dnaDistances = pairwiseTaxonDistanceMatrix(net2)
distanceMismatch = pairwiseDistanceLSscore(net, divergenceTimes, dnaDistances)
@test distanceMismatch==130.7669683062202


# on a network:
net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
net2 = readTopology("(((D,(C)#H1:::0.9,(B,#H1:::0.1))),A);")
divergenceTimes = Dict(1 => 0.0, 2 => 0.0, 4 => 0.0, 5 => 0.0, -4 => 10.0, -3 => 15.0, -2 => 25.0, -6 => 7.0, 3 => 5.0)  
setBLfromDivergenceTimes!(net2, divergenceTimes)
dnaDistances = pairwiseTaxonDistanceMatrix(net2)
distanceMismatch = pairwiseDistanceLSscore(net, divergenceTimes, dnaDistances)
@test distanceMismatch==63.75829357816911

end
