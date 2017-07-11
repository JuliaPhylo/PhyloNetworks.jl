@testset "Testing parsimony score & reconstruction" begin

# on a tree:
net = readTopology("(A,(B,(C,D)));")
tips = Dict("A" => 0, "B" => 0, "C" => 1, "D" => 1)
score, states = parsimonyDiscrete(net, tips)
@test score==1
@test states==Dict(4=>Set([1]),-4=>Set([1]),-3=>Set([0]),
         2=>Set([0]),3=>Set([1]),-2=>Set([0]),1=>Set([0]))

# on a network:
net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
score, states = parsimonyDiscrete(net, tips)
@test score==1
@test states==Dict(4=>Set([1]),-4=>Set([0]),-3=>Set([1]),
         2=>Set([0]),-2=>Set([1]),5=>Set([1]),1=>Set([0]))

tips = Dict("A" => 0, "B" => 1, "C" => 0, "D" => 1)
score, states = parsimonyDiscrete(net, tips)
@test score==2
@test states==Dict(4=>Set([0]),-6=>Set([0]),-4=>Set([0]),
  -3=>Set([0]),2=>Set([1]),-2=>Set([0,1]),5=>Set([1]),1=>Set([0]))

# from a data frame and with missing data:
dat = DataFrame(taxon=["A","B","C","D"], trait=[0,0,1,1])
dat[1,2]=NA
score, states = parsimonyDiscrete(net, dat)
@test score==1
@test states==Dict(4=>Set([1]),-6=>Set([1]),-4=>Set([0]),
  -3=>Set([1]),2=>Set([0]),-2=>Set([1]),5=>Set([1]))

end
