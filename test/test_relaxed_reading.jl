@testset "test: newick parsing" begin
global net
@testset "readTopology White Symbol Tests" begin
    global n1, n2, n3
	#Test newlines, spaces, tabs, and carriage returns
	n1 = readTopology("(A,((B,#H1),(C,(D)#H1)));")
	n2 = readTopology("(A\n,\t (\n\r (B , #H1 ),( C ,(D )#H1 ) \n\n\n\n\n) );")
	@test writeTopology(n1) == writeTopology(n2)

	#Test spaces in troublesome locations
	n3 = readTopology("( A, ( ( B , #H 1 ) ,( C , ( D )#H1 ) ) ) ;")
	@test writeTopology(n1) == writeTopology(n3)

	#Test spaces in names
	resultant = readTopology("('Homo sapiens', ((B,#H1),(C,(D)#H1)));")
	expected = readTopology("('Homosapiens',((B,#H1),(C,(D)#H1)));")
	@test writeTopology(resultant) == writeTopology(expected)

	#Test other white symbols in names
	resultant = readTopology("('H\tom\ro\nsa   p \n\n\n\n\n\n\n\ni\t\t\t\t\t\t\t\ten s   ', ((B,#H1),(C,(D)#H1)));")
	expected = readTopology("('Homosapiens',((B,#H1),(C,(D)#H1)));")
	@test writeTopology(resultant) == writeTopology(expected)
end
@testset "isMajor and gamma consistency" begin
	net = readTopology("((((B)#H1)#H2,((D,C,#H2:::0.8),(#H1,A))));");
	@test writeTopology(net, round=true, digits=8) == "(#H2:::0.2,((D,C,((B)#H1)#H2:::0.8),(#H1,A)));"
	net = readTopologyLevel1("(E,((B)#H1:::.5,((D,C),(#H1:::.5,A))));");
	@test writeTopology(net) == "(D:1.0,C:1.0,((#H1:1.0::0.5,A:1.0):1.0,((B:1.0)#H1:1.0::0.5,E:1.0):1.0):1.0);"
    originalstdout = stdout
	redirect_stdout(open("/dev/null", "w")) # not portable to Windows
	@test_logs PhyloNetworks.printEverything(net)
	redirect_stdout(originalstdout)
end
@testset "internal nodes" begin
	@test writeTopology(readTopology("(a,b):0.5;")) == "(a,b);"
	@test writeTopology(readTopology("((a,(b)#H1)i1,(#H1,c))r;")) == "((a,(b)#H1)i1,(#H1,c))r;"
	@test writeTopology(readTopology("((a,(b)#H1)i1,(#H1,c))r;"), internallabel=false) == "((a,(b)#H1),(#H1,c));"
	@test_logs readTopology("((a,(b)#H1)i1,(#H1,c)i2)root:0.5;");
	@test_logs readTopology("(((a,(b)#H1)i1,(#H1,c)i2)root:0.5);"); # root edge was deleted
	# writeTopology(net) == "((a,(b)#H1)i1,(#H1,c)i2);"
	# readTopology("((((a,(b)#H1)i1,(#H1,c)i2)root:0.5));"); still has 1 (of the 2) root edges
end
end
