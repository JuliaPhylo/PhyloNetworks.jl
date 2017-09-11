@testset "readTopology White Symbol Tests" begin
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
