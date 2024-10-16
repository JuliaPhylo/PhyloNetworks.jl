@testset "test: newick parsing" begin
global net
@testset "readTopology: spaces and comments" begin
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

    # nexus-style comments and negative edge or gamma values
    n1 = (@test_logs (:error, r"expecting non-negative value") (:error, r"expecting non-negative value") readTopology("(E,((B)#H1[&com:1],((D:-1[&theta=2][&prob=0.3],C:[&length:4]1.2[&com 5])[&internalname=G],(#H1:::-0.1,A))))[&com=6];"))
    @test writeTopology(n1) == "(E,((B)#H1:::1.0,((D:0.0,C:1.2),(#H1:::0.0,A))));"

    # edge cases that should error
    @test_throws Exception readTopology("(E,((B)#3,((D,C),(#3,A))));") # name not starting with letter
    @test_logs (:warn,r"^Expected H") (:warn,r"^Expected H") readTopology("(E,((B)#P1,((D,C),(#P1,A))));")
    @test_logs (:warn,r"received h") (
    @test_throws Exception readTopology("(E,((B)#H#h,((D,C),(#H#h,A))));") # 2 # sign in the name
    )
    @test_throws Exception readTopology("(E,((B)#H1,((D,C),((F)#H1,A))));") # both H1 internal
    @test_throws Exception readTopology("(E,((B)#H1") # doesn't end with ;
    @test_throws Exception readTopology(IOBuffer("E;")) # Expected beginning of tree with (
end
@testset "isMajor and gamma consistency" begin
    net = readTopology("((((B)#H1)#H2,((D,C,#H2:::0.8),(#H1,A))));");
    @test writeTopology(net, round=true, digits=8) == "(#H2:::0.2,((D,C,((B)#H1)#H2:::0.8),(#H1,A)));"
    net = readTopology("(E,((B)#H1:::.5,((D,C),(#H1:::.5,A))));");
    @test writeTopology(net) == "(E,((B)#H1:::0.5,((D,C),(#H1:::0.5,A))));"

end
@testset "internal nodes, writemulti" begin
    @test writeTopology(readTopology("(a,b):0.5;")) == "(a,b);"
    @test writeTopology(readTopology("((a,(b)#H1)i1,(#H1,c))r;")) == "((a,(b)#H1)i1,(#H1,c))r;"
    @test writeTopology(readTopology("((a,(b)#H1)i1,(#H1,c))r;"), internallabel=false) == "((a,(b)#H1),(#H1,c));"
    n11 = (@test_logs readTopology("((a,(b)#H1)i1,(#H1,c)i2)root:0.5;"));
    n12 = (@test_logs readTopology("(((a,(b)#H2)i1,(#H2,c)i2)root:0.5);")); # root edge was deleted
    # writeTopology(net) == "((a,(b)#H1)i1,(#H1,c)i2);"
    # readTopology("((((a,(b)#H1)i1,(#H1,c)i2)root:0.5));"); still has 1 (of the 2) root edges
    # writeMultiTopology([n1,n2], stdout)
    n11.root = 2 # below the hybrid node: will trigger RootMismatch and message below
    originalstdout = stdout
    redirect_stdout(devnull)
    writeMultiTopology([n11,n12], "test_relaxedreading.net")
    rm("test_relaxedreading.net")
    redirect_stdout(originalstdout)
end
@testset "write for hybrid-lambda" begin
    net = readTopology("((a:1,(bH:1)#H1:1::0.8)H2:5,(#H1:0::0.2,c:1):1);");
    @test hybridlambdaformat(net) == "((a:1.0,(bH:1.0)H1#0.8:1.0)I1:5.0,(H1#0.8:0.0,c:1.0)I2:1.0)I3;"
    net = readTopology("((#H1:::0.2,c),(I3,(b)#H1:::0.8):5);")
    @test hybridlambdaformat(net) == "((H1#0.2,c)I4,(I3,(b)H1#0.2)I5:5.0)I6;"
    net = readTopology("((((B)#H1:::0.7,D)100:4,(#H1:::0.3,E)98:6.2):2,O);") # bootstrap values
    @test hybridlambdaformat(net) == "((((B)H1#0.7,D)I1:4.0,(H1#0.7,E)I2:6.2)I3:2.0,O)I4;"
    net = readTopology("((((D)#H1:::0.7,D)1:4,(#H1:::0.3,E)1:6.2):2,O);") # 2 tips named D
    @test_throws Exception hybridlambdaformat(net)
    net = readTopology("((((B)#H1,D)1:4,(#H1,E)1:6.2):2,O);") # missing gamma
    @test_logs (:error, r"^edge.*0.5$") (@test hybridlambdaformat(net) == "((((B)H1#0.5,D)I1:4.0,(H1#0.5,E)I2:6.2)I3:2.0,O)I4;")
end
end
