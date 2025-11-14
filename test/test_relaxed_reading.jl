@testset "test: newick parsing" begin
global net
@testset "readnewick: spaces and comments" begin
    global n1, n2, n3
    #Test newlines, spaces, tabs, and carriage returns
    n1 = readnewick("(A,((B,#H1),(C,(D)#H1)));")
    n2 = readnewick("(A\n,\t (\n\r (B , #H1 ),( C ,(D )#H1 ) \n\n\n\n\n) );")
    @test writenewick(n1) == writenewick(n2)

    #Test spaces in troublesome locations
    n3 = readnewick("( A, ( ( B , #H 1 ) ,( C , ( D )#H1 ) ) ) ;")
    @test writenewick(n1) == writenewick(n3)

    #Test spaces in names
    resultant = readnewick("('Homo sapiens', ((B,#H1),(C,(D)#H1)));")
    expected = readnewick("('Homosapiens',((B,#H1),(C,(D)#H1)));")
    @test writenewick(resultant) == writenewick(expected)

    #Test other white symbols in names
    resultant = readnewick("('H\tom\ro\nsa   p \n\n\n\n\n\n\n\ni\t\t\t\t\t\t\t\ten s   ', ((B,#H1),(C,(D)#H1)));")
    expected = readnewick("('Homosapiens',((B,#H1),(C,(D)#H1)));")
    @test writenewick(resultant) == writenewick(expected)

    # nexus-style comments and negative edge or gamma values
    n1 = (@test_logs (:error, r"expecting non-negative value") (:error, r"expecting non-negative value") readnewick("(E,((B)#H1[&com:1],((D:-1[&theta=2][&prob=0.3],C:[&length:4]1.2[&com 5])[&internalname=G],(#H1:::-0.1,A))))[&com=6];"))
    @test writenewick(n1) == "(E,((B)#H1:::1.0,((D:0.0,C:1.2),(#H1:::0.0,A))));"

    # edge cases that should error
    @test_throws Exception readnewick("(E,((B)#3,((D,C),(#3,A))));") # name not starting with letter
    @test_logs (:warn,r"^Expected H") (:warn,r"^Expected H") readnewick("(E,((B)#P1,((D,C),(#P1,A))));")
    @test_logs (:warn,r"received h") (
    @test_throws Exception readnewick("(E,((B)#H#h,((D,C),(#H#h,A))));") # 2 # sign in the name
    )
    @test_throws Exception readnewick("(E,((B)#H1,((D,C),((F)#H1,A))));") # both H1 internal
    @test_throws Exception readnewick("(E,((B)#H1") # doesn't end with ;
    @test_throws Exception readnewick(IOBuffer("E;")) # Expected beginning of tree with (
end
@testset "edge parameters" begin
 net = readnewick("(A:2e-03,((B:1.2e1,#H1:1e+1):2.2E-02,(D:1.1E+2)#H1:4.4E+2::5.5E-10));")
 @test [e.length for e in net.edge] == [0.002,12,10,0.022,110,440,-1]
 @test [e.gamma for e in net.edge if e.hybrid] == [0.99999999945, 5.5e-10]
end
@testset "ismajor & gamma consistency, and miscellaneous" begin
    net = readnewick(" ((((B)#H1)#H2,((D,C,#H2:::0.8),(#H1,A))));"); # leading white space in string
    @test writenewick(net, round=true, digits=8) == "(#H2:::0.2,((D,C,((B)#H1)#H2:::0.8),(#H1,A)));"
    net = readnewick("(E,((B)#H1:::.5,((D,C),(#H1:::.5,A:))));"); # no edge length after colon
    @test writenewick(net) == "(E,((B)#H1:::0.5,((D,C),(#H1:::0.5,A))));"
    @test !PhyloNetworks.shrinkedge!(net,net.edge[6]) # parent to (D,C)
    @test_throws "at a hybrid" PhyloNetworks.resolvetreepolytomy!(net, net.node[3])
    PhyloNetworks.resolvetreepolytomy!(net, net.node[8])
    @test all(e.y == -1 for e in net.edge) # missing (-1) support values, even after resolving polytomy
    @test all(e.z == -1 for e in net.edge)
    @test net.node[11].number == 7 # new node
    @test length(net.edge) == 11
    @test net.edge[11].number == 6 # new edge, recycled old unused number
    @test net.edge[11].length == 0 # from old node -5 towards new node 7
    @test (getparent(net.edge[11]).number,getchild(net.edge[11]).number)==(-5,7)
    @test (getparent(net.edge[4]).number, getchild(net.edge[4]).number) == (7,4)
    @test (getparent(net.edge[5]).number, getchild(net.edge[5]).number) == (7,5)
    @test_logs (:warn, r"^first colon") (
        @test_throws r"digit after ':' but found I." readnewick( "(((B::50)#H1:.1:.2:.6,(C:1:100,(#H1::.3,A))):Inf:25);"))
    net = readnewick("((((B::50)#H1:.1::.6,(C:.7:100,(#H1::.3,A))):1e30:25));"); # root node pruned. new root node still degree 1
    @test [e.gamma for e in net.edge] == [1,.6,1,.4,1, 1, 1, 1]
    @test [e.y for e in net.edge] == [50,-1,100,.3,-1,-1,-1,25]
    @test all(e.z == -1 for e in net.edge)
    @test writenewick(net) == "(((B)#H1:0.1::0.6,(C:0.7,(#H1:::0.4,A))):1.0e30);"
    @test writenewick(net, support=true) ==
      "(((B::50.0)#H1:0.1::0.6,(C:0.7:100.0,(#H1::0.3:0.4,A))):1.0e30:25.0);"
    net = readnewick("(((a,(b:.5:50)#H2:0.1)i1::70,(#H2::2,c:1)::90)root:0.5);")
    @test writenewick(net, support=true, internallabel=true) ==
      "((a,(b:0.5:50.0)#H2:0.1)i1::70.0,(#H2::2.0,c:1.0)::90.0)root;"
end
@testset "internal nodes, writemulti" begin
    @test writenewick(readnewick("(a,b):0.5;")) == "(a,b);"
    @test writenewick(readnewick("((a,(b)#H1)i1,(#H1,c))r;")) == "((a,(b)#H1)i1,(#H1,c))r;"
    @test writenewick(readnewick("((a,(b)#H1)i1,(#H1,c))r;"), internallabel=false) == "((a,(b)#H1),(#H1,c));"
    n11 = (@test_logs readnewick("((a,(b)#H1)i1,(#H1,c)i2)root:0.5;"));
    n12 = (@test_logs readnewick("(((a,(b)#H2)i1,(#H2,c)i2)root:0.5);")); # root edge was deleted
    # writenewick(net) == "((a,(b)#H1)i1,(#H1,c)i2);"
    # readnewick("((((a,(b)#H1)i1,(#H1,c)i2)root:0.5));"); still has 1 (of the 2) root edges
    # writemultinewick([n1,n2], stdout)
    n11.rooti = 2 # below the hybrid node: will trigger RootMismatch and message below
    originalstdout = stdout
    redirect_stdout(devnull)
    writemultinewick([n11,n12], "test_relaxedreading.net")
    rm("test_relaxedreading.net")
    redirect_stdout(originalstdout)
end
@testset "write for hybrid-lambda" begin
    net = readnewick("((a:1,(bH:1)#H1:1::0.8)H2:5,(#H1:0::0.2,c:1):1);");
    @test hybridlambdaformat(net) == "((a:1.0,(bH:1.0)H1#0.8:1.0)I1:5.0,(H1#0.8:0.0,c:1.0)I2:1.0)I3;"
    net = readnewick("((#H1:::0.2,c),(I3,(b)#H1:::0.8):5);")
    @test hybridlambdaformat(net) == "((H1#0.2,c)I4,(I3,(b)H1#0.2)I5:5.0)I6;"
    net = readnewick("((((B)#H1:::0.7,D)100:4,(#H1:::0.3,E)98:6.2):2,O);") # bootstrap values
    @test hybridlambdaformat(net) == "((((B)H1#0.7,D)I1:4.0,(H1#0.7,E)I2:6.2)I3:2.0,O)I4;"
    net = readnewick("((((D)#H1:::0.7,D)1:4,(#H1:::0.3,E)1:6.2):2,O);") # 2 tips named D
    @test_throws Exception hybridlambdaformat(net)
    net = readnewick("((((B)#H1,D)1:4,(#H1,E)1:6.2):2,O);") # missing gamma
    @test_logs (:error, r"^edge.*0.5$") (@test hybridlambdaformat(net) == "((((B)H1#0.5,D)I1:4.0,(H1#0.5,E)I2:6.2)I3:2.0,O)I4;")
end
end
