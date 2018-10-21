include("test_functions_5taxon_read.jl")

tests = ["F","G","H","J","I"];
wrong = AbstractString[];

function whichtree(t::String)
    if(t == "tree")
        tree = "(((6:0.1,4:1.5)1:0.2,7:0.2)5:0.1,8:0.1,10:0.1);" # normal tree
    elseif(t == "C")
        tree = "((((6:0.1,4:1.5),(7:0.2)11#H1),11#H1),8:0.1,10:0.1);" # Case C: bad triangle II
    elseif(t == "F")
        tree = "(((6:0.1,(4)11#H1)1:0.2,(11#H1,7))5:0.1,8:0.1,10:0.1);" # Case F: bad diamond I
    elseif(t == "G")
        tree = "((((6:0.1,4:1.5)1:0.2,(7)11#H1)5:0.1,(11#H1,8)),10:0.1);" # Case G
    elseif(t == "H")
        tree = "((((6,4),#H1),7),(8)#H1,10);" # Case H
    elseif(t == "J")
        tree = "((((6)#H1,4),7),8,(#H1,10));" # Case J
    elseif(t == "D")
        tree = "((((6,4))#H1,(#H1,7)),8,10);" # Case D Bad triangle I
    elseif(t == "E")
        tree = "(((((8,10))#H1,7),#H1),6,4);" # Case E Bad triangle I
    elseif(t == "I")
        tree = "((((8,10))#H1,7),6,(4,#H1));" # Case I Bad diamond II
    else
        error("not a known 5 taxon network case")
    end
    return tree
end


for t in tests
    #println("running $(t)")
    net = nothing;
    tree = whichtree(t)
    net = readTopologyLevel1(tree);
    if(t == "tree")
        try
            testTree(net)
        catch
            println("error in $(t)")
            push!(wrong,t);
        end
    elseif(t == "C")
        try
            testCaseC(net)
        catch
            println("error in $(t)")
            push!(wrong,t);
        end
    elseif(t == "F")
        try
            testCaseF(net)
        catch
            println("error in $(t)")
            push!(wrong,t);
        end
    elseif(t == "G")
        try
            testCaseG(net)
        catch
            println("error in $(t)")
            push!(wrong,t);
        end
    elseif(t == "H")
        try
            testCaseH(net)
        catch
            println("error in $(t)")
            push!(wrong,t);
        end
    elseif(t == "J")
        try
            testCaseJ(net)
        catch
            println("error in $(t)")
            push!(wrong,t);
        end
    elseif(t == "D")
        try
            testCaseD(net)
        catch
            println("error in $(t)")
            push!(wrong,t);
        end
    elseif(t == "E")
        try
            testCaseE(net)
        catch
            println("error in $(t)")
            push!(wrong,t);
        end
    elseif(t == "I")
        try
            testCaseI(net)
        catch
            println("error in $(t)")
            push!(wrong,t);
        end
    else
        error("not a known 5 taxon network case")
    end
end

## if(!isempty(wrong))
##     for t in wrong
##         println("running $(t)")
##         net = nothing;
##         tree = whichtree(t)
##         f = open("prueba_tree.txt","w")
##         write(f,tree)
##         close(f)
##         net = readTopologyUpdate("prueba_tree.txt");
##         if(t == "tree")
##             testTree(net)
##         elseif(t == "C")
##             testCaseC(net)
##         elseif(t == "F")
##             testCaseF(net)
##         elseif(t == "G")
##             testCaseG(net)
##         elseif(t == "H")
##             testCaseH(net)
##         elseif(t == "J")
##             testCaseJ(net)
##         elseif(t == "D")
##             testCaseD(net)
##         elseif(t == "E")
##             testCaseE(net)
##         elseif(t == "I")
##             testCaseI(net)
##         else
##             error("not a known 5 taxon network case")
##         end
##     end
## else
##     println("----------NO ERRORS!----------");
## end

@test isempty(wrong)
