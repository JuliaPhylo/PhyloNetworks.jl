# tests with the 5 taxon networks write to parenthetical format
# Claudia March 2015
# need to draw in dendroscope to compare, still not automatic


# types in "types.jl"
include("../src/types.jl")
include("../src/functions.jl")

# needed modules:
using Base.Collections # for updateInCycle with priority queue

tests = ["C","F","G","H","J","D","E","I"];
wrong = String[];

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

tests = ["F","G","H","J","I"];
t="I"
for t in tests
    println("running $(t)")
    net = nothing;
    tree = whichtree(t)
    f = open("prueba_tree.txt","w")
    write(f,tree)
    close(f)
    net = readTopologyUpdate("prueba_tree.txt");
    printEdges(net)
    printNodes(net)
    written = writeTopologyLevel1(net)
    written2 = writeTopologyLevel1(net,true)
    tree
end
