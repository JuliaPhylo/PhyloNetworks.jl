# tests with the 5 taxon networks parameters: net.ht, net.numht
# Claudia January 2015

# test functions

#tests = ["C","F","G","H","J","D","E","I"];
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


function whichtest(t::String,net::HybridNetwork)
    if(t == "C")
        all(map(approxEq,net.ht,[exp(-1.)*(1-0.1*(1-exp(-1.))),exp(-1.)*(1-0.9*(1-exp(-1.))),-exp(-2.)*(0.1*(1-exp(-1.))*(1-exp(-1.))*0.9)])) || error("net.ht wrong in case $(t)")
        net.numht == [41,42,43] || error("net.numht wrong in case $(t)")
    elseif(t == "F")
        all(map(approxEq,net.ht,[0.1,0.9*(1-exp(-.2)),0.1*(1-exp(-1.))])) || error("net.ht wrong in case $(t)")
        net.numht == [8,31,32] || error("net.numht wrong in case $(t)")
        net.index == [8,4,6] || error("net.index wrong in case $(t)")
    elseif(t == "G")
        all(map(approxEq,net.ht,[0.1,0.2,0.1,1.0])) || error("net.ht wrong in case $(t)")
        net.numht == [7,3,6,9] || error("net.numht wrong in case $(t)")
        net.index == [7,3,6,9] || error("net.index wrong in case $(t)")
    elseif(t == "H")
        all(map(approxEq,net.ht,[0.1,1.0,1.0,1.0])) || error("net.ht wrong in case $(t)")
        net.numht == [4,3,5,7] || error("net.numht wrong in case $(t)")
        net.index == [4,3,5,7] || error("net.index wrong in case $(t)")
    elseif(t == "J")
        all(map(approxEq,net.ht,[0.1,1.0,1.0,1.0])) || error("net.ht wrong in case $(t)")
        net.numht == [8,4,6,10] || error("net.numht wrong in case $(t)")
        net.index == [8,4,6,10] || error("net.index wrong in case $(t)")
    elseif(t == "D")
        all(map(approxEq,net.ht,[1.0,0.1*(1-exp(-1.)),0.1*0.1*(1-exp(-1.))+0.9*0.9*(1-exp(-2.))])) || error("net.ht wrong in case $(t)")
        net.numht == [8,31,32] || error("net.numht wrong in case $(t)")
    elseif(t == "E")
        all(map(approxEq,net.ht,[1.0,0.9*(1-exp(-1.)),0.1*0.1*(1-exp(-2.))+0.9*0.9*(1-exp(-1.))])) || error("net.ht wrong in case $(t)")
        net.numht == [8,31,32] || error("net.numht wrong in case $(t)")
    elseif(t == "I")
        all(map(approxEq,net.ht,[0.1,2.0,1.0,1.0,1.0])) || error("net.ht wrong in case $(t)")
        net.numht == [9,4,6,9,10] || error("net.numht wrong in case $(t)")
        net.index == [9,4,6,9,10] || error("net.index wrong in case $(t)")
    else
        error("not a known 5 taxon network case")
    end
end

for t in tests
    #println("running $(t)")
    global net = readTopologyLevel1(whichtree(t));
    try
        whichtest(t,net)
    catch
        println("error in $(t)")
        push!(wrong,t);
    end
end

@test isempty(wrong)
