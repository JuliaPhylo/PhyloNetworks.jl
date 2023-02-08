# test for hasEdge of a QuartetNetwork
# Claudia January 2015
# Also, tests for net.ht, net.numht, qnet.indexht

#println("----- Case G ------")
include("../examples/case_g_example.jl");
# include(joinpath(dirname(pathof(PhyloNetworks)),  "..","examples","case_g_example.jl"))

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
#parameters!(net)

extractQuartet!(net,d)
error1 = false
try
    q1.qnet.hasEdge == [true,true,true,true] ? nothing : error("q1 wrong hasEdge")
    q2.qnet.hasEdge == [true,false,true,true] ? nothing : error("q2 wrong hasEdge")
    q3.qnet.hasEdge == [true,false,true,true] ? nothing : error("q3 wrong hasEdge")
    q4.qnet.hasEdge == [false,true,true,false] ? nothing : error("q4 wrong hasEdge")
    q5.qnet.hasEdge == [true, true,true,false] ? nothing : error("q5 wrong hasEdge")


    net.ht == [0.1,0.2,0.1,1.0] ? nothing : error("net.ht not correct")
    net.numht == [7,3,6,9] ? nothing : error("net.numth not correct")

    q1.qnet.indexht == [1,2,3,4] ? nothing : error("q1.qnet.indexht not correct")
    q2.qnet.indexht == [1,3,4] ? nothing : error("q2.qnet.indexht not correct")
    q3.qnet.indexht == [1,3,4] ? nothing : error("q3.qnet.indexht not correct")
    q4.qnet.indexht == [2,3] ? nothing : error("q4.qnet.indexht not correct")
    q5.qnet.indexht == [1,2,3] ? nothing : error("q5.qnet.indexht not correct")

    q1.qnet.index == [7,3,6,9] ? nothing : error("q1.qnet.index not correct")
    q2.qnet.index == [5,4,7] ? nothing : error("q2.qnet.index not correct")
    q3.qnet.index == [5,4,7] ? nothing : error("q3.qnet.index not correct")
    q4.qnet.index == [3,4] ? nothing : error("q4.qnet.index not correct")
    q5.qnet.index == [7,3,6] ? nothing : error("q5.qnet.index not correct")
catch
    println("---- error in case G -----")
    global error1 = true
end


#println("----- Case F: bad diamond ------")
include("../examples/case_f_example.jl");
# include(joinpath(dirname(pathof(PhyloNetworks)),  "..","examples","case_f_example.jl"))

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
parameters!(net)

extractQuartet!(net,d)
error1 = false
try
    q1.qnet.hasEdge == [false,true,true] ? nothing : error("q1 wrong hasEdge")
    q2.qnet.hasEdge == [true,false,false] ? nothing : error("q2 wrong hasEdge")
    q3.qnet.hasEdge == [true,false,true] ? nothing : error("q3 wrong hasEdge")
    q4.qnet.hasEdge == [true,true,false] ? nothing : error("q4 wrong hasEdge")
    q5.qnet.hasEdge == [false, true,true] ? nothing : error("q5 wrong hasEdge")


    all(map(approxEq,net.ht,[0.1,0.7*(1-exp(-0.2)),0.3*(1-exp(-0.1))])) ? nothing : error("net.ht not correct")
    net.numht == [9,21,22] ? nothing : error("net.numth not correct")

    q1.qnet.indexht == [2,3] ? nothing : error("q1.qnet.indexht not correct")
    q2.qnet.indexht == [1] ? nothing : error("q2.qnet.indexht not correct")
    q3.qnet.indexht == [1,3] ? nothing : error("q3.qnet.indexht not correct")
    q4.qnet.indexht == [1,2] ? nothing : error("q4.qnet.indexht not correct")
    q5.qnet.indexht == [2,3] ? nothing : error("q5.qnet.indexht not correct")

    q1.qnet.index == [1,3] ? nothing : error("q1.qnet.index not correct")
    q2.qnet.index == [4] ? nothing : error("q2.qnet.index not correct")
    q3.qnet.index == [5,3] ? nothing : error("q3.qnet.index not correct")
    q4.qnet.index == [5,3] ? nothing : error("q4.qnet.index not correct")
    q5.qnet.index == [1,3] ? nothing : error("q5.qnet.index not correct")
catch
    println("---- error in case F -----")
    global error1 = true
end

#println("----- Case I: bad diamondII ------")
include("../examples/case_i_example.jl");

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);

extractQuartet!(net,d)
error1 = false
try
    q1.qnet.hasEdge == [true,false,true,false,true] ? nothing : error("q1 wrong hasEdge")
    q2.qnet.hasEdge == [true,true,true,true,true] ? nothing : error("q2 wrong hasEdge")
    q3.qnet.hasEdge == [true,true,true,true,true] ? nothing : error("q3 wrong hasEdge")
    q4.qnet.hasEdge == [true,true,true,true,true] ? nothing : error("q4 wrong hasEdge")
    q5.qnet.hasEdge == [true,false,true,false,true] ? nothing : error("q5 wrong hasEdge")


    net.ht == [0.1,2.,1.,1.,1.] ? nothing : error("net.ht not correct")
    net.numht == [9,4,6,9,10] ? nothing : error("net.numth not correct")

    q1.qnet.indexht == [1,3,5] ? nothing : error("q1.qnet.indexht not correct")
    q2.qnet.indexht == [1,2,3,4,5] ? nothing : error("q2.qnet.indexht not correct")
    q3.qnet.indexht == [1,2,3,4,5] ? nothing : error("q3.qnet.indexht not correct")
    q4.qnet.indexht == [1,2,3,4,5] ? nothing : error("q4.qnet.indexht not correct")
    q5.qnet.indexht == [1,3,5] ? nothing : error("q5.qnet.indexht not correct")

    q1.qnet.index == [7,4,8] ? nothing : error("q1.qnet.index not correct")
    q2.qnet.index == [8,4,6,8,9] ? nothing : error("q2.qnet.index not correct")
    q3.qnet.index == [8,4,6,8,9] ? nothing : error("q3.qnet.index not correct")
    q4.qnet.index == [8,4,5,8,9] ? nothing : error("q4.qnet.index not correct")
    q5.qnet.index == [7,4,8] ? nothing : error("q5.qnet.index not correct")

catch
    println("---- error in case I -----")
    global error1 = true
end

if error1
    throw("errors in has edge")
end


