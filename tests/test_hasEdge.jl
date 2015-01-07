# test for hasEdge of a QuartetNetwork
# Claudia January 2015

println("----- Case G ------")
include("../case_g_example.jl")

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(net,d)


q1.qnet.hasEdge == vcat([d==1 for d in ones(9)],false) ? nothing : error("q1 wrong hasEdge")
q2.qnet.hasEdge == vcat([true, false, false],[d==1 for d in ones(7)]) ? nothing : error("q2 wrong hasEdge")
q3.qnet.hasEdge == vcat([false,true, false],[d==1 for d in ones(7)]) ? nothing : error("q3 wrong hasEdge")
q4.qnet.hasEdge == [true, true,true,false,false,true,false,true,false,true] ? nothing : error("q4 wrong hasEdge")
q5.qnet.hasEdge == vcat([d==1 for d in ones(7)],[false,false,true]) ? nothing : error("q5 wrong hasEdge")

# fixit: need to do bad triangle/bad diamond
