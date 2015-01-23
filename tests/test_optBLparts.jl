# test the components in optBL separately
# Claudia January 2015

println("--------- Case G --------------")
include("../case_g_example.jl");

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(net,d)

oldht = net.ht

x = [0.3,1.0,1.5,2.0]
error = false

println("x is $(x) all changed-----------")

try
    update!(q1.qnet,x,net)
    q1.qnet.edge[3].length !=x[2] || q1.qnet.edge[6].length !=x[3] || q1.qnet.edge[9].length !=x[4] ? error("qnet edges lengths not correct") : nothing
    q1.qnet.edge[5].gamma !=1-x[1] || q1.qnet.edge[7].gamma !=x[1] ? error("qnet edges gammas not correct") : nothing

    update!(q2.qnet,x,net)
    q2.qnet.edge[4].length !=x[3] || q2.qnet.edge[7].length !=x[4] ? error("qnet edges lengths not correct") : nothing
    q2.qnet.edge[3].gamma !=1-x[1] || q2.qnet.edge[5].gamma !=x[1] ? error("qnet edges gammas not correct") : nothing

    update!(q3.qnet,x,net)
    q3.qnet.edge[4].length !=x[3] || q3.qnet.edge[7].length !=x[4] ? error("qnet edges lengths not correct") : nothing
    q3.qnet.edge[3].gamma !=1-x[1] || q3.qnet.edge[5].gamma !=x[1] ? error("qnet edges gammas not correct") : nothing

    update!(q4.qnet,x,net)
    q4.qnet.edge[3].length !=x[2] || q4.qnet.edge[4].length !=x[3] ? error("qnet edges lengths not correct") : nothing

    update!(q5.qnet,x,net)
    q5.qnet.edge[3].length !=x[2] || q5.qnet.edge[6].length !=x[3] ? error("qnet edges lengths not correct") : nothing
    q5.qnet.edge[5].gamma !=1-x[1] || q5.qnet.edge[7].gamma !=x[1] ? error("qnet edges gammas not correct") : nothing

    all([q.qnet.changed for q in d.quartet]) || error("all qnet should be changed")

    update!(net,x)
    net.ht == x || ("net.ht not correctly changed to x with update")
catch
    error = true
end

for q in d.quartet
    update!(q.qnet,oldht,net)
end
update!(net,oldht)


x = [0.1,0.2,0.1,2.0] # changing t9 only
println("x is $(x) changing t9 only-----------")

try
    update!(q1.qnet,x,net)
    q1.qnet.edge[3].length !=x[2] || q1.qnet.edge[6].length !=x[3] || q1.qnet.edge[9].length !=x[4] ? error("qnet edges lengths not correct") : nothing
    q1.qnet.edge[5].gamma !=1-x[1] || q1.qnet.edge[7].gamma !=x[1] ? error("qnet edges gammas not correct") : nothing

    update!(q2.qnet,x,net)
    q2.qnet.edge[4].length !=x[3] || q2.qnet.edge[7].length !=x[4] ? error("qnet edges lengths not correct") : nothing
    q2.qnet.edge[3].gamma !=1-x[1] || q2.qnet.edge[5].gamma !=x[1] ? error("qnet edges gammas not correct") : nothing

    update!(q3.qnet,x,net)
    q3.qnet.edge[4].length !=x[3] || q3.qnet.edge[7].length !=x[4] ? error("qnet edges lengths not correct") : nothing
    q3.qnet.edge[3].gamma !=1-x[1] || q3.qnet.edge[5].gamma !=x[1] ? error("qnet edges gammas not correct") : nothing

    update!(q4.qnet,x,net)
    q4.qnet.edge[3].length !=x[2] || q4.qnet.edge[4].length !=x[3] ? error("qnet edges lengths not correct") : nothing

    update!(q5.qnet,x,net)
    q5.qnet.edge[3].length !=x[2] || q5.qnet.edge[6].length !=x[3] ? error("qnet edges lengths not correct") : nothing
    q5.qnet.edge[5].gamma !=1-x[1] || q5.qnet.edge[7].gamma !=x[1] ? error("qnet edges gammas not correct") : nothing

    [q.qnet.changed for q in d.quartet] == [true,true,true,false,false] || error("q.qnet.changed not correct for all quartets")

    update!(net,x)
    net.ht == x || error("net.ht not correctly changed to x with update")
catch
    error = true
end

for q in d.quartet
    update!(q.qnet,oldht,net)
end
update!(net,oldht)

x = [0.3,0.2,0.1,2.0] # changing gamma and t9 only
println("x is $(x) changing gamma and t9 only-----------")

try
    update!(q1.qnet,x,net)
    q1.qnet.edge[3].length !=x[2] || q1.qnet.edge[6].length !=x[3] || q1.qnet.edge[9].length !=x[4] ? error("qnet edges lengths not correct") : nothing
    q1.qnet.edge[5].gamma !=1-x[1] || q1.qnet.edge[7].gamma !=x[1] ? error("qnet edges gammas not correct") : nothing

    update!(q2.qnet,x,net)
    q2.qnet.edge[4].length !=x[3] || q2.qnet.edge[7].length !=x[4] ? error("qnet edges lengths not correct") : nothing
    q2.qnet.edge[3].gamma !=1-x[1] || q2.qnet.edge[5].gamma !=x[1] ? error("qnet edges gammas not correct") : nothing

    update!(q3.qnet,x,net)
    q3.qnet.edge[4].length !=x[3] || q3.qnet.edge[7].length !=x[4] ? error("qnet edges lengths not correct") : nothing
    q3.qnet.edge[3].gamma !=1-x[1] || q3.qnet.edge[5].gamma !=x[1] ? error("qnet edges gammas not correct") : nothing

    update!(q4.qnet,x,net)
    q4.qnet.edge[3].length !=x[2] || q4.qnet.edge[4].length !=x[3] ? error("qnet edges lengths not correct") : nothing

    update!(q5.qnet,x,net)
    q5.qnet.edge[3].length !=x[2] || q5.qnet.edge[6].length !=x[3] ? error("qnet edges lengths not correct") : nothing
    q5.qnet.edge[5].gamma !=1-x[1] || q5.qnet.edge[7].gamma !=x[1] ? error("qnet edges gammas not correct") : nothing

    [q.qnet.changed for q in d.quartet] == [true,true,true,false,true] || error("q.qnet.changed not correct for all quartets")

    update!(net,x)
    net.ht == x || error("net.ht not correctly changed to x with update")
catch
    error = true
end

for q in d.quartet
    update!(q.qnet,oldht,net)
end
update!(net,oldht)


# ---- calculateExpCF
x = [0.3,0.2,0.1,2.0] # changing gamma and t9 only
println("---- calculate expCF for $(x)")
try
    calculateExpCFAll!(d,x,net)
    all(map(approxEq,q1.qnet.expCF,[(1-x[1])/3*exp(-x[2])+x[1]/3*exp(-x[2]-x[3]-x[4]),(1-x[1])*(1-2/3*exp(-x[2]))+x[1]*(1-2/3*exp(-x[2]-x[3]-x[4])),
                                    (1-x[1])/3*exp(-x[2])+x[1]/3*exp(-x[2]-x[3]-x[4])])) || error("q1 expCF wrong")
    all(map(approxEq,q2.qnet.expCF, [(1-x[1])*(1-2/3*exp(-x[3]))+x[1]*(1/3*exp(-x[4])),(1-x[1])*(1/3*exp(-x[3]))+x[1]*(1-2/3*exp(-x[4])),
                                     (1-x[1])/3*exp(-x[3])+x[1]/3*exp(-x[4])])) || error("q2 expCF wrong")
    all(map(approxEq,q3.qnet.expCF, [(1-x[1])/3*exp(-x[3])+x[1]/3*exp(-x[4]),(1-x[1])*(1/3*exp(-x[3]))+x[1]*(1-2/3*exp(-x[4])),(1-x[1])*(1-2/3*exp(-x[3]))+x[1]*(1/3*exp(-x[4]))])) || error("q3 expCF wrong")
    all(map(approxEq,q4.qnet.expCF,[1/3*exp(-x[2]-x[3]),1-2/3*exp(-x[2]-x[3]),1/3*exp(-x[2]-x[3])])) || error("q4 expCF wrong")
    all(map(approxEq,q5.qnet.expCF,[(1-x[1])/3*exp(-x[2])+x[1]/3*exp(-x[2]-x[3]),(1-x[1])*(1-2/3*exp(-x[2]))+x[1]*(1-2/3*exp(-x[2]-x[3])),
                                    (1-x[1])/3*exp(-x[2])+x[1]/3*exp(-x[2]-x[3])])) || error("q5 expCF wrong")
catch
    error = true
end

#logPseudoLik(d)

if(!error)
    println("-------------Case G: NO ERRORS!------------")
end

println("--------- Case F Bad Diamond I --------------")
include("../case_f_example.jl");
parameters!(net)

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(net,d)

oldht = net.ht

x = [0.3,0.2,0.1]
error = false

println("x is $(x) all changed-----------")

try
    update!(q1.qnet,x,net)
    q1.qnet.node[1].gammaz !=x[2] || q1.qnet.node[3].gammaz !=x[3]  ? error("qnet gammaz not correct") : nothing

    update!(q2.qnet,x,net)
    q2.qnet.edge[4].length !=x[1]  ? error("qnet edges lengths not correct") : nothing

    update!(q3.qnet,x,net)
    q3.qnet.edge[5].length !=x[1] ? error("qnet edges lengths not correct") : nothing
    q3.qnet.edge[3].length !=-log(1-x[3]) ? error("qnet edges gammaz not correct") : nothing

    update!(q4.qnet,x,net)
    q4.qnet.edge[5].length !=x[1] ? error("qnet edges lengths not correct") : nothing
    q4.qnet.edge[3].length !=-log(1-x[2]) ? error("qnet edges gammaz not correct") : nothing

    update!(q5.qnet,x,net)
    q1.qnet.node[1].gammaz !=x[2] || q1.qnet.node[3].gammaz !=x[3]  ? error("qnet gammaz not correct") : nothing

    all([q.qnet.changed for q in d.quartet]) || error("all qnet should be changed")

    update!(net,x)
    net.ht == x || ("net.ht not correctly changed to x with update")
catch
    error = true
end

for q in d.quartet
    update!(q.qnet,oldht,net)
end
update!(net,oldht)


x = [0.1,0.2,0.1] # changing gammaz1, gammaz2 only
println("x is $(x) changing gammaz1, gammaz2 only-----------")

try
    update!(q1.qnet,x,net)
    q1.qnet.node[1].gammaz !=x[2] || q1.qnet.node[3].gammaz !=x[3]  ? error("qnet gammaz not correct") : nothing

    update!(q2.qnet,x,net)
    q2.qnet.edge[4].length !=x[1]  ? error("qnet edges lengths not correct") : nothing

    update!(q3.qnet,x,net)
    q3.qnet.edge[5].length !=x[1] ? error("qnet edges lengths not correct") : nothing
    q3.qnet.edge[3].length !=-log(1-x[3]) ? error("qnet edges gammaz not correct") : nothing

    update!(q4.qnet,x,net)
    q4.qnet.edge[5].length !=x[1] ? error("qnet edges lengths not correct") : nothing
    q4.qnet.edge[3].length !=-log(1-x[2]) ? error("qnet edges gammaz not correct") : nothing

    update!(q5.qnet,x,net)
    q1.qnet.node[1].gammaz !=x[2] || q1.qnet.node[3].gammaz !=x[3]  ? error("qnet gammaz not correct") : nothing

    [q.qnet.changed for q in d.quartet] == [true,false,true,true,true] || error("q.qnet.changed not correct for all quartets")

    update!(net,x)
    net.ht == x || error("net.ht not correctly changed to x with update")
catch
    error = true
end

for q in d.quartet
    update!(q.qnet,oldht,net)
end
update!(net,oldht)


# ---- calculateExpCF
x = [0.3,0.2,0.1]
println("---- calculate expCF for $(x)")
try
    calculateExpCFAll!(d,x,net)
    all(map(approxEq,q1.qnet.expCF,[(1-x[2]-x[3])/3,(1+2*x[2]-x[3])/3,(1-x[2]+2*x[3])/3])) || error("q1 expCF wrong")
    all(map(approxEq,q2.qnet.expCF, [1-2/3*exp(-x[1]),1/3*exp(-x[1]),1/3*exp(-x[1])])) || error("q2 expCF wrong")
    all(map(approxEq,q3.qnet.expCF, [1/3*exp(-x[1]+log(1-x[3])),1/3*exp(-x[1]+log(1-x[3])),1-2/3*exp(-x[1]+log(1-x[3]))])) || error("q3 expCF wrong")
    all(map(approxEq,q4.qnet.expCF,[1/3*exp(-x[1]+log(1-x[2])),1-2/3*exp(-x[1]+log(1-x[2])),1/3*exp(-x[1]+log(1-x[2]))])) || error("q4 expCF wrong")
    all(map(approxEq,q5.qnet.expCF,[(1-x[2]-x[3])/3,(1+2*x[2]-x[3])/3,(1-x[2]+2*x[3])/3])) || error("q5 expCF wrong")

catch
    error = true
end

#logPseudoLik(d)

if(!error)
    println("-------------Case F: NO ERRORS!------------")
end
