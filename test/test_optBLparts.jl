# test the components in optBL separately
# Claudia January 2015

globalerror = false
#println("--------- Case G --------------")
include("../examples/case_g_example.jl");

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(net,d)

oldht = net.ht

x = [0.3,1.0,1.5,2.0]
err = false

#println("x is $(x) all changed-----------")

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

    reduce(&,[q.qnet.changed for q in d.quartet]) || error("all qnet should be changed")

    update!(net,x)
    net.ht == x || ("net.ht not correctly changed to x with update")
catch
    global err = true
end

for q in d.quartet
    update!(q.qnet,oldht,net)
end
update!(net,oldht)


x = [0.1,0.2,0.1,2.0] # changing t9 only
#println("x is $(x) changing t9 only-----------")

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
    global err = true
end

for q in d.quartet
    update!(q.qnet,oldht,net)
end
update!(net,oldht)

x = [0.3,0.2,0.1,2.0] # changing gamma and t9 only
#println("x is $(x) changing gamma and t9 only-----------")

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
    global err = true
end

for q in d.quartet
    update!(q.qnet,oldht,net)
end
update!(net,oldht)


# ---- calculateExpCF
x = [0.3,0.2,0.1,2.0] # changing gamma and t9 only
#println("---- calculate expCF for $(x)")
try
    calculateExpCFAll!(d,x,net)
    reduce(&,map(approxEq,q1.qnet.expCF,[(1-x[1])/3*exp(-x[2])+x[1]/3*exp(-x[2]-x[3]-x[4]),(1-x[1])*(1-2/3*exp(-x[2]))+x[1]*(1-2/3*exp(-x[2]-x[3]-x[4])),
                                    (1-x[1])/3*exp(-x[2])+x[1]/3*exp(-x[2]-x[3]-x[4])])) || error("q1 expCF wrong")
    reduce(&,map(approxEq,q2.qnet.expCF, [(1-x[1])*(1-2/3*exp(-x[3]))+x[1]*(1/3*exp(-x[4])),(1-x[1])*(1/3*exp(-x[3]))+x[1]*(1-2/3*exp(-x[4])),
                                     (1-x[1])/3*exp(-x[3])+x[1]/3*exp(-x[4])])) || error("q2 expCF wrong")
    reduce(&,map(approxEq,q3.qnet.expCF, [(1-x[1])/3*exp(-x[3])+x[1]/3*exp(-x[4]),(1-x[1])*(1/3*exp(-x[3]))+x[1]*(1-2/3*exp(-x[4])),(1-x[1])*(1-2/3*exp(-x[3]))+x[1]*(1/3*exp(-x[4]))])) || error("q3 expCF wrong")
    reduce(&,map(approxEq,q4.qnet.expCF,[1/3*exp(-x[2]-x[3]),1-2/3*exp(-x[2]-x[3]),1/3*exp(-x[2]-x[3])])) || error("q4 expCF wrong")
    reduce(&,map(approxEq,q5.qnet.expCF,[(1-x[1])/3*exp(-x[2])+x[1]/3*exp(-x[2]-x[3]),(1-x[1])*(1-2/3*exp(-x[2]))+x[1]*(1-2/3*exp(-x[2]-x[3])),
                                    (1-x[1])/3*exp(-x[2])+x[1]/3*exp(-x[2]-x[3])])) || error("q5 expCF wrong")
catch
    global err = true
end

#logPseudoLik(d)

if !err
    #println("-------------Case G: NO ERRORS!------------")
else
    globalerror = true
end

#println("--------- Case F Bad Diamond I --------------")
include("../examples/case_f_example.jl");
parameters!(net)

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(net,d)

oldht = net.ht

x = [0.4,0.2,0.1]
err = false

#println("x is $(x) all changed-----------")

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

    reduce(&,[q.qnet.changed for q in d.quartet]) || error("all qnet should be changed")

    update!(net,x)
    net.ht == x || ("net.ht not correctly changed to x with update")
catch
    global err = true
end

for q in d.quartet
    update!(q.qnet,oldht,net)
end
update!(net,oldht)


x = [0.1,0.2,0.1] # changing gammaz1, gammaz2 only
#println("x is $(x) changing gammaz1, gammaz2 only-----------")

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
    global err = true
end

for q in d.quartet
    update!(q.qnet,oldht,net)
end
update!(net,oldht)


# ---- calculateExpCF
x = [0.3,0.2,0.1]
#println("---- calculate expCF for $(x)")
try
    calculateExpCFAll!(d,x,net)
    reduce(&,map(approxEq,q1.qnet.expCF,[(1-x[2]-x[3])/3,(1+2*x[2]-x[3])/3,(1-x[2]+2*x[3])/3])) || error("q1 expCF wrong")
    reduce(&,map(approxEq,q2.qnet.expCF, [1-2/3*exp(-x[1]),1/3*exp(-x[1]),1/3*exp(-x[1])])) || error("q2 expCF wrong")
    reduce(&,map(approxEq,q3.qnet.expCF, [1/3*exp(-x[1]+log(1-x[3])),1/3*exp(-x[1]+log(1-x[3])),1-2/3*exp(-x[1]+log(1-x[3]))])) || error("q3 expCF wrong")
    reduce(&,map(approxEq,q4.qnet.expCF,[1/3*exp(-x[1]+log(1-x[2])),1-2/3*exp(-x[1]+log(1-x[2])),1/3*exp(-x[1]+log(1-x[2]))])) || error("q4 expCF wrong")
    reduce(&,map(approxEq,q5.qnet.expCF,[(1-x[2]-x[3])/3,(1+2*x[2]-x[3])/3,(1-x[2]+2*x[3])/3])) || error("q5 expCF wrong")

catch
    global err = true
end

#logPseudoLik(d)

if !err
    #println("-------------Case F: NO ERRORS!------------")
else
    globalerror = true
end


#println("--------- Case I Bad Diamond II --------------")
include("../examples/case_i_example.jl");

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);

d = DataCF([q1,q2,q3,q4,q5]);
extractQuartet!(net,d)

oldht = net.ht

x = [0.2,0.2,0.1,0.1,0.1]
err = false

#println("x is $(x) all changed-----------")

try
    update!(q1.qnet,x,net)
    q1.qnet.edge[7].gamma != x[1] || q1.qnet.edge[2].gamma != 1-x[1]  ? error("qnet gamma not correct") : nothing
    q1.qnet.edge[4].length != x[3] || q1.qnet.edge[8].length != x[4] ? error("qnet edge lengths not correct") : nothing

    update!(q2.qnet,x,net)
    q2.qnet.edge[8].gamma != x[1] || q2.qnet.edge[4].gamma != 1-x[1]  ? error("qnet gamma not correct") : nothing
    (q2.qnet.edge[4].length !=x[2] || q2.qnet.edge[6].length !=x[3] || q2.qnet.edge[8].length !=x[4] || q2.qnet.edge[9].length !=x[5])  ? error("qnet edges lengths not correct") : nothing

    update!(q3.qnet,x,net)
    q3.qnet.edge[8].gamma != x[1] || q3.qnet.edge[4].gamma != 1-x[1]  ? error("qnet gamma not correct") : nothing
    (q2.qnet.edge[4].length !=x[2] || q2.qnet.edge[6].length !=x[3] || q2.qnet.edge[8].length !=x[4] || q2.qnet.edge[9].length !=x[5])  ? error("qnet edges lengths not correct") : nothing

    update!(q4.qnet,x,net)
    q4.qnet.edge[8].gamma != x[1] || q4.qnet.edge[4].gamma != 1-x[1]  ? error("qnet gamma not correct") : nothing
    (q2.qnet.edge[4].length !=x[2] || q2.qnet.edge[6].length !=x[3] || q2.qnet.edge[8].length !=x[4] || q2.qnet.edge[9].length !=x[5])  ? error("qnet edges lengths not correct") : nothing

    update!(q5.qnet,x,net)
    q5.qnet.edge[7].gamma != x[1] || q5.qnet.edge[2].gamma != 1-x[1]  ? error("qnet gamma not correct") : nothing
    q5.qnet.edge[4].length != x[3] || q5.qnet.edge[8].length != x[4] ? error("qnet edge lengths not correct") : nothing

    reduce(&,[q.qnet.changed for q in d.quartet]) || error("all qnet should be changed")

    update!(net,x)
    net.ht == x || ("net.ht not correctly changed to x with update")
catch
    global err = true
end

for q in d.quartet
    update!(q.qnet,oldht,net)
end
update!(net,oldht)

x = [0.1,0.2,1.,1.,1.] # changing t3 only
err = false

#println("x is $(x) changing t3 only-----------")

try
    update!(q1.qnet,x,net)
    q1.qnet.edge[7].gamma != x[1] || q1.qnet.edge[2].gamma != 1-x[1]  ? error("qnet gamma not correct") : nothing
    q1.qnet.edge[4].length != x[3] || q1.qnet.edge[8].length != x[4] ? error("qnet edge lengths not correct") : nothing

    update!(q2.qnet,x,net)
    q2.qnet.edge[8].gamma != x[1] || q2.qnet.edge[4].gamma != 1-x[1]  ? error("qnet gamma not correct") : nothing
    (q2.qnet.edge[4].length !=x[2] || q2.qnet.edge[6].length !=x[3] || q2.qnet.edge[8].length !=x[4] || q2.qnet.edge[9].length !=x[5])  ? error("qnet edges lengths not correct") : nothing

    update!(q3.qnet,x,net)
    q3.qnet.edge[8].gamma != x[1] || q3.qnet.edge[4].gamma != 1-x[1]  ? error("qnet gamma not correct") : nothing
    (q2.qnet.edge[4].length !=x[2] || q2.qnet.edge[6].length !=x[3] || q2.qnet.edge[8].length !=x[4] || q2.qnet.edge[9].length !=x[5])  ? error("qnet edges lengths not correct") : nothing

    update!(q4.qnet,x,net)
    q4.qnet.edge[8].gamma != x[1] || q4.qnet.edge[4].gamma != 1-x[1]  ? error("qnet gamma not correct") : nothing
    (q2.qnet.edge[4].length !=x[2] || q2.qnet.edge[6].length !=x[3] || q2.qnet.edge[8].length !=x[4] || q2.qnet.edge[9].length !=x[5])  ? error("qnet edges lengths not correct") : nothing

    update!(q5.qnet,x,net)
    q5.qnet.edge[7].gamma != x[1] || q5.qnet.edge[2].gamma != 1-x[1]  ? error("qnet gamma not correct") : nothing
    q5.qnet.edge[4].length != x[3] || q5.qnet.edge[8].length != x[4] ? error("qnet edge lengths not correct") : nothing

    [q.qnet.changed for q in d.quartet] == [false, true, true, true, false] || error("not all qnet should be changed")

    update!(net,x)
    net.ht == x || ("net.ht not correctly changed to x with update")
catch
    global err = true
end

for q in d.quartet
    update!(q.qnet,oldht,net)
end
update!(net,oldht)



# ---- calculateExpCF
x = [0.2,0.2,0.1,0.1,0.1]
#println("---- calculate expCF for $(x)")
try
    calculateExpCFAll!(d,x,net)
    reduce(&,map(approxEq,q1.qnet.expCF,[(1-x[1])*(1/3*exp(-x[3]))+x[1]*(1-2/3*exp(-x[4])),(1-x[1])*(1-2/3*exp(-x[3]))+x[1]*(1/3*exp(-x[4])),(1-x[1])*(1/3*exp(-x[3]))+x[1]*(1/3*exp(-x[4]))])) || error("q1 expCF wrong")
    t=-log(1+x[1]*(1-exp(-x[3]))-x[1]*x[1]*(1-exp(-x[5]-x[4]))-x[1]*x[1]*(1-exp(-x[3]))-(1-x[1])*(1-x[1])*(1-exp(-x[2])))
    reduce(&,map(approxEq,q2.qnet.expCF, [1-2/3*exp(-t),1/3*exp(-t),1/3*exp(-t)])) || error("q2 expCF wrong")
    t=-log(1+x[1]*(1-exp(-x[3]-x[5]))-x[1]*x[1]*(1-exp(-x[4]))-x[1]*x[1]*(1-exp(-x[3]-x[5]))-(1-x[1])*(1-x[1])*(1-exp(-x[2])))
    reduce(&,map(approxEq,q3.qnet.expCF, [1/3*exp(-t),1/3*exp(-t),1-2/3*exp(-t)])) || error("q3 expCF wrong")
    t=-log(1+x[1]*(1-exp(-x[5]))-x[1]*x[1]*(1-exp(-x[5]))-x[1]*x[1]*(1-exp(-x[4]))-(1-x[1])*(1-x[1])*(1-exp(-x[3]-x[2])))
    reduce(&,map(approxEq,q4.qnet.expCF, [1/3*exp(-t),1-2/3*exp(-t),1/3*exp(-t)])) || error("q4 expCF wrong")
    reduce(&,map(approxEq,q5.qnet.expCF,[(1-x[1])*(1/3*exp(-x[3]))+x[1]*(1-2/3*exp(-x[4])),(1-x[1])*(1-2/3*exp(-x[3]))+x[1]*(1/3*exp(-x[4])),(1-x[1])*(1/3*exp(-x[3]))+x[1]*(1/3*exp(-x[4]))])) || error("q5 expCF wrong")

catch
    global err = true
end

#logPseudoLik(d)

if !err
    #println("-------------Case I: NO ERRORS!------------")
else
    globalerror = true
end

@test !globalerror
