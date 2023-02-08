# test of calculateExpCF function
# Claudia December 2014

# Case G -----------------

#println("------ Case G ----------")
include("../examples/case_g_example.jl")
# include(joinpath(dirname(pathof(PhyloNetworks)),  "..","examples","case_g_example.jl"))
error1 = false
ind = 0

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q1);
try
    identifyQuartet!(qnet)
    qnet.which != 1 ? error("qnet which not correctly assigned") : nothing
    qnet.hybrid[1].k != 3 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
    qnet.hybrid[1].typeHyb != 2 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
    qnet.hybrid[1].prev.number != -7 ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 0 || qnet.numHybrids != 0 ? error("qnet should not have hybrid nodes anymore") : nothing
    qnet.t1 != 0.2-log(1-0.1*(1-exp(-1.1))) ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [1,1,2,2] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [2,1,2] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [1/3*exp(-qnet.t1),1-2/3*exp(-qnet.t1),1/3*exp(-qnet.t1)] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 1")
    global error1 |= true
    global ind = 1
end


q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q2);
try
    identifyQuartet!(qnet)
    qnet.which != 2 ? error("qnet which not correctly assigned") : nothing
    qnet.hybrid[1].k != 4 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
    qnet.hybrid[1].typeHyb != 5 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
    !isa(qnet.hybrid[1].prev,Nothing) ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 1 || qnet.numHybrids != 1 ? error("qnet should have 1 hybrid nodes") : nothing
    qnet.t1 != -1 ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [-1,-1,-1,-1] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [-1,-1,-1] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [0.9*(1-2/3*exp(-0.1))+0.1*1/3*exp(-1.0), 0.9*(1/3*exp(-0.1))+0.1*(1-2/3*exp(-1.0)), 0.9*(1/3*exp(-0.1))+0.1*(1/3*exp(-1.0))] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 2")
    global error1 |= true
    global ind = 2
end

q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q3);
try
    identifyQuartet!(qnet)
    qnet.which != 2 ? error("qnet which not correctly assigned") : nothing
    qnet.hybrid[1].k != 4 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
    qnet.hybrid[1].typeHyb != 5 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
    !isa(qnet.hybrid[1].prev,Nothing) ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 1 || qnet.numHybrids != 1 ? error("qnet should have 1 hybrid nodes") : nothing
    qnet.t1 != -1 ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [-1,-1,-1,-1] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [-1,-1,-1] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [0.9*(1/3*exp(-0.1))+0.1*1/3*exp(-1.0), 0.9*(1/3*exp(-0.1))+0.1*(1-2/3*exp(-1.0)), 0.9*(1-2/3*exp(-0.1))+0.1*(1/3*exp(-1.0))] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 3")
    global error1 |= true
    global ind = 3
end

q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q4);
try
    identifyQuartet!(qnet)
    qnet.which != 1 ? error("qnet which not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 0 || qnet.numHybrids != 0 ? error("qnet should not have hybrid nodes anymore") : nothing
    qnet.t1 > 0.30001 || qnet.t1 < 0.2999999999 ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [1,1,2,2] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [2,1,2] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [1/3*exp(-qnet.t1),1-2/3*exp(-qnet.t1),1/3*exp(-qnet.t1)] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 4")
    global error1 |= true
    global ind = 4
end


q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q5);
try
    identifyQuartet!(qnet)
    qnet.which != 1 ? error("qnet which not correctly assigned") : nothing
    qnet.hybrid[1].k != 3 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
    qnet.hybrid[1].typeHyb != 2 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
    qnet.hybrid[1].prev.number != -3 ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 0 || qnet.numHybrids != 0 ? error("qnet should not have hybrid nodes anymore") : nothing
    qnet.t1 != 0.2-log(1-0.1*(1-exp(-0.1))) ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [1,1,2,2] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [2,1,2] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [1/3*exp(-qnet.t1),1-2/3*exp(-qnet.t1),1/3*exp(-qnet.t1)] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 5")
    global error1 |= true
    global ind = 5
end

if error1
    throw("error Case G in quartet $(ind)")
end



# Case F Bad Diamond I -----------------

#println("------ Case F Bad diamond I ----------")
include("../examples/case_f_example.jl");
# include(joinpath(dirname(pathof(PhyloNetworks)),  "..","examples","case_f_example.jl"))
error1 = false
ind = 0
parameters!(net)

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q1);
try
    identifyQuartet!(qnet)
    qnet.which != 2 ? error("qnet which not correctly assigned") : nothing
    qnet.hybrid[1].k != 4 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
    qnet.hybrid[1].typeHyb != 5 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
    !isa(qnet.hybrid[1].prev,Nothing) ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 1 || qnet.numHybrids != 1 ? error("qnet should not have hybrid nodes anymore") : nothing
    qnet.t1 != -1 ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [-1,-1,-1,-1] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [-1,-1,-1] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [(1-(0.7*(1-exp(-0.2)))-(0.3*(1-exp(-0.1))))/3, (1+2*(0.7*(1-exp(-0.2)))-(0.3*(1-exp(-0.1))))/3,(1-(0.7*(1-exp(-0.2)))+2*(0.3*(1-exp(-0.1))))/3] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 1")
    global error1 |= true
    global ind = 1
end


q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q2);
try
    identifyQuartet!(qnet)
    qnet.which != 1 ? error("qnet which not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 0 || qnet.numHybrids != 0 ? error("qnet should have 0 hybrid nodes") : nothing
    qnet.t1 != 0.1 ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [1,1,2,2] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [1,2,2] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [1-2/3*(exp(-0.1)),1/3*(exp(-0.1)),1/3*(exp(-0.1))] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 2")
    global error1 |= true
    global ind = 2
end

q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q3);
try
    identifyQuartet!(qnet)
    qnet.which != 1 ? error("qnet which not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 0 || qnet.numHybrids != 0 ? error("qnet should have 0 hybrid nodes") : nothing
    qnet.t1 != 0.1-log(1-(0.3*(1-exp(-0.1)))) ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [1,1,2,2] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [2,2,1] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [1/3*exp(-qnet.t1),1/3*exp(-qnet.t1),1-2/3*exp(-qnet.t1)] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 3")
    global error1 |= true
    global ind = 3
end


q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q4);

try
    identifyQuartet!(qnet)
    qnet.which != 1 ? error("qnet which not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 0 || qnet.numHybrids != 0 ? error("qnet should have 0 hybrid nodes") : nothing
    qnet.t1 != 0.1-log(1-(0.7*(1-exp(-0.2)))) ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [1,1,2,2] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [2,1,2] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [1/3*exp(-qnet.t1),1-2/3*exp(-qnet.t1),1/3*exp(-qnet.t1)] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 4")
    global error1 |= true
    global ind = 4
end


q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q5);
try
    identifyQuartet!(qnet)
    qnet.which != 2 ? error("qnet which not correctly assigned") : nothing
    qnet.hybrid[1].k != 4 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
    qnet.hybrid[1].typeHyb != 5 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
    !isa(qnet.hybrid[1].prev,Nothing) ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 1 || qnet.numHybrids != 1 ? error("qnet should not have hybrid nodes anymore") : nothing
    qnet.t1 != -1 ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [-1,-1,-1,-1] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [-1,-1,-1] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [(1-0.7*(1-exp(-0.2))-0.3*(1-exp(-0.1)))/3,(1+2*0.7*(1-exp(-0.2))-0.3*(1-exp(-0.1)))/3,(1-0.7*(1-exp(-0.2))+2*0.3*(1-exp(-0.1)))/3] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 5")
    global error1 |= true
    global ind = 5
end

if error1
    throw("error Case F in quartet $(ind)")
end


# Case I Bad Diamond II -----------------

#println("------ Case I Bad diamond II ----------")
include("../examples/case_i_example.jl");
# include(joinpath(dirname(pathof(PhyloNetworks)),  "..","examples","case_i_example.jl"))
error1 = false
ind = 0

q1 = Quartet(1,["6","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q1);
try
    identifyQuartet!(qnet)
    qnet.which != 2 ? error("qnet which not correctly assigned") : nothing
    qnet.hybrid[1].k != 4 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
    qnet.hybrid[1].typeHyb != 5 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
    !isa(qnet.hybrid[1].prev,Nothing) ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 1 || qnet.numHybrids != 1 ? error("qnet should not have hybrid nodes anymore") : nothing
    qnet.t1 != -1 ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [-1,-1,-1,-1] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [-1,-1,-1] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [0.9*(1/3*exp(-1))+0.1*(1-2/3*exp(-1)),0.9*(1-2/3*exp(-1))+0.1*1/3*exp(-1),0.9*(1/3*exp(-1))+0.1*(1/3*exp(-1))] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 1")
    global error1 |= true
    global ind = 1
end


q2 = Quartet(2,["6","7","10","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q2);
try
    identifyQuartet!(qnet)
    qnet.which != 1 ? error("qnet which not correctly assigned") : nothing
    qnet.hybrid[1].k != 3 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
    qnet.hybrid[1].typeHyb != 4 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
    (qnet.hybrid[1].prev.number != -2 && qnet.hybrid[1].prev.number != -3) ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 0 || qnet.numHybrids != 0 ? error("qnet should have 0 hybrid nodes") : nothing
    !approxEq(qnet.t1,-log(1+0.1*(1-exp(-1))-0.1*0.1*(1-exp(-1-1))-0.1*0.1*(1-exp(-1))-0.9*0.9*(1-exp(-2.)))) ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [1,1,2,2] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [1,2,2] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [1-2/3*(exp(-qnet.t1)),1/3*(exp(-qnet.t1)),1/3*(exp(-qnet.t1))] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 2")
    global error1 |= true
    global ind = 2
end

q3 = Quartet(3,["10","7","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q3);
try
    identifyQuartet!(qnet)
    qnet.which != 1 ? error("qnet which not correctly assigned") : nothing
    qnet.hybrid[1].k != 3 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
    qnet.hybrid[1].typeHyb != 4 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
    (qnet.hybrid[1].prev.number != -6 && qnet.hybrid[1].prev.number != -3) ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 0 || qnet.numHybrids != 0 ? error("qnet should have 0 hybrid nodes") : nothing
    !approxEq(qnet.t1,-log(1+0.1*(1-exp(-2))-0.1*0.1*(1-exp(-1))-0.1*0.1*(1-exp(-2))-0.9*0.9*(1-exp(-2.)))) ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [1,1,2,2] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [2,2,1] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [1/3*exp(-qnet.t1),1/3*exp(-qnet.t1),1-2/3*exp(-qnet.t1)] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 3")
    global error1 |= true
    global ind = 3
end


q4 = Quartet(4,["6","10","4","8"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q4);

try
    identifyQuartet!(qnet)
    qnet.which != 1 ? error("qnet which not correctly assigned") : nothing
    qnet.hybrid[1].k != 3 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
    qnet.hybrid[1].typeHyb != 4 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
    (qnet.hybrid[1].prev.number != -6 && qnet.hybrid[1].prev.number != -3) ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 0 || qnet.numHybrids != 0 ? error("qnet should have 0 hybrid nodes") : nothing
    !approxEq(qnet.t1,-log(1+0.1*(1-exp(-1))-0.1*0.1*(1-exp(-1))-0.1*0.1*(1-exp(-1))-0.9*0.9*(1-exp(-3)))) ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [1,1,2,2] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [2,1,2] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [1/3*exp(-qnet.t1),1-2/3*exp(-qnet.t1),1/3*exp(-qnet.t1)] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 4")
    global error1 |= true
    global ind = 4
end


q5 = Quartet(5,["6","7","4","10"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net,q5);
try
    identifyQuartet!(qnet)
    qnet.which != 2 ? error("qnet which not correctly assigned") : nothing
    qnet.hybrid[1].k != 4 ? error("qnet.hybrid[1].k not correctly assigned") : nothing
    qnet.hybrid[1].typeHyb != 5 ? error("qnet.hybrid[1].typeHyb not correctly assigned") : nothing
    !isa(qnet.hybrid[1].prev,Nothing) ? error("qnet.hybrid[1].prev not correctly assigned") : nothing

    eliminateHybridization!(qnet)
    size(qnet.hybrid,1) != 1 || qnet.numHybrids != 1 ? error("qnet should not have hybrid nodes anymore") : nothing
    qnet.t1 != -1 ? error("internal edge length not correctly updated") : nothing

    updateSplit!(qnet)
    qnet.split != [-1,-1,-1,-1] ? error("qnet.split not correctly assigned") : nothing

    updateFormula!(qnet)
    qnet.formula != [-1,-1,-1] ? error("qnet.formula not correctly assigned") : nothing

    calculateExpCF!(qnet)
    qnet.expCF != [0.9*(1/3*exp(-1))+0.1*(1-2/3*exp(-1)),0.9*(1-2/3*exp(-1))+0.1*1/3*exp(-1),0.9*(1/3*exp(-1))+0.1*(1/3*exp(-1))] ? error("qnet.expCF wrongly calculated") : nothing
catch
    println("errors in quartet 5")
    global error1 |= true
    global ind = 5
end

if error1
    throw("error Case I in quartet $(ind)")
end

