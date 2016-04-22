# test the function to do bootstrap on snaq
# Claudia October 2015

include("../src/types.jl")
include("../src/functions.jl")

if !isdefined(:doalltests) doalltests = false; end
# below: to run the tests from anywhere using "include(path/test_bootsnaq.jl")
scriptfile = @__FILE__
exdir = (scriptfile==nothing ? joinpath("..","examples") :
        joinpath(dirname(dirname(@__FILE__)), "examples"))

info("Testing that bootsnaq runs with no error")
T=readStartTop(joinpath(exdir,"startTree.txt"),d);
df=readtable(joinpath(exdir,"tableCFCI.txt"),separator=';')
if (doalltests)
  d=readTableCF(joinpath(exdir,"tableCF.txt"));
  net1=snaq!(T,d,filename="net1_snaq");
  bootnet = bootsnaq(T,df,nrep=2,bestNet=net1);
end
bootnet = bootsnaq(T,df,nrep=2);
