# things to do:


1. finished debugging

2. do a checkQnet to check an extracted quartet that it does not have
redundant cycles, nor internal nodes flying

3. improve tests: add runtest.jl that runs all tests, using checkNet
and checkQnet

3. add necessary files for a package, read how to convert to package (link)

4. include documentation: by Docile package, pdf and iJulia notebook maybe

###IMPROVEMENTS 

I had two ideas in mind: 

1) help develop more functions for the julia package: extract the
trees compatible with a network, extract the main tree, distance
between networks, are two networks the same, etc.  This will be
related to what John is doing, but it would be a matter of organizing
who will do what.

2) help improve the performance of the existing code of estimation: in
particular, there are three areas:

2.1) help implement a local optimization (subset of parameters) after
a move instead of a global optimization (all the parameters). also
help with bad diamond I, existing algorithm not efficient

2.2) help to implement more network moves: I still need the "SPR"
move, delete/add at random in strange place, and there are other
interesting moves that we could consider in the paper I just sent
you. I have not read carefully, but they claim that the moves the
propose would help search the space of networks more efficiently.

2.3) the paralellization


###VERSIONS

Pkg.tag(...)
v0.0.1: September 22,2015