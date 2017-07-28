# fixit: it seems that we need to do this to get plots to work,
#        and this exact order.
# using RCall
# using PhyloNetworks

info("checking RCall-based plots")
# checking that for no errors, not testing for correctness

#net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
net = readTopology("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);");
plot(net,:RCall);
plot(net,:RCall, useEdgeLength=true);
println("(A message is normal here: about missing BL being given length 1.0)")
plot(net,:RCall, showTipLabel=false);
plot(net,:RCall, showNodeNumber=true, showIntNodeLabel=true);
plot(net,:RCall, tipOffset=1, showGamma=true);
plot(net,:RCall, showEdgeLength=true, showEdgeNumber=true);
plot(net,:RCall, edgeColor="tomato4",minorHybridEdgeColor="skyblue",
          majorHybridEdgeColor="tan");
dat = DataFrame(node=[-5,-10,-1],bs=["90","95","100"],edge=[11,22,26]);
plot(net,:RCall, nodeLabel=dat);
plot(net,:RCall, edgeLabel=dat[:,[:edge,:bs]]);

# plot based on RCall and ape:
tre = readTopology("(((((((1,2),3),4),5),(6,7)),(8,9)),10);");
# fixit: plot(tre, :ape)
