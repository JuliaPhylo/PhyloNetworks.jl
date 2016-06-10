# Fixed Network

## Optimizing branch lengths and inheritance probabilities for a given network

For a given network topology, we can optimize the branch lengths and
inheritance probabilities (γ) with the pseudolikelihood.
This is useful if we have a few candidate networks to compare.
Each network can be optimized individually, and the network with the best
pseudolikelihood can be chosen.

The score being optimized is the pseudo-deviance, i.e.
the negative log pseudo-likelihood up to an additive constant,
such that a perfect fit corresponds to a deviance of 0.0 (the lower the better).
```julia
net1topo = readTopology("(2,(4,(3,(5,(6,#H1)))),(1)#H1);");
net1par = topologyMaxQPseudolik!(net1topo,d)
net1par.loglik # pseudo deviance, actually
```
For a more thorough optimization, we may increase the requirements before
the search stops:
```julia
net1par = topologyMaxQPseudolik!(net1topo,d, xtolRel=1e-10, xtolAbs=1e-10)
net1par.loglik
```
## Network Score with no optimization

For a network with given branch lengths and γ heritabilies,
we can compute the pseudolikelihood with:
```julia
net1withBL = readTopology("(2,(4,(3,(5,(6,#H6:1.0::0.3):5.006):0.518):0.491):1.533,(1)#H6:1.0::0.7);");
topologyQPseudolik!(net1withBL,d)
net1withBL.loglik
```
This function is not maximizing the pseudolikelihood, it is simply computing the
pseudolikelihood (or deviance) for the given branch lengths and probabilities of
inheritance. At the moment, both of these functions require that the
given network is of level 1 (cycles don't overlap).

