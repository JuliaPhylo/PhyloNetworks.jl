# Candidate Networks

## Optimizing parameters for a given network

For a given network topology, we can optimize the branch lengths and
inheritance probabilities (γ) with the pseudolikelihood.
This is useful if we have a few candidate networks to compare.
Each network can be optimized individually, and the network with the best
pseudolikelihood can be chosen.

The score being optimized is the pseudo-deviance, i.e.
the negative log pseudo-likelihood up to an additive constant
(the lower the better).

Following our example in [Getting a Network](@ref),
we can optimize parameters on the true network
(the one originally used to simulate the data):

```@setup fixednetworkoptim
using PhyloNetworks
using Logging # to suppress info messages below
baselogger = global_logger()
mkpath("../assets/figures")
raxmltrees = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","raxmltrees.tre")
raxmlCF = readTrees2CF(raxmltrees, writeTab=false, writeSummary=false)
```

```@repl fixednetworkoptim
truenet = readTopology("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
net1alt = topologyMaxQPseudolik!(truenet, raxmlCF);
writeTopology(net1alt, round=true)
net1alt.loglik # pseudo deviance, actually
```
```@example fixednetworkoptim
using PhyloPlots, RCall
R"name <- function(x) file.path('..', 'assets', 'figures', x)" 
R"svg(name('truenet_opt.svg'), width=4, height=4)" 
R"par"(mar=[0,0,0,0])
plot(net1alt, :R, showGamma=true);
R"dev.off()" 
nothing # hide
```
![truenet_opt](../assets/figures/truenet_opt.svg)

We get a score of 29.941,
which is comparable to the score of the SNaQ network (net1: 28.315),
especially compared to the score of the best tree (net0: 53.532).
This begs the question: is the true network within the "range" of uncertainty?
We can run a [Bootstrap](@ref) analysis to measure uncertainty
in our network inference.

For a more thorough optimization, we may increase the requirements before
the search stops (but the optimization will take longer).
It makes no difference on this small data set.
```julia
net1par = topologyMaxQPseudolik!(truenet, raxmlCF, ftolRel=1e-10, xtolAbs=1e-10)
net1par.loglik
```

## Network Score with no optimization

For a network with given branch lengths and γ heritabilies,
we can compute the pseudolikelihood with:
```@repl fixednetworkoptim
topologyQPseudolik!(truenet,raxmlCF);
truenet.loglik
```
This function is not maximizing the pseudolikelihood, it is simply computing the
pseudolikelihood (or deviance) for the given branch lengths and probabilities of
inheritance. At the moment, both of these functions require that the
given network is of level 1 (cycles don't overlap).

## Candidate networks compatible with a known outgroup

If the network was estimated via `snaq!`, it might turn out to be impossible
to root our estimated network with a known outgroup (see section
[What if the root conflicts with the direction of a reticulation?](@ref).)
At this time, `snaq!` does not impose any rooting constraint on the network:
the search for the lowest score considers all level-1 networks, including those
that are incompatible with a known outgroup.
(The monophyly of outgroups is not imposed either, like in many other methods.)

If the estimated network cannot be rooted with the known outgroup,
we can check the `.networks` output file.
It has a list of networks that are slight modifications of the best network,
where the modifications changed the direction of one reticulation at a time.
For each modified network, the score was calculated. So if we find in this list
a modified network that has a score close to that of the best network,
and that can be re-rooted with our known root position, then this modified network
is a better candidate than the network with the best score.

Below is what the `net1.networks` file looks like, after performing
the analysis in the section [Network Estimation](@ref).
Scroll to the right to see the scores.

    (C,D,((O,(E,#H7:::0.19558838614943078):0.31352437658618976):0.6640664399202987,(B,(A)#H7:::0.8044116138505693):10.0):10.0);, with -loglik 28.31506721890958 (best network found, remaining sorted by log-pseudolik; the smaller, the better)
    (C,D,((O,(E)#H7:::0.8150784689693145):0.9336405757682176,(B,(A,#H7:::0.18492153103068557):0.25386142779877724):1.8758156446611114):10.0);, with -loglik 31.535560380783814
    (B,#H7:9.90999345612101::0.2555404440833535,(A,(E,(O,((C,D):10.0)#H7:0.3419231810962026::0.7444595559166465):0.19994859441332047):2.5014911511063644):0.7957621793330066);, with -loglik 56.64548310161462
    (C,D,((O,(E,((B)#H7:::0.7957543284159452,A):4.786202415937916):0.004527712280136759):1.7952610454570868,#H7:::0.20424567158405482):10.0);, with -loglik 67.17775727492258
    (C,D,(#H7:::0.32947301811471164,(B,(A,(E,(O)#H7:::0.6705269818852884):1.371799259141243):0.0):6.397073999864152):7.677245926003807);, with -loglik 199.11401961057143

We can read this file and look at its list of networks like this:

```@repl fixednetworkoptim
file = "net1.networks";
# or use the example file available with the package:
file = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","net1.networks");
netlist = readMultiTopology(file) # read the full list of networks in that file
```
Next, we would like to extract the network scores from the file.
Below is a one-liner to do this
(we make Julia send a `sed` command to the shell --sorry, Mac or Linux for this.)
```@repl fixednetworkoptim
scoresInString = read(`sed -E 's/.+with -loglik ([0-9]+.[0-9]+).+/\1/' $file`, String)
scores = parse.(Float64, split(scoresInString))
# next: update the "loglik" of each network with the score read from the file
for i in eachindex(netlist)
   netlist[i].loglik = scores[i]
   println("net $i in the list: score = ",scores[i])
end
```
The first network in the list is the best network returned by `snaq!`.
We see that the second network has a score that's not too far, but the other networks
have worse scores. The best network and its best modification (second network in the
list) are shown below. We chose to show edge numbers, to use them later
to re-root the networks.

```@example fixednetworkoptim
R"svg(name('fixednetworkoptim_othernets1.svg'), width=7, height=4)" # hide
R"layout(matrix(1:2,1,2))"; # hide
R"par"(mar=[0,0,0,0]) # hide
plot(netlist[1], :R, showGamma=true, showEdgeNumber=true, tipOffset=0.1);
R"mtext"("best net, score=28.3", line=-1);
plot(netlist[2], :R, showGamma=true, showEdgeNumber=true, tipOffset=0.1);
R"mtext"("direction modified, score=31.5", line=-1);
R"dev.off()"; # hide
```
![othernets before reroot](../assets/figures/fixednetworkoptim_othernets1.svg)

Now imagine that our outgroup is taxon A.
- best network: we would get a "RootMismatch" error if we tried to set
  the root on the external edge 9 to A, with `rootatnode!(netlist[1], "A")`
  (see section
  [What if the root conflicts with the direction of a reticulation?](@ref)).
  But we could root the best network on the major parent edge to A, edge 10
  (rooted network on the left below).
- For the second best network in our list, there are 2 ways to root it
  with A: on the external edge 8 to A (top right), or on its parent edge 10
  (bottom right). These 2 options give quite different rooted versions
  of the network, one of which requires the existence of an unsampled taxon,
  sister to BOECD, that would have contributed to introgression into
  an ancestor of E. The second rooted version says that an ancestor of
  (or sister to) A contributed to the introgression into the ancestor of E.
  A is an outgroup in both cases, but the second case is more parsimonious,
  in the sense that it does not require the existence of an unsampled taxon.

```@example fixednetworkoptim
R"svg(name('fixednetworkoptim_othernets2.svg'), width=7, height=7)" # hide
R"layout(matrix(c(1,4,2,3),2,2))"; # hide
R"par"(mar=[0,0,0.5,0]) # hide
rootonedge!(netlist[1], 10); # root best net to make A outgroup
rotate!(netlist[1], -4); # to 'un-cross' edges
rotate!(netlist[1], -6);
plot(netlist[1], :R, showGamma=true, tipOffset=0.1);
R"mtext"("best net, score=28.3", line=-1);
global_logger(NullLogger()); # hide
rootatnode!(netlist[2], "A"); # net with modified direction: first way to make A outgroup
global_logger(baselogger);   # hide
plot(netlist[2], :R, showGamma=true, tipOffset=0.1);
R"mtext"("second best in list, score=31.5\nrequires unsampled population", line=-2);
rootonedge!(netlist[2], 10) # net with modified direction: second way to make A outgroup
plot(netlist[2], :R, showGamma=true, tipOffset=0.1);
R"mtext"("second best in list, score=31.5\ndifferent root position", line=-2);
R"dev.off()"; # hide
```
![othernets after reroot](../assets/figures/fixednetworkoptim_othernets2.svg)
