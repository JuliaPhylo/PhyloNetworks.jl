# Tests for function simulate
using PhyloNetworks
using Base.Test

## Get a network and make it ultrametric
net = readTopology("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);");

for i = 1:27
	  net.edge[i].length = 1;
end
net.edge[27].length = 7;
net.edge[26].length = 3;
net.edge[25].length = 4;
net.edge[24].length = 4;
net.edge[21].length = 5;
net.edge[19].length = 4;
net.edge[16].length = 2;
net.edge[8].length = 2;
net.edge[3].length = 2;
net.edge[1].length = 5;

## Simulate a BM
srand(17920921); # fix the seed
pars = paramsBM(1, 0.1); # params of a BM
@show pars

sim = simulate(net, pars); # simulate according to a BM
@show sim

# Extract simulated values
traitsTips = sim[:Tips];
traitsNodes = sim[:InternalNodes];

# Expected values
traitsTipsExp = [0.6455995230091043 -0.22588106270381064 0.05703904710270408 -0.692650796714688 1.578622599565194 1.4106438068675058 1.9166557600811194 1.0579005662214953 1.2340762902144904 1.4130757789427886 0.7115737497673081 2.201943319276716];

traitsNodesExp = [-0.3481603206484607 -0.6698437934551933 -0.018135478212541654 -0.33844527112230455 -0.0717742134084467 0.19417331380691694 1.3919535151447147 1.5106942025265466 1.2526948727806593 1.1552248152172964 1.224823113083187 1.0617270280846993 1.0436547766241817 1.0];

@test_approx_eq traitsTips traitsTipsExp
@test_approx_eq traitsNodes traitsNodesExp


###############################################################################
## Test of distibution
###############################################################################

## Generate some values
srand(18480224); # fix the seed
pars = paramsBM(1, 0.1); # params of a BM
N = 50000
S = length(tipLabels(net));
values = zeros(Float64, (S, N));
for i = 1:N
	values[:,i] = simulate(net, pars)[:Tips]
end

## Check that each tip has same mean (1)
for s in 1:S
	@test_approx_eq_eps mean(values[s, :]) pars.mu 1e-2
end

## Check for variances
V = sharedPathMatrix(net);
Sig = V[:Tips] * pars.sigma2;
for s in 1:S
	for t in s:S
		@test_approx_eq_eps cov(values[s, :], values[t,:]) Sig[s, t] 1e-2 
	end
end


