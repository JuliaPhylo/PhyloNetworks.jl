# functions to test the functions for multiple alleles per species
# Claudia October 2015 - Cecile July 2016

if !isdefined(:individualtest) individualtest = false; end

if(individualtest)
    include("../src/types.jl")
    include("../src/functions.jl")
end

#--------------------------------------------------------------#
# old code from October 2015: please adapt as you wish Claudia #
#--------------------------------------------------------------#

if (false)
treefile = "(6,(5,(7,(3,4))));"
tree = readTopologyUpdate(treefile);
printEdges(tree)
printNodes(tree)
#repSpecies=["7"]
#expandLeaves!(repSpecies,tree)
#writeTopology(tree)

df = readtable("CFtable1.csv")
alleleDF=DataFrame(allele=["1","2"],species=["7","7"])
newdf = mapAllelesCFtable!(alleleDF,df,true,"CFmapped.csv")
mapD= readTableCF("CFmapped.csv")
end

#----------------------------------------------------------#
#   testing sorting of taxa and CFs                        #
#----------------------------------------------------------#
info("testing sorttaxa!")

letters = ["a","b","c","d"]; cfvalues = [0.6, 0.39, 0.01] # for ab_cd, ac_bd, ad_bc
d = DataFrame(t1=Array{ASCIIString}(24),t2=Array{ASCIIString}(24),t3=Array{ASCIIString}(24),t4=Array{ASCIIString}(24),
              CF12_34=Array{Float64}(24), CF13_24=Array{Float64}(24), CF14_23=Array{Float64}(24));
irow=1        # d will contain 6!=24 rows: for all permutations on 4 letters
for i1 in 1:4
  ind234 = deleteat!(collect(1:4),i1)
  for i2 in ind234
    ind34 = deepcopy(ind234)
    deleteat!(ind34,findfirst(ind34, i2))
    for j in 1:2
      i3=ind34[j]; i4=ind34[3-j]
      d[:t1][irow]=letters[i1]; d[:t2][irow]=letters[i2]; d[:t3][irow]=letters[i3]; d[:t4][irow]=letters[i4]
      # CF12_34 corresponds to CFi1i2_i3i4
      if     (i1,i2)∈[(1,2),(2,1),(3,4),(4,3)] d[:CF12_34][irow] = cfvalues[1]
      elseif (i1,i2)∈[(1,3),(3,1),(2,4),(4,2)] d[:CF12_34][irow] = cfvalues[2]
      elseif (i1,i2)∈[(1,4),(4,1),(2,3),(3,2)] d[:CF12_34][irow] = cfvalues[3]
      end # next: set CF13_24
      if     (i1,i3)∈[(1,2),(2,1),(3,4),(4,3)] d[:CF13_24][irow] = cfvalues[1]
      elseif (i1,i3)∈[(1,3),(3,1),(2,4),(4,2)] d[:CF13_24][irow] = cfvalues[2]
      elseif (i1,i3)∈[(1,4),(4,1),(2,3),(3,2)] d[:CF13_24][irow] = cfvalues[3]
      end # nest: set CF14_23
      if     (i1,i4)∈[(1,2),(2,1),(3,4),(4,3)] d[:CF14_23][irow] = cfvalues[1]
      elseif (i1,i4)∈[(1,3),(3,1),(2,4),(4,2)] d[:CF14_23][irow] = cfvalues[2]
      elseif (i1,i4)∈[(1,4),(4,1),(2,3),(3,2)] d[:CF14_23][irow] = cfvalues[3]
      end
      irow += 1
    end
  end
end
# d
d2 = deepcopy(d);
sorttaxa!(d2);
d3 = DataFrame(t1=repeat([letters[1]],outer=[24]),t2=repeat([letters[2]],outer=[24]),
               t3=repeat([letters[3]],outer=[24]),t4=repeat([letters[4]],outer=[24]),
               CF12_34=repeat([cfvalues[1]],outer=[24]),CF13_24=repeat([cfvalues[2]],outer=[24]),CF14_23=repeat([cfvalues[3]],outer=[24]));
@test d2==d3

dat = readTableCF(d);
net = readTopologyLevel1("(a,((b)#H1,((#H1,c),d)));");
topologyQPseudolik!(net, dat);
println("these 4 warnings are normal: only 4 taxa in network")
sorttaxa!(dat)

@test [q.obsCF for q in dat.quartet] == [[0.6,0.39,0.01] for i in 1:24]
@test [q.qnet.expCF for q in dat.quartet] == [[0.6915349833361827,0.12262648039048075,0.1858385362733365] for i in 1:24]
@test [q.taxon for q in dat.quartet] == [letters for i in 1:24]
@test [q.qnet.quartetTaxon for q in dat.quartet] == [letters for i in 1:24]
