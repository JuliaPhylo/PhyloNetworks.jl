@testset "multiple alleles" begin
global tree, df, d, net, currT

@testset "test: map alleles to species" begin
    tree = readTopology("(6,(5,(7,(3,4))));");
    PhyloNetworks.expandLeaves!(["7"],tree)
    @test writeTopology(tree) == "(6,(5,((7:0.0,7__2:0.0):1.0,(3,4))));"
    PhyloNetworks.mergeLeaves!(tree)
    @test writeTopology(tree) == "(6,(5,(7:1.0,(3,4))));"
    alleleDF=DataFrame(allele=["1","2"], species=["7","7"])
    # df = CSV.read(joinpath(@__DIR__, "..", "examples", "tableCFCI.csv"))
    # PhyloNetworks.mapAllelesCFtable!(df,alleleDF,[1,2,3,4],true,"CFmapped.csv")
    CSV.write("tmp.csv", alleleDF);
    df = (@test_logs (:warn, r"^not all alleles were mapped") mapAllelesCFtable("tmp.csv",
      joinpath(@__DIR__, "..", "examples", "tableCFCI.csv"),
      # joinpath(dirname(pathof(PhyloNetworks)), "..", "examples", "tableCFCI.csv"),
      filename="CFmapped.csv"))
    rm("CFmapped.csv")
    rm("tmp.csv")
    @test df[!,:t4] == ["4","7","3","7","3","3","7","3","3","3","7","3","3","3","3"]
end

#----------------------------------------------------------#
#   testing sorting of taxa and CFs                        #
#----------------------------------------------------------#
@testset "sorttaxa!" begin

letters = ["a","b","c","d"]; cfvalues = [0.6, 0.39, 0.01] # for ab_cd, ac_bd, ad_bc
d = DataFrame(t1=Array{String}(undef,24),t2=Array{String}(undef,24),t3=Array{String}(undef,24),t4=Array{String}(undef,24),
              CF12_34=Array{Float64}(undef,24), CF13_24=Array{Float64}(undef,24), CF14_23=Array{Float64}(undef,24));
irow=1        # d will contain 6!=24 rows: for all permutations on 4 letters
for i1 in 1:4
  ind234 = deleteat!(collect(1:4),i1)
  for i2 in ind234
    ind34 = deepcopy(ind234)
    deleteat!(ind34, findfirst(isequal(i2), ind34))
    for j in 1:2
      i3=ind34[j]; i4=ind34[3-j]
      d[irow,:t1]=letters[i1]; d[irow,:t2]=letters[i2]; d[irow,:t3]=letters[i3]; d[irow,:t4]=letters[i4]
      # CF12_34 corresponds to CFi1i2_i3i4
      if     (i1,i2)∈[(1,2),(2,1),(3,4),(4,3)] d[irow,:CF12_34] = cfvalues[1]
      elseif (i1,i2)∈[(1,3),(3,1),(2,4),(4,2)] d[irow,:CF12_34] = cfvalues[2]
      elseif (i1,i2)∈[(1,4),(4,1),(2,3),(3,2)] d[irow,:CF12_34] = cfvalues[3]
      end # next: set CF13_24
      if     (i1,i3)∈[(1,2),(2,1),(3,4),(4,3)] d[irow,:CF13_24] = cfvalues[1]
      elseif (i1,i3)∈[(1,3),(3,1),(2,4),(4,2)] d[irow,:CF13_24] = cfvalues[2]
      elseif (i1,i3)∈[(1,4),(4,1),(2,3),(3,2)] d[irow,:CF13_24] = cfvalues[3]
      end # nest: set CF14_23
      if     (i1,i4)∈[(1,2),(2,1),(3,4),(4,3)] d[irow,:CF14_23] = cfvalues[1]
      elseif (i1,i4)∈[(1,3),(3,1),(2,4),(4,2)] d[irow,:CF14_23] = cfvalues[2]
      elseif (i1,i4)∈[(1,4),(4,1),(2,3),(3,2)] d[irow,:CF14_23] = cfvalues[3]
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
net = (@test_logs readTopologyLevel1("(a,((b)#H1,((#H1,c),d)));"));
# earlier warning: "net does not have identifiable branch lengths"
@test_logs topologyQPseudolik!(net, dat);
sorttaxa!(dat)

@test [q.obsCF for q in dat.quartet] == [[0.6,0.39,0.01] for i in 1:24]
@test [q.qnet.expCF for q in dat.quartet] == [[0.6915349833361827,0.12262648039048075,0.1858385362733365] for i in 1:24]
@test [q.taxon for q in dat.quartet] == [letters for i in 1:24]
@test [q.qnet.quartetTaxon for q in dat.quartet] == [letters for i in 1:24]

end # of testset: sorttaxa!

@testset "snaq on multiple alleles" begin
df=DataFrame(t1=["6","6","10","6","6","7","7","7","7","7","7"],
             t2=["7","7","7","10","7","7","7","7","7","7","7"],
             t3=["4","10","4","4","4","8","8","8","10","10","6"],
             t4=["8","8","8","8","10","10","4","6","4","6","4"],
             CF1234=[0.2729102510259939, 0.3967750546426937, 0.30161247267865315, 0.24693940689390592, 0.2729102510259939, 0.155181,  0.792153,  0.486702,  0.962734,  0.202531,  0.486886],
             CF1324=[0.45417949794801216, 0.30161247267865315, 0.30161247267865315, 0.5061211862121882, 0.45417949794801216, 0.673426 ,0.145408,  0.391103,  0.023078,  0.714826,  0.419015],
             CF1423=[0.2729102510259939, 0.30161247267865315, 0.3967750546426937, 0.24693940689390592, 0.2729102510259939, 0.171393,  0.062439,  0.122195,  0.014188,  0.082643,  0.094099])
d = readTableCF(df)
@test !isempty(d.repSpecies)
@test d.repSpecies == ["7"]

tree = "((6,4),(7,8),10);"
currT = readTopology(tree);

originalstdout = stdout
redirect_stdout(open("/dev/null", "w")) # not portable to Windows
estNet = snaq!(currT,d,hmax=1,seed=1010, runs=1, filename="", Nfail=10)
redirect_stdout(originalstdout)
@test 185.27 < estNet.loglik < 185.29 # or: wrong loglik
@test estNet.hybrid[1].k == 4 # or: wrong k
@test estNet.numTaxa == 5 # or: wrong number of taxa

redirect_stdout(open("/dev/null", "w")) # not portable to Windows
estNet = snaq!(currT,d,hmax=1,seed=8378, runs=1, filename="", Nfail=10,
               ftolAbs=1e-6,ftolRel=1e-5,xtolAbs=1e-4,xtolRel=1e-3)
redirect_stdout(originalstdout)
@test 174.58 < estNet.loglik < 174.59 # or: loglik wrong
@test estNet.hybrid[1].k == 5 # or: wrong k in hybrid
@test estNet.numTaxa == 5 # or: wrong # taxa

# net = snaq!(currT,d,hmax=1,seed=8378,filename="")
net = readTopology("(((4,#H7:::0.47411636966376686):0.6360197250223204,10):0.09464128563363322,(7:0.0,(6)#H7:::0.5258836303362331):0.36355727108454877,8);")
@test topologyQPseudolik!(net, d) ≈ 174.58674796123705
@test net.loglik ≈ 174.58674796123705
net = readTopology("(((4,#H1),10),(7,(6)#H1),8);")
net = topologyMaxQPseudolik!(net,d,  # loose tolerance for faster test
        ftolRel=1e-2,ftolAbs=1e-2,xtolAbs=1e-2,xtolRel=1e-2)
@test net.loglik > 174.5
end # test of snaq on multiple alleles

end # overall multiple allele sets of testests
