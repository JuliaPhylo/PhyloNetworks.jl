"""
    readmultinewick_files(listfile; relative2listfile=true)

Read the list of file names in `listfile`, then read all the trees in each of
these files. Output: vector of vectors of trees (networks with h>0 allowed).

`listfile` should be the name of a file containing the path/name to multiple
bootstrap files, one on each line (no header). Each named bootstrap file should
contain multiple trees, one per line (such as bootstrap trees from a single gene).

The path/name to each bootstrap file should be relative to `listfile`.
Otherwise, use option `relative2listfile=false`, in which case the file names
are interpreted as usual: relative to the user's current directory
if not given as absolute paths.
"""
function readmultinewick_files(filelist::AbstractString; relative2listfile::Bool=true)
    filelistdir = dirname(filelist)
    bootfiles = DataFrame(CSV.File(filelist, header=false, types=Dict(1=>String));
        copycols=false)
    size(bootfiles)[2] > 0 ||
        error("there should be a column in file $filelist: with a single bootstrap file name on each row (no header)")
    ngenes = size(bootfiles)[1]
    bf = (relative2listfile ? joinpath.(filelistdir, bootfiles[!,1]) : bootfiles[!,1])
    treelists = Array{Vector{HybridNetwork}}(undef, ngenes)
    for igene in 1:ngenes
        treelists[igene] = readmultinewick(bf[igene])
        print("read $igene/$ngenes bootstrap tree files\r") # using \r for better progress display
    end
    return treelists
end

"""
    samplebootstrap_multiloci(vector of tree lists; seed=0, generesampling=false, row=0)
    samplebootstrap_multiloci!(tree list, vector of tree lists; seed=0, generesampling=false, row=0)

Sample bootstrap gene trees, 1 tree per gene.
Set the seed with keyword argument `seed`, which is 0 by default.
When `seed=0`, the actual seed is set using the clock.
Assumes a vector of vectors of networks (see `readmultinewick_files`),
each one of length 1 or more (error if one vector is empty, tested in `bootsnaq`).

- site resampling: always, from sampling one bootstrap tree from each given list.
  This tree is sampled at **random** unless `row>0` (see below).
- gene resampling: if `generesampling=true` (default is false),
  genes (i.e. lists) are sampled with replacement.
- `row=i`: samples the ith bootstrap tree for each gene.
  `row` is turned back to 0 if gene resampling is true.

output: one vector of trees. the modifying function (!) modifies the input tree list and returns it.
"""
function samplebootstrap_multiloci(
    trees::Vector{Vector{HybridNetwork}};
    seed::Integer=0,
    generesampling::Bool=false,
    row::Integer=0
)
    bootTrees = Array{HybridNetwork}(undef, length(trees))
    samplebootstrap_multiloci!(bootTrees, trees, seed=seed, generesampling=generesampling, row=row)
end

function samplebootstrap_multiloci!(
    bootTrees::Vector{HybridNetwork},
    trees::Vector{Vector{HybridNetwork}};
    seed::Integer=0,
    generesampling::Bool=false,
    row::Integer=0
)
    numgen = length(trees) ## number of genes
    numgen>0 || error("needs at least 1 array of trees")
    numgen <= length(bootTrees) || error("the input tree list needs to be of length $numgen at least")
    if (generesampling) row=0; end
    if row==0
      if seed == 0
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) #better seed based on clock
        println("using seed $(seed) for bootstrap trees")
      end
      Random.seed!(seed)
      if generesampling
        indxg = sample(1:numgen, numgen) # default is with replacement. good!
      end
    end
    for g in 1:numgen
        ig = (generesampling ? indxg[g] : g )
        if row==0
          indxt = sample(1:length(trees[ig]),1)[1]
        else
          indxt = row
          length(trees[ig]) >= row || error("gene $g has fewer than $row bootstrap trees.")
        end
        bootTrees[g] = trees[ig][indxt]
    end
    return bootTrees
end


"""
    treeedges_support(sample_net::Vector{HybridNetwork}, ref_net::HybridNetwork)`

Read a sample of networks `sample_net` (such as a bootstrap sample or a
posterior sample) and a reference network (`ref_net`),
and calculate the support for the tree edges in the reference network.
All minor hybrid edges (γ<0.5) are removed to extract the major tree from
each network. All remaining edges are tree edges, each associated with a bipartition.

output:

- a data frame with one row per tree edge and two columns: edge number, bootstrap support
  (as a percentage)
- the major tree from the reference network, where minor hybrid edges (with γ<0.5)
  have been removed.
"""
function treeedges_support(net::Vector{HybridNetwork}, net0::HybridNetwork)
    # estimated network, major tree and matrix
    S = tiplabels(net0)
    tree0 = majortree(net0, unroot=true)
    M0 = tree2Matrix(tree0,S, rooted=false)

    M = Matrix[]
    tree = HybridNetwork[]
    for n in net
        t = majortree(n, unroot=true)
        push!(tree,t)
        mm = tree2Matrix(t,S, rooted=false)
        push!(M,mm)
    end

    df = DataFrame(edgeNumber=Int[], proportion=Float64[])

    for i in axes(M0, 1) # rows in M0: internal edges
        cnt = 0 # count
        for j in eachindex(M) # for every M
            for k in axes(M[j],1) # check every row in M
                if M0[i,2:end] == M[j][k,2:end] || M0[i,2:end] == map(x->(x+1)%2,M[j][k,2:end]) #same row
                    #println("found same row: $(M0[i,2:end]) and $(M[j][k,2:end])")
                    cnt += 1
                    break
                end
            end
        end
        push!(df,[M0[i,1] cnt*100/length(M)])
        # if edges had an attribute for bootstrap/posterior support, we could
        # set it to cnt/length(M) for the edge numbered M0[i,1].
    end
    @info """edge numbers in the data frame correspond to the current edge numbers in the network.
       If the network is modified, the edge numbers in the (modified) network might not correspond
       to those in the bootstrap table. Plot the bootstrap values onto the current network with
       plot(network_name, edgelabel=bootstrap_table_name)"""
    return df, tree0
end


"""
    hybridclades_support(boot_net::Vector{HybridNetwork}, ref_net::HybridNetwork; rooted=false)

Match hybrid nodes in a reference network with those in an array of networks,
like bootstrap networks or in a posterior sample of networks.
All networks must be fully resolved, and on the same taxon set.
If `rooted=true`, all networks are assumed to have been properly rooted beforehand.
Otherwise, the origin of each hybrid edge is considered as an unrooted bipartition (default).

Two hybrid edges in two networks are said to match if they share the same "hybrid" clade
(or recipient) and the same "donor clade", which is a sister to the hybrid clade in the network.
Since a hybrid clade has 2 parent edges, it is sister to two clades simultaneously: one is its
major sister (following the major hybrid edge with γ>0.5) and one is its minor sister (following
the major hybrid edge with γ<0.5).

To calculate these hybrid and sister clades at a given hybrid node, all other hybrid edges are
first removed from the network. Then, the hybrid clade is the hardwired cluster (descendants) of
either hybrid edge and major/minor clade is the hardwired cluster of the sibling edge of the
major/minor hybrid parent.
If `rooted=false`, sister clades are considered as bipartitions.

Output:

1. a "node" data frame (see below)
2. an "edge" data frame (see below)
3. a "clade" data frame to describe the make up of all clades found as hybrids or sisters,
   starting with a column `taxa` that lists all taxa. All other columns correspond to a given
   clade and contain true/false values. `true` means that a given taxon belongs in a given clade.
   For a clade named `H1`, for instance, and if the data frame was named `cla`, the
   list of taxa in this clade can be obtained with `cla[:taxa][cla[:H1]]`.
4. an array of gamma values, with one row for each network in the sample and
   two columns (major/minor) for each hybrid edge in the reference network.
   If this hybrid edge was found in the sampled network
   (i.e. same hybrid and sister clades, after removal of all other hybrid nodes),
   its inheritcance gamma value is recorded here. Otherwise, the gamma entry is 0.0.
5. a vector with the number of each hybrid edge in the reference network,
   in the same order as for the columns in the array of gamma values above.

The "node" data frame has one row per clade and 9 columns giving:

   - `:clade`: the clade's name, like the taxon name (if a hybrid is a single taxon) or
     the hybrid tag (like 'H1') in the reference network
   - `:node`: the node number in the reference network. missing if the clade is not in this network.
   - `:hybridnode`: typically the same node number as above, except for hybrid clades in the
     reference network. For those, the hybrid node number is listed here.
   - `:edge`: number of the parent edge, parent to the node in column 2,
     if found in the ref network. missing otherwise.
   - `:BS_hybrid`: percentage of sample networks in which the clade is found to
     be a hybrid clade.
   - `:BS_sister`: percentage of sample networks in which the clade is found to
     be sister to some hybrid clade (sum of the next 2 columns)
   - `:BS_major_sister`: percentage of sample networks in which the clade is
     found to be the major sister to some hybrid clade
   - `:BS_minor_sister`: same as previous, but minor
   - `:BS_hybrid_samesisters`: percentage of sample networks in which the clade
     is found to be a hybrid and with the same set of sister clades as in the
     reference network.
     Applies to hybrid clades found in the reference network only, missing for all other clades.

The "edge" data frame has one row for each pair of clades, and 8 columns:

  - `:edge`: hybrid edge number, if the edge appears in the reference network. missing otherwise.
  - `:hybrid_clade`: name of the clade found to be a hybrid, descendent of 'edge'
  - `:hybrid`: node number of that clade, if it appears in the reference network. missing otherwise.
  - `:sister_clade`: name of the clade that is sister to 'edge', i.e. be sister to a hybrid
  - `:sister`: node number of that clade, if in the ref network.
  - `:BS_hybrid_edge`: percentage of sample networks in which 'edge' is found to be a hybrid
     edge, i.e. when the clade in the 'hybrid' column is found to be a hybrid and the clade in
     the 'sister' column is one of its sisters.
  - `:BS_major`: percentage of sample networks in which 'edge' is found to be a major hybrid
     edge, i.e. when 'hybrid' is found to be a hybrid clade and 'sister' is found to be its
     major sister.
  - `:BS_minor`: same as previous, but minor
"""
function hybridclades_support(
    nets::Vector{HybridNetwork},
    refnet::HybridNetwork;
    rooted::Bool=false
)
    numNets = length(nets)
    numNets>0 || error("there aren't any networks in the sample")
    numHybs = refnet.numhybrids
    numHybs>0 || error("there aren't any hybrid in reference network")
    try directedges!(refnet)
    catch err
        if isa(err, RootMismatch)
            err.msg *= "\nPlease change the root in reference network (see rootatnode! or rootonedge!)"
        end
        rethrow(err)
    end
    taxa = tiplabels(refnet)
    ntax = length(taxa)

    # extract hardwired clusters of each tree edge in major tree of reference net,
    # and of each hybrid edge after all other hybrids are removed following the major edge.
    clade = Vector{Bool}[] # list all clades in reference and in bootstrap networks
    treenode = Int[]     # node number for each clade: of tree node if found in reference net
    treeedge = Int[]     # number of tree edge corresponding to the clade (parent of treenode)
    leafname = AbstractString[] # "" if internal node, leaf name if leaf
    # clade, treenode, treeedge: same size. leafname and hybparent: same size initially only
    hybind = Int[]       # indices in 'clade', that appear as hybrid clades in reference
    hybnode = Int[]      # for those hybrid clades: number of hybrid node
    majsisedge = Int[]   # for those hybrid clades: number of major sister edge (in ref net)
    minsisedge = Int[]   #                                    minor
    majsisind = Int[]    # for those hybrid clades: index of major sister clade in 'clade'
    minsisind = Int[]    #                                   minor
    # hybind, hybnode, majsis*, minsis*: same size = number of hybrid nodes in reference network

    if !rooted && length(refnet.node[refnet.rooti].edge)==2 && any(e -> e.hybrid, refnet.node[refnet.rooti].edge)
        refnet = deepcopy(refnet) # new binding inside function
        fuseedgesat!(refnet.rooti, refnet)
    end # issues otherwise: correct tree edge for root bipartitions, find sister clades, ...

    reftre = majortree(refnet, unroot=true)
    skipone = (!rooted && length(reftre.node[reftre.rooti].edge)<3) # not count same bipartition twice
    for pe in reftre.edge
        hwc = hardwiredcluster(pe,taxa) # not very efficient, but human readable
        if skipone && refnet.node[refnet.rooti] ≡ getparent(pe) && sum(hwc)>1
            skipone = false             # wrong algo for trivial 2-taxon rooted tree (A,B);
            println("skip edge $(pe.number)")
        else
            push!(clade, hwc)
            push!(treeedge, pe.number)
            cn = getchild(pe) # child node of pe
            push!(treenode, cn.number)
            push!(leafname, (cn.leaf ? cn.name : ""))
        end
    end
    hybparent = zeros(Int,length(clade)) # 0 if has no hybrid parent node. Index in hybnode otherwise.
    for trueh = 1:numHybs
        net0 = deepcopy(refnet)
        displayednetworkat!(net0, net0.hybrid[trueh]) # removes all minor hybrid edges but one
        hn = net0.hybrid[1]
        hemaj, hemin, ce = hybridEdges(hn) # assumes no polytomy at hybrid nodes and correct node.hybrid
        (hemin.hybrid && !hemin.ismajor) || error("edge should be hybrid and minor")
        (hemaj.hybrid &&  hemaj.ismajor) || error("edge should be hybrid and major")
        ic = findfirst(isequal(ce.number), treeedge)
        ic !== nothing || error("hybrid node $(hn.number): child edge not found in major tree")
        hybparent[ic] = trueh
        push!(hybind,ic)
        push!(hybnode, hn.number)
        push!(majsisedge, hemaj.number)
        push!(minsisedge, hemin.number)
        for sis in ["min","maj"]
          he = (sis=="min" ? hemin : hemaj)
          pn = getparent(he) # parent node of sister (origin of gene flow if minor)
          atroot = (!rooted && pn ≡ net0.node[net0.rooti]) # polytomy at root pn of degree 3: will exclude one child edge
          hwc = zeros(Bool,ntax) # new binding each time. pushed to clade below.
          for ce in pn.edge    # important if polytomy
            if ce ≢ he && pn ≡ getparent(ce)
                hw = hardwiredcluster(ce,taxa)
                if atroot && any(hw .& clade[ic]) # sister clade intersects child clade
                    (hw .& clade[ic]) == clade[ic] ||
                        @warn "weird clusters at the root in reference, hybrid node $(hn.number)"
                else
                    hwc .|= hw
                end
            end
          end
          i = findfirst(isequal(hwc), clade)
          if (!rooted && i===nothing) i = findfirst(isequal(.!hwc), clade) end
          i !== nothing || error(string("hyb node $(hn.number): ",sis,"or clade not found in main tree"))
          if (sis=="min") push!(minsisind, i)
          else            push!(majsisind, i)
          end
          if sis=="min"  # need to get clade not in main tree: hybrid + minor sister
            pe = nothing # looking for the (or one) parent edge of pn
            for ce in pn.edge
              if pn ≡ getchild(ce)
                  pe=ce
                  break
              end
            end
            # pe == nothing if minor hybrid is at root and rooted=true. Will just miss edge and node number
            # for that clade, but in that case, the minor sister might have been assigned that edge anyway...
            if pe != nothing
              hwc = hardwiredcluster(pe,taxa)
              i = findfirst(isequal(hwc), clade) # i>0: (hybrid + minor sister) can be in main tree if
              # hybrid origin is ancestral, i.e. hybrid clade is nested within minor sister.
              if i===nothing
                push!(clade, hwc)
                push!(treenode, pn.number)
                push!(treeedge, pe.number)
                push!(hybparent, 0)
                push!(leafname, (pn.leaf ? pn.name : ""))
              end
            end
          end
        end
    end
    # for cl in clade @show taxa[cl]; end; @show treenode; @show treeedge; @show hybparent; @show leafname
    # @show hybind; @show hybnode; @show majsisedge; @show minsisedge; @show majsisind; @show minsisind

    # node-associated summaries:
    nclades = length(clade)
    nh = length(hybnode)
    nedges = 2*nh
    BShyb        = zeros(Float64, nclades) # clade = hybrid
    BSmajsis     = zeros(Float64, nclades) # clade = major sister of some hybrid
    BSminsis     = zeros(Float64, nclades) # clade = minor sister of some hybrid
    # edge-associated summaries, i.e. associated to a pair of clades (hyb,sis):
    hybcladei   = repeat(hybind, inner=[2]) # indices in 'clade'
    siscladei   = Array{Int}(undef, nedges) # edge order: (major then minor) for all hybrids
    edgenum     = Array{Int}(undef, nedges)
    for i=1:nh
        siscladei[2*i-1] = majsisind[i]
        siscladei[2*i]   = minsisind[i]
        edgenum[2*i-1] = majsisedge[i]
        edgenum[2*i]   = minsisedge[i]
    end # less desirable order: [majsisind; minsisind], corresponds to outer=[2] for hybcladei
    BShybmajsis = zeros(Float64, nedges) # one clade = hybrid, one clade = major sister
    BShybminsis = zeros(Float64, nedges) # one clade = hybrid, one clade = minor sister
    # network*edge-associated values: only keep those for edges in ref net
    gamma = zeros(Float64,numNets,nedges)
    # node-associate 3-way partition summary:
    BShybsamesis = zeros(Float64, nh) # same hybrid & sister set as in ref net

    nextnum = max(maximum([n.number for n in refnet.edge]),
                  maximum([n.number for n in refnet.node]) ) + 1

    for i = 1:numNets
        net = nets[i]
        length(tiplabels(net))==ntax || error("networks have non-matching taxon sets")
        try directedges!(net) # make sure the root is admissible
        catch err
          if isa(err, RootMismatch)
            err.msg *= "\nPlease change the root in test network (see rootatnode! or rootatedge!)"
          end
          rethrow(err)
        end
        for esth = 1:net.numhybrids   # try to match estimated hybrid edge
            hwcPar = zeros(Bool,ntax)   # minor sister clade
            hwcChi = zeros(Bool,ntax)   #       child  clade
            hwcSib = zeros(Bool,ntax)   # major sister clade
            net1 = deepcopy(net)
            displayednetworkat!(net1, net1.hybrid[esth])
            hn = net1.hybrid[1]
            if !rooted && length(net1.node[net1.rooti].edge)==2 && any(e -> e.hybrid, net1.node[net1.rooti].edge)
                fuseedgesat!(net1.rooti, net1)
            end
            hemaj, hemin, ce = hybridEdges(hn) # assumes no polytomy at hybrid node
            (hemin.hybrid && !hemin.ismajor) || error("edge should be hybrid and minor")
            (hemaj.hybrid &&  hemaj.ismajor) || error("edge should be hybrid and major")
            hardwiredcluster!(hwcChi,hemin,taxa)
            for sis in ["min","maj"]
                he = (sis=="min" ? hemin : hemaj)
                pn = getparent(he) # parent of hybrid edge
                atroot = (!rooted && pn ≡ net1.node[net1.rooti])
                # if at root: exclude the child edge in the same cycle as he.
                # its cluster includes hwcChi. all other child edges do not interest hwcChi.
                # if (atroot) @show i; @warn "$(sis)or edge is at the root!"; end
                for ce in pn.edge
                  if ce ≢ he && pn ≡ getparent(ce)
                    hwc = hardwiredcluster(ce,taxa)
                    if !atroot || sum(hwc .& hwcChi) == 0 # empty intersection
                      if (sis=="maj") hwcSib .|= hwc;
                      else            hwcPar .|= hwc; end
                    elseif (hwc .& hwcChi) != hwcChi
                        @warn "weird clusters at the root. net i=$i in the sample, hybrid $(net.hybrid[esth].name)"
                    end
                  end
                end # will use complement too: test network may be rooted differently
            end
            # @show taxa[hwcChi]; @show taxa[hwcPar]
            if all(hwcPar) || all(hwcSib) || all(.!hwcPar) || all(.!hwcSib)
                @warn "parent or sibling cluster is full or empty. net i=$i in the sample, hybrid $(net.hybrid[esth].name)"
            end

            ihyb = findfirst(isequal(hwcChi), clade)
            newhyb = (ihyb===nothing) # hwcChi not found in clade list
            if newhyb
              push!(clade, hwcChi)
              push!(treenode, nextnum)
              push!(treeedge, nextnum)
              push!(BShyb,    1.0)
              push!(BSmajsis, 0.0)
              push!(BSminsis, 0.0)
              ihyb = length(clade) # was nothing; converted to nex available integer
              nextnum += 1
            else
              BShyb[ihyb] += 1.0
            end

            iSmaj = findfirst(isequal(hwcSib), clade)
            iSmin = findfirst(isequal(hwcPar), clade)
            if (!rooted && iSmaj===nothing) iSmaj = findfirst(isequal(.!hwcSib), clade) end
            if (!rooted && iSmin===nothing) iSmin = findfirst(isequal(.!hwcPar), clade) end
            newmaj = (iSmaj===nothing)
            if newmaj
              push!(clade, hwcSib)
              push!(treenode, nextnum)
              push!(treeedge, nextnum)
              push!(BShyb,    0.0)
              push!(BSmajsis, 1.0)
              push!(BSminsis, 0.0)
              iSmaj = length(clade) # was nothing; now integer
              nextnum += 1
            else
              BSmajsis[iSmaj] += 1.0
            end
            newmin = (iSmin===nothing)
            if newmin
              push!(clade, hwcPar)
              push!(treenode, nextnum)
              push!(treeedge, nextnum)
              push!(BShyb,    0.0)
              push!(BSmajsis, 0.0)
              push!(BSminsis, 1.0)
              iSmin = length(clade) # was nothing; now integer
              nextnum += 1
            else
              BSminsis[iSmin] += 1.0
            end
            samehyb = (newhyb ? Int[] : findall(isequal(ihyb), hybcladei))
            for sis in ["maj","min"]
                newpair = newhyb || (sis=="min" ? newmin : newmaj)
                iSsis = (sis=="min" ? iSmin : iSmaj)
                if !newpair # hyb and sis clades were already found, but not sure if together
                    iish = findfirst(isequal(iSsis), siscladei[samehyb])
                    if iish !== nothing # pair was indeed found
                        ipair = samehyb[iish]
                        if (sis=="min") BShybminsis[ipair] += 1.0
                        else            BShybmajsis[ipair] += 1.0 end
                        if ipair <= nedges # pair = edge in reference network
                            gamma[i,ipair] = (sis=="min" ? hemin.gamma : hemaj.gamma)
                        end
                    else newpair = true; end
                end
                if newpair
                    push!(hybcladei, ihyb)
                    push!(siscladei, iSsis)
                    push!(BShybminsis, (sis=="min" ? 1.0 : 0.0))
                    push!(BShybmajsis, (sis=="min" ? 0.0 : 1.0))
                end
            end
            if ihyb<=nclades && hybparent[ihyb]>0 # hyb clade match in ref net
                th = hybparent[ihyb]
                if ((iSmaj==majsisind[th] && iSmin==minsisind[th]) ||
                    (iSmin==majsisind[th] && iSmaj==minsisind[th]) )
                    BShybsamesis[th] += 1
                end
            end
        end # loop over hybrid nodes
    end # loop over bootstrap nets

    fac = 100/numNets
    BShyb    *= fac # size: length(clade)
    BSmajsis *= fac
    BSminsis *= fac
    BShybmajsis  *= fac # size: length of hybcladei or siscladei
    BShybminsis  *= fac
    BShybsamesis *= fac # size: length(hybnode) = nedges/2

    # combine results into 2 dataframes: node & edge summaries
    # and 1 dataframe to describe clades
    # first: detect nodes with all BS=0
    keepc = ones(Bool,length(clade)) # keep clade h in output tables?
    for h=nclades:-1:1
        if BShyb[h]==0.0 && BSmajsis[h]==0.0 && BSminsis[h]==0.0
            keepc[h] = false
            deleteat!(BShyb,h); deleteat!(BSmajsis,h); deleteat!(BSminsis,h)
        end
    end
    nkeepc = sum(keepc)
    # clade descriptions
    resCluster = DataFrame(taxa=taxa)
    cladestr = Array{String}(undef, length(clade))
    rowh = 1
    for h=1:length(clade)
        nn = treenode[h] # node number in ref net
        cladestr[h] = string("c_", (nn<0 ? "minus" : ""), abs(nn))
        # column name for a clade at node 5: "c_5". At node -5: "c_minus5"
        # because symbol :c_-5 causes an error.
        if h <= nclades &&  leafname[h] != "" # replace "c5" by leaf name, e.g. "taxon1"
            cladestr[h] = leafname[h]
        end
        if h <= nclades &&  hybparent[h]>0 # replace "c5" or leaf name by hybrid name, e.g. "H1"
            na = refnet.hybrid[hybparent[h]].name
            cladestr[h] = (na=="" ? string("H", refnet.hybrid[hybparent[h]].number) : replace(na, r"^#" => ""))
        end
        if keepc[h]
            rowh += 1
            insertcols!(resCluster, rowh, Symbol(cladestr[h]) => clade[h])
        end
    end
    # node summaries
    resNode = DataFrame(clade=cladestr[keepc],
                        node=allowmissing(treenode[keepc]),
                        hybridnode=allowmissing(treenode[keepc]),
                        edge=allowmissing(treeedge[keepc]),
                        BS_hybrid=BShyb, BS_sister = BSmajsis + BSminsis,
                        BS_major_sister=BSmajsis, BS_minor_sister=BSminsis,
                        BS_hybrid_samesisters=Vector{Union{Missing,Float64}}(undef, nkeepc))
    rowh = 1
    for h=1:length(clade)
        if h <= nclades && keepc[h] && hybparent[h]>0
            resNode[rowh,:hybridnode]            =      hybnode[hybparent[h]]
            resNode[rowh,:BS_hybrid_samesisters] = BShybsamesis[hybparent[h]]
        elseif keepc[h]
            resNode[rowh,:BS_hybrid_samesisters] = missing
        end
        if h>nclades # clade *not* in the reference network
            resNode[rowh,:node] = missing
            resNode[rowh,:hybridnode] = missing
            resNode[rowh,:edge] = missing
        end
        if keepc[h]  rowh += 1; end
    end
    insertcols!(resNode, 10, :BS_all => resNode[!,:BS_hybrid]+resNode[!,:BS_sister])
    sort!(resNode, [:BS_all,:BS_hybrid]; rev=true)
    select!(resNode, Not(:BS_all))
    # edge summaries
    resEdge = DataFrame(edge = Vector{Union{Int, Missing}}(undef, length(hybcladei)),
                        hybrid_clade=cladestr[hybcladei],
                        hybrid=Vector{Union{Int, Missing}}(treenode[hybcladei]),
                        sister_clade=cladestr[siscladei],
                        sister=Vector{Union{Int, Missing}}(treenode[siscladei]),
                        BS_hybrid_edge = BShybmajsis+BShybminsis,
                        BS_major=BShybmajsis, BS_minor=BShybminsis)
    for i=1:length(hybcladei)
        h = hybcladei[i]
        if h <= nclades && hybparent[h]>0
            resEdge[i,:hybrid] = hybnode[hybparent[h]]
        end
        if h>nclades            resEdge[i,:hybrid]=missing; end
        if siscladei[i]>nclades resEdge[i,:sister]=missing; end
        if i <= nedges
             resEdge[i,:edge] = edgenum[i]
        else resEdge[i,:edge] = missing
        end
    end
    o = [1:nedges; sortperm(resEdge[nedges+1:length(hybcladei),:BS_hybrid_edge],rev=true) .+ nedges]
    return resNode, resEdge[o,:], resCluster, gamma, edgenum
end
