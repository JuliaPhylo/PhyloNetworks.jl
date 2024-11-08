__precompile__()

module PhyloNetworks

    # stdlib (standard libraries)
    using Dates
    using Distributed # used for parsimony search, currently broken
    using LinearAlgebra: diag, I, logdet, norm, LowerTriangular, mul!, lmul!, rmul!,
            Diagonal, cholesky, qr, BLAS
    # alternative: drop support for julia v1.4, because LinearAlgebra.rotate! requires julia v1.5
    # using LinearAlgebra # bring all of LinearAlgebra into scope
    # import LinearAlgebra.rotate! # allow re-definition of rotate!
    using Printf: @printf, @sprintf
    using Random
    using Statistics: mean, quantile, median

    # other libraries, indicate compatible version in Project.toml
    using BioSequences
    using BioSymbols
    using Combinatorics: combinations
    using CSV
    using DataFrames # innerjoin new in v0.21
    using DataStructures # for updateInCycle with priority queue
    using Distributions #for RateVariationAcrossSites
    using FASTX
    using Functors: fmap
    using NLopt # for branch lengths optimization
    using StatsBase # cov2cor! for vcv
    #using StatsFuns # logsumexp, logaddexp, log2π, various cdf

    import Base: show

    const DEBUGC = false # even more debug messages
    global CHECKNET = false # for debugging only

    export
        ## Network Definition
        HybridNetwork,
        readnewick,
        writenewick,
        writeSubTree!,
        hybridlambdaformat,
        deleteleaf!,
        deleteaboveLSA!,
        removedegree2nodes!,
        shrink2cycles!,
        shrink3cycles!,
        printEdges,
        printNodes,
        readInputTrees,
        readnexus_treeblock,

        ## getters
        # fixit: add ancestors? getsibling? getdescendants (currently descendants)?
        getroot,
        istimeconsistent,
        getnodeheights,
        getnodeheights!,
        getnodeheights_average,
        getnodeheights_majortree,
        isrootof,
        isleaf,
        isexternal,
        isparentof,
        ischildof,
        hassinglechild,
        getchild,
        getchildren,
        getchildedge,
        getparent,
        getparents,
        getparentminor,
        getparentedge,
        getparentedgeminor,
        getpartneredge,
        tipLabels,
        ## Network Manipulation
        rootatnode!,
        rootonedge!,
        directEdges!,
        preorder!,
        cladewiseorder!,
        rotate!,
        setlength!,
        setlengths!,
        setGamma!,
        deleteHybridThreshold!,
        displayedTrees,
        majorTree,
        minorTreeAt,
        displayedNetworkAt!,
        hardwiredClusters,
        hardwiredCluster,
        hardwiredCluster!,
        hardwiredClusterDistance,
        biconnectedComponents,
        blobDecomposition!,
        blobDecomposition,
        nni!,
        checkroot!,
        treeedgecomponents,
        ## Network support
        treeedges_support,
        hybridclades_support,
        readmultinewick_files,
        writemultinewick,
        readmultinewick,
        ## Network Calibration
        getNodeAges,
        pairwiseTaxonDistanceMatrix,
        calibrateFromPairwiseDistances!,
        # recursion
        sharedPathMatrix,
        descendenceMatrix,
        vcv,
        ## parsimony
        parsimonySoftwired,
        parsimonyGF,
        #maxParsimonyNet # broken after v0.17 refactoring: fix network search
        readfastatodna,
        # neighbor joining
        nj

    ##Constants
    const fAbsBL = 1e-10
    const fRelBL = 1e-12
    const xAbsBL = 1e-10
    const xRelBL = 1e-10

    include("types.jl")
    include("auxiliary.jl")
    include("generate_topology.jl")
    include("moves_semidirected.jl")
    include("addHybrid.jl")
    include("readwrite.jl")
    include("descriptive.jl")
    include("manipulateNet.jl")
    include("bootstrap.jl")
    include("compareNetworks.jl")
    include("recursion_routines.jl")
    include("recursion_matrices.jl")
    include("parsimony.jl")
    include("pairwiseDistanceLS.jl")
    include("interop.jl")
    include("graph_components.jl")
    include("deprecated.jl")
    include("nj.jl")

end #module
