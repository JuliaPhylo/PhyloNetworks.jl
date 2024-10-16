__precompile__()

module PhyloNetworks

    # stdlib (standard libraries)
    using Dates
    using Distributed
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
    using GLM # for the lm function
    using NLopt # for branch lengths optimization
    using StaticArrays
    using StatsBase # sample, coef etc.
    using StatsFuns # logsumexp, logaddexp, log2π, various cdf
    using StatsModels # re-exported by GLM. for ModelFrame ModelMatrix Formula etc

    import Base: show
    import GLM: ftest
    import StatsModels: coefnames

    const DEBUGC = false # even more debug messages
    global CHECKNET = false # for debugging only

    export ftest
    export
        ## Network Definition
        HybridNetwork,
        readTopology,
        writeTopology,
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
        getGammas,
        getHeights,
        getHeights!,
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
        setGammas!,
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
        ## Network Bootstrap
        treeEdgesBootstrap,
        hybridDetection,
        summarizeHFdf,
        hybridBootstrapSupport,
        readBootstrapTrees,
        writeMultiTopology,
        readMultiTopology,
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
        maxParsimonyNet,
        readfastatodna,
        # neighbor joining
        nj

    include("types.jl")
    include("auxiliary.jl")
    include("generate_topology.jl")
    include("addHybrid.jl")
    include("moves_semidirected.jl")
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
