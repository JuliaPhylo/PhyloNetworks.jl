__precompile__()

module PhyloNetworks

    using DataStructures # for updateInCycle with priority queue
    using DataFrames
    using CSV
    using StatsBase # sample, coef etc.
    using StatsModels
    using GLM # for the lm function
    using NLopt # for branch lengths optimization
    using Combinatorics.combinations
    using StaticArrays
    using BioSequences
    using BioSymbols
    using StaticArrays
    using StatsFuns

    import Base.show
    import GLM.ftest

    global DEBUG = false #for debugging only
    const DEBUGC = false #more detailed prints
    global CHECKNET = false #for debugging only
    global REDIRECT = false # changed for debugging to a file

    export ftest
    export
        ## Network Definition
        HybridNetwork,
        DataCF,
        Quartet,
        readTopology,
        readTopologyLevel1,
        tipLabels,
        writeTopology,
        writeSubTree!,
        deleteleaf!,
        printEdges,
        printNodes,
        sorttaxa!,
        ## SNAQ
        readTrees2CF,
        readTableCF,
        readTableCF!,
        readInputTrees,
        summarizeDataCF,
        snaq!,
        readSnaqNetwork,
        snaqDebug,
        topologyMaxQPseudolik!,
        topologyQPseudolik!,
        ## Network Manipulation
        # getParent, getParents, getMajorParentEdge, getMinorParentEdge, getChildren,
        # functions above: first rename them throughout to be consistent with other packages, like:
        # parent child parents children parentmajor parentminor ancestor sibling offspring
        rootatnode!,
        rootonedge!,
        directEdges!,
        preorder!,
        cladewiseorder!,
        fittedQuartetCF,
        rotate!,
        setLength!,
        setGamma!,
        mapAllelesCFtable,
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
        ## Network Bootstrap
        treeEdgesBootstrap,
        hybridDetection,
        summarizeHFdf,
        hybridBootstrapSupport,
        bootsnaq,
        readBootstrapTrees,
        writeMultiTopology,
        readMultiTopologyLevel1,
        readMultiTopology,
        hybridatnode!,
        undirectedOtherNetworks,
        ## Network Calibration
        getNodeAges,
        pairwiseTaxonDistanceMatrix,
        calibrateFromPairwiseDistances!,
        ## Network PCM
        phyloNetworklm,
        PhyloNetworkLinearModel,
        simulate,
        TraitSimulation,
        ParamsBM,
        ShiftNet,
        shiftHybrid,
        getShiftEdgeNumber,
        getShiftValue,
        sharedPathMatrix,
        descendenceMatrix,
        regressorShift,
        regressorHybrid,
        ancestralStateReconstruction,
        ReconstructedStates,
        sigma2_estim,
        mu_estim,
        lambda_estim,
        expectations,
        expectationsPlot,
        predint,
        predintPlot,
        ## Discrete Trait PCM
        parsimonySoftwired,
        parsimonyGF,
        TraitSubstitutionModel,
        EqualRatesSubstitutionModel,
        BinaryTraitSubstitutionModel,
        TwoBinaryTraitSubstitutionModel,
        nStates,
        Q, P,
        randomTrait,
        randomTrait!

    include("types.jl")
    include("auxiliary.jl")
    include("update.jl")
    include("undo.jl")
    include("addHybrid.jl")
    include("deleteHybrid.jl")
    include("moves.jl")
    include("readwrite.jl")
    include("readData.jl")
    include("optimization.jl")
    include("pseudolik.jl")
    include("descriptive.jl")
    include("manipulateNet.jl")
    include("bootstrap.jl")
    include("multipleAlleles.jl")
    include("compareNetworks.jl")
    include("traits.jl")
    include("parsimony.jl")
    include("pairwiseDistanceLS.jl")
    include("interop.jl")
    include("substitutionModels.jl")
    include("biconnectedComponents.jl")
    include("traitsLikDiscrete.jl")

end #module
