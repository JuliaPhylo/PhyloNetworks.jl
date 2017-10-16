__precompile__()

module PhyloNetworks

    using DataStructures # for updateInCycle with priority queue
    using DataFrames # for functions to read/write tables, and names()
    using GLM # for the lm function
    using NLopt # for branch lengths optimization
    using Gadfly # for plots
    using ColorTypes # used by Gadfly already. To resolve data type names (Colorant)
    using StatsBase: sample
    using Combinatorics.combinations
    using RCall
    using RCall: protect, unprotect, RClass
    using NullableArrays

    import Base.show
    import Gadfly.plot
    import GLM.ftest
    import RCall.sexp
    export plot, ftest, sexp

    global DEBUG = false #for debugging only
    const DEBUGC = false #more detailed prints
    global CHECKNET = false #for debugging only
    global REDIRECT = false # changed for debugging to a file

    export
        ## Network Definition
        HybridNetwork,
        DataCF,
        Quartet,
        readTopology,
        readTopologyLevel1,
        tipLabels,
        writeTopology,
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
        parsimonyDiscrete,
        EqualRatesSubstitutionModel,
        BinaryTraitSubstitutionModel,
        randomTrait,
        randomTrait!,
        #
        apeRExport

    # export part

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
    include("plotsGadfly.jl")
    include("plotsRCall.jl")
    include("compareNetworks.jl")
    include("traits.jl")
    include("parsimony.jl")
    include("pairwiseDistanceLS.jl")
    include("apeRExport.jl")

end #module
