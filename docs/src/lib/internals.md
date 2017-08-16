```@meta
CurrentModule = PhyloNetworks
DocTestSetup = quote
  using PhyloNetworks
end
```

# Internal Documentation

## Contents

```@contents
Pages = ["internals.md"]
```

## Index

```@index
Pages = ["internals.md"]
```

## Types

```@docs
ANode
MatrixTopologicalOrder
```

## Functions and methods

```@docs
assignhybridnames!
deleteNode!
sampleBootstrapTrees
sampleCFfromCI
setNonIdBL!
recursionPreOrder
recursionPreOrder!
updatePreOrder!
sharedPathMatrix
recursionPostOrder
updatePostOrder!
getindex(::MatrixTopologicalOrder, ::Symbol)
getMajorParent
getMinorParent
getMajorParentEdge
getMinorParentEdge
```

## Main functions in SNaQ

```@docs
optTopRuns!
optTopRun1!
optTopLevel!
proposedTop!
optBL!
afterOptBL!
gammaZero!
moveHybrid!
afterOptBLRepeat!
afterOptBLAll!
```

```@meta
DocTestSetup = nothing
```

