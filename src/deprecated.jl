@deprecate phyloNetworklm phylolm
@deprecate sigma2_estim sigma2_phylo
@deprecate mu_estim mu_phylo
@deprecate getChild getchild false
@deprecate getChildren getchildren false
@deprecate getChildEdge getchildedge false
@deprecate getParents     getparents false
@deprecate getParent      getparent  false
@deprecate getMajorParent getparent  false
@deprecate getMinorParent getparentminor false
@deprecate getMajorParentEdge getparentedge false
@deprecate getMinorParentEdge getparentedgeminor false
@deprecate getPartner getpartneredge false

function Base.getproperty(mm::PhyloNetworkLinearModel, f::Symbol)
    if f === :model
        Base.depwarn("""accessing the 'model' field of PhyloNetworkLinearModel's is
            deprecated, as they are no longer wrapped in a TableRegressionModel.
            They can be used directly now.""",
            :getproperty) # force=true to show warning. currently silent by default
        return mm
    elseif f === :mf
        Base.depwarn("""accessing the 'mf' field of PhyloNetworkLinearModel's is
            deprecated, as they are no longer wrapped in a TableRegressionModel.
            Use `formula(m)` to access the model formula.""", :getproperty)
        form = formula(mm)
        return ModelFrame{Nothing, typeof(mm)}(form, nothing, nothing, typeof(mm))
    elseif f === :mm
        Base.depwarn("""accessing the 'mm' field of PhyloNetworkLinearModel's is
            deprecated, as they are no longer wrapped in a TableRegressionModel.
            Use `modelmatrix(m)` to access the model matrix.""", :getproperty)
        modmatr = modelmatrix(mm)
        return ModelMatrix{typeof(modmatr)}(modmatr, Int[])
    else
        return getfield(mm, f)
    end
end