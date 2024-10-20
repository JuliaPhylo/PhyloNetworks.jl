"""
    startree_newick(n, l)

String for the Newick parenthetical description of the star tree
with `n` tips, and all branch lengths equal to `l`.
"""
function startree_newick(n::Integer, ell=1.0)
    return "(" * join(["t$i:$ell" for i=1:n], ",") * ");"
end

"""
    symmetrictree_newick(n::Int, ell::Real, i=1)

String for the Newick parenthetical description of a symmetric tree with
2^n tips, numbered from i to i+2^n-1. All branch lengths are set equal to `ell`.
The tree can be created later by reading the string with [`readTopology`](@ref).
"""
function symmetrictree_newick(n::Int, ell::Real, i::Int=1)
    tree = "A$(i-1+2^n):$(ell)"
    if n==0 return("("*"A$(i-1+2^n):0"*");") end
    for k in 1:(n-1)
        tree = "(" * tree * "," * tree * "):$(ell)"
    end
    tree = "(" * tree * "," * tree * ");"
    # Rename tips
    for k in (2^n-1):-1:1
        tree = replace(tree, "A$(i-1+k+1):$(ell)"=>"A$(i-1+k):$(ell)", count=k)
    end
    return(tree)
end

"""
    symmetricnet_newick(n::Int, i::Int, j::Int, γ::Real, l::Real)

Newick string for a network with a symmetric major tree with 2^n tips,
numbered from 1 to 2^n. All the branch lengths of the major tree are set to `l`.
One hybrid branch, going from level i to level j is added, cutting in half
each initial edge in the tree. The new edge has length `l` and inheritance `γ`.
"""
function symmetricnet_newick(n::Int, i::Int, j::Int, gamma::Real, ell::Real)
    (n < i || i < j || j < 1) && error("must have n >= i >= j > 0")
    # Underlying tree
    tree = symmetrictree_newick(n, ell)
    ## start hyb
    op = "(A1:$(ell)"
    clo = "A$(2^(i-1)):$(ell)):$(ell)"
    for k in 3:i
        op = "("*op
        clo = clo*"):$(ell)"
    end
    clobis = clo[1:(length(clo)-length("$(ell)"))]*"$(0.5*ell)"
    tree = replace(tree,  op => "(#H:$(ell)::$(gamma),"*op)
    tree = replace(tree, clo => clobis*"):$(0.5*ell)")
    ## end hyb
    op = "A$(2^(i-1)+1):$(ell)"
    clo = "A$(2^(i-1) + 2^(j-1)):$(ell)"
    for k in 2:j
        op = "("*op
        clo = clo*"):$(ell)"
    end
    clobis = clo[1:(length(clo)-length("$(ell)"))]*"$(0.5*ell)"
    tree = replace(tree, op  => "("*op)
    tree = replace(tree, clo => clobis*")#H:$(0.5*ell)::$(1-gamma)")
    return(tree)
end


"""
    symmetricnet_newick(n::Int, h::Int, γ::Real)

Create a string for a symmetric network with 2^n tips, numbered from 1 to 2^n,
with a symmetric major tree, whose branch lengths are all equal.
2^(n-h) hybrids are added from level h to h-1 "symmetrically".
The network is time-consistent and ultrametric, with a total height of 1.
"""
function symmetricnet_newick(n::Int, h::Int, gamma::Real, i::Int=1)
    (n < h || h < 2) && error("must have n >= h > 1.")
    # length of branch
    ell = 1.0/n
    # Element net
    net = symmetricnet_newick(h, h, h-1, gamma, ell)
    # Iterate
    if n==h return(net) end
    net = chop(net)
    for k in 1:(n-h)
        net = "(" * net * ":$(ell)," * net * ":$(ell))"
    end
    net = net * ";"
    # Rename hybrids
    net = replace(net, "H" => "H$(2^(n-h))")
    for k in (2^(n-h)-1):-1:1
        net = replace(net, "H$(k+1)" => "H$(k)", count=2*k)
    end
    # Rename tips
    for k in (2^h):-1:1
        net = replace(net, "A$(k):" => "A$(i-1+2^n):")
    end
    for k in (2^n-1):-1:1
        net = replace(net, "A$(i-1+k+1):" => "A$(i-1+k):", count=k)
    end
    return(net)
end
