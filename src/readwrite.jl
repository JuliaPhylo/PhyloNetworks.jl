# functions to read/write networks topologies

#= peek the next non-white-space char
   removing white spaces from the IOStream/IOBuffer.
code similar to [skipchars](https://docs.julialang.org/en/v1/base/io-network/#Base.skipchars)
in julia/base/io.jl, except that
- the peeked character is returned (instead of io)
- 'missing' is returned if io is empty (i.e. eof(io))
=#
function peekskip(io::IO, linecomment=nothing)
    for c in readeach(io, Char) # empty iterable if eof(io)
        if c === linecomment
            readline(io)
        elseif !isspace(c)
            skip(io, -ncodeunits(c))
            return c # break
        end
    end
    return missing
end

# aux function to read the next non-white-symbol char in s, advances s
# input: s IOStream/IOBuffer
function readskip!(io::IO)
    c = missing
    while !eof(io)
        c = read(io, Char)
        if !isspace(c)
            break
        end
    end
    return c
end


# advance stream in readSubtree
function advance!(s::IO, numLeft::Array{Int,1})
    c = readskip!(s)
    if ismissing(c)
        error("Tree ends prematurely while reading subtree after left parenthesis $(numLeft[1]-1).")
    end
    return c
end

"""
    readnexus_comment(s::IO, c::Char)

Read (and do nothing with) nexus-style comments: `[& ... ]`  
Assumption: 'c' is the next character to be read from s.  
Output: nothing.

Comments can appear after (or instead of) a node or leaf name,
before or after an edge length, and after another comment.
"""
function readnexus_comment(s::IO, c::Char)
    while c == '[' # a comment could be followed by another
        # [ should be followed by &, otherwise bad newick string
        # read [ and next character: don't skip spaces
        read(s, Char) == '[' || error("I was supposed to read '['")
        read(s, Char) == '&' || error("read '[' but not followed by &")
        skipchars(!isequal(']'), s)
        eof(s) && error("comment without ] to end it")
        read(s, Char) # to read ] and advance s
        c = peekskip(s)
    end
    return
end

"""
    readnewick_nodename(s::IO, c::Char, net, numLeft)

Auxiliary function to read a taxon name during newick parsing.
output: tuple (number, name, pound_boolean)

Names may have numbers: numbers are treated as strings.
Accepts `#` as part of the name (but excludes it from the name), in which
case `pound_boolean` is true. `#` is used in extended newick to flag
hybrid nodes.

Nexus-style comments following the node name, if any, are read and ignored.
"""
function readnewick_nodename(s::IO, c::Char, net::HybridNetwork, numLeft::Array{Int,1})
    pound = 0
    name = ""
    while isvalidnewick_namesymbol(c)
        readskip!(s) # advance s past c (discard c, already read)
        if c == '#'
            pound += 1
            c = readskip!(s) # read the character after the #
            if !isletter(c)
                a = read(s, String);
                error("Expected name after # but received $(c) in left parenthesis $(numLeft[1]-1). remaining: $(a).")
            end
            if c != 'L' && c != 'H' && c != 'R'
                @warn "Expected H, R or LGT after # but received $(c) in left parenthesis $(numLeft[1]-1)."
            end
            name = c * name # put the H, R or LGT first
        else
            name *= c # original c: put it last
        end
        c = peekskip(s);
    end
    if pound >1
        a = read(s, String);
        error("strange node name with $(pound) # signs: $name. remaining: $(a).")
    end
    readnexus_comment(s,c)
    return size(net.names,1)+1, name, pound == 1
end


# aux function to read floats like length or gamma values, to be read after a colon
function readnewick_float(s::IO, c::Char)
    if !(isdigit(c) || c in ['.','e','-','E', '+'])
        a = read(s, String);
        error("Expected float digit after ':' but found $(c). remaining is $(a).");
    end
    num = ""
    while isdigit(c) || c in ['.','e','-', 'E', '+']
        d = read(s, Char) # reads c and advances IO
        num = string(num,d);
        c = peekskip(s);
    end
    f = try
        parse(Float64, num)
    catch
        error("problem with number read $(num), not a float number")
    end
    return f
end

# aux function to identify if a symbol in taxon name is valid
# allowed: letters, numbers, underscores, # (for hybrid name)
# NOT allowed: white space () [] : ; ' ,
# according to richnewick.pdf
# actually: ' is allowed. coded weirdly imo!
function isvalidnewick_namesymbol(c::Char)
    return isletter(c) || isnumeric(c) || c=='_' || c=='#' ||
      (!isspace(c) && c != '(' && c != ')' && c != '[' && c != ']' && c != ':' && c != ';' && c != ',')
end

"""
    parsenewick_remainingsubtree!(s::IO, numLeft, net, hybrids)

Create internal node. Helper for [`readnewick_subtree!`](@ref),
which creates the parent edge of the node created by `parsenewick_remainingsubtree!`:
`readnewick_subtree!` calls `parsenewick_remainingsubtree!`, and vice versa.
Called once a `(` has been read in a tree topology and reads until the corresponding `)` has been found.
This function performs the recursive step for `readnewick_subtree!`.
Advances `s` past the subtree, adds discovered nodes and edges to `net`, and `hybrids`.

Does *not* read the node name and the edge information of the subtree root:
this is done by [`readnewick_subtree!`](@ref)
"""
@inline function parsenewick_remainingsubtree!(
    s::IO,
    numLeft::Array{Int,1},
    net::HybridNetwork,
    hybrids::Vector{String}
)
    numLeft[1] += 1
    n = Node(-1*numLeft[1],false);
    # @debug "creating node $(n.number)"
    keepon = true;
    c = readskip!(s)
    while (keepon)
        bl = readnewick_subtree!(s,n,numLeft,net,hybrids)
        c = advance!(s,numLeft)
        if c == ')'
            keepon = false
        elseif c != ','
            a = read(s, String);
            error("Expected right parenthesis after left parenthesis $(numLeft[1]-1) but read $(c). The remainder of line is $(a).")
        end
    end
    return n
end

"""
    parsenewick_hybridnode!(node, parentNode, hybridName, net, hybrids)

Helper for [`readnewick_subtree!`](@ref). Create the parent edge for `node`.
Return this edge, and the hybrid node retained (`node` or its clone in the newick string).
Insert new edge and appropriate node into `net` and `hybrids` accordingly.
Handles any type of given hybrid node.
Called after a `#` has been found in a tree topology.
"""
@inline function parsenewick_hybridnode!(
    n::Node,
    parent::Node,
    name::String,
    net::HybridNetwork,
    hybrids::Vector{String}
)
    n.hybrid = true;
    ind = findfirst(isequal(name), hybrids) # index of 'name' in the list 'hybrid'. nothing if not found
    e = Edge(net.numedges+1) # ismajor = true by default
    if n.leaf e.ismajor = false; end
    e.hybrid = true
    e.gamma = -1.0
    if ind !== nothing # the hybrid name was seen before
        ni = findfirst(isequal(name), [no.name for no in net.node])
        ni !== nothing || error("hybrid name $name was supposed to be in the network, but not found")
        other = net.node[ni]
        if !n.leaf && !other.leaf
            error("both hybrid nodes are internal nodes: successors of the hybrid node must only be included in the node list of a single occurrence of the hybrid node.")
        elseif n.leaf
            pushEdge!(net,e);
            setNode!(e,[other,parent]); # ischild1 = true by default constructor
            setEdge!(other,e);
            setEdge!(parent,e);
            n = other # original 'n' dropped, 'other' retained: 'n' output to modify 'n' outside
        else # !n.leaf : delete 'other' from the network
            # @debug "n is not leaf, other is leaf"
            size(other.edge,1) == 1 || # other should be a leaf
               error("strange: node $(other.number) is a leaf hybrid node. should have only 1 edge but has $(size(other.edge,1))")
            # @debug "other is $(other.number), n is $(n.number), edge of other is $(other.edge[1].number)"
            otheredge = other.edge[1];
            removeNode!(other,otheredge);
            deleteNode!(net,other);
            setNode!(otheredge,n);
            setEdge!(n,otheredge);
            # @debug "creating hybrid edge $(e.number) between n $(n.number) and parent $(parent.number)"
            setNode!(e,[n,parent]);
            setEdge!(n,e);
            setEdge!(parent,e);
            pushNode!(net,n);
            pushEdge!(net,e);
            n.number = other.number; # modifies original negative node number, to positive node #
            n.name = other.name;
            # @debug "edge $(e.number) boole1? $(e.boole1)"
            # @debug "otheredge $(otheredge.number) boole1? $(otheredge.boole1)"
        end
    else # ind==nothing: hybrid name not seen before
        # @debug "$(name) not found in hybrids list"
        # @debug "$(name) is leaf? $(n.leaf)"
        n.hybrid = true;
        nam = string(name)
        push!(net.names, nam);
        n.name = nam;
        push!(hybrids, nam);
        pushNode!(net,n);
        # @debug "creating hybrid edge $(e.number)"
        pushEdge!(net,e);
        setNode!(e,[n,parent]);
        setEdge!(n,e);
        setEdge!(parent,e);
        # @debug "edge $(e.number) boole1? $(e.boole1)"
    end
    e.containroot = !e.hybrid # not good: but necessay for SNaQ functions
    return (e,n)
end

"""
    parsenewick_treenode!(node, parentNode, net)

Helper for [`readnewick_subtree!`](@ref).
Insert the input tree node and associated edge (created here) into `net`.
"""
@inline function parsenewick_treenode!(n::Node, parent::Node, net::HybridNetwork)
    pushNode!(net,n);
    e = Edge(net.numedges+1);
    pushEdge!(net,e);
    setNode!(e,[n,parent]);
    setEdge!(n,e);
    setEdge!(parent,e);
    return e
end

"""
    parsenewick_getfloat!(s::IO, int, numLeft::Array{Int,1})

Helper for [`parsenewick_edgedata!`](@ref).
Read a single floating point edge data value in a tree topology.
Ignore (and skip) nexus-style comments before & after the value
(see [`readnexus_comment`](@ref)) then
advance `s` past the next colon character.

output:
- -1.0 if no value exists before the next colon, closing parenthesis or comma
- the numerical value as a float otherwise, read by `readnewick_float`, which
  may error if what comes next in `s` cannot be parsed as a numerical value.
"""
@inline function parsenewick_getfloat!(s::IO, call::Int, numLeft::Array{Int,1})
    errors = ["first colon : without double after in left parenthesis $(numLeft[1]-1), ignored.",
              "second colon : without double after in left parenthesis $(numLeft[1]-1), ignored.",
              "third colon : without gamma value after in $(numLeft[1]-1) left parenthesis, ignored"]
    c = peekskip(s)
    if c == '[' # e.g. comments only, no value, but : after
        readnexus_comment(s,c)
        c = peekskip(s)
    end
    if c == ':' || c == ')' || c == ',' # no value after the previous colon after all
        return -1.0
    else # value is present: read it, and any following comment(s)
    # elseif isdigit(c) || c == '.' || c == '-'
        val = try
            readnewick_float(s, c)
        catch e
            @warn errors[call]
            rethrow(e)
        end
        if val < 0.0
            @error "expecting non-negative value but read '-', left parenthesis $(numLeft[1]-1). will set to 0."
            val = 0.0
        end
        readnexus_comment(s,peekskip(s))
        return val
    end
end

"""
    parsenewick_edgedata!(s::IO, edge, numberOfLeftParentheses::Array{Int,1})

Helper for [`readnewick_subtree!`](@ref). Modify `edge` according to the
specified length, support and gamma values in the newick string, then
advance the stream `s` past any existing edge data.
Edges in a topology may optionally be followed by ":length:support:gamma"
where length, support (such as bootstrap support), and gamma are numerical values.
Nexus-style comments `[&...]`, if any, are ignored.
"""
@inline function parsenewick_edgedata!(s::IO, e::Edge, numLeft::Array{Int,1})
    read(s, Char); # to read the first ":"
    e.length = parsenewick_getfloat!(s, 1, numLeft)
    if peekskip(s) == ':'
        readskip!(s)
        support = parsenewick_getfloat!(s, 2, numLeft)
        e.y = support # e.y originally -1 by default
    end
    # e.gamma = -1.0 by default when e is created by parsenewick_hybridnode!
    if peekskip(s) == ':'
        readskip!(s)
        e.gamma = parsenewick_getfloat!(s, 3, numLeft)
    end
    if e.gamma != -1.0 && !e.hybrid && e.gamma != 1.0
        @warn "γ read for edge $(e.number) but it is not hybrid, so γ=$(e.gamma) ignored"
        e.gamma = 1.0
    end
end

"""
    synchronizepartnersdata!(e::Edge, n::Node)

Synchronize γ and ismajor for edges `e` and its partner,
both hybrid edges with the same child `n`:

- if one γ is missing and the other is not: set the missing γ to 1 - the other
- γ's should sum up to 1.0
- update `ismajor` to match the γ information: the major edge is the one with γ > 0.5.

**Warnings**: does not check that `e` is a hybrid edge,
nor that `n` is the child of `e`.
"""
@inline function synchronizepartnersdata!(e::Edge, n::Node)
    partners = Edge[] # The edges having n as a child, other than e
    for e2 in n.edge
        if e2.hybrid && e2!=e && n==getchild(e2)
            push!(partners, e2)
        end
    end
    numPartners = length(partners)
    if numPartners == 0
        return nothing
    end
    if numPartners > 1 # 3 or more hybrid parents
        error("γ value parsed but hybrid edge has $numPartners partners (should be 0 or 1)")
    end
    # numPartners == 1 then. partner edge has been read, may have gamma set
    partner = partners[1]
    if e.gamma == -1. && partner.gamma == -1. # exact -1.0 means missing
        return nothing
    end
    # γ non-missing for e and/or partner, then
    # update γ and ismajor of both edges, to be consistent with each other
    if e.gamma == -1.
        if partner.gamma < 0.0
            @warn "partners: $(e.number) with no γ, $(partner.number) with γ<0. will turn to 1 and 0"
            partner.gamma = 0.0
            e.gamma = 1.0
        elseif partner.gamma > 1.0
            @warn "partners: $(e.number) with no γ, $(partner.number) with γ>1. will turn to 0 and 1"
            partner.gamma = 1.0
            e.gamma = 0.0
        else # partner γ in [0,1]
            e.gamma = 1. - partner.gamma
        end
    elseif partner.gamma == -1.
        if e.gamma < 0.0
            @warn "partners: $(partner.number) with no γ, $(e.number) with γ<0. will turn to 1 and 0"
            e.gamma = 0.0
            partner.gamma = 1.0
        elseif e.gamma > 1.0
            @warn "partners: $(partner.number) with no γ, $(e.number) with γ>1. will turn to 0 and 1"
            e.gamma = 1.0
            partner.gamma = 0.0
        else # γ in [0,1]
            partner.gamma = 1. - e.gamma
        end
    else # both γ's are non-missing. won't check if negative
        gammasum = e.gamma + partner.gamma
        if !isapprox(gammasum, 1.0)
            e.gamma = e.gamma/gammasum
            partner.gamma = 1. - e.gamma
        end
    end
    # next: update ismajor, originally based on which node was a leaf and which was not
    # if both γ are 0.5: keep ismajor as is. Otherwise: γ's take precedence.
    if !isapprox(e.gamma, 0.5)
        emajor = e.gamma > 0.5
        if e.ismajor != emajor # leaf status was inconsistent with γ info
            e.ismajor = emajor
            partner.ismajor = !emajor
        end
    end
end

"""
    readnewick_subtree!(s::IO, parentNode, numLeft, net, hybrids)

Recursive helper method for `readnewick`:
read a subtree from an extended Newick topology.
input `s`: IOStream/IOBuffer.

Reads additional info formatted as: `:length:bootstrap:gamma`.
Allows for name of internal nodes without # after closing parenthesis: (1,2)A.
Warning if hybrid edge without γ, or if γ (ignored) without hybrid edge
"""
function readnewick_subtree!(s::IO, parent::Node, numLeft::Array{Int,1}, net::HybridNetwork, hybrids::Vector{String})
    c = peekskip(s)
    e = nothing;
    hasname = false; # to know if the current node has name
    pound = false;
    if c == '('
        # read the rest of the subtree (perform the recursive step!)
        n = parsenewick_remainingsubtree!(s, numLeft, net, hybrids)
        c = peekskip(s);
        # read potential internal node name (and skip comments)
        num, name, pound = readnewick_nodename(s, c, net, numLeft);
        if name != ""
            hasname = true;
            n.number = num; # n was given <0 number by parsenewick_remainingsubtree!, now >0
        end
    else # leaf, it should have a name
        hasname = true;
        num, name, pound = readnewick_nodename(s, c, net, numLeft)
        if name == ""
            a = read(s, String);
            error("Expected digit, alphanum or # at the start of taxon name, but received $(c). remaining: $(a).");
        end
        n = Node(num, true); # positive node number to leaves in the newick-tree description
        # @debug "creating node $(n.number)"
    end
    if pound # found pound sign in name
        # merge the 2 nodes corresponding the hybrid: make n the node that is retained
        e,n = parsenewick_hybridnode!(n, parent, name, net, hybrids) # makes hybrid node number >0
    else
        if hasname
            push!(net.names, name);
            n.name = name
        end
        e = parsenewick_treenode!(n, parent, net)
    end
    c = peekskip(s);
    e.length = -1.0
    if c == ':'
        parsenewick_edgedata!(s, e, numLeft)
    end
    if e.hybrid
        # if hybrid edge: 'e' might have no info, but its partner may have had info
        synchronizepartnersdata!(e, n) # update γ and ismajor of e and/or its partner
    end
    return true
end


# function to read topology from parenthetical format
# input: file name or tree in parenthetical format
# calls readnewick(s::IO)
# warning: crashes if file name starts with (
function readnewick(input::AbstractString)
    input = lstrip(input) # strips leading whitespaces
    if input[1] == '(' # input = parenthetical description
       s = IOBuffer(input)
    else # input = file name
        s = try
            open(input)
        catch
            error("Could not find or open $(input) file");
        end
    end
    net = readnewick(s)
    return net
end

"""
    readnewick(file name)
    readnewick(parenthetical description)
    readnewick(IO)

Read a tree or network topology from parenthetical format (extended Newick).

Input: text file or parenthetical format directly. Leading white spaces are ignored.
The input string may start with white spaces followed by the file name.
The file name may *not* start with spaces. Also, it may *not* start with a left parenthesis,
otherwise the file name itself would be interpreted as the parenthetical description.
Nexus-style comments (`[&...]`) are ignored, and may be placed
after (or instead) of a node name, and before/after an edge length.

A root edge, not enclosed within a pair a parentheses, is ignored.
If the root node ρ0 has a single child edge ρ0 → ρ1,
then this initial root ρ0 and its child edge are deleted from the topology.
Its initial child ρ1 becomes the new root.

See also: [`readnexus_treeblock`](@ref)
"""
function readnewick(s::IO)
    net = HybridNetwork()
    line = readuntil(s,";", keep=true);
    if(line[end] != ';')
        error("file does not end in ;")
    end
    seekstart(s)
    c = peekskip(s)
    numLeft = [1]; # made Array to make it mutable; start at 1 to avoid node -1 which breaks undirectedOtherNetworks
    hybrids = String[];
    if(c == '(')
        numLeft[1] += 1;
        #println(numLeft)
        n = Node(-1*numLeft[1],false);
        c = readskip!(s)
        b = false;
        while(c != ';')
            b |= readnewick_subtree!(s,n,numLeft,net,hybrids)
            c = readskip!(s);
            if eof(s)
                error("Tree ended while reading in subtree beginning with left parenthesis number $(numLeft[1]-1).")
            elseif c == ','
                continue;
            elseif c == ')'
                # read potential root name (or comments)
                c = peekskip(s);
                num, name, pound = readnewick_nodename(s, c, net, numLeft)
                if name != ""
                    n.name = name
                    # log warning or error if pound > 0?
                end
                c = peekskip(s)
                if(c == ':') # skip information on the root edge, if it exists
                    # @warn "root edge ignored"
                    while c != ';'
                        c = readskip!(s)
                    end
                end
            end
        end
        # @debug "after readsubtree:"
        # @debug begin printedges(net); "printed edges" end
        # delete the root edge, if present
        if size(n.edge,1) == 1 # root node n has only one edge
            edge = n.edge[1]
            child = getOtherNode(edge,n);
            removeEdge!(child,edge);
            net.rooti = getIndex(child,net);
            deleteEdge!(net,edge);
        else
            pushNode!(net,n);
            net.rooti = getIndex(n,net);
        end
    else
        a = read(s, String)
        error("Expected beginning of tree with ( but received $(c) instead, rest is $(a)")
    end
    storeHybrids!(net)
    checkNumHybEdges!(net)
    directedges!(net; checkMajor=true) # to update edges containroot: true until hybrid, false below hybrid
    net.isrooted = true
    return net
end

"""
    checkNumHybEdges!(net)

Check for consistency between hybrid-related attributes in the network:
- for each hybrid node: 2 or more hybrid edges
- exception: allows for a leaf to be attached to a single hybrid edge
- exactly 2 incoming parent hybrid edges
Run after `storeHybrids!`.
"""
function checkNumHybEdges!(net::HybridNetwork)
    if isTree(net) return nothing; end
    !isempty(net.hybrid) || error("net.hybrid should not be empty for this network")
    for n in net.hybrid
        hyb = sum([e.hybrid for e in n.edge]); # number of hybrid edges attached to node
        if hyb == 1
            if net.numhybrids == 1
                error("only one hybrid node $(n.number) named $(n.name) found with one hybrid edge attached")
            else
                error("hybrid node $(n.number) named $(n.name) has only one hybrid edge attached. there are $(net.numhybrids-1) other hybrids out there but this one remained unmatched")
            end
        elseif hyb == 0
            if length(n.edge) == 0
                error("strange hybrid node $(n.number) attached to 0 edges")
            elseif length(n.edge) == 1
                n.leaf || error("hybrid node $(n.number) has only one tree edge attached and it is not a leaf")
            elseif length(n.edge) >= 2
                @warn "hybrid node $(n.number) not connected to any hybrid edges. Now transformed to tree node"
                n.hybrid = false;
            end
        elseif hyb >=2 # check: exactly 2 incoming, no more.
            nhybparents = 0
            for e in n.edge
                if n == getchild(e)
                    if e.hybrid
                        nhybparents += 1
                    else @error "node $(n.number) has parent tree edge $(e.number): wrong ischild1 for this edge?"
                    end
                end
            end
            if nhybparents < 2
                error("hybrid node $(n.number) with $nhybparents = fewer than 2 hybrid parents")
            elseif nhybparents >2
                error("hybrid node $(n.number) with $nhybparents: we don't allow such polytomies")
            end
        end
    end
    return nothing
end

"""
    resolvetreepolytomy!(net::HybridNetwork, n::Node, e1::Edge, e2::Edge)

Create new node `v` and new edge `n --> v`, detach edges `e1` and `e2` from `n`
to connect them to `v` instead. Return the new edge, whose child is `v`.

Warning: assumes that `n` is incident to both `e1` and `e2`. It is best if
these edges are originally children of `n` because the new edge is directed
from `n` to the newly created node.
"""
function resolvetreepolytomy!(
    net::HybridNetwork,
    n::Node,
    ce1::Edge,
    ce2::Edge,
    b1::Bool=getparent(ce1)===n,
    b2::Bool=getparent(ce2)===n
)
    removeEdge!(n,ce1)
    removeEdge!(n,ce2)
    removeNode!(n,ce1)
    removeNode!(n,ce2)
    unused_enum = setdiff(1:net.numedges, e.number for e in net.edge)
    new_enum = isempty(unused_enum) ? (net.numedges+1) : unused_enum[1]
    e_new = Edge(new_enum, 0.0)
    max_node = maximum(no.number for no in net.node)
    n_new = Node(max_node+1,false,false,[ce1,ce2,e_new])
    setEdge!(n, e_new)
    setNode!(ce1, n_new) # n_new is second in ce1.node
    ce1.ischild1 = b1    # n_new is parent of ce1 if n was originally parent
    setNode!(ce2, n_new)
    ce2.ischild1 = b2
    setNode!(e_new, [n_new,n]) # ischild1 true, so `n` is the parent
    pushNode!(net, n_new)
    pushEdge!(net, e_new)
    return e_new
end

"""
    resolvetreepolytomy!(net::HybridNetwork, n::Node)

Create new nodes and edges from a tree node `n`, until it has degree 3.
It is best if the network's tree edges are correctly directed beforehand.
"""
function resolvetreepolytomy!(net::HybridNetwork, n::Node)
    !n.hybrid || error("cannot solve polytomy at a hybrid node $(n.number).")
    while size(n.edge,1) > 3
        ce1, ce2, b1, b2 = find_twochildedges(n)
        resolvetreepolytomy!(net, n, ce1, ce2, b1, b2)
    end
end

# tuple of 2 edges incident to n,
# children of n according to ischild1 *if possible*
function find_twochildedges(n::Node)
    ee = n.edge
    ischildofn = [getparent(e)===n for e in ee]
    i2 = findlast(ischildofn)
    if isnothing(i2) # pick edges 1 and 2
        return(ee[1], ee[2], false, false) # false: not children of n
    end
    i1 = findprev(ischildofn, i2-1)
    i1_notfound = isnothing(i1)
    if i1_notfound # set i1 to 1 (or 2 if i2=1)
        i1 = (i2==1 ? 2 : 1)
    end
    return ee[i1], ee[i2], !i1_notfound, true
end

# aux function to add a child to a leaf hybrid
function addChild!(net::HybridNetwork, n::Node)
    n.hybrid || error("cannot add child to tree node $(n.number).")
    ed1 = Edge(net.numedges+1,0.0)
    ed1.containroot = false
    n1num = size(net.names,1)+1
    n1 = Node(n1num,true,false,[ed1])
    n1.name = "$(n.name)_$n1num"
    push!(net.names, n1.name)
    setEdge!(n,ed1)
    setNode!(ed1,[n1,n])
    pushNode!(net,n1)
    pushEdge!(net,ed1)
    if n.leaf
        n.leaf = false
        i = findfirst(isequal(n), net.leaf)
        deleteat!(net.leaf, i)
    end
    return n1
end
# expand the children of a hybrid node: n -newedge-> newnode - -> previous children
function expandChild!(net::HybridNetwork, n::Node)
    n.hybrid || error("cannot expand children of a tree node.")
    suma = count(!e.hybrid for e in n.edge)
    ed1 = Edge(net.numedges+1, 0.0)
    ed1.containroot = false
    n1num = size(net.names,1)+1
    n1 = Node(n1num, false, false, [ed1])
    n1.name = "$(n.name)_$n1num"
    for e in reverse(n.edge) # n.edge is modified inside
        e.hybrid && continue # skip hybrid edges: to skip parents of n
        # now e is a tree edge: *assumes no hybrid ladder*
        # println("se va a borrar a $(e.number)")
        removeEdge!(n,e)
        removeNode!(n,e)
        setEdge!(n1,e)
        setNode!(e,n1)
    end
    #println("now node $(n1.number) has the edges $([e.number for e in n1.edge])")
    setEdge!(n,ed1)
    setNode!(ed1,[n1,n]) # child is n1
    pushNode!(net,n1)
    pushEdge!(net,ed1)
    if size(n1.edge,1) > 3
        resolvetreepolytomy!(net,n1)
    end
end


# function to search for the hybrid nodes in a read network after cleaning it
# and store this information as a network's attribute
function storeHybrids!(net::HybridNetwork)
    flag = true;
    hybrid = nothing
    try
        hybrid = searchHybridNode(net)
    catch
        #@warn "topology read is a tree as it has no hybrid nodes"
        flag = false;
    end
    if(flag)
        net.hybrid = hybrid;
        net.numhybrids = size(hybrid,1);
    end
    return nothing
end


"""
    writesubtree!(IO, network, dendroscope::Bool, namelabel::Bool,
                  round_branch_lengths::Bool, digits::Integer,
                  internallabel::Bool, support=nothing)

Write to IO the extended newick format (parenthetical description)
of a network.
If written for dendroscope, inheritance γ's are not written.
If `namelabel` is true, taxa are labelled by their names;
otherwise taxa are labelled by their number IDs.
If unspecified, branch lengths and γ's are rounded to 3 digits.
Use `internallabel=false` to suppress the labels of internal nodes.
"""
function writesubtree!(
    s::IO,
    net::HybridNetwork,
    di::Bool,
    namelabel::Bool,
    roundBL::Bool,
    digits::Integer,
    internallabel::Bool,
    support=nothing, # default for backward compatibility, in SNaQ.writenewick_level1
)
    rootnode = getroot(net)
    if net.numnodes > 1
        print(s,"(")
        degree = length(rootnode.edge)
        for e in rootnode.edge
            writesubtree!(s,getOtherNode(e,rootnode),e,di,namelabel,roundBL,digits,internallabel,support)
            degree -= 1
            degree == 0 || print(s,",")
        end
        print(s,")")
    end
    if internallabel || net.numnodes == 1
        print(s, (namelabel ? rootnode.name : rootnode.number))
    end
    print(s,";")
    return nothing
end


"""
    writesubtree!(IO, node, edge, dendroscope::Bool, namelabel::Bool,
                  round_branch_lengths::Bool, digits::Integer,
                  internallabel::Bool, support)

Write the extended newick format of the sub-network rooted at
`node` and assuming that `edge` is a parent of `node`.

If the parent `edge` is `nothing`, the edge attribute `ischild1` is used
and assumed to be correct to write the subtree rooted at `node`.
This is useful to write a subtree starting at a non-root node.

In the example below, we run `directedges!` just to be safe
(edges are directed correctly right after `readnewick`), as a reminder
that `writesubtree!` assumes tree edges to be correctly directed away
from the current root (while `writenewick` checks).

```jldoctest
julia> net = readnewick("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);") |> directedges!;

julia> getroot(net).number
-2

julia> crownABC = net.node[7].number # not the root: crown of ABC clade
-3

julia> s = IOBuffer(); writesubtree!(s, crownABC, nothing, false,true);

julia> String(take!(s))
"((A,(B)#H1:::0.9),(C,#H1:::0.1));"

julia> ancestorBC = net.node[6]; # ancestor of B and C, but not stable

julia> writesubtree!(s, ancestorBC, nothing, false,true);

julia> String(take!(s)) # includes only 1 of the 2 partner hybrid edges
"(C,#H1:::0.1);"
```

Used by [`writenewick`](@ref).
"""
writesubtree!(s::IO,n::Node,parent,di,namelabel) = # form never used in PN, but tested
    writesubtree!(s,n,parent,di,namelabel, true,3,true,nothing)

# "parent' is assumed to be adjancent to "node". not checked.
# algorithm comes from "parent": do not traverse again.
function writesubtree!(
    s::IO,
    n::Node,
    parent::Union{Edge,Nothing},
    di::Bool,
    namelabel::Bool,
    roundBL::Bool,
    digits::Integer,
    internallabel::Bool,
    support,
)
    # subtree below node n:
    if !n.leaf && (parent == nothing || parent.ismajor) # do not descent below a minor hybrid edge
        print(s,"(")
        firstchild = true
        for e in n.edge
            e !== parent || continue # skip parent edge where we come from
            if isnothing(parent)     # skip if n = child of e
                n != getchild(e) || continue
            end
            (e.hybrid && getchild(e)==n) && continue # no going up minor hybrid
            firstchild || print(s, ",")
            firstchild = false
            child = getOtherNode(e,n)
            writesubtree!(s,child,e, di,namelabel, roundBL, digits, internallabel, support)
        end
        print(s,")")
    end
    # node label:
    if parent != nothing && parent.hybrid
        print(s, "#")
        print(s, (namelabel ? n.name : string("H", n.number)))
        n.name != "" || parent.ismajor || @warn "hybrid node $(n.number) has no name"
    elseif internallabel || n.leaf
        print(s, (namelabel ? n.name : n.number))
    end
    # branch lengths, support value and γ, if desired and available:
    printBL = false
    if parent != nothing && parent.length != -1.0 # -1.0 means missing
        print(s,string(":",(roundBL ? round(parent.length, digits=digits) : parent.length)))
        printBL = true
    end
    # fixit: print support here
    if !isnothing(parent) && parent.hybrid && !di # && (!printID || !n.booln2))
        if(parent.gamma != -1.0)
            if(!printBL) print(s,":"); end
            print(s,string("::",(roundBL ? round(parent.gamma, digits=digits) : parent.gamma)))
        end
    end
    if isnothing(parent)
        print(s, ";")
    end
end


"""
    readmultinewick(filename::AbstractString, fast=true)
    readmultinewick(newicktrees_list::Vector{<:AbstractString})


Read a list of networks in parenthetical format, either from a file
(one network per line) if the input is a string giving the path
to the file, or from a vector of strings with each string corresponding to
a newick-formatted topology.
By default (`fast=true`), `Functors.fmap` is used for repeatedly
reading the newick trees into of HybridNetwork-type objects.
The option `fast=false` corresponds to the behavior up until v0.14.3:
with a file name as input, it prints a message (without failing) when a
phylogeny cannot be parsed, and allows for empty lines.
Each network is read with [`readnewick`](@ref).

Return an array of HybridNetwork objects.

# Examples

```julia
julia> multitreepath = joinpath(dirname(Base.find_package("PhyloNetworks")), "..", "examples", "multitrees.newick");
julia> multitree = readmultinewick(multitreepath) # vector of 25 HybridNetworks
julia> multitree = readmultinewick(multitreepath, false) # same but slower & safer
julia> treestrings = readlines(multitreepath) # vector of 25 strings
julia> multitree = readmultinewick(treestrings)
julia> readmultinewick(treestrings, false) # same, but slower
```

"""
function readmultinewick(topologies::Vector{<:AbstractString}, fast::Bool=true)
    return (fast ? fmap(readnewick, topologies) : map(readnewick, topologies))
end
function readmultinewick(file::AbstractString, fast::Bool=true)
    if fast
        return readmultinewick(readlines(file), true)
    end
    s = open(file)
    numl = 1
    vnet = HybridNetwork[];
    for line in eachline(s)
        line = strip(line) # remove spaces
        c = isempty(line) ? "" : line[1]
        if(c == '(')
           try
               push!(vnet, readnewick(line))
           catch err
               print("skipped phylogeny on line $(numl) of file $file: ")
               if :msg in fieldnames(typeof(err)) println(err.msg); else println(typeof(err)); end
           end
        end
        numl += 1
    end
    close(s)
    return vnet
end


@doc raw"""
    readnexus_treeblock(filename, treereader=readnewick, args...;
                        reticulate=true, stringmodifier=[r"#(\d+)" => s"#H\1"])

Read the *first* "trees" block of a nexus-formatted file, using the translate
table if present, and return a vector of `HybridNetwork`s.
Information inside `[&...]` are interpreted as comments and are discarded by the
default tree reader. Optional arguments `args` are passed to the tree reader.

For the nexus format, see
[Maddison, Swofford & Maddison (1997)](https://doi.org/10.1093/sysbio/46.4.590).

Unless `reticulate` is false, the following is done to read networks with reticulations.

Prior to reading each phylogeny, each instance of `#number` is replaced by
`#Hnumber` to fit the standard extended Newick format at hybrid nodes.
This behavior can be changed with option `stringmodifier`, which should be a
vector of pairs accepted by `replace`.

Inheritance γ values are assumed to be given within "comment" blocks at *minor*
hybrid edges (cut as tips to form the extended Newick) like this for example,
as output by bacter ([Vaughan et al. 2017](http://dx.doi.org/10.1534/genetics.116.193425)):

    #11[&conv=0, relSize=0.08, ...

or like this, as output by SpeciesNetwork
([Zhang et al. 2018](https://doi.org/10.1093/molbev/msx307)):

    #H11[&gamma=0.08]

In this example, the corresponding edge to hybrid H11 has γ=0.08.
"""
function readnexus_treeblock(
    file::AbstractString,
    treereader::Function=readnewick,
    args...;
    reticulate=true,
    stringmodifier=[r"#(\d+)\b" => s"#H\1"], # add H
)
    vnet = HybridNetwork[]
    rx_start = r"^\s*begin\s+trees\s*;"i
    rx_end = r"^\s*end\s*;"i
    rx_tree = r"^\s*tree\s+[^(]+(\([^;]*;)"i
    treeblock = false
    translate = false
    id2name = nothing
    open(file) do s
        numl = 0
        for line in eachline(s)
            numl += 1
            if treeblock
                occursin(rx_end, line) && break # exit if end of tree block
            elseif occursin(rx_start, line)     # start reading tree block
                    treeblock = true
                    line, translate, id2name = readnexus_translatetable(s)
            else continue
            end
            # now we are inside the treeblock
            m = match(rx_tree, line)
            isnothing(m) && continue
            phy = m.captures[1]
            if reticulate # fix #Hn and extract γ from #Hn[&conv=n, relSize=γ]
                phy = replace(phy, stringmodifier...)
                id2gamma = readnexus_extractgamma(phy)
            end
            net = nothing
            try
                net = treereader(phy, args...)
            catch err
                warnmsg = "skipped phylogeny on line $(numl) of file\n$file\n" *
                    (:msg in fieldnames(typeof(err)) ? err.msg : string(typeof(err)))
                @warn warnmsg
                continue # don't push to vnet
            end
            reticulate && readnexus_assigngammas!(net, id2gamma)
            if translate
                for tip in net.leaf
                    id = parse(Int, tip.name)
                    tip.name = id2name[id]
                end
            end
            push!(vnet, net)
        end
    end
    return vnet
end

"""
    readnexus_translatetable(io)

Read translate table from IO object `io`, whose first non-empty line should contain
"translate". Then each line should have "number name" and the end of the table
is indicated by a ;. Output tuple:
- line that was last read, and is not part of the translate table, taken from `io`
- translate: boolean, whether a table was successfully read
- id2name: dictionary mapping number to name.
"""
function readnexus_translatetable(io)
    rx_translate = r"^\s*translate"i
    rx_emptyline = r"^\s*$"
    line = readline(io)
    translate = false
    id2name = Dict{Int,String}()
    while true
        if occursin(rx_translate, line)
            translate = true
            break
        elseif occursin(rx_emptyline, line)
            line = readline(io)
        else
            translate = false
            break
        end
    end
    if translate # then read the table
        rx_end = r"^\s*;"
        rx_idname = r"\s*(\d+)\s+(\w+)\s*([,;]?)"
        while true
            line = readline(io)
            occursin(rx_emptyline, line) && continue
            if occursin(rx_end, line)
                line = readline(io)
                break
            end
            m = match(rx_idname, line)
            if isnothing(m)
                @warn "problem reading the translate table at line $line.\nnumbers won't be translated to names"
                translate = false
                break
            end
            push!(id2name, parse(Int,m.captures[1]) => String(m.captures[2]))
            if m.captures[3] == ";"
                line = readline(io)
                break
            end
        end
    end
    return line, translate, id2name
end

"""
    readnexus_extractgamma(nexus_string)

Extract γ from comments and return a dictionary hybrid number ID => γ, from
one single phylogeny given as a string.
The output from BEAST2 uses this format for reticulations at *minor* edges,
as output by bacter ([Vaughan et al. 2017](http://dx.doi.org/10.1534/genetics.116.193425)):

    #11[&conv=0, relSize=0.08, ...

or as output by SpeciesNetwork ([Zhang et al. 2018](https://doi.org/10.1093/molbev/msx307)):

    #H11[&gamma=0.08]

The function below assumes that the "H" was already added back if not present
already (from bacter), like this:

    #H11[&conv=0, relSize=0.19, ...

The bacter format is tried first. If this format doesn't give any match,
then the SpeciesNetwork format is tried next.  
See [`readnexus_assigngammas!`](@ref).
"""
function readnexus_extractgamma(nexstring)
    rx_gamma_v1 = r"#H(\d+)\[&conv=\d+,\s*relSize=(\d+\.\d+)"
    rx_gamma_v2 = r"#H(\d+)\[&gamma=(\d+\.\d+)"
    id2gamma = Dict{Int,Float64}()
    # first: try format v1
    for m in eachmatch(rx_gamma_v1, nexstring)
        push!(id2gamma, parse(Int, m.captures[1]) => parse(Float64,m.captures[2]))
    end
    if isempty(id2gamma) # then try format v2
      for m in eachmatch(rx_gamma_v2, nexstring)
        push!(id2gamma, parse(Int, m.captures[1]) => parse(Float64,m.captures[2]))
      end
    end
    return id2gamma
end

"""
    readnexus_assigngammas!(net, d::Dict)

Assign d[i] as the `.gamma` value of the minor parent edge of hybrid "Hi",
if this hybrid node name is found, and if its minor parent doesn't already
have a non-missing γ. See [`readnexus_extractgamma`](@ref)
"""
function readnexus_assigngammas!(net::HybridNetwork, id2gamma::Dict)
    for (i,gam) in id2gamma
        nam = "H$i"
        j = findfirst(n -> n.name == nam, net.hybrid)
        if isnothing(j)
            @warn "didn't find any hybrid node named $nam."
            continue
        end
        hn = net.hybrid[j]
        he = getparentedgeminor(hn)
        if he.gamma == -1.0
            setgamma!(he, gam)
        else
            @warn "hybrid edge number $(he.number) has γ=$(he.gamma). won't erase with $gam."
        end
    end
    return net
end

"""
    writemultinewick(nets, file_name; append=false)
    writemultinewick(nets, IO)

Write an array of networks in parenthetical extended Newick format, one network per line.
Use the option append=true to append to the file. Otherwise, the default is to create a new
file or overwrite it, if it already existed.
Each network is written with `writenewick`.

# Examples
```julia
julia> net = [readnewick("(D,((A,(B)#H7:::0.864):2.069,(F,E):3.423):0.265,(C,#H7:::0.1361111):10);"),
              readnewick("(A,(B,C));"),readnewick("(E,F);"),readnewick("(G,H,F);")];

julia> writemultinewick(net, "fournets.net") # to (over)write to file "fournets.net"
julia> writemultinewick(net, "fournets.net", append=true) # to append to this file
julia> writemultinewick(net, stdout)         # to write to the screen (standard out)
(D,((A,(B)#H7:::0.864):2.069,(F,E):3.423):0.265,(C,#H7:::0.1361111):10.0);
(A,(B,C));
(E,F);
(G,H,F);
```
"""
function writemultinewick(n::Vector{HybridNetwork},file::AbstractString; append::Bool=false)
    mode = (append ? "a" : "w")
    open(file, mode) do s
    writemultinewick(n,s)
    end # closes file safely
end

function writemultinewick(net::Vector{HybridNetwork},s::IO)
    for i in 1:length(net)
      try
        # writeTopologyLevel1(net[i],s,false,true,"none",false,false,3)
        writenewick(net[i],s) # no rounding, not for dendroscope
        write(s,"\n")
      catch err
        if isa(err, RootMismatch) # continue writing other networks in list
            @error "\nError with topology $i:\n" * err.msg
        else rethrow(err); end
      end
    end
end


"""
    writenewick(net; kwargs...)
    writenewick(net, filename; kwargs...)
    writenewick(net, IO, args...)

Write the parenthetical extended Newick format of a network,
as a string, to a file or to an IO buffer / stream.
Optional arguments (default values):

- append (false): if true, appends to the file
- round (false): rounds branch lengths and heritabilities γ
- digits (3): digits after the decimal place for rounding
- di (false): write in format for Dendroscope
- internallabel (true): if true, writes internal node labels
- support (nothing): value to write in the "support" slot for edges

If the current root placement is not admissible, other placements are tried.
The network is updated with this new root placement, if successful.

Uses lower-level function [`writesubtree!`](@ref).
"""
function writenewick(
    n::HybridNetwork,
    file::AbstractString;
    append::Bool=false,
    round::Bool=false,
    digits::Integer=3,
    di::Bool=false,
    internallabel::Bool=true,
    support=nothing,
)
    mode = (append ? "a" : "w")
    s = open(file, mode)
    writenewick(n,s,round,digits,di,internallabel,support)
    write(s,"\n")
    close(s)
end

function writenewick(
    n::HybridNetwork;
    round::Bool=false,
    digits::Integer=3,
    di::Bool=false,
    internallabel::Bool=true,
    support=nothing,
)
    s = IOBuffer()
    writenewick(n,s,round,digits,di,internallabel,support)
    return String(take!(s))
end

function writenewick(
    net::HybridNetwork,
    s::IO,
    round::Bool=false,
    digits::Integer=3,
    di::Bool=false,
    internallabel::Bool=true,
    support=nothing,
)
    # check/find admissible root: otherwise could be trapped in infinite loop
    rootsaved = net.rooti
    changeroot = false
    msg = ""
    try
        directedges!(net)
    catch err
        if isa(err, RootMismatch)
            println(err.msg * "\nCannot write topology with current root.")
            changeroot = true
        else rethrow(err); end
    end
    while changeroot
        for e in net.edge
          # parents of hybrid edges should be sufficient, but gives weird look
          #if e.hybrid
            i = getIndex(getparent(e), net)
            net.rooti = i
            try
                directedges!(net)
                print("Setting root at node $(net.node[i].number) (net.rooti = $i)\n\n")
                print(msg)
                changeroot = false
                break # stop loop over edges
            catch err
                if !isa(err, RootMismatch) rethrow(err); end
            end
          #end
        end
        if changeroot # none of hybrid edges worked
            net.rooti = rootsaved
            throw(RootMismatch("Could not find admissible root. Cannot write topology."))
            changeroot=false # safety exit of while (but useless)
        end
    end
    if net.node[net.rooti].leaf
        @warn """Root is placed at a leaf node, so the parenthetical format will look strange.
                 Use rootatnode! or rootonedge! to change the root position
              """
    end
    # finally, write parenthetical format
    writesubtree!(s,net,di,true,round,digits,internallabel,support)
    # namelabel = true: to print leaf & node names (labels), not numbers
end

"""
    hybridlambdaformat(net::HybridNetwork; prefix="I")

Output `net` as a string in the format that the
[Hybrid-Lambda](https://github.com/hybridLambda/hybrid-Lambda)
simulator expects, namely:
- all internal nodes are named, including the root, with names
  that are unique and start with a letter.
- hybrid nodes are written as `H6#γ1:length1` and `H6#γ1:length2`
  instead of `#H6:length1::γ1` and `#H6:length2::γ2`
  (note the samme γ value expected by Hybrid-Lambda)

This is a modified version of the
[extended Newick](https://doi.org/10.1186/1471-2105-9-532) format.

Optional keyword argument `prefix`: must start with a letter, other than "H".
Internal nodes are given names like "I1", "I2", etc. Existing internal non-hybrid
node names are **replaced**, which is crucial if some of them don't start with a
letter (e.g. in case node names are bootstrap values).
See [`nameinternalnodes!`](@ref) to add node names.

# examples

```jldoctest
julia> net = readnewick("((a:1,(b:1)#H1:1::0.8):5,(#H1:0::0.2,c:1):1);");

julia> hybridlambdaformat(net) # net is unchanged here
"((a:1.0,(b:1.0)H1#0.8:1.0)I1:5.0,(H1#0.8:0.0,c:1.0)I2:1.0)I3;"

julia> # using PhyloPlots; plot(net, shownodenumber=true) # shows that node -2 is the root

julia> rotate!(net, -2)

julia> writenewick(net) # now the minor edge with γ=0.2 appears first
"((#H1:0.0::0.2,c:1.0):1.0,(a:1.0,(b:1.0)#H1:1.0::0.8):5.0);"

julia> hybridlambdaformat(net)
"((H1#0.2:0.0,c:1.0)I2:1.0,(a:1.0,(b:1.0)H1#0.2:1.0)I1:5.0)I3;"

julia> net = readnewick("((((B)#H1:::.6)#H2,((D,C,#H2:::0.8),(#H1,A))));"); # 2 reticulations, no branch lengths

julia> writenewick(net, round=true)
"(#H2:::0.2,((D,C,((B)#H1:::0.6)#H2:::0.8),(#H1:::0.4,A)));"

julia> hybridlambdaformat(net; prefix="int")
"(H2#0.2,((D,C,((B)H1#0.6)H2#0.2)int1,(H1#0.6,A)int2)int3)int4;"
```
"""
function hybridlambdaformat(net::HybridNetwork; prefix="I")
  startswith(prefix, r"[a-zA-GI-Z]") || error("unsafe prefix $prefix: please start with a letter, but not H")
  leafnames = tiplabels(net)
  length(Set(leafnames)) == length(leafnames) || error("taxon names must be unique: $(sort(leafnames))")
  net = deepcopy(net) # binding to new object
  for e in net.edge
    if e.hybrid && e.ismajor && e.gamma == -1.0
      @error("edge number $(e.number) is missing gamma: will use 0.5")
      setgamma!(e, 0.5)
    end
  end
  for no in net.node
    (no.leaf || no.hybrid) && continue # skip leaves & hybrid nodes
    no.name = "" # erase any exisiting name: especially bootstrap values
  end
  nameinternalnodes!(net, prefix)
  str1 = writenewick(net, round=true, digits=15) # internallabels=true by default
  rx_noBL = r"#(H[\w\d]+)::\d*\.?\d*(?:e[+-]?\d+)?:(\d*\.?\d*(?:e[+-]?\d+)?)"
  subst_noBL = s"\1#\2"
  rx_withBL = r"#(H[\w\d]+):(\d*\.?\d*(?:e[+-]?\d+)?):\d*\.?\d*(?:e[+-]?\d+)?:(\d*\.?\d*(?:e[+-]?\d+)?)"
  subst_withBL = s"\1#\3:\2"
  str2 = replace(replace(str1, rx_noBL => subst_noBL), rx_withBL => subst_withBL)
  ## next: replace the γ2 of the second occurrence by γ1 from the first occurrence:
  ## this is what Hybrid-Lambda wants...
  nh = length(net.hybrid)
  m = eachmatch(r"[)(,](H[^#)(,]*#)", str2)
  hboth = collect(h.captures[1] for h in m)
  hone = unique(hboth)
  length(hboth) == 2length(hone) || error("did not find Hname# twice for some one (or more) of the hybrid names.")
  str3 = str2
  for hname in hone
    rx = Regex( hname * raw"(?<gamma1>\d*\.?\d*)(?<middle>.*)" * hname * raw"(?<gamma2>\d*\.?\d*)")
    subst = SubstitutionString(hname * raw"\g<gamma1>\g<middle>" * hname * raw"\g<gamma1>")
    str3 = replace(str3, rx => subst)
  end
  return str3
end

"""
    nameinternalnodes!(net::HybridNetwork, prefix="i")

Add names to nodes in `net` that don't already have a name.
Leaves already have names; but if not, they will be given names as well.
New node names will be of the form "prefixk" where `k` is an integer.
So by default, new node names will be of the form "i1", "i2", etc.

# examples
```jldoctest
julia> net = readnewick("((a:1,(b:1)#H1:1::0.8):5,(#H1:0::0.2,c:1):1);");

julia> nameinternalnodes!(net, "I") # by default, shown without internal node names
HybridNetwork, Rooted Network
7 edges
7 nodes: 3 tips, 1 hybrid nodes, 3 internal tree nodes.
tip labels: a, b, c
((a:1.0,(b:1.0)#H1:1.0::0.8)I1:5.0,(#H1:0.0::0.2,c:1.0)I2:1.0)I3;

julia> writenewick(net; internallabel=false) # by default, writenewick shows internal names if they exist
"((a:1.0,(b:1.0)#H1:1.0::0.8):5.0,(#H1:0.0::0.2,c:1.0):1.0);"

julia> net = readnewick("((int5:1,(b:1)#H1:1::0.8):5,(#H1:0::0.2,c:1):1);"); # one taxon name starts with "int"

julia> nameinternalnodes!(net, "int");

julia> writenewick(net)
"((int5:1.0,(b:1.0)#H1:1.0::0.8)int6:5.0,(#H1:0.0::0.2,c:1.0)int7:1.0)int8;"
```
"""
function nameinternalnodes!(net::HybridNetwork, prefix="i")
  # get maximum index I of nodes whose names are already like: prefixI
  rx = Regex("^$(prefix)(\\d+)\$")
  nexti = 1
  for node in net.node
    node.name != "" || continue # skip nodes with empty names
    m = match(rx, node.name)
    m !== nothing || continue
    nexti = max(nexti, parse(Int, m.captures[1])+1)
  end
  # assign names like: prefixI for I = nexti, nexti+1, etc.
  for node in net.node
    node.name == "" || continue # skip nodes with non-empty names
    node.name = prefix * string(nexti)
    nexti += 1
  end
  return net
end


"""
    readphylip(file::AbstractString)

Read a PHYLIP-formatted DNA sequence alignment (see for example
[wikipedia](https://en.wikipedia.org/wiki/PHYLIP#File_format)).
Currently, only sequential PHYLIP alignments are supported,
i.e., each sequence appears on one line.
The reader makes no assumption about the length of the sequence identifier
so IDs (or taxon names) may be longer than 10 characters. However:
- at least one space is required between the sequence ID and the sequence itself
- spaces are not allowed inside the sequence.

Output: tuple `(species, sequences)`
where `species` is a vector of Strings with identifier names, and
`sequences` is a vector of BioSequences, each of type `sequencetype`
(DNA by default)

See also: [`PhyloNetworks.readfastatoarray`](@ref)
"""
function readphylip(
    file::AbstractString,
    sequencetype=BioSequences.LongDNA{4}
)
    sequences = Array{BioSequences.BioSequence}(undef, 0)
    species = String[]
    open(file) do f
        readline(f)  #skip header information
        while !eof(f)
            line = readline(f)
            line = strip(line)
            isempty(line) && continue
            ind = split(line)
            record=FASTARecord(ind[1],ind[2])
            push!(sequences, FASTX.sequence(sequencetype, record))
            push!(species, FASTX.identifier(record))
        end
    end
    return((species,sequences))
end
