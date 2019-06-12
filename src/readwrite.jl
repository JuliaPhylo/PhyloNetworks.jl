# functions to read/write networks topologies

# peek the next non-white-space char
# removes white spaces from the IOStream/IOBuffer
# see skipchars(predicate, io::IO; linecomment=nothing) in io.jl
# https://github.com/JuliaLang/julia/blob/3b02991983dd47313776091720871201f75f644a/base/io.jl#L971
ncodeunits(c::Char) = write(devnull, c) # for Julia v0.7. Remove for Julia v1.0.3+
function peekskip(io::IO, linecomment=nothing)
    c = missing
    while !eof(io)
        c = read(io, Char)
        if c === linecomment
            readline(io)
        elseif !isspace(c)
            skip(io, -ncodeunits(c))
            break
        end
    end
    return c # skipchar returns io instead
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


# auxiliary function to read a taxon name.
# allows names with letters and numbers: treats numbers as strings
# it also reads # as part of the name and returns pound=true
# it returns the node name as string as well to check if it exists already (as hybrid)
function readnodename(s::IO, c::Char, net::HybridNetwork, numLeft::Array{Int,1})
    if !isValidSymbol(c)
        a = read(s, String);
        error("Expected digit, alphanum or # at the start of taxon name, but received $(c). remaining: $(a).");
    end
    pound = 0
    name = ""
    while isValidSymbol(c)
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
    return size(net.names,1)+1, name, pound == 1
end


# aux function to read floats like length or gamma values, to be read after a colon
function readFloat(s::IO, c::Char)
    if !(isdigit(c) || c in ['.','e','-','E'])
        a = read(s, String);
        error("Expected float digit after ':' but found $(c). remaining is $(a).");
    end
    num = ""
    while isdigit(c) || c in ['.','e','-', 'E']
        d = read(s, Char) # reads c and advances IO
        num = string(num,d);
        c = peekskip(s);
    end
    f = 0.0
    try
        f = parse(Float64, num)
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
function isValidSymbol(c::Char)
    return isletter(c) || isnumeric(c) || c=='_' || c=='#' ||
      (!isspace(c) && c != '(' && c != ')' && c != '[' && c != ']' && c != ':' && c != ';' && c != ',')
end

"""
    parseRemainingSubtree!(s::IO, numLeft, net, hybrids)

Create internal node. Helper for [`readSubtree!`](@ref),
which creates the parent edge of the node created by `parseRemainingSubtree!`:
`readSubtree!` calls `parseRemainingSubtree!`, and vice versa.
Called once a `(` has been read in a tree topology and reads until the corresponding `)` has been found.
This function performs the recursive step for `readSubtree!`.
Advances `s` past the subtree, adds discovered nodes and edges to `net`, and `hybrids`.

Does *not* read the node name and the edge information of the subtree root:
this is done by [`readSubtree!`](@ref)
"""
@inline function parseRemainingSubtree!(s::IO, numLeft::Array{Int,1}, net::HybridNetwork, hybrids::Vector{String})
    numLeft[1] += 1
    DEBUGC && @debug "" numLeft
    n = Node(-1*numLeft[1],false);
    @debug "creating node $(n.number)"
    keepon = true;
    c = readskip!(s)
    while (keepon)
        bl = readSubtree!(s,n,numLeft,net,hybrids)
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
    parseHybridNode!(node, parentNode, hybridName, net, hybrids)

Helper function for `readSubtree!`. Create the parent edge for `node`.
Return this edge, and the hybrid node retained (`node` or its clone in the newick string).
Insert new edge and appropriate node into `net` and `hybrids` accordingly.
Handles any type of given hybrid node.
Called after a `#` has been found in a tree topology.
"""
@inline function parseHybridNode!(n::Node, parent::Node, name::String, net::HybridNetwork, hybrids::Vector{String})
    @debug "found pound in $(name)"
    n.hybrid = true;
    DEBUGC && @debug "got hybrid $(name)"
    DEBUGC && @debug "hybrids list has length $(length(hybrids))"
    ind = findfirst(isequal(name), hybrids) # index of 'name' in the list 'hybrid'. nothing if not found
    e = Edge(net.numEdges+1) # isMajor = true by default
    if n.leaf e.isMajor = false; end
    e.hybrid = true
    e.gamma = -1.0
    if ind !== nothing # the hybrid name was seen before
        @debug "$(name) was found in hybrids list"
        ni = findfirst(isequal(name), [no.name for no in net.node])
        ni !== nothing || error("hybrid name $name was supposed to be in the network, but not found")
        other = net.node[ni]
        @debug "other is $(other.number)"
        DEBUGC && @debug "other is leaf? $(other.leaf), n is leaf? $(n.leaf)"
        if !n.leaf && !other.leaf
            error("both hybrid nodes are internal nodes: successors of the hybrid node must only be included in the node list of a single occurrence of the hybrid node.")
        elseif n.leaf
            @debug "n is leaf"
            @debug "creating hybrid edge $(e.number) attached to other $(other.number) and parent $(parent.number)"
            pushEdge!(net,e);
            setNode!(e,[other,parent]); # isChild1 = true by default constructor
            setEdge!(other,e);
            setEdge!(parent,e);
            n = other # original 'n' dropped, 'other' retained: 'n' output to modify 'n' outside
            @debug "e $(e.number )istIdentifiable? $(e.istIdentifiable)"
        else # !n.leaf : delete 'other' from the network
            @debug "n is not leaf, other is leaf"
            size(other.edge,1) == 1 || # other should be a leaf
               error("strange: node $(other.number) is a leaf hybrid node. should have only 1 edge but has $(size(other.edge,1))")
            DEBUGC && @debug "other is $(other.number), n is $(n.number), edge of other is $(other.edge[1].number)"
            otheredge = other.edge[1];
            otherparent = getOtherNode(otheredge,other);
            @debug "otheredge is $(otheredge.number)"
            @debug "parent of other is $(otherparent.number)"
            removeNode!(other,otheredge);
            deleteNode!(net,other);
            setNode!(otheredge,n);
            setEdge!(n,otheredge);
            ## otheredge.istIdentifiable = true ## setNode should catch this, but when fixed, causes a lot of problems
            @debug "setting otheredge to n $(n.number)"
            @debug "creating hybrid edge $(e.number) between n $(n.number) and parent $(parent.number)"
            setNode!(e,[n,parent]);
            setEdge!(n,e);
            setEdge!(parent,e);
            pushNode!(net,n);
            pushEdge!(net,e);
            n.number = other.number; # modifies original negative node number, to positive node #
            n.name = other.name;
            @debug "edge $(e.number) istIdentifiable? $(e.istIdentifiable)"
            @debug "otheredge $(otheredge.number) istIdentifiable? $(otheredge.istIdentifiable)"
        end
    else # ind==nothing: hybrid name not seen before
        @debug "$(name) not found in hybrids list"
        @debug "$(name) is leaf? $(n.leaf)"
        n.hybrid = true;
        nam = string(name)
        push!(net.names, nam);
        n.name = nam;
        DEBUGC && @debug "put $(nam) in hybrids name list"
        push!(hybrids, nam);
        pushNode!(net,n);
        @debug "creating hybrid edge $(e.number)"
        pushEdge!(net,e);
        setNode!(e,[n,parent]);
        setEdge!(n,e);
        setEdge!(parent,e);
        @debug "edge $(e.number) istIdentifiable? $(e.istIdentifiable)"
    end
    e.containRoot = !e.hybrid # not good: but necessay for SNaQ functions
    return (e,n)
end

"""
    parseTreeNode!(node, parentNode, net)

Helper function for `readSubtree!`.
Insert the input tree node and associated edge (created here) into `net`.
"""
@inline function parseTreeNode!(n::Node, parent::Node, net::HybridNetwork)
    pushNode!(net,n);
    e = Edge(net.numEdges+1);
    pushEdge!(net,e);
    setNode!(e,[n,parent]);
    setEdge!(n,e);
    setEdge!(parent,e);
    return e
end

"""
    getdataValue!(s::IO, int, numLeft::Array{Int,1})

Helper function for `parseEdgeData!`.
Read a single floating point edge data value in a tree topology.
Return -1.0 if no value exists before the next colon, return the value as a float otherwise.
Modifies s by advancing past the next colon character.
Only call this function to read a value when you know a numerical value exists!
"""
@inline function getDataValue!(s::IO, call::Int, numLeft::Array{Int,1})
    errors = ["one colon read without double in left parenthesis $(numLeft[1]-1), ignored.",
              "second colon : read without any double in left parenthesis $(numLeft[1]-1), ignored.",
              "third colon : without gamma value after in $(numLeft[1]-1) left parenthesis, ignored"]
    c = peekskip(s)

    # Value is present
    if isdigit(c) || c == '.'
        val = readFloat(s, c)
        return val
    # No value
    elseif c == ':'
        return -1.0
    else
        @warn errors[call]
        return -1.0
    end
end

"""
    parseEdgeData!(s::IO, edge, node, numberOfLeftParentheses::Array{Int,1})

Helper function for readSubtree!, fixes a bug from using setGamma
Modifies `e` according to the specified edge length and gamma values in the tree topology.
Advances the stream `s` past any existing edge data.
Edges in a topology may optionally be followed by ":edgeLen:bootstrap:gamma"
where edgeLen, bootstrap, and gamma are decimal values.
"""
@inline function parseEdgeData!(s::IO, e::Edge, n::Node, numLeft::Array{Int,1})
    read(s, Char); # to read the first ":"
    e.length = getDataValue!(s, 1, numLeft)
    bootstrap = nothing;
    if peekskip(s) == ':'
        readskip!(s)
        bootstrap = getDataValue!(s, 2, numLeft)
    end
    # e.gamma = -1.0 by default when e is created by parseHybridNode!
    if peekskip(s) == ':'
        readskip!(s)
        e.gamma = getDataValue!(s, 3, numLeft)
    end
    if e.gamma != -1.0 && !e.hybrid && e.gamma != 1.0
        @warn "γ read for edge $(e.number) but it is not hybrid, so γ=$(gamma) ignored"
        e.gamma = 1.0
    end
end

"""
    synchronizePartnersData!(e::Edge, n::Node)

Synchronize γ and isMajor for edges `e` and its partner,
both hybrid edges with the same child `n`:

- if one γ is missing and the other is not: set the missing γ to 1 - the other
- γ's should sum up to 1.0
- update `isMajor` to match the γ information: the major edge is the one with γ > 0.5.

**Warnings**: does not check that `e` is a hybrid edge,
nor that `n` is the child of `e`.
"""
@inline function synchronizePartnersData!(e::Edge, n::Node)
    partners = Edge[] # The edges having n as a child, other than e
    for e2 in n.edge
        if e2.hybrid && e2!=e && n==getChild(e2)
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
    # update γ and isMajor of both edges, to be consistent with each other
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
    # next: update isMajor, originally based on which node was a leaf and which was not
    # if both γ are 0.5: keep isMajor as is. Otherwise: γ's take precedence.
    if !isapprox(e.gamma, 0.5)
        emajor = e.gamma > 0.5
        if e.isMajor != emajor # leaf status was inconsistent with γ info
            e.isMajor = emajor
            partner.isMajor = !emajor
        end
    end
end

"""
    readSubtree!(s::IO, parentNode, numLeft, net, hybrids)

Recursive helper method for `readTopology`:
read a subtree from an extended Newick topology.
input `s`: IOStream/IOBuffer.

Reads additional info formatted as: `:length:bootstrap:gamma`.
Allows for name of internal nodes without # after closing parenthesis: (1,2)A.
Warning if hybrid edge without γ, or if γ (ignored) without hybrid edge
"""
function readSubtree!(s::IO, parent::Node, numLeft::Array{Int,1}, net::HybridNetwork, hybrids::Vector{String})
    c = peekskip(s)
    e = nothing;
    hasname = false; # to know if the current node has name
    pound = false;
    if c == '('
        # read the rest of the subtree (perform the recursive step!)
        n = parseRemainingSubtree!(s, numLeft, net, hybrids)
        c = peekskip(s);
        if isValidSymbol(c) # internal node has name
            hasname = true;
            num, name, pound = readnodename(s, c, net, numLeft);
            n.number = num; # n was given <0 number by parseRemainingSubtree!, now >0
            c = peekskip(s);
        end
    else # leaf, it should have a name
        hasname = true;
        num, name, pound = readnodename(s, c, net, numLeft)
        n = Node(num, true); # positive node number to leaves in the newick-tree description
        @debug "creating node $(n.number)"
    end
    if pound # found pound sign in name
        # merge the 2 nodes corresponding the hybrid: make n the node that is retained
        e,n = parseHybridNode!(n, parent, name, net, hybrids) # makes hybrid node number >0
    else
        if hasname
            push!(net.names, name);
            n.name = name
        end
        e = parseTreeNode!(n, parent, net)
    end
    c = peekskip(s);
    e.length = -1.0
    if c == ':'
        parseEdgeData!(s, e, n, numLeft)
    end
    if e.hybrid
        # if hybrid edge: 'e' might have no info, but its partner may have had info
        synchronizePartnersData!(e, n) # update γ and isMajor of e and/or its partner
    end
    return true
end


# function to read topology from parenthetical format
# input: file name or tree in parenthetical format
# calls readTopology(s::IO)
# warning: crashes if file name starts with (
function readTopology(input::AbstractString,verbose::Bool)
    if(input[1] == '(') # input = parenthetical description
       s = IOBuffer(input)
    else # input = file name
        try
            s = open(input)
        catch
            error("Could not find or open $(input) file");
        end
       s = open(input)
    end
    net = readTopology(s,verbose)
    return net
end

"""
    readTopology(file name)
    readTopology(parenthetical description)

Read tree or network topology from parenthetical format (extended Newick).
If the root node has a single child: ignore (i.e. delete from the topology)
the root node and its child edge.

Input: text file or parenthetical format directly.
The file name may not start with a left parenthesis, otherwise the file
name itself would be interpreted as the parenthetical description.

A root edge, not enclosed within a pair a parentheses, is ignored.
If the root node has a single edge, this one edge is removed.
"""
readTopology(input::AbstractString) = readTopology(input,true)

function readTopology(s::IO,verbose::Bool)
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
            b |= readSubtree!(s,n,numLeft,net,hybrids)
            c = readskip!(s);
            if eof(s)
                error("Tree ended while reading in subtree beginning with left parenthesis number $(numLeft[1]-1).")
            elseif c == ','
                continue;
            elseif c == ')'
                c = peekskip(s);
                if isValidSymbol(c) # the root has a name
                    num, name, pound = readnodename(s, c, net, numLeft);
                    n.name = name
                    # log warning or error if pound > 0?
                    c = peekskip(s);
                end
                if(c == ':') # skip information on the root edge, if it exists
                    # @warn "root edge ignored"
                    while c != ';'
                        c = readskip!(s)
                    end
                end
            end
        end
        @debug "after readsubtree:"
        @debug begin printEdges(net); "printed edges" end
        # delete the root edge, if present
        if size(n.edge,1) == 1 # root node n has only one edge
            edge = n.edge[1]
            child = getOtherNode(edge,n);
            removeEdge!(child,edge);
            net.root = getIndex(child,net);
            deleteEdge!(net,edge);
        else
            pushNode!(net,n);
            net.root = getIndex(n,net);
        end
    else
		a = read(s, String)
        error("Expected beginning of tree with ( but received $(c) instead, rest is $(a)")
    end
    storeHybrids!(net)
    checkNumHybEdges!(net)
       ## if(verbose)
       ##     any([(e.length == -1.0 && e.istIdentifiable) for e in net.edge]) && println("edges lengths missing, so assigned default value of -1.0") #fixit: not best approach, better to add a flag inside readSubTree, careful bool not modified inside function, need bool array
       ## end
    net.isRooted = true
    return net
end

readTopology(s::IO) = readTopology(s,true)

"""
    `checkNumHybEdges!(net)`

Check for consistency between hybrid-related attributes in the network:
- for each hybrid node: 2 or more hybrid edges
- exception: allows for a leaf to be attached to a single hybrid edge
- exactly 2 incoming parent hybrid edges
Run after `storeHybrids!`. See also `check2HybEdges`.
"""
function checkNumHybEdges!(net::HybridNetwork)
    if isTree(net) return nothing; end
    !isempty(net.hybrid) || error("net.hybrid should not be empty for this network")
    for n in net.hybrid
        hyb = sum([e.hybrid for e in n.edge]); # number of hybrid edges attached to node
        if hyb == 1
            if net.numHybrids == 1
                error("only one hybrid node $(n.number) named $(n.name) found with one hybrid edge attached")
            else
                error("hybrid node $(n.number) named $(n.name) has only one hybrid edge attached. there are $(net.numHybrids-1) other hybrids out there but this one remained unmatched")
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
                if n == getChild(e)
                    if e.hybrid
                        nhybparents += 1
                    else @error "node $(n.number) has parent tree edge $(e.number): wrong isChild1 for this edge?"
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

# aux function to send an error if the number of hybrid attached to every
# hybrid node is >2
function check2HybEdges(net::HybridNetwork)
    for n in net.hybrid
        hyb = sum([e.hybrid for e in n.edge]);
        if hyb > 2
            error("hybrid node $(n.number) has more than two hybrid edges attached to it: polytomy that cannot be resolved without intersecting cycles.")
        end
    end
end


# aux function to solve a polytomy
# warning: chooses one resolution at random
function solvePolytomyRecursive!(net::HybridNetwork, n::Node)
    if(size(n.edge,1) == 4)
        edge1 = n.edge[1];
        edge2 = n.edge[2];
        edge3 = n.edge[3];
        edge4 = n.edge[4];
        removeEdge!(n,edge3);
        removeEdge!(n,edge4);
        removeNode!(n,edge3);
        removeNode!(n,edge4);
        ednew = Edge(net.numEdges+1,0.0);
        max_node = maximum([e.number for e in net.node]);
        n1 = Node(max_node+1,false,false,[edge3,edge4,ednew]);
        setEdge!(n,ednew);
        setNode!(edge3,n1);
        setNode!(edge4,n1);
        setNode!(ednew,[n,n1]);
        pushNode!(net,n1);
        pushEdge!(net,ednew);
    else
        edge1 = n.edge[1];
        removeEdge!(n,edge1);
        solvePolytomyRecursive!(net,n);
        setEdge!(n,edge1);
    end
end

# function to solve a polytomy among tree edges recursively
function solvePolytomy!(net::HybridNetwork, n::Node)
    !n.hybrid || error("cannot solve polytomy in a hybrid node $(n.number).")
    while(size(n.edge,1) > 3)
        solvePolytomyRecursive!(net,n);
    end
end

# aux function to add a child to a leaf hybrid
function addChild!(net::HybridNetwork, n::Node)
    n.hybrid || error("cannot add child to tree node $(n.number).")
    ed1 = Edge(net.numEdges+1,0.0);
    n1 = Node(size(net.names,1)+1,true,false,[ed1]);
    setEdge!(n,ed1);
    setNode!(ed1,[n,n1]);
    pushNode!(net,n1);
    pushEdge!(net,ed1);
end
# aux function to expand the children of a hybrid node
function expandChild!(net::HybridNetwork, n::Node)
    if n.hybrid
        suma = count([!e.hybrid for e in n.edge]);
        #println("create edge $(net.numEdges+1)")
        ed1 = Edge(net.numEdges+1,0.0);
        n1 = Node(size(net.names,1)+1,false,false,[ed1]);
        #println("create node $(n1.number)")
        hyb = Edge[];
        for i in 1:size(n.edge,1)
            if !n.edge[i].hybrid push!(hyb,n.edge[i]); end
        end
        #println("hyb tiene $([e.number for e in hyb])")
        for e in hyb
            #println("se va a borrar a $(e.number)")
            removeEdge!(n,e);
            removeNode!(n,e);
            setEdge!(n1,e);
            setNode!(e,n1);
        end
        #println("now node $(n1.number) has the edges $([e.number for e in n1.edge])")
        setEdge!(n,ed1);
        setNode!(ed1,[n,n1]);
        pushNode!(net,n1);
        pushEdge!(net,ed1);
        if size(n1.edge,1) > 3
            solvePolytomy!(net,n1);
        end
    else
        error("cannot expand children of a tree node.")
    end
end

# function to clean topology after readTopology
# looks for:
# TREE:
# - all tree edges must have gamma=1. fixit: cannot point out which doesn't,
#   only shows error.
# - internal nodes with only 2 edges and solves this case.
# - polytomies and choose one resolution at random, issuing a warning
# NETWORK:
# - number of hybrid edges per hybrid node:
#   if 0,1: error (with warning in old functions)
#   if >2: error of hybrid polytomy
#   if 2: check number of tree edges
# - number of tree edges per hybrid node:
#   if 0: leaf hybrid, add child
#   if >1: expand child
#   if 1: check values of gamma:
# - gammas: need to sum to one and be present.
#   error if they do not sum up to one
#   default values of 0.1,0.9 if not present
# leaveRoot=true: leaves the root even if it has only 2 edges (for plotting), default=false
function cleanAfterRead!(net::HybridNetwork, leaveRoot::Bool)
    nodes = copy(net.node)
    for n in nodes
        if isNodeNumIn(n,net.node) # very important to check
            if size(n.edge,1) == 2
                if !n.hybrid
                    if !leaveRoot || !isEqual(net.node[net.root],n) #if n is the root
                        deleteIntNode!(net,n);
                    end
                else
                    hyb = count([e.hybrid for e in n.edge]);
                    if hyb == 1
                        deleteIntNode!(net,n);
                    end
                end
            end
            if !n.hybrid
                if size(n.edge,1) > 3
                    @debug "warning: polytomy found in node $(n.number), random resolution chosen"
                    solvePolytomy!(net,n);
                end
                hyb = count([e.hybrid for e in n.edge])
                if hyb == 1
                    n.hasHybEdge == true;
                elseif hyb > 1
                    @warn "strange tree node $(n.number) with more than one hybrid edge, intersecting cycles maybe"
                end
            else
                hyb = count([e.hybrid for e in n.edge]);
                tre = length(n.edge) - hyb
                if hyb > 2
                    error("hybrid node $(n.number) has more than two hybrid edges attached to it: polytomy that cannot be resolved without intersecting cycles.")
                elseif hyb == 1
                    hybnodes = count([n.hybrid for n in net.node]);
                    if hybnodes == 1
                        error("only one hybrid node number $(n.number) with name $(net.names[n.number]) found with one hybrid edge attached")
                    else
                        error("current hybrid node $(n.number) with name $(net.names[n.number]) has only one hybrid edge attached. there are other $(hybnodes-1) hybrids out there but this one remained unmatched")
                    end
                elseif hyb == 0
                    @warn "hybrid node $(n.number) is not connected to any hybrid edges, it was transformed to tree node"
                    n.hybrid = false;
                else # 2 hybrid edges
                    if tre == 0 #hybrid leaf
                        @warn "hybrid node $(n.number) is a leaf, so we add an extra child"
                        addChild!(net,n);
                    elseif tre > 1
                        @warn "hybrid node $(n.number) has more than one child so we need to expand with another node"
                        expandChild!(net,n);
                    end
                    suma = sum([e.hybrid ? e.gamma : 0.0 for e in n.edge]);
                    # synchronizePartnersData! already made suma ≈ 1.0, when non-missing,
                    # and already synchronized isMajor, even when γ's ≈ 0.5
                    if suma == -2.0 # hybrid edges have no gammas in newick description
                        println("hybrid edges for hybrid node $(n.number) have missing gamma's, set default: 0.9,0.1")
                        for e in n.edge
                            if e.hybrid
                                e.gamma = (e.isMajor ? 0.9 : 0.1)
                            end
                        end
                    end
                end
            end
        end
    end
    for e in net.edge
        if e.hybrid
          if e.gamma < 0.0 || e.gamma > 1.0 # no -1.0 (missing) γ's by now
            error("hybrid edge $(e.number) with γ = $(e.gamma) outside of [0,1]")
          end
        else
          e.gamma == 1.0 ||
            error("tree edge $(e.number) with γ = $(e.gamma) instead of 1.0")
        end
    end
end

cleanAfterRead!(net::HybridNetwork) = cleanAfterRead!(net,false)

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
        net.numHybrids = size(hybrid,1);
    end
    return nothing
end

# function to update the read topology after reading
# it will go over the net.hybrid array and check each
# of the hybridization events defined to update:
# - in cycle
# - contain root
# - gammaz
# it uses updateAllNewHybrid! function that
# returns: success (bool), hybrid, flag, nocycle, flag2, flag3
# if tree read, also check that contain root is true for all, ishybrid and hashybedge is false
# warning: needs to have run storeHybrids! before
# warning: it will stop when finding one conflicting hybrid
function updateAllReadTopology!(net::HybridNetwork)
    if(isTree(net))
        #@warn "not a network read, but a tree as it does not have hybrid nodes"
        all((e->e.containRoot), net.edge) ? nothing : error("some tree edge has contain root as false")
        all((e->!e.hybrid), net.edge) ? nothing : error("some edge is hybrid and should be all tree edges in a tree")
        all((n->!n.hasHybEdge), net.node) ? nothing : error("some tree node has hybrid edge true, but it is a tree, there are no hybrid edges")
    else
        if(!net.cleaned)
            for n in net.hybrid
                success,hyb,flag,nocycle,flag2,flag3 = updateAllNewHybrid!(n,net,false,true,false)
                if(!success)
                    @warn "current hybrid $(n.number) conflicts with previous hybrid by intersecting cycles: $(!flag), nonidentifiable topology: $(!flag2), empty space for contain root: $(!flag3), or does not create a cycle (probably problem with the root placement): $(nocycle)."
                    #net.cleaned = false
                end
            end
            @debug "before update partition"
            @debug begin printPartitions(net); "printed partitions" end
            for n in net.hybrid #need to updatePartition after all inCycle
                nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,n);
                updatePartition!(net,nodesInCycle)
                @debug begin printPartitions(net);
                    "partitions after updating partition for hybrid node $(n.number)"
                end
            end
        end
    end
end

# cleanAfterReadAll includes all the step to clean a network after read
function cleanAfterReadAll!(net::HybridNetwork, leaveRoot::Bool)
    @debug "check for 2 hybrid edges at each hybrid node -----"
    check2HybEdges(net)
    @debug "cleanBL -----"
    cleanBL!(net)
    @debug "cleanAfterRead -----"
    cleanAfterRead!(net,leaveRoot)
    @debug "updateAllReadTopology -----"
    updateAllReadTopology!(net) #fixit: it could break if leaveRoot = true (have not checked it), but we need to updateContainRoot
    if(!leaveRoot)
        @debug "parameters -----"
        parameters!(net)
    end
    @debug "check root placement -----"
    checkRootPlace!(net)
    net.node[net.root].leaf && @warn "root node $(net.node[net.root].number) is a leaf, so when plotting net, it can look weird"
    net.cleaned = true #fixit: set to false inside updateAllReadTopology if problem encountered
    net.isRooted = false
end

cleanAfterReadAll!(net::HybridNetwork) = cleanAfterReadAll!(net,false)

# function to read a topology from file name/tree directly and update it
# by calling updateAllReadTopology after
# leaveRoot=true if the root will not be deleted even if it has only 2 edges
# used for plotting (default=false)
# warning: if leaveRoot=true, net should not be used outside plotting, things will crash
function readTopologyUpdate(file::AbstractString, leaveRoot::Bool,verbose::Bool)
    @debug "readTopology -----"
    net = readTopology(file,verbose)
    cleanAfterReadAll!(net,leaveRoot)
    return net
end

readTopologyUpdate(file::AbstractString) = readTopologyUpdate(file, false, true)
readTopologyUpdate(file::AbstractString,verbose::Bool) = readTopologyUpdate(file, false, verbose)

"""
    readTopologyLevel1(filename)
    readTopologyLevel1(parenthetical format)

same as readTopology, reads a tree or network from parenthetical
format, but this function enforces the necessary conditions for any
starting topology in SNaQ: non-intersecting cycles, no polytomies,
unrooted. It sets any missing branch length to 1.0.

If the network has a bad diamond II (in which edge lengths are γ's are not identifiable)
and if the edge below this diamond has a length `t` different from 0, then this length is
set back to 0 and the major parent hybrid edge is lengthened by `t`.
"""
readTopologyLevel1(file::AbstractString) = readTopologyUpdate(file, false, true)


# aux function to check if the root is placed correctly, and re root if not
# warning: it needs updateContainRoot set
function checkRootPlace!(net::HybridNetwork; verbose=false::Bool, outgroup="none"::AbstractString)
    if(outgroup == "none")
        if(!canBeRoot(net.node[net.root]))
            verbose && println("root node $(net.node[net.root].number) placement is not ok, we will change it to the first found node that agrees with the direction of the hybrid edges")
            for i in 1:length(net.node)
                if(canBeRoot(net.node[i]))
                    net.root = i
                    break
                end
            end
        end
    else # outgroup
        tmp = findall(n -> n.name == outgroup, net.leaf)
        if length(tmp)==0
            error("leaf named $(outgroup) was not found in the network.")
        elseif length(tmp)>1
            error("several leaves were found with name $(outgroup).")
        end
        leaf = net.leaf[tmp[1]]
        leaf.leaf || error("found outgroup not a leaf: $(leaf.number), $(outgroup)")
        length(leaf.edge) == 1 || error("found leaf with more than 1 edge: $(leaf.number)")
        other = getOtherNode(leaf.edge[1],leaf);
        if(canBeRoot(other))
            net.root = getIndexNode(other.number,net)
        else
            throw(RootMismatch("outgroup $(outgroup) contradicts direction of hybrid edges"))
        end
    end
    canBeRoot(net.node[net.root]) || error("tried to place root, but couldn't. root is node $(net.node[net.root])")
end

"""
    writeSubTree!(IO, network, dendroscope::Bool, namelabel::Bool,
                  round_branch_lengths::Bool, digits::Integer,
                  internallabel::Bool)

Write to IO the extended newick format (parenthetical description)
of a network.
If written for dendroscope, inheritance γ's are not written.
If `namelabel` is true, taxa are labelled by their names;
otherwise taxa are labelled by their number IDs.
If unspecified, branch lengths and γ's are rounded to 3 digits.
Use `internallabel=false` to suppress the labels of internal nodes.
"""
function writeSubTree!(s::IO, net::HybridNetwork, di::Bool, namelabel::Bool,
                       roundBL::Bool, digits::Integer, internallabel::Bool)
    rootnode = net.node[net.root]
    if net.numNodes > 1
        print(s,"(")
        degree = length(rootnode.edge)
        for e in rootnode.edge
            writeSubTree!(s,getOtherNode(e,rootnode),e,di,namelabel,roundBL,digits,internallabel)
            degree -= 1
            degree == 0 || print(s,",")
        end
        print(s,")")
    end
    if internallabel || net.numNodes == 1
        print(s, (namelabel ? rootnode.name : rootnode.number))
    end
    print(s,";")
    return nothing
end


"""
    writeSubTree!(IO, node, edge, dendroscope::Bool, namelabel::Bool,
                  round_branch_lengths::Bool, digits::Integer, internallabel::Bool)

Write the extended newick format of the sub-network rooted at
`node` and assuming that `edge` is a parent of `node`.

If the parent `edge` is `nothing`, the edge attribute `isChild1` is used
and assumed to be correct to write the subtree rooted at `node`.
This is useful to write a subtree starting at a non-root node.
Example:

```julia
net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
directEdges!(net)
s = IOBuffer()
writeSubTree!(s, net.node[7], nothing, false, true)
String(take!(s))
```

Used by [`writeTopology`](@ref).
"""
writeSubTree!(s,n,parent,di,namelabel) =
    writeSubTree!(s,n,parent,di,namelabel, true,3,true)

# "parent' is assumed to be adjancent to "node". not checked.
# algorithm comes from "parent": do not traverse again.
function writeSubTree!(s::IO, n::Node, parent::Union{Edge,Nothing},
    di::Bool, namelabel::Bool, roundBL::Bool, digits::Integer, internallabel::Bool)
    # subtree below node n:
    if !n.leaf && (parent == nothing || parent.isMajor) # do not descent below a minor hybrid edge
        print(s,"(")
        firstchild = true
        for e in n.edge
            e != parent || continue # skip parent edge where we come from
            if parent == nothing    # skip if n = child of e
                n != getChild(e) || continue
            end
            (e.hybrid && getChild(e)==n) && continue # no going up minor hybrid
            firstchild || print(s, ",")
            firstchild = false
            child = getOtherNode(e,n)
            writeSubTree!(s,child,e, di,namelabel, roundBL, digits, internallabel)
        end
        print(s,")")
    end
    # node label:
    if parent != nothing && parent.hybrid
        print(s, "#")
        print(s, (namelabel ? n.name : string("H", n.number)))
        n.name != "" || parent.isMajor || @warn "hybrid node $(n.number) has no name"
    elseif internallabel || n.leaf
        print(s, (namelabel ? n.name : n.number))
    end
    # branch lengths and γ, if available:
    printBL = false
    if parent != nothing && parent.length != -1.0 # -1.0 means missing
        print(s,string(":",(roundBL ? round(parent.length, digits=digits) : parent.length)))
        printBL = true
    end
    if parent != nothing && parent.hybrid && !di # && (!printID || !n.isBadDiamondI))
        if(parent.gamma != -1.0)
            if(!printBL) print(s,":"); end
            print(s,string("::",(roundBL ? round(parent.gamma, digits=digits) : parent.gamma)))
        end
    end
    if parent == nothing
        print(s, ";")
    end
end

# see full docstring below
# Need HybridNetwork input, since QuartetNetwork does not have root.
function writeTopologyLevel1(net0::HybridNetwork, di::Bool, str::Bool, namelabel::Bool,outgroup::AbstractString, printID::Bool, roundBL::Bool, digits::Integer, multall::Bool)
    s = IOBuffer()
    writeTopologyLevel1(net0,s,di,namelabel,outgroup,printID,roundBL,digits, multall)
    if str
        return String(take!(s))
    else
        return s
    end
end

# warning: I do not want writeTopologyLevel1 to modify the network if outgroup is given! thus, we have updateRoot, and undoRoot
function writeTopologyLevel1(net0::HybridNetwork, s::IO, di::Bool, namelabel::Bool,
           outgroup::AbstractString, printID::Bool, roundBL::Bool, digits::Integer, multall::Bool)
    global CHECKNET
    net = deepcopy(net0) #writeTopologyLevel1 needs containRoot, but should not alter net0
    # net.numBad == 0 || println("network with $(net.numBad) bad diamond I. Some γ and edge lengths t are not identifiable, although their γ * (1-exp(-t)) are.")
    if printID
        setNonIdBL!(net) # changes non identifiable BL to -1.0, except those in/below bad diamonds/triangles.
    end
    assignhybridnames!(net)
    if(net.numNodes == 1)
        print(s,string(net.node[net.root].number,";")) # error if 'string' is an argument name.
    else
        if(!isTree(net) && !net.cleaned)
            @debug "net not cleaned inside writeTopologyLevel1, need to run updateContainRoot"
            for n in net.hybrid
                flag,edges = updateContainRoot!(net,n)
                flag || error("hybrid node $(n.hybrid) has conflicting containRoot")
            end
        end
        updateRoot!(net,outgroup)
        #@debug begin printEverything(net); "printed everything" end
        CHECKNET && canBeRoot(net.node[net.root])
        if(multall)
            mergeLeaves!(net)
        end
        writeSubTree!(s, net, di,namelabel, roundBL,digits,false)
    end
    # outgroup != "none" && undoRoot!(net) # not needed because net is deepcopy of net0
    # to delete 2-degree node, for snaq.
end

writeTopologyLevel1(net::HybridNetwork,di::Bool,str::Bool,namelabel::Bool,outgroup::AbstractString,printID::Bool) = writeTopologyLevel1(net,di,str,namelabel,outgroup,printID, false,3, false)
# above: default roundBL=false (at unused digits=3 decimal places)
writeTopologyLevel1(net::HybridNetwork,printID::Bool) = writeTopologyLevel1(net,false, true,true,"none",printID, false, 3, false)
writeTopologyLevel1(net::HybridNetwork,outgroup::AbstractString) = writeTopologyLevel1(net,false, true,true,outgroup,true, false, 3, false)
writeTopologyLevel1(net::HybridNetwork,di::Bool,outgroup::AbstractString) = writeTopologyLevel1(net,di, true,true,outgroup,true, false, 3, false)

"""
`writeTopologyLevel1(net::HybridNetwork)`

Write the extended Newick parenthetical format of a
level-1 network object with many optional arguments (see below).
Makes a deep copy of net: does *not* modify `net`.

- di=true: write in format for Dendroscope (default false)
- namelabel=true: If `namelabel` is true, taxa are labelled by their names;
otherwise taxa are labelled by their numbers (unique identifiers).
- outgroup (string): name of outgroup to root the tree/network.
  if "none" is given, the root is placed wherever possible.
- printID=true, only print branch lengths for identifiable egdes
  according to the snaq estimation procedure (default false)
  (true inside of `snaq!`.)
- round: rounds branch lengths and heritabilities γ (default: true)
- digits: digits after the decimal place for rounding (defult: 3)
- string: if true (default), returns a string,
  otherwise returns an IOBuffer object.
- multall: (default false). set to true when there are multiple
  alleles per population.

The topology may be written using a root different than net.root,
if net.root is incompatible with one of more hybrid node.
Missing hybrid names are written as "#Hi" where "i" is the hybrid node number if possible.
""" #"
writeTopologyLevel1(net::HybridNetwork; di=false::Bool, string=true::Bool, namelabel=true::Bool,
    outgroup="none"::AbstractString, printID=false::Bool, round=false::Bool, digits=3::Integer,
    multall=false::Bool) =
writeTopologyLevel1(net, di, string, namelabel, outgroup, printID, round, digits, multall)

# function to check if root is well-placed
# and look for a better place if not
# searches on net.node because net.root is the index in net.node
# if we search in net.edge, we then need to search in net.node
function updateRoot!(net::HybridNetwork, outgroup::AbstractString)
    checkroot = false
    if(outgroup == "none")
        @debug "no outgroup defined"
        checkroot = true
    else
        println("outgroup defined $(outgroup)")
        try
            index = getIndex(true,[isequal(outgroup,n.name) for n in net.node])
        catch
            error("outgroup $(outgroup) not in net.names $(net.names)")
        end
        index = getIndex(true,[isequal(outgroup,n.name) for n in net.node])
        node = net.node[index]
        node.leaf || error("outgroup $(outgroup) is not a leaf in net")
        length(net.node[index].edge) == 1 || error("strange leaf $(outgroup), node number $(net.node[index].number) with $(length(net.node[index].edge)) edges instead of 1")
        edge = net.node[index].edge[1]
        if(edge.containRoot)
            DEBUGC && @debug "creating new node in the middle of the external edge $(edge.number) leading to outgroup $(node.number)"
            othernode = getOtherNode(edge,node)
            removeEdge!(othernode,edge)
            removeNode!(othernode,edge)
            max_edge = maximum([e.number for e in net.edge]);
            max_node = maximum([e.number for e in net.node]);
            newedge = Edge(max_edge+1) #fixit: maybe this edge not identifiable, need to add that check
            newnode = Node(max_node+1,false,false,[edge,newedge])
            if(net.cleaned && !isTree(net) && !isempty(net.partition)) # fixit: this will crash if network estimated with snaq, and then manipulated
                part = whichPartition(net,edge)
                push!(net.partition[part].edges,newedge)
            end
            setNode!(edge,newnode)
            setNode!(newedge,newnode)
            setEdge!(othernode,newedge)
            setNode!(newedge,othernode)
            pushEdge!(net,newedge)
            pushNode!(net,newnode)
            t = edge.length
                setLength!(edge,t/2)
                setLength!(newedge,t/2)
            net.root = length(net.node) #last node is root
       else
            @warn "external edge $(net.node[index].edge[1].number) leading to outgroup $(outgroup) cannot contain root, root placed wherever"
            checkroot = true
       end
    end
    if(checkroot && !isTree(net))
        checkRootPlace!(net)
    end
 end

# function to check if a node could be root
# by the containRoot attribute of edges around it
function canBeRoot(n::Node)
    !n.hybrid || return false
    #!n.hasHybEdge || return false #need to allow for some reason, check ipad notes
    !n.leaf || return false
    return any([e.containRoot for e in n.edge])
end

# function to delete the extra node created in updateRoot
# this extra node is needed to be able to compare networks with the distance function
# but if left in the network, everything crashes (as everything assumes three edges per node)
# fromUpdateRoot=true if called after updateRoot (in which case leaf has to be a leaf), ow it is used in readTopUpdate
function undoRoot!(net::HybridNetwork, fromUpdateRoot::Bool)
    if(length(net.node[net.root].edge) == 2)
        root = net.node[net.root]
        leaf = getOtherNode(root.edge[1],root).leaf ? getOtherNode(root.edge[1],root) : getOtherNode(root.edge[2],root)
        (fromUpdateRoot && leaf.leaf) || error("root should have one neighbor leaf which has to be the outgroup defined")
        deleteIntLeafWhile!(net,root,leaf);
    end
end

undoRoot!(net::HybridNetwork) = undoRoot!(net, true)

# function to read the .out file from snaq (optTopRuns) function
function readOutfile(file::AbstractString)
    try
        s = open(file)
    catch
        error("Could not find or open $(file) file");
    end
    s = open(file)
    line = readline(s)
    @debug "line read $(line)"
    c = line[1]
    if(c == '(')
       println("Estimated network from file $(file): $(line)")
       net = readTopologyUpdate(line)
       vec = split(line,"-Ploglik = ")
       net.loglik = parse(Float64,vec[2])
    else
        error("output file $(filename).out does not contain a tree in the first line, instead it has $(line); or we had trouble finding ploglik.")
    end
    return net
end

"""
`readSnaqNetwork(output file)`

function to read the estimated network from an .out file generated by the snaq function
"""
readSnaqNetwork(file::AbstractString) = readOutfile(file)


# function to change negative branch lengths to 1.0 for starting topology
# and to change big branch lengths to 10.0
# also uses setLength for all edges
function cleanBL!(net::HybridNetwork)
    ##println("missing branch lengths will be set to 1.0")
    for e in net.edge
        if(e.length < 0)
            setLength!(e,1.0)
        elseif(e.length > 10.0)
            setLength!(e,10.0)
        else
            setLength!(e,e.length)
        end
    end
end


# function to read multiple topologies
# - calls readInputTrees in readData.jl, which
#   calls readTopologyUpdate here, for level 1 networks.
# - read a file and create one object per line read
# (each line starting with "(" will be considered a topology)
# the file can have extra lines that are ignored
# returns an array of HybridNetwork objects (that can be trees)
function readMultiTopologyLevel1(file::AbstractString)
    readInputTrees(file)
end

"""
`readMultiTopology(file)`

Read a text file with a list of networks in parenthetical format (one per line).
Each network is read with [`readTopology`](@ref).
Return an array of HybridNetwork object.
"""
function readMultiTopology(file::AbstractString)
    s = open(file)
    numl = 1
    vnet = HybridNetwork[];
    for line in eachline(s)
        line = strip(line) # remove spaces
        c = isempty(line) ? "" : line[1]
        if(c == '(')
           try
               push!(vnet, readTopology(line,false)) # false for non-verbose
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

"""
    writeMultiTopology(nets, file_name; append=false)
    writeMultiTopology(nets, IO)

Write an array of networks in parenthetical extended Newick format, one network per line.
Use the option append=true to append to the file. Otherwise, the default is to create a new
file or overwrite it, if it already existed.
Each network is written with `writeTopology`.

# Examples #"
```
julia> net = [readTopology("(D,((A,(B)#H7:::0.864):2.069,(F,E):3.423):0.265,(C,#H7:::0.1361111):10);"),
              readTopology("(A,(B,C));"),readTopology("(E,F);"),readTopology("(G,H,F);")];

julia> writeMultiTopology(net, "fournets.net") # to (over)write to file "fournets.net"
julia> writeMultiTopology(net, "fournets.net", append=true) # to append to this file
julia> writeMultiTopology(net, stdout)         # to write to the screen (standard out)
(D,((A,(B)#H7:::0.864):2.069,(F,E):3.423):0.265,(C,#H7:::0.1361111):10.0);
(A,(B,C));
(E,F);
(G,H,F);
```
""" #"
function writeMultiTopology(n::Vector{HybridNetwork},file::AbstractString; append::Bool=false)
    mode = (append ? "a" : "w")
    s = open(file, mode)
    writeMultiTopology(n,s)
    close(s)
end

function writeMultiTopology(net::Vector{HybridNetwork},s::IO)
    for i in 1:length(net)
      try
        # writeTopologyLevel1(net[i],s,false,true,"none",false,false,3)
        writeTopology(net[i],s) # no rounding, not for dendroscope
        write(s,"\n")
      catch err
        if isa(err, RootMismatch) # continue writing other networks in list
            @error "\nError with topology $i:\n" * err.msg
        else rethrow(err); end
      end
    end
end


"""
    writeTopology(net)
    writeTopology(net, filename)
    writeTopology(net, IO)

Write the parenthetical extended Newick format of a network,
as a string, to a file or to an IO buffer / stream.
Optional arguments (default values):

- di (false): write in format for Dendroscope
- round (false): rounds branch lengths and heritabilities γ
- digits (3): digits after the decimal place for rounding
- append (false): if true, appends to the file
- internallabel (true): if true, writes internal node labels

If the current root placement is not admissible, other placements are tried.
The network is updated with this new root placement, if successful.

Uses lower-level function [`writeSubTree!`](@ref).
"""
function writeTopology(n::HybridNetwork, file::AbstractString; append::Bool=false,
        round=false::Bool, digits=3::Integer, di=false::Bool, # keyword arguments
        internallabel=true::Bool)
    mode = (append ? "a" : "w")
    s = open(file, mode)
    writeTopology(n,s,round,digits,di,internallabel)
    write(s,"\n")
    close(s)
end

function writeTopology(n::HybridNetwork;
        round=false::Bool, digits=3::Integer, di=false::Bool, # keyword arguments
        internallabel=true::Bool)
    s = IOBuffer()
    writeTopology(n,s,round,digits,di,internallabel)
    return String(take!(s))
end

function writeTopology(net::HybridNetwork, s::IO,
        round=false::Bool, digits=3::Integer, di=false::Bool, # optional arguments
        internallabel=true::Bool)
    # check/find admissible root: otherwise could be trapped in infinite loop
    rootsaved = net.root
    changeroot = false
    msg = ""
    try
        directEdges!(net)
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
            i = getIndex(getParent(e), net)
            net.root = i
            try
                directEdges!(net)
                print("Setting root at node $(net.node[i].number) (net.root = $i)\n\n")
                print(msg)
                changeroot = false
                break # stop loop over edges
            catch err
                if !isa(err, RootMismatch) rethrow(err); end
            end
          #end
        end
        if changeroot # none of hybrid edges worked
            net.root = rootsaved
            throw(RootMismatch("Could not find admissible root. Cannot write topology."))
            changeroot=false # safety exit of while (but useless)
        end
    end
    # finally, write parenthetical format
    writeSubTree!(s,net,di,true,round,digits,internallabel)
    # namelabel = true: to print leaf & node names (labels), not numbers
end


###############################################################################
## Generate symetric tree
###############################################################################

"""
    symmetricTree(n, i=1)

Create a string with a symmetric tree with 2^n tips, numbered from i to i+2^n-1.
All the branch length are set equal to 1.
The tree can be created with function readTopology.
"""
function symmetricTree(n::Int, ell::Real, i=1::Int)
    # Build tree
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
    symmetricNet(n, i, j, gamma)

Create a string with a symmetric net with 2^n tips, numbered from 1 to 2^n
All the branch length are set equal to 1.
One hybrid branch, going from level i to level j is added, with weigth gamma.
The tree can be created with function readTopology.
"""
function symmetricNet(n::Int, i::Int, j::Int, gamma::Real, ell::Real)
    # Checks
    if (n < i || i < j || j < 1) error("Must be n > i > j > 0") end
    # Underlying tree
    tree = symmetricTree(n, ell)
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
    symmetricNet(n, h, gamma)

Create a string with a symmetric net with 2^n tips, numbered from 1 to 2^n
The total height of the network is set to 1.
Hybrids are added from level h to h-1 symmetrically.
"""
function symmetricNet(n::Int, h::Int, gamma::Real, i=1::Int)
    # Checks
    if (n < h || h < 2) error("Must be n >= h > 1.") end
    # length of branch
    ell = 1.0/n
    # Element net
    net = symmetricNet(h, h, h-1, gamma, ell)
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
