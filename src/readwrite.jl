# functions to read/write networks topologies
# originally in functions.jl
# Claudia March 2015

# ------------------------------- Read topology ------------------------------------------

function peekchar(s::IOStream)
    Base.peekchar(s)
end

function peekchar(s::IOBuffer)
    mark(s)
    c = read(s,Char)
    reset(s)
    return c
end

# aux function to detect characters which should be ignored by readTopology
function iswhitesymbol(c::Char)
    whitesymbols = [' ', '\n', '\r', '\t']
    return c in whitesymbols
end

# aux function to peek the next non-white-symbol char in s, does not modify s
# input: s IOStream/IOBuffer
function peekskip(s::IO)
    c = peekchar(s)
    mark(s)
    while iswhitesymbol(c)
        c = read(s, Char)
    end
    reset(s)
	return c
end

# aux function to read the next non-white-symbol char in s, advances s
# input: s IOStream/IOBuffer
function readskip!(s::IO)
    c = read(s, Char)
    while iswhitesymbol(c)
        c = read(s, Char)
    end
    return c
end


# aux function to advance stream in readSubtree
# input: s IOStream/IOBuffer
function advance!(s::IO, c::Char, numLeft::Array{Int,1})
    c = peekskip(s)
    if(Base.eof(s))
        error("Tree ends prematurely while reading subtree after left parenthesis $(numLeft[1]-1).")
    end
    return readskip!(s)
end


# aux function to read all digits of taxon name
# it allows names with letters and numbers
# it also reads # as part of the name and returns pound=true
# it returns the node name as string as well to check if it exists already (as hybrid)
# warning: treats digit taxon numbers as strings to avoid repeated node numbers
function readNum(s::IO, c::Char, net::HybridNetwork, numLeft::Array{Int,1})
    pound = 0;
    if(isalnum(c) || isValidSymbol(c) || c == '#')
        pound += (c == '#') ? 1 : 0
        name = readskip!(s)
        c = peekskip(s)
        while(isalnum(c) || isValidSymbol(c) || c == '#')
            d = readskip!(s)
            name = string(name,d)
            if(d == '#')
                pound += 1;
               c = peekskip(s);
               if(isalnum(c))
                   if(c != 'L' && c != 'H' && c != 'R')
                       warn("Expected H, R or LGT after # but received $(c) in left parenthesis $(numLeft[1]-1).")
                   end
               else
                   a = readstring(s);
                   error("Expected name after # but received $(c) in left parenthesis $(numLeft[1]-1). Remaining is $(a).")
               end
            end
            c = peekskip(s);
        end
        if(pound == 0)
            return size(net.names,1)+1, name, false
        elseif(pound == 1)
            return size(net.names,1)+1, name, true
        else
            a = readstring(s);
            error("strange node name with $(pound) # signs. remaining is $(a).")
        end
    else
        a = readstring(s);
        error("Expected int digit, alphanum or # but received $(c). remaining is $(a).");
    end
end


# aux function to read floats like length
function readFloat(s::IO, c::Char)
    if(isdigit(c) || in(c, ['.','e','-']))
        num = string(readskip!(s));
        c = peekskip(s);
        while(isdigit(c) || in(c, ['.','e','-']))
            d = readskip!(s);
            num = string(num,d);
            c = peekskip(s);
        end
        f = 0.0
        try
            f = float(num)
        catch
            error("problem with number read $(num), not a float number")
        end
        return f
    else
        a = readstring(s);
        error("Expected float digit after : but found $(c). remaining is $(a).");
    end
end

# aux function to identify if a symbol in taxon name is valid
# symbol cannot be: () [] : ; ' , . space \t \r \n
# according to richnewick.pdf
function isValidSymbol(c::Char)
    return !isspace(c) && c != '(' && c != ')' && c != '[' && c != ']' && c != ':' && c != ';' && c != ',' #&& c != '.'
end

"""
    parseRemainingSubtree!(s::IO, numLeft, net, hybrids, index)

Helper function for readSubtree!
Called once a `(` has been read in a tree topology and reads until the corresponding `)` has been found.
This function performs the recursive step for readSubtree!
Advances `s` past the subtree, adds discovered nodes and edges to `net`, `hybrids`, and `index`.
"""
@inline function parseRemainingSubtree!(s::IO, numLeft::Array{Int,1}, net::HybridNetwork, hybrids::Array{String, 1}, index::Array{Int, 1})
    numLeft[1] += 1
    DEBUGC && println(numLeft)
    n = Node(-1*numLeft[1],false);
    DEBUG && println("creating node $(n.number)")
    keepon = true;
    c = readskip!(s)
    while (keepon)
        bl = readSubtree!(s,n,numLeft,net,hybrids,index)
        c = advance!(s,c,numLeft)
        if (c == ')')
            keepon = false
        elseif (c != ',')
            a = readstring(s);
            error("Expected right parenthesis after left parenthesis $(numLeft[1]-1) but read $(c). The remainder of line is $(a).")
        end
    end
    return n
end

"""
    parseHybridNode!(node, parentNode, hybridName, net, hybrids, index)

Helper function for readSubtree!
Handles any type of given hybrid node.
To be called once a `#` has been found in a tree topology.
Creates a hybrid node and its edges, then inserts those into `net`, `hybrids` and `index` accordingly.
"""
@inline function parseHybridNode!(n::Node, parent::Node, name::String, net::HybridNetwork, hybrids::Array{String, 1}, index::Array{Int, 1})
    DEBUG && println("found pound in $(name)")
    n.hybrid = true;
    DEBUGC && println("encontro un hybrid $(name).")
    DEBUGC && println("hybrids list tiene size $(size(hybrids,1))")
    if(in(name, hybrids))
        DEBUG && println("dice que $(name) esta en hybrids")
        ind = getIndex(name,hybrids);
        other = net.node[index[ind]];
        DEBUG && println("other is $(other.number)")
        DEBUGC && println("other is leaf? $(other.leaf), n is leaf? $(n.leaf)")
        if(!n.leaf && !other.leaf)
            error("both hybrid nodes are internal nodes: successors of the hybrid node must only be included in the node list of a single occurrence of the hybrid node.")
        elseif(n.leaf)
            DEBUG && println("n is leaf")
            e = Edge(net.numEdges+1);
            DEBUG && println("creating hybrid edge $(e.number) attached to other $(other.number) and parent $(parent.number)")
            e.hybrid = true
            e.isMajor = false;
            pushEdge!(net,e);
            setNode!(e,[other,parent]);
            setEdge!(other,e);
            setEdge!(parent,e);
            DEBUG && println("e $(e.number )istIdentifiable? $(e.istIdentifiable)")
        else # !n.leaf
            DEBUG && println("n is not leaf, other is leaf")
            if(size(other.edge,1) == 1) #other should be a leaf
                DEBUGC && println("other is $(other.number), n is $(n.number), edge of other is $(other.edge[1].number)")
                otheredge = other.edge[1];
                otherparent = getOtherNode(otheredge,other);
                DEBUG && println("otheredge is $(otheredge.number)")
                DEBUG && println("parent of other is $(otherparent.number)")
                removeNode!(other,otheredge);
                deleteNode!(net,other);
                setNode!(otheredge,n);
                setEdge!(n,otheredge);
##                    otheredge.istIdentifiable = true ## setNode should catch this, but when fixed, causes a lot of problems
                DEBUG && println("setting otheredge to n $(n.number)")
                e = Edge(net.numEdges+1);
                DEBUG && println("creating hybrid edge $(e.number) between n $(n.number) and parent $(parent.number)")
                e.hybrid = true
                setNode!(e,[n,parent]);
                setEdge!(n,e);
                setEdge!(parent,e);
                pushNode!(net,n);
                pushEdge!(net,e);
                n.number = other.number;
                n.name = other.name;
                DEBUG && println("edge $(e.number) istIdentifiable? $(e.istIdentifiable)")
                DEBUG && println("otheredge $(otheredge.number) istIdentifiable? $(otheredge.istIdentifiable)")
            else
                error("strange: node $(other.number) is a leaf hybrid node so it should have only one edge and it has $(size(other.edge,1))")
            end
        end
    else
        DEBUG && println("dice que $(name) no esta en hybrids")
        DEBUG && println("$(name) es leaf? $(n.leaf)")
        n.hybrid = true;
        push!(net.names,string(name));
        n.name = string(name);
        DEBUGC && println("aqui vamos a meter a $(name) en hybrids")
        push!(hybrids,string(name));
        pushNode!(net,n);
        push!(index,size(net.node,1));
        e = Edge(net.numEdges+1);
        DEBUG && println("creating hybrid edge $(e.number)")
        e.hybrid = true
        n.leaf ? e.isMajor = false : e.isMajor = true
        pushEdge!(net,e);
        setNode!(e,[n,parent]);
        setEdge!(n,e);
        setEdge!(parent,e);
        DEBUG && println("thus edge $(e.number) istIdentifiable? $(e.istIdentifiable)")
    end
    e.containRoot = !e.hybrid
    return e
end

"""
    parseNonhybridNode!(node, parentNode, net)

Helper function for readSubtree!
Inserts a nonhybrid node and associated edge into `net`.
"""
@inline function parseNonhybridNode!(n::Node, parent::Node, net::HybridNetwork)
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

Helper function for parseEdgeData!
Reads a single floating point edge data value in a tree topology.
Returns -1 if no value exists before the next colon, returns a float if a value is present.
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
        warn(errors[call])
        return -1.0
    end
end

#TODO: Check the windows version

"""
    parseEdgeData!(s::IO, edge, node, numberOfLeftParentheses::Array{Int,1})

Helper function for readSubtree!, fixes a bug from using setGamma
Modifies `e` according to the specified edge length and gamma values in the tree topology.
Advances the stream `s` past any existing edge data.
Edges in a topology may optionally be followed by ":edgeLen:bootstrap:gamma"
where edgeLen, bootstrap, and gamma are decimal values.
"""
@inline function parseEdgeData!(s::IO, e::PhyloNetworks.Edge, n::PhyloNetworks.Node, numLeft::Array{Int,1})
    bootstrap = nothing;
    gamma = -1.0;
    read(s, Char);
    edgeLen = getDataValue!(s, 1, numLeft)
    if peekskip(s) == ':'
        readskip!(s)
        bootstrap = getDataValue!(s, 2, numLeft)
    end
    if peekskip(s) == ':'
        readskip!(s)
        gamma = getDataValue!(s, 3, numLeft)
    end

    edgeLenPresent = (edgeLen != -1.0);
    gammaPresent   = (gamma != -1.0 && gamma != nothing);

    e.length = edgeLen
    if gammaPresent
        if(!e.hybrid)
            warn("gamma read for current edge $(e.number) but it is not hybrid, so gamma=$(gamma) ignored")
        else
            parentEdges = Edge[] # The edges which have n as a child
            for e in n.edge
                if e.hybrid
                    nIsChild = (e.isChild1 && e.node[1] == n || !e.isChild1 && e.node[2] == n)
                    if nIsChild
                        push!(parentEdges, e)
                    end
                end
            end
            numParents = size(parentEdges)[1]
            if (numParents == 0) || (numParents == 1)
                e.gamma = gamma
                # isMajor will be set when the other edge is found
            elseif (numParents == 2) #other edge has been read, may have gamma set
                local otheredge::PhyloNetworks.Edge
                parentEdges[1] == e ? otheredge = parentaledges[2] :
                                        otheredge = parentEdges[1]
                # Note: only need to correct isMajor if the gamma values are not equal
                if (!approxEq(gamma + otheredge.gamma, 1)) # gammas do not sum to 1
                    # some gamma < 0
                    if (gamma <= 0 && otheredge.gamma > 0)
                        gamma = 1 - otheredge.gamma
                    elseif (otheredge.gamma <= 0 && gamma > 0)
                        otheredge.gamma = 1 - gamma
                    end
                    # rescale so that they sum to 1
                    e.gamma, otheredge.gamma = ((gamma)/(gamma + otheredge.gamma),
                                                (otheredge.gamma)/(gamma + otheredge.gamma))
                else # both gamma values sum to 1
                    e.gamma = gamma
                end
                if (!approxEq(e.gamma, 0.5))
                        if e.gamma > 0.5
                            e.isMajor = true
                            otheredge.isMajor = false
                        else
                            e.isMajor = false
                            otheredge.isMajor = true
                        end
                    end
            else
                error("hybrid edge gamma value parsed but hybrid node has $(size(parentEdges)[1]) parent edges! (should be between 0 and 2)")
            end
        end
    end
end

"""
    readSubtree!(s::IO, parentNode, numLeft, net, hybrids, index)

A recursive helper method for readTopology
Reads a subtree of a Extended Newick tree topology
input: s IOStream/IOBuffer
warning: reads additional info :length:bootstrap:gamma
warning: allows for name of internal nodes without # after: (1,2)A,...
warning: warning if hybrid edge without gamma value, warning if gamma value (ignored) without hybrid edge
modified from original Cecile c++ code to allow polytomies
"""

function readSubtree!(s::IO, parent::Node, numLeft::Array{Int,1}, net::HybridNetwork, hybrids::Array{String,1}, index::Array{Int,1})
    c = peekskip(s)
    e = nothing;
    hasname = false; # to know if the current node has name
    pound = false;
    if(c == '(')
        # read the rest of the subtree (perform the recursive step!)
        n = parseRemainingSubtree!(s, numLeft, net, hybrids, index)
        c = peekskip(s);
        if(isalnum(c) || isValidSymbol(c) || c == '#') # internal node has name
            hasname = true;
            num, name, pound = readNum(s, c, net, numLeft);
            n.number = num;
            c = peekskip(s);
            #if(!pound)
            #    warn("internal node with name without it being a hybrid node. node name might be meaningless after tree modifications.")
            #end
        end
    elseif(isalnum(c) || isValidSymbol(c) || c == '#')
        hasname = true;
        bl = true;
        num, name, pound = readNum(s, c, net, numLeft)
        n = Node(num, true);
        DEBUG && println("creating node $(n.number)")
    else
        a = readstring(s);
        error("Expected beginning of subtree but read $(c) after left parenthesis $(numLeft[1]-1), remaining is $(a).");
    end
    if(pound) # found pound sign in name
        e = parseHybridNode!(n, parent, name, net, hybrids, index)
    else
        if(hasname)
            push!(net.names,string(name));
            n.name = string(name)
        end
        e = parseNonhybridNode!(n, parent, net)
    end
    c = peekskip(s);
    e.length = -1.0
    e.gamma = e.hybrid? -1.0 : 1.0
    if(c == ':')
        parseEdgeData!(s, e, n, numLeft)
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

Input: text file or parenthetical format directly.
The file name may not start with a left parenthesis, otherwise the file
name itself would be interpreted as the parenthetical description.
"""
readTopology(input::AbstractString) = readTopology(input,true)

function readTopology(s::IO,verbose::Bool)
    net = HybridNetwork()
    line = readuntil(s,";");
    if(line[end] != ';')
        error("file does not end in ;")
    end
    seekstart(s)
    c = peekskip(s)
    numLeft = [1]; # made Array to make it mutable; start at 1 to avoid node -1 which breaks undirectedOtherNetworks
    hybrids = String[];
    index = Int[];
    if(c == '(')
        numLeft[1] += 1;
        #println(numLeft)
        n = Node(-1*numLeft[1],false);
        c = readskip!(s)
        b = false;
        while(c != ';')
            b |= readSubtree!(s,n,numLeft,net,hybrids,index)
            c = readskip!(s);
            if(eof(s))
                error("Tree ended while reading in subtree beginning with left parenthesis number $(numLeft[1]-1).")
            elseif(c == ',')
                continue;
            elseif(c == ')')
                c = peekskip(s);
                if(c == ':')
                    while(c != ';')
                        c = readskip!(s)
                    end
                end
            end
        end
        DEBUG && println("after readsubtree:")
        DEBUG && printEdges(net)
        if(size(n.edge,1) == 1) # root has only one child
            edge = n.edge[1]; # assume it has only one edge
            child = getOtherNode(edge,n);
            removeEdge!(child,edge);
            net.root = getIndex(child,net);
            deleteEdge!(net,edge);
        else
            pushNode!(net,n);
            net.root = getIndex(n,net);
        end
    else
		a = readstring(s)
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

# aux function to send an error if the number of hybrid attached to every
# hybrid node are 0,1 or >2
# to be used in readTopology
# need to be run after storeHybrids
function checkNumHybEdges!(net::HybridNetwork)
    if(!isTree(net))
        isempty(net.hybrid) && error("net.hybrid should not be empty for this network")
        for n in net.hybrid
            hyb = sum([e.hybrid?1:0 for e in n.edge]);
            if(hyb > 2)
                error("hybrid node $(n.number) has more than two hybrid edges attached to it: polytomy that cannot be resolved without intersecting cycles.")
            elseif(hyb == 1)
                if(net.numHybrids == 1)
                    error("only one hybrid node number $(n.number) with name $(net.names[n.number]) found with one hybrid edge attached")
                else
                    error("current hybrid node $(n.number) with name S(net.names[n.number]) has only one hybrid edge attached. there are other $(net.numHybrids-1) hybrids out there but this one remained unmatched")
                end
            elseif(hyb == 0)
                if(length(n.edge) == 0)
                    error("strange hybrid node $(n.number) attached to 0 edges")
                elseif(length(n.edge) == 1)
                    n.leaf || error("hybrid node $(n.number) has only one tree edge attached and it is not a leaf")
                elseif(length(n.edge) >= 2)
                    warn("hybrid node $(n.number) is not connected to any hybrid edges, it was transformed to tree node")
                    n.hybrid = false;
                end
            end
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
    if(n.hybrid)
        suma = sum([!e.hybrid?1:0 for e in n.edge]);
        #println("create edge $(net.numEdges+1)")
        ed1 = Edge(net.numEdges+1,0.0);
        n1 = Node(size(net.names,1)+1,false,false,[ed1]);
        #println("create node $(n1.number)")
        hyb = Edge[];
        for i in 1:size(n.edge,1)
            !n.edge[i].hybrid ? push!(hyb,n.edge[i]) : nothing
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
        if(size(n1.edge,1) > 3)
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
    mod(sum([!e.hybrid?e.gamma:0 for e in net.edge]),1) == 0 ? nothing : error("tree (not network) read and some tree edge has gamma different than 1")
    nodes = copy(net.node)
    for n in nodes
        if(isNodeNumIn(n,net.node)) # very important to check
            if(size(n.edge,1) == 2)
                if(!n.hybrid)
                    if(!leaveRoot || !isEqual(net.node[net.root],n)) #if n is the root
                        deleteIntNode!(net,n);
                    end
                else
                    hyb = sum([e.hybrid?1:0 for e in n.edge]);
                    if(hyb == 1)
                        deleteIntNode!(net,n);
                    end
                end
            end
            if(!n.hybrid)
                if(size(n.edge,1) > 3)
                    DEBUG && warn("polytomy found in node $(n.number), random resolution chosen")
                    solvePolytomy!(net,n);
                end
                hyb = sum([e.hybrid?1:0 for e in n.edge]);
                if(hyb == 1)
                    n.hasHybEdge == true;
                elseif(hyb > 1)
                    warn("strange tree node $(n.number) with more than one hybrid edge, intersecting cycles maybe")
                end
            else
                hyb = sum([e.hybrid?1:0 for e in n.edge]);
                tre = sum([!e.hybrid?1:0 for e in n.edge]);
                if(hyb > 2)
                    error("hybrid node $(n.number) has more than two hybrid edges attached to it: polytomy that cannot be resolved without intersecting cycles.")
                elseif(hyb == 1)
                    hybnodes = sum([n.hybrid?1:0 for n in net.node]);
                    if(hybnodes == 1)
                        error("only one hybrid node number $(n.number) with name $(net.names[n.number]) found with one hybrid edge attached")
                    else
                        error("current hybrid node $(n.number) with name $(net.names[n.number]) has only one hybrid edge attached. there are other $(hybnodes-1) hybrids out there but this one remained unmatched")
                    end
                elseif(hyb == 0)
                    warn("hybrid node $(n.number) is not connected to any hybrid edges, it was transformed to tree node")
                    n.hybrid = false;
                else # 2 hybrid edges
                    if(tre == 0) #hybrid leaf
                        warn("hybrid node $(n.number) is a leaf, so we add an extra child")
                        addChild!(net,n);
                    elseif(tre > 1)
                        warn("hybrid node $(n.number) has more than one child so we need to expand with another node")
                        expandChild!(net,n);
                    end
                    suma = sum([e.hybrid?e.gamma:0 for e in n.edge]);
                    if(suma == -2)
                        #warn("hybrid edges in read network without gammas")
                        println("hybrid edges for hybrid node $(n.number) have missing gamma's, set default: 0.9,0.1")
                        for e in n.edge
                            if(e.hybrid)
                                if (e.isMajor)
                                    e.gamma = 0.9;
                                    e.isMajor = true
                                else
                                    e.gamma = 0.1;
                                    e.isMajor = false
                                end
                            end
                        end
                    elseif(suma != 1)
                        ed1 = nothing
                        ed2 = nothing
                        for e in n.edge
                            if(e.hybrid)
                                isa(ed1,Void) ? ed1=e : ed2=e
                            end
                        end
                        if(ed1.gamma > 0 && ed2.gamma > 0 && ed1.gamma < 1 && ed2.gamma < 1) #both gammas were set, but contradictory
                            error("hybrid edges for hybrid node $(n.number) have gammas that do not sum up to one: $(ed1.gamma),$(ed2.gamma)")
                        elseif(ed1.gamma != -1.0)
                            warn("only one hybrid edge of hybrid node $(n.number) has gamma value $(ed1.gamma) set, the other edge will be assigned $(1-ed1.gamma).")
                            setGamma!(ed2,1-ed1.gamma, false);
                        else
                            warn("only one hybrid edge of hybrid node $(n.number) has gamma value $(ed2.gamma) set, the other edge will be assigned $(1-ed2.gamma).")
                            setGamma!(ed1,1-ed2.gamma, false);
                        end
                    elseif(suma == 1)
                        for e in n.edge
                            if(e.hybrid)
                                if(approxEq(e.gamma,0.5))
                                    e.isMajor = false
                                    break
                                else
                                    break
                                end
                            end
                        end
                    end
                end
            end
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
        #warn("topology read is a tree as it has no hybrid nodes")
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
        #warn("not a network read, but a tree as it does not have hybrid nodes")
        all((e->e.containRoot), net.edge) ? nothing : error("some tree edge has contain root as false")
        all((e->!e.hybrid), net.edge) ? nothing : error("some edge is hybrid and should be all tree edges in a tree")
        all((n->!n.hasHybEdge), net.node) ? nothing : error("some tree node has hybrid edge true, but it is a tree, there are no hybrid edges")
    else
        if(!net.cleaned)
            for n in net.hybrid
                success,hyb,flag,nocycle,flag2,flag3 = updateAllNewHybrid!(n,net,false,true,false)
                if(!success)
                    warn("current hybrid $(n.number) conflicts with previous hybrid by intersecting cycles: $(!flag), nonidentifiable topology: $(!flag2), empty space for contain root: $(!flag3), or does not create a cycle (probably problem with the root placement): $(nocycle).")
                    #net.cleaned = false
                end
            end
            DEBUG && println("before update partition")
            DEBUG && printPartitions(net)
            for n in net.hybrid #need to updatePartition after all inCycle
                nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,n);
                updatePartition!(net,nodesInCycle)
                DEBUG && println("after updating partition for hybrid node $(n.number)")
                DEBUG && printPartitions(net)
            end
        end
    end
end

# cleanAfterReadAll includes all the step to clean a network after read
function cleanAfterReadAll!(net::HybridNetwork, leaveRoot::Bool)
    DEBUG && println("cleanBL -----")
    cleanBL!(net)
    DEBUG && println("cleanAfterRead -----")
    cleanAfterRead!(net,leaveRoot)
    DEBUG && println("updateAllReadTopology -----")
    updateAllReadTopology!(net) #fixit: it could break if leaveRoot = true (have not checked it), but we need to updateContainRoot
    if(!leaveRoot)
        DEBUG && println("parameters -----")
        parameters!(net)
    end
    DEBUG && println("check root placement -----")
    checkRootPlace!(net)
    net.node[net.root].leaf && warn("root node $(net.node[net.root].number) is a leaf, so when plotting net, it can look weird")
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
    DEBUG && println("readTopology -----")
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
        tmp = findin([n.name for n in net.leaf], [outgroup])
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




    # --------------------------- write topology -------------------------------------
# function to write a node and its descendants in parenthetical format
# di=true in densdroscope format, names=true, prints names
writeSubTree!(s::IO, n::Node, parent::Edge, di::Bool, names::Bool) =
    writeSubTree!(s,n,parent,di,names, true,3)

# method to start at the root and write the whole tree/network
function writeSubTree!(s::IO, net::HybridNetwork, di::Bool, names::Bool,
                       roundBL::Bool, digits::Integer)
    if net.numNodes == 1
        print(s, (names? net.node[net.root].name : string(net.node[net.root].number)))
    elseif net.numNodes > 1
        print(s,"(")
        degree = length(net.node[net.root].edge)
        for e in net.node[net.root].edge
            writeSubTree!(s,getOtherNode(e,net.node[net.root]),e,di,names,roundBL,digits)
            degree -= 1
            degree == 0 || print(s,",")
        end
        print(s,")")
    end
    print(s,";")
    return nothing
end

# method to start at a node coming from an adjacent (parent) edge
function writeSubTree!(s::IO, n::Node, parent::Edge,di::Bool,names::Bool,
                       roundBL::Bool, digits::Integer)
    # subtree below node n:
    if (parent.isMajor && !n.leaf) # do not descent below a minor hybrid edge
        print(s,"(")
        firstchild = true
        for e in n.edge
            e == parent && continue # skip parent edge where we come from
            (e.hybrid && e.node[(e.isChild1 ? 1 : 2)]==n) && continue # no going up minor hybrid
            firstchild || print(s, ",")
            firstchild = false
            child = getOtherNode(e,n)
            writeSubTree!(s,child,e, di,names, roundBL, digits)
        end
        print(s,")")
    end
    # node label:
    if (parent.hybrid)
        print(s, (names ? n.name : string("#H",n.number)))
        n.name != "" || parent.isMajor || warn("hybrid node $(n.number) has no name")
    elseif (n.leaf)
        print(s, (names ? n.name : n.number))
    end
    # branch lengths and γ, if available:
    printBL = false
    if(parent.length != -1.0) # -1.0 means missing
        print(s,string(":",(roundBL ? round(parent.length,digits) : parent.length)))
        printBL = true
    end
    if(parent.hybrid && !di) # && (!printID || !n.isBadDiamondI))
        if(parent.gamma != -1.0)
            if(!printBL) print(s,":"); end
            print(s,string("::",(roundBL ? round(parent.gamma,digits) : parent.gamma)))
        end
    end
end

# function to writeTopology for level 1 networks. bad net.root okay
#                           deepcopies net -> does *not* modify it
# if string=true, returns a string with network in parenthetical format
#                 ow returns the IOBuffer object
# need as input HybridNetwork, since QuartetNetwork does not have root
# input di=true if written for Dendroscope (without gammas)
# names=true, writes the names instead of node numbers, default true
# outgroup: place the root in the external edge of this taxon if possible,
# if none given, placed the root wherever possible
# printID=true, only print identifiable BL as determined by setNonIdBL!.
#         false by default. true inside snaq.
# multall = true, multiple alleles case
function writeTopologyLevel1(net0::HybridNetwork, di::Bool, str::Bool, names::Bool,outgroup::AbstractString, printID::Bool, roundBL::Bool, digits::Integer, multall::Bool)
    s = IOBuffer()
    writeTopologyLevel1(net0,s,di,names,outgroup,printID,roundBL,digits, multall)
    if(str)
        return String(s)
    else
        return s
    end
end

# warning: I do not want writeTopologyLevel1 to modify the network if outgroup is given! thus, we have updateRoot, and undoRoot
function writeTopologyLevel1(net0::HybridNetwork, s::IO, di::Bool, names::Bool,
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
            DEBUG && println("net not cleaned inside writeTopologyLevel1, need to run updateContainRoot")
            for n in net.hybrid
                flag,edges = updateContainRoot!(net,n)
                flag || error("hybrid node $(n.hybrid) has conflicting containRoot")
            end
        end
        updateRoot!(net,outgroup)
        #DEBUG && printEverything(net)
        CHECKNET && canBeRoot(net.node[net.root])
        if(multall)
            mergeLeaves!(net)
        end
        writeSubTree!(s, net ,di,names, roundBL,digits)
    end
    # outgroup != "none" && undoRoot!(net) # not needed because net is deepcopy of net0
    # to delete 2-degree node, for snaq.
end

writeTopologyLevel1(net::HybridNetwork,di::Bool,str::Bool,names::Bool,outgroup::AbstractString,printID::Bool) = writeTopologyLevel1(net,di,str,names,outgroup,printID, false,3, false)
# above: default roundBL=false (at unused digits=3 decimal places)
writeTopologyLevel1(net::HybridNetwork,printID::Bool) = writeTopologyLevel1(net,false, true,true,"none",printID, false, 3, false)
writeTopologyLevel1(net::HybridNetwork,outgroup::AbstractString) = writeTopologyLevel1(net,false, true,true,outgroup,true, false, 3, false)
writeTopologyLevel1(net::HybridNetwork,di::Bool,outgroup::AbstractString) = writeTopologyLevel1(net,di, true,true,outgroup,true, false, 3, false)

"""
`writeTopologyLevel1(net::HybridNetwork)`

writes the parenthetical format of a HybridNetwork object with many optional arguments:

- di=true: write in format for Dendroscope (default false)
- names=false: write the leaf nodes numbers instead of taxon names (default true)
- outgroup (string): name of outgroup to root the tree/network
- printID=true, only print branch lengths for identifiable egdes according to the snaq estimation procedure (default false)
- round: rounds branch lengths and heritabilities γ (default: true)
- digits: digits after the decimal place for rounding (defult: 3)

The topology may be written using a root different than net.root,
if net.root is incompatible with one of more hybrid node.
Missing hybrid names are written as "#Hi" where "i" is the hybrid node number if possible.
The network object is *not* modified.
""" #"
writeTopologyLevel1(net::HybridNetwork; di=false::Bool, string=true::Bool, names=true::Bool,outgroup="none"::AbstractString, printID=false::Bool, round=false::Bool, digits=3::Integer, multall=false::Bool) = writeTopologyLevel1(net, di, string, names,outgroup,printID, round,digits, multall)

# function to check if root is well-placed
# and look for a better place if not
# searches on net.node because net.root is the index in net.node
# if we search in net.edge, we then need to search in net.node
function updateRoot!(net::HybridNetwork, outgroup::AbstractString)
    checkroot = false
    if(outgroup == "none")
        DEBUG && println("no outgroup defined")
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
            DEBUGC && println("creating new node in the middle of the external edge $(edge.number) leading to outgroup $(node.number)")
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
            warn("external edge $(net.node[index].edge[1].number) leading to outgroup $(outgroup) cannot contain root, root placed wherever")
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
    DEBUG && println("line read $(line)")
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
Each network is read with `readTopology`.
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
           catch(err)
               println("could not read tree on line $(numl) of file $file. error was this:")
               rethrow(err)
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
julia> writeMultiTopology(net, STDOUT)         # to write to the screen (standard out)
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
            warn("\nError with topology $i:\n" * err.msg)
        else rethrow(err); end
      end
    end
end


"""
    writeTopology(net)
    writeTopology(net, filename)

write the parenthetical extended Newick format of a HybridNetwork object, as a string or to a file.
Optional arguments (default values):

- di (false): write in format for Dendroscope
- round (false): rounds branch lengths and heritabilities γ
- digits (3): digits after the decimal place for rounding

If the current root placement is not admissible, other placements are tried.
The network is updated with this new root placement, if successful.
"""
function writeTopology(n::HybridNetwork, file::AbstractString; append::Bool=false,
        round=false::Bool, digits=3::Integer, di=false::Bool) # keyword arguments
    mode = (append ? "a" : "w")
    s = open(file, mode)
    writeTopology(n,s,round,digits,di)
    write(s,"\n")
    close(s)
end

function writeTopology(n::HybridNetwork;
        round=false::Bool, digits=3::Integer, di=false::Bool) # keyword arguments
    s = IOBuffer()
    writeTopology(n,s,round,digits,di)
    return String(s)
end

function writeTopology(net::HybridNetwork, s::IO,
        round=false::Bool, digits=3::Integer, di=false::Bool) # optional arguments
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
            i = getIndex(e.node[e.isChild1? 2 : 1], net)
            net.root = i
            try
                directEdges!(net)
                print("Setting root at node $(net.node[i].number) (net.root = $i)\n\n")
                print(msg)
                changeroot = false
                break # stop loop over edges
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
    writeSubTree!(s,net,di,true,round,digits)
    # names = true: to print leaf names (labels), not numbers
    ## printID = false: print all branch lengths, not just identifiable ones
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
        tree = replace(tree, "A$(i-1+k+1):$(ell)", "A$(i-1+k):$(ell)", k)
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
    tree = replace(tree, op, "(#H:$(ell)::$(gamma),"*op)
    tree = replace(tree, clo, clobis*"):$(0.5*ell)")
    ## end hyb
    op = "A$(2^(i-1)+1):$(ell)"
    clo = "A$(2^(i-1) + 2^(j-1)):$(ell)"
    for k in 2:j
        op = "("*op
        clo = clo*"):$(ell)"
    end
    clobis = clo[1:(length(clo)-length("$(ell)"))]*"$(0.5*ell)"
    tree = replace(tree, op, "("*op)
    tree = replace(tree, clo, clobis*")#H:$(0.5*ell)::$(1-gamma)")
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
    net = replace(net, "H", "H$(2^(n-h))")
    for k in (2^(n-h)-1):-1:1
        net = replace(net, "H$(k+1)", "H$(k)", 2*k)
    end
    # Rename tips
    for k in (2^h):-1:1
        net = replace(net, "A$(k):", "A$(i-1+2^n):")
    end
    for k in (2^n-1):-1:1
        net = replace(net, "A$(i-1+k+1):", "A$(i-1+k):", k)
    end
    return(net)
end
