# Stage2 pseudolikelihood implementation: sketch of classes (types) in Julia
# Claudia (June 2014)
# Classes based on "ane/public/quartetNetwork/notes" and "mbsumtree.h"
# Methods (functions) for each class are coded outside the class and they can
# be used by any class
#
# Working code: after tests, code goes to types.jl or functions.jl
#
# Ctrl S: fixit, todo, check
######################################################################

# types in "types.jl"
include("types.jl")

# functions in "functions.jl"
include("functions.jl")

# needed modules:
using DataStructures # for updateInCycle with queue
using Base.Collections # for updateInCycle with priority queue

# examples
include("case_f_example.jl");
include("bad_triangle_example.jl");

include("tree_example.jl");

# -------------- NETWORK ----------------------- #

# aux function to advance stream in readSubtree
# warning: s IOStream needs to be a file, not a stream converted
#          from string
function advance!(s::IOStream, c::Char, numLeft::Int64)
    c = Base.peekchar(s);
    if(Base.eof(s))
        error("Tree ends prematurely while reading subtree after left parenthesis $(numLeft).")
    end
    c = read(s,Char);
end


# aux function to read all digits of taxon name
function readNum(s::IOStream, c::Char)
    if(isdigit(c))
        num = read(s,Char);
        while(isdigit(c))
            d = read(s,Char);
            num = string(num,d);
            c = Base.peekchar(s);
        end
        return int(num)
    else
        error("Expected digit but received $(c)");
    end
end

# aux function to read floats like length
function readFloat(s::IOStream, c::Char)
    if(isdigit(c))
        num = read(s,Char);
        while(isdigit(c) || c == '.')
            d = read(s,Char);
            num = string(num,d);
            c = Base.peekchar(s);
        end
        return float(num)
    else
        error("Expected digit after : but received $(c)");
    end
end


# aux function to read subtree
# warning: s IOStream needs to be a file, not a stream converted
#          from string
# warning: assumes taxon numbers, not taxon names
# fixit: allow taxon names also
function readSubtree(s::IOStream, parent::Node, numLeft::Int64,net::HybridNetwork)
    c = Base.peekchar(s);
    e = nothing;
    if(c =='(')
       n = Node(numLeft,false);
       c = read(s,Char);
       numLeft += 1;
       bl = readSubtree(s,n,numLeft,net);
       advance!(s,c,numLeft);
       br = false;
       if(c == ',')
           br = readSubtree(s,n,numLeft,net);
       else
           a = readall(s);
           error("Expected comma after left parenthesis $(numLeft) but read $(c). The remainder of line is $(a).")
       end
       advance!(s,c,numLeft);
       if(c != ')')
           a = readall(s);
           error("Expected right parenthesis after left parenthesis $(numLeft) but read $(c). The remainder of line is $(a).")
       end
       if(bl && br)
           pushNode!(net,n);
           e = Edge(net.numEdges+1,1.0);
           pushEdge!(net,e);
           setNode!(e,[n,parent]);
           setEdge!(n,e);
           setEdge!(parent,e);
           n.leaf = false;
       else
           if(size(n.edge,1) == 1) # root only has one child
               edge = n.edge[1]; # assume it has only one edge
               child = getOtherNode(edge,n);
               setNode!(edge,[child,parent]);
               setEdge!(parent,edge);
           end
           return true
       end
    elseif(isdigit(c))
        num = readNum(s,c);
        println("creating node $(num)")
        n = Node(num,true);
        pushNode!(net,n);
        e = Edge(net.numEdges+1,1.0);
        pushEdges!(net,e);
        setNode!(e,[n,parent]);
        setEdge!(n,e);
        setEdge!(parent,e);
    else
        error("Expected beginning of subtree but read $(c)");
    end
    c = Base.peekchar(s);
    if(c == ':')
        c = read(s,Char);
        c = Base.peekchar(s);
        length = readFloat(s,c);
    end
    if(isa(e,Nothing))
        return false
    else
        setLength!(e,length);
    end
end

# function to read topology from parenthetical format
# input: file name
# warning: at the moment, assumes a tree
# warning: assumes taxon numbers, not taxon names
function readTopology(file::String)
    net = HybridNetwork();
    try
        s = open(file);
    catch
        error("Could not find or open $(file) file");
    end
    s = open(file);
    c = Base.peekchar(s);
    numLeft = 0;
    if(c == '(')
       numLeft += 1;
       n = Node(numLeft,false);
       c = read(s,Char);
       b = false;
       while(c != ';')
           b |= readSubtree(s,n,numLeft,net);
           c = read(s,Char);
           if(eof(s))
               error("Tree ended while reading in subtree beginning with left parenthesis number $(numLeft).")
           elseif(c == ',')
               continue;
           elseif(c == ')')
               c = Base.peekchar(s);
           end
       end
       if(size(n.edge,1) == 1) # root has only one child
           edge = n.edge[1]; # assume it has only one edge
           child = getOtherNode(edge,n);
           removeEdge!(child,edge);
           net.root = getIndex(child,net);
       else
           pushNode!(net,n);
           net.root = getIndex(n,net);
       end
    else
       error("Expected beginning of tree with ( but received $(c) instead")
    end
    return net
end


# cecile: check updategammaz function, maybe we need two functions,
# one to update when changing length one to update when changing
# gamma? what i like about updategammaz is that you use that directly
# at the beginning of network, so maybe we should consider doing
# things ourselves inside setLength and setGamma, instead of calling
# update gamma



# function to traverse the network
# simply prints the traversal path, can be modified to do other things
# needs:
visited  =  [false for i  =  1:size(net.node,1)];

function traverse(net::HybridNetwork, node::Node, visited::Array{Bool,1})
    println("estamos en $(node.number)");
    visited[getIndex(node,net)]  =  true;
    if(node.leaf)
        println("llegamos a leaf $(node.number)");
    else
        for(i in 1:size(node.edge,1))
            other  =  getOtherNode(node.edge[i],node);
            if(!visited[getIndex(other,net)])
                println("vamos a ir a $(other.number)");
                traverse(net,other,visited);
            end
        end
    end
end

# need function to check if after updateContainRoot! there is no place for the root
# careful because updateContainRoot changes things, so maybe we want to be careful and only change
# if the new hybridization is going to stay

# think of the process of adding a hybrid edge:
# updateInCycle: what happens if cycle intersects, can we go back?
# updateContainRoot: what happens if containRoot is empty, can we go back?

# todo: function to create an hybrid edge:
# - make sure the hybridization is "identifiable": not between the same edge, or in a cherry
# - detect whether the new cycle would overlap with another cycle already in the network.
#   just check that the 2 edges to be connected are not already marked as
#   being on a cycle: updateInCycle! returns false
# - detect where the cycle is: i think it always starts in the hybrid node, so simply use searchHybridNode, or use the hybrid node just created
# - mark edges along the cycle with the number of the hybrid edge/node: updateInCycle!
# - create the new nodes and edges, with correct hybrid labels
# - mark which edges can contain the root, check that the set of edges that
#   can contain the root is non-empty: updateContainRoot, still need function to check if empty
# - check cycle configuration (value of k, and clade sizes ni)
#   if bad triangle: set gammaz and gamma2z for appropriate nodes
#   if bad diamond: gammaz for the "other" node (one just created) of each hybrid edge
#   if some parameters need to be set to 0:
# - identify the second hybrid edge, mark it as hybrid
# - depending on gamma, mark one of the 2 edges as the major "tree" edge

# todo: functions to propose a new network
# example: pick 2 edges and add a hybrid edge to link the 2
# todo: function to change direction of hybrid edge (hybrid edge  =  hybrid&&!isMajor),
#                    source or recipient of either hybrid edge, to propose new network


# todo: function readNetwork!(network::Network, string) # check string as parameter
# C function to read in tree (recursive) in mbsum*,
# maybe start reading a tree, and then add the hybrid edge(s)
# string will contain the parenthetical format, maybe not needed as parameter, but as return


# todo: function printTopology!(string, network::Network) # parameters

# todo: function network2Tree(network::Network) function to remove a hybrid edge and transform the network in tree?

# todo: function to reduce network to quartet: think of rules of how to remove hybrid edges, and when do we need to keep them and when not.

# todo: function to check that everything in network makes sense (gamma, t, gammaz, hybrid edges pointing at hybrid nodes, 2 hybrid edges: one major, one minor)

# todo: function to identify bad diamond/triangle in a network?
