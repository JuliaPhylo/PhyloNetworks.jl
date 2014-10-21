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

# aux function to read subtree
# warning: s IOStream needs to be a file, not a stream converted
#          from string
# warning: reads additional info :length:bootstrap:gamma
#          if gamma is read, assumes hybrid edge
# warning: does not allow for name of internal nodes without # after: (1,2)A,...
# fixit: it would be better to assure that all tree edges have gamma 1.0
#        two stages: clean network to verify all this: gammas, sum of gamma =1, etc
function readSubtree!(s::IOStream, parent::Node, numLeft::Array{Int64,1},net::HybridNetwork)
    c = Base.peekchar(s)
    e = nothing;
    if(c =='(')
       numLeft[1] += 1
       #println(numLeft)
       n = Node(-1*numLeft[1],false);
       c = read(s,Char)
       bl = readSubtree!(s,n,numLeft,net)
       c = advance!(s,c,numLeft)
       br = false;
       if(c == ',')
           br = readSubtree!(s,n,numLeft,net);
           c = advance!(s,c,numLeft)
       end
       if(c != ')')
           a = readall(s);
           error("Expected right parenthesis after left parenthesis $(numLeft[1]) but read $(c). The remainder of line is $(a).")
       end
        c = Base.peekchar(s);
        if(isdigit(c) || isalpha(c))
            num = readNum(s,c,net,true);
            n.number = num;
            c = Base.peekchar(s);
            if(c == '#')
                c = read(s,Char);
                c = Base.peekchar(s);
               if(isdigit(c) || isalpha(c))
                   if(c != 'L' && c != 'R' && c != 'H')
                       warn("Expected H, R or LGT after # but received $(c) in left parenthesis $(numLeft[1]).")
                   end
                   n.hybrid = true;
                   num = readNum(s,c,net,false); #ignored because there is node name
               else
                   a = readall(s);
                   error("Expected name after # but received $(c) in left parenthesis $(numLeft[1]). Remaining is $(a).")
               end
            else
                a = readall(s);
                error("Expected # but received $(c): internal node $(num) with name but no hybrid #. Remaining is $(a).")
            end
        elseif(c == '#')
            c = read(s,Char);
            c = Base.peekchar(s);
            if(isdigit(c) || isalpha(c))
                if(c != 'L' && c != 'R' && c != 'H')
                    warn("Expected H, R or LGT after # but received $(c) in left parenthesis $(numLeft[1])")
                end
                println("making $(n.number) hybrid")
                n.hybrid = true;
                num = readNum(s,c,net,true);
                n.number = num;
            else
                a = readall(s);
                error("Expected name after # but received $(c) in left parenthesis $(numLeft[1]). Remaining is $(a).")
            end
        end
        println("aqui n.hybrid es $(n.hybrid), bl is $(bl), br is $(br)")
       if((bl && br) || n.hybrid)
           pushNode!(net,n);
           e = Edge(net.numEdges+1);
           pushEdge!(net,e);
           setNode!(e,[n,parent]);
           setEdge!(n,e);
           setEdge!(parent,e);
           n.leaf = false;
           n.hybrid ? e.hybrid = true : e.hybrid = false
       else
           println("vamos a borrar el node $(n.number) xq solo tiene un child y es hybrid=$(n.hybrid)")
           if(size(n.edge,1) == 1) # internal node only has one child but is not hybrid
               edge = n.edge[1]; # assume it has only one edge
               child = getOtherNode(edge,n);
               setNode!(edge,[child,parent]);
               setEdge!(parent,edge);
           end
           println("not sure if this return should be here")
           return true # fixit: do we want this return here? what about (1,(2)3:0.9)?
       end
    elseif(isdigit(c) || isalpha(c) || c == '#')
        hybrid = false;
        if(c == '#') # maybe leaf has no name, only #H1
            c = read(s,Char);
            c = Base.peekchar(s);
            hybrid = true;
        end
        num = readNum(s,c,net,true)
        #println("creating node $(num)")
        n = Node(num,true);
        n.hybrid = hybrid;
        pushNode!(net,n);
        e = Edge(net.numEdges+1);
        pushEdge!(net,e);
        setNode!(e,[n,parent]);
        setEdge!(n,e);
        setEdge!(parent,e);
        n.hybrid ? e.hybrid = true : e.hybrid =false
    else
        error("Expected beginning of subtree but read $(c)");
    end
    c = Base.peekchar(s)
    if(isa(e,Nothing))
        return false
    end
    if(c == ':')
        c = read(s,Char);
        c = Base.peekchar(s);
        if(isdigit(c))
            length = readFloat(s,c);
            setLength!(e,length);
            c = Base.peekchar(s);
            if(c == ':')
                c = read(s,Char);
                c = Base.peekchar(s);
                if(isdigit(c))
                    length = readFloat(s,c); #bootstrap value
                    c = Base.peekchar(s);
                    if(c == ':')
                        c = read(s, Char);
                        c = Base.peekchar(s);
                        if(isdigit(c))
                            length = readFloat(s,c); #gamma
                            if(!e.hybrid)
                                warn("gamma read for current edge $(e.number) but it is not hybrid, so gamma=$(length) ignored")
                            else
                                setGamma!(e,length);
                            end
                        else
                            error("third colon : without gamma value after in $(numLeft[1]) left parenthesis")
                        end
                    else
                        e.hybrid ? error("hybrid edge $(e.number) read but without gamma value in left parenthesis $(numLeft[1])") : nothing
                    end
                elseif(c == ':')
                    c = read(s, Char);
                    c = Base.peekchar(s);
                    if(isdigit(c))
                        length = readFloat(s,c); #gamma
                        if(!e.hybrid)
                            warn("gamma read for current edge $(e.number) but it is not hybrid, so gamma=$(length) ignored")
                        else
                            setGamma!(e,length);
                        end
                    else
                        error("third colon : without gamma value after in $(numLeft[1]) left parenthesis")
                    end
                end
            end
        elseif(c == ':')
            c = read(s,Char);
            c = Base.peekchar(s);
            if(isdigit(c))
                length = readFloat(s,c); #bootstrap value
                c = Base.peekchar(s);
                if(c == ':')
                    c = read(s, Char);
                    c = Base.peekchar(s);
                    if(isdigit(c))
                        length = readFloat(s,c); #gamma
                        if(!e.hybrid)
                            warn("gamma read for current edge $(e.number) but it is not hybrid, so gamma=$(length) ignored")
                        else
                            setGamma!(e,length);
                        end
                    else
                        error("third colon : without gamma value after in $(numLeft[1]) left parenthesis")
                    end
                else
                    e.hybrid ? error("hybrid edge $(e.number) read but without gamma value in left parenthesis $(numLeft[1])") : nothing
                end
            elseif(c == ':')
                c = read(s, Char);
                c = Base.peekchar(s);
                if(isdigit(c))
                    length = readFloat(s,c); #gamma
                    if(!e.hybrid)
                        warn("gamma read for current edge $(e.number) but it is not hybrid, so gamma=$(length) ignored")
                    else
                        setGamma!(e,length);
                    end
                else
                    error("third colon : without gamma value after in left parenthesis number $(numLeft[1])")
                end
            end
        end
    end
    return true
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
