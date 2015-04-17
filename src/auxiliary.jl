# auxiliary functions for all the other methods
# originally in functions.jl
# Claudia February 2015
#####################

# ----- aux general functions ---------------

#based in coupon's collector: E+sqrt(V)
function coupon(n::Number)
    return n*log(n) + n
end

function binom(n::Number,k::Number)
    n >= k || return 0
    n == 1 && return 1
    k == 0 && return 1
    binom(n-1,k-1) + binom(n-1,k) #recursive call
end

function approxEq(a::Number,b::Number,absTol::Number,relTol::Number)
    if(a<eps() || b<eps())
        abs(a-b) < absTol
    else
        abs(a-b) < relTol*eps(abs(a)+abs(b))
    end
end

approxEq(a::Number,b::Number) = approxEq(a,b,1e-5,100)


function isEqual(n1::Node,n2::Node)
    return (n1.number == n2.number && approxEq(n1.gammaz,n2.gammaz) && n1.inCycle == n2.inCycle)
end

function isEqual(n1::Edge,n2::Edge)
    return (n1.number == n2.number && approxEq(n1.length,n2.length))
end

function isEqual(net1::HybridNetwork, net2::HybridNetwork)
    result = true
    result &= (net1.numTaxa == net2.numTaxa)
    result &= (net1.numNodes == net2.numNodes)
    result &= (net1.numEdges == net2.numEdges)
    ## result &= (net1.node == net2.node)
    ## result &= (net1.edge == net2.edge)
    result &= (net1.root == net2.root)
    result &= (net1.names == net2.names)
##    result &= (net1.hybrid == net2.hybrid)
    result &= (net1.numHybrids == net2.numHybrids)
##    result &= (net1.leaf == net2.leaf)
    result &= (net1.ht == net2.ht)
    result &= (net1.numht == net2.numht)
    result &= (net1.numBad == net2.numBad)
    result &= (net1.hasVeryBadTriangle == net2.hasVeryBadTriangle)
    result &= (net1.index == net2.index)
    result &= (net1.loglik == net2.loglik)
    return result
end


#------------- EDGE functions --------------------#

# warning: node needs to be defined as hybrid before adding to a
#          hybrid edge. First, an edge is defined as hybrid, and then
#          the nodes are added to it. If the node added is leaf, the
#          edge length is set unidentifiable (as it is external edge)
function setNode!(edge::Edge, node::Node)
    size(edge.node,1)  !=  2 || error("vector of nodes already has 2 values");
    push!(edge.node,node);
    if(size(edge.node,1) == 1)
        if(edge.hybrid)
            if(node.hybrid)
                edge.isChild1 = true;
            else
                edge.isChild1 = false;
	    end
        end
        if(node.leaf)
            edge.istIdentifiable = false;
        else
            edge.istIdentifiable = true
        end
    else
        if(node.leaf)
            !edge.node[1].leaf || error("edge $(edge.number) has two leaves")
            edge.istIdentifiable = false;
        else
            if(edge.hybrid)
	        if(node.hybrid)
                    !edge.node[1].hybrid || error("hybrid edge $(edge.number) has two hybrid nodes");
                    edge.isChild1 = false;
	        else
	            edge.node[1].hybrid || error("hybrid edge $(edge.number) has no hybrid nodes");
	            edge.isChild1 = true;
	        end
            else #edge is tree
                if(!edge.node[1].leaf)
                    if(!node.hybrid && !edge.node[1].hybrid)
                        if(edge.fromBadDiamondI)
                            edge.istIdentifiable = false
                        else
                            edge.istIdentifiable = true
                        end
                    else
                        if(node.hybrid && (node.isBadDiamondI || node.isBadDiamondII || node.isBadTriangle))
                            edge.istIdentifiable = false
                        elseif(edge.node[1].hybrid && (edge.node[1].isBadDiamondI ||edge.node[1].isBadDiamondII || edge.node[1].isBadTriangle))
                            edge.istIdentifiable = false
                        else
                            edge.istIdentifiable = true
                        end
                    end
                else
                    edge.istIdentifiable = false
                end
            end
        end
    end
end

# warning: node needs to be defined as hybrid before adding to a hybrid edge.
#          First, an edge is defined as hybrid, and then the nodes are added to it.
#          If there is a leaf in node, the edge.istIdentifiable=false
function setNode!(edge::Edge,node::Array{Node,1})
    size(node,1) ==  2 || error("vector of nodes must have exactly 2 values")
    edge.node = node;
    if(edge.hybrid)
      if(node[1].hybrid)
          edge.isChild1 = true;
      else
          node[2].hybrid || error("hybrid edge without hybrid node");
          edge.isChild1 = false;
      end
    end
    if(edge.node[1].leaf || edge.node[2].leaf)
        edge.istIdentifiable = false;
    else
        edge.istIdentifiable = true;
    end
end


# -------------- NODE -------------------------#

function setEdge!(node::Node,edge::Edge)
   push!(node.edge,edge);
   all([!e.hybrid for e in node.edge]) ? node.hasHybEdge = false : node.hasHybEdge = true;
end

function getOtherNode(edge::Edge,node::Node)
  isequal(edge.node[1],node) ? edge.node[2] : edge.node[1]
end
# -------------- NETWORK ----------------------- #

function getIndex(node::Node, net::Network)
    i = 1;
    while(i<= size(net.node,1) && !isEqual(node,net.node[i]))
        i = i+1;
    end
    i>size(net.node,1)?error("node not in network"):return i;
end

function getIndex(edge::Edge, net::Network)
    i = 1;
    while(i<= size(net.edge,1) && !isEqual(edge,net.edge[i]))
        i = i+1;
    end
    i>size(net.edge,1)?error("edge not in network"):return i;
end

function getIndex(bool::Bool, array::Array{Bool,1})
    i = 1;
    while(i<= size(array,1) && !isequal(bool,array[i]))
        i = i+1;
    end
    i>size(array,1)?error("$(bool) not in array"):return i;
end

function getIndex(bool::Bool, array::Array{Bool,1})
    i = 1;
    while(i<= size(array,1) && !isequal(bool,array[i]))
        i = i+1;
    end
    i>size(array,1)?error("$(bool) not in array"):return i;
end

function getIndex(bool::Bool, array::Array{Any,1})
    i = 1;
    while(i<= size(array,1) && !isequal(bool,array[i]))
        i = i+1;
    end
    i>size(array,1)?error("$(bool) not in array"):return i;
end


# aux function to find the index of a string in a
# string array
function getIndex(name::ASCIIString, array::Array{ASCIIString,1})
    i = 1;
    while(i<= size(array,1) && !isequal(name,array[i]))
        i = i+1;
    end
    i>size(array,1)?error("$(name) not in array"):return i;
end

# aux function to find the index of a int64 in a
# int64 array
function getIndex(name::Int64, array::Array{Int64,1})
    i = 1;
    while(i<= size(array,1) && !isequal(name,array[i]))
        i = i+1;
    end
    i>size(array,1)?error("$(name) not in array"):return i;
end


# aux function to find the index of a node in a
# node array
function getIndex(name::Node, array::Array{Node,1})
    i = 1;
    while(i<= size(array,1) && !isequal(name,array[i]))
        i = i+1;
    end
    i>size(array,1)?error("$(name.number) not in array"):return i;
end


function getIndexNode(number::Int64,net::Network)
    try
        getIndex(true,[number==n.number for n in net.node])
    catch
        error("node number not in net.node")
    end
    return getIndex(true,[number==n.number for n in net.node])
end

function getIndexEdge(number::Int64,net::Network)
    try
        getIndex(true,[number==n.number for n in net.edge])
    catch
        error("edge number not in net.edge")
    end
    return getIndex(true,[number==n.number for n in net.edge])
end

# find the index of an edge in node.edge
function getIndexEdge(edge::Edge,node::Node)
    getIndex(true,[isequal(edge,e) for e in node.edge])
end

# find the index of an edge with given number in node.edge
function getIndexEdge(number::Int64,node::Node)
    getIndex(true,[isequal(edge,e) for e in node.edge])
end

# find the index of a node in edge.node
function getIndexNode(edge::Edge,node::Node)
    size(edge.node,1) == 2 || warn("this edge $(edge.number) has more or less than 2 nodes: $([n.number for n in edge.node])")
    if(isequal(node,edge.node[1]))
        return 1
    elseif(isequal(node,edge.node[2]))
        return 2
    else
        error("node not in edge.node")
    end
end

# function to find hybrid index in net.hybrid
function getIndexHybrid(node::Node, net::Network)
    node.hybrid || error("node $(node.number) is not hybrid so it cannot be in net.hybrid")
    i = 1;
    while(i<= size(net.hybrid,1) && !isEqual(node,net.hybrid[i]))
        i = i+1;
    end
    i>size(net.hybrid,1)?error("hybrid node not in network"):return i;
end

# function to find leaf index in qnet.leaf
function getIndexLeaf(node::Node, net::Network)
    if(node.leaf)
        i = 1;
        while(i<= size(net.leaf,1) && !isEqual(node,net.leaf[i]))
            i = i+1;
        end
        i>size(net.leaf,1)?error("leaf node not in network"):return i;
    else
        error("node $(node.number) is not leaf so it cannot be in net.leaf")
    end
end


# function that given a hybrid node, it gives you the minor hybrid edge
function getHybridEdge(node::Node)
    node.hybrid || error("node $(node.number) is not hybrid node, cannot get hybrid edges")
    a = nothing;
    for(e in node.edge)
        (e.hybrid && !e.isMajor) ? a = e : nothing;
    end
    isa(a,Nothing) ? error("hybrid node $(node.number) does not have minor hybrid edge, edges: $([e.number for e in node.edge])") : return a
end


# function that given two nodes, it gives you the edge that connects them
# returns error if they are not connected by an edge
function getConnectingEdge(node1::Node,node2::Node)
    found = false;
    i = 1;
    while(i<= size(node1.edge,1) && !found)
        if(isequal(getOtherNode(node1.edge[i],node1),node2))
            found = true;
        end
        i = i+1;
    end
    if(found)
        return node1.edge[i-1]
    else
        error("nodes not connected")
    end
end

# function to check in an edge is in an array by comparing
# the edges numbers (uses isEqual)
# needed for updateHasEdge
function isEdgeNumIn(edge::Edge,array::Array{Edge,1})
    return all([!isEqual(edge,e) for e in array]) ? false : true
end

# function to check in a leaf is in an array by comparing
# the numbers (uses isEqual)
# needed for updateHasEdge
function isNodeNumIn(node::Node,array::Array{Node,1})
    return all([!isEqual(node,e) for e in array]) ? false : true
end

# function to push a Node in net.node and
# update numNodes and numTaxa
function pushNode!(net::Network, n::Node)
    push!(net.node,n);
    net.numNodes += 1;
    if(n.leaf)
        net.numTaxa += 1
        push!(net.leaf,n);
    end
    if(n.hybrid)
        pushHybrid!(net,n)
    end
end

# function to push an Edge in net.edge and
# update numEdges
function pushEdge!(net::Network, e::Edge)
    push!(net.edge,e);
    net.numEdges += 1;
end


# function to push a hybrid Node in net.hybrid and
# update numHybrids
function pushHybrid!(net::Network, n::Node)
    if(n.hybrid)
        push!(net.hybrid,n);
        net.numHybrids += 1;
    else
        error("node $(n.number) is not hybrid, so cannot be pushed in net.hybrid")
    end
end


# function to delete a Node in net.node and
# update numNodes and numTaxa for HybridNetwork
# if hybrid node, it deletes also from net.hybrid
# and updates numHybrids
# note that net.names is never updated to keep it
# accurate
function deleteNode!(net::HybridNetwork, n::Node)
    try
        index = getIndex(n,net);
    catch
        error("Node $(n.number) not in network");
    end
    index = getIndex(n,net);
    deleteat!(net.node,index);
    net.numNodes -= 1;
    if(net.root == index)
        warn("Root node deleted, first node in net.node set arbitrarily as new root")
        net.root = 1 #arbitrarily chosen, not used until writeTopology
    end
    if(n.hybrid)
       removeHybrid!(net,n)
    end
    if(n.leaf)
        removeLeaf!(net,n)
    end
end

# function to delete a Node in net.node and
# update numNodes and numTaxa for QuartetNetwork
# if hybrid node, it deletes also from net.hybrid
# and updates numHybrids
# note that net.names is never updated to keep it
# accurate
# if n is leaf, we delete from qnet.leaf
function deleteNode!(net::QuartetNetwork, n::Node)
    try
        index = getIndex(n,net);
    catch
        error("Node $(n.number) not in network");
    end
    index = getIndex(n,net);
    deleteat!(net.node,index);
    net.numNodes -= 1;
    net.numTaxa -= n.leaf ? 1 : 0;
    if(n.hybrid)
       removeHybrid!(net,n)
    end
    if(n.leaf)
        index = getIndexLeaf(n,net)
        deleteat!(net.leaf,index)
    end
end

# function to delete an Edge in net.edge and
# update numEdges from a Network
function deleteEdge!(net::Network, e::Edge)
    try
        index = getIndex(e,net);
    catch
        error("Edge not in network");
    end
    index = getIndex(e,net);
    deleteat!(net.edge,index);
    net.numEdges -= 1;
end


# function to delete a hybrid Node in net.hybrid and
# update numHybrid
# used when you do not want to delete the actual node
# only remove it from net.hybrid
function removeHybrid!(net::Network, n::Node)
    n.hybrid || error("cannot delete node $(n.number) from net.hybrid because it is not hybrid")
    try
        index = getIndexHybrid(n,net);
    catch
        error("Hybrid Node $(n.number) not in network");
    end
    index = getIndexHybrid(n,net);
    deleteat!(net.hybrid,index);
    net.numHybrids -= 1;
end

# function to delete a leaf node in net.leaf
# and update numTaxa
function removeLeaf!(net::Network,n::Node)
    if(n.leaf)
        try
            index = getIndexLeaf(n,net)
        catch
            error("Leaf node $(n.number) not in network")
        end
        index = getIndexLeaf(n,net)
        deleteat!(net.leaf,index)
        net.numTaxa -= 1
    else
        error("cannot delete node $(n.number) from net.leaf because it is not leaf")
    end
end

# function to delete an internal node with only 2 edges
function deleteIntNode!(net::Network, n::Node)
    size(n.edge,1) == 2 || error("node $(n.number) does not have only two edges")
    index = n.edge[1].number < n.edge[2].number ? 1 : 2;
    edge1 = n.edge[index];
    edge2 = n.edge[index==1?2:1];
    (!edge1.hybrid && !edge2.hybrid) || error("the two edges $([edge1.number,edge2.number]) attached to node $(node.number) must be tree edges to delete node")
    node1 = getOtherNode(edge1,n);
    node2 = getOtherNode(edge2,n);
    removeEdge!(node2,edge2);
    removeNode!(n,edge1);
    setEdge!(node2,edge1);
    setNode!(edge1,node2);
    deleteNode!(net,n);
    deleteEdge!(net,edge2);
end


# search the hybrid node(s) in network: returns the hybrid node(s)
# in an array
# throws error if no hybrid in network
function searchHybridNode(net::Network)
    suma = sum([net.node[i].hybrid?1:0 for i = 1:size(net.node,1)]);
    if(suma == 0)
        error("network has no hybrid node");
    end
    k = getIndex(true,[net.node[i].hybrid for i = 1:size(net.node,1)]);
    if(suma>1)
        a = [net.node[k]];
        count = suma-1;
        index = k;
        vect = [net.node[i].hybrid for i = 1:size(net.node,1)];
        while(count>0 && count<size(net.node,1))
            index == 1 ? vect = [false,vect[2:size(net.node,1)]] : vect = [vect[1:(index-1)],false,vect[(index+1):size(net.node,1)]]
            index = getIndex(true,vect);
            push!(a,net.node[index]);
            count = count-1;
        end
        return a
    else
        return [net.node[k]]
    end
end

# search the hybrid edges in network: returns the hybrid edges
# hybrid edges come in pairs, both edges are returned
# throws error if no hybrid in network
# check: change to return only the minor edge?
function searchHybridEdge(net::Network)
    suma = sum([net.edge[i].hybrid?1:0 for i = 1:size(net.edge,1)]);
    if(suma == 0)
        error("network has no hybrid edge");
    end
    k = getIndex(true,[net.edge[i].hybrid for i = 1:size(net.edge,1)]);
    if(suma>1)
        a = [net.edge[k]];
        count = suma-1;
        index = k;
        vect = [net.edge[i].hybrid for i = 1:size(net.edge,1)];
        while(count>0 && count<size(net.edge,1))
            index == 1 ? vect = [false,vect[2:size(net.node,1)]] : vect = [vect[1:(index-1)],false,vect[(index+1):size(net.node,1)]]
            index = getIndex(true,vect);
            push!(a,net.edge[index]);
            count = count-1;
        end
        return a
    else
        return net.edge[k]
    end
end

# print for every edge, nodes, inCycle, containRoot, istIdentifiable
function printEdges(net::QuartetNetwork)
    println("Edge\tNode1\tNode2\tInCycle\tcontainRoot\tistIdentitiable\tLength\tisHybrid\tGamma")
    for e in net.edge
        println("$(e.number)\t$(e.node[1].number)\t$(e.node[2].number)\t$(e.inCycle)\t$(e.containRoot)\t\t$(e.istIdentifiable)\t\t$(round(e.length,2))\t$(e.hybrid)\t$(round(e.gamma,4))")
    end
end

function printEdges(net::HybridNetwork)
    if(net.numBad > 0)
        warn("net has $(net.numBad) bad diamond I, gammas and some branch lengths are not identifiable, and therefore, meaningless")
    end
    println("Edge\tNode1\tNode2\tInCycle\tcontainRoot\tistIdentitiable\tLength\tisHybrid\tGamma")
    for e in net.edge
        println("$(e.number)\t$(e.node[1].number)\t$(e.node[2].number)\t$(e.inCycle)\t$(e.containRoot)\t\t$(e.istIdentifiable)\t\t$(round(e.length,2))\t$(e.hybrid)\t$(round(e.gamma,4))")
    end
end

# print for every node, inCycle and edges
function printNodes(net::Network)
    println("Node\tIn Cycle\tisHybrid\thasHybEdge\tEdges numbers")
    for n in net.node
        print("$(n.number)\t$(n.inCycle)\t\t$(n.hybrid)\t$(n.hasHybEdge)\t")
        for e in n.edge
            print("\t$(e.number)")
        end
        print("\n")
    end
end

# find the edges for a given hybrid node
# in the order: hybrid major, hybrid minor, tree edge
# if node is tree node with hybrid edges, it returns
# hybrid edge, tree edge in cycle, tree edge not in cycle
# warning: assumes any tree node with hybrid edge has two tree edges
#          one in cycle, the other not in cycle
function hybridEdges(node::Node)
    size(node.edge,1) == 3 || error("node $(node.number) has $(size(node.edge,1)) edges instead of 3");
    if(node.hybrid)
        hybmajor = nothing;
        hybminor = nothing;
        tree = nothing;
        for(e in node.edge)
            (e.hybrid && e.isMajor) ? hybmajor = e : nothing
            (e.hybrid && !e.isMajor) ? hybminor = e : nothing
            !e.hybrid ? tree = e : nothing
        end
        return hybmajor, hybminor, tree
    elseif(node.hasHybEdge)
        hybrid = nothing;
        treecycle = nothing;
        tree = nothing;
        for(e in node.edge)
            (e.hybrid) ? hybrid = e : nothing
            (!e.hybrid && e.inCycle != -1) ? treecycle = e : nothing
            (!e.hybrid && e.inCycle == -1) ? tree = e : nothing
        end
        return hybrid, treecycle, tree
    else
        #warn("node $(node.number) is not hybrid $(node.hybrid) nor tree with hybrid edges (hasHybEdge) $(node.hasHybEdge), return the node.edge in order, unless a leaf is attached, then the edge attached to leaf is last");
        edge1 = nothing
        edge2 = nothing
        edge3 = nothing
        leaffound = false
        ind = 1
        for(i in 1:3)
            if(getOtherNode(node.edge[i],node).leaf)
                leaffound = true
                edge3 = node.edge[i]
                ind = i
                break
            end
        end
        if(leaffound)
            if(ind == 1)
                return node.edge[2], node.edge[3], edge3
            elseif(ind == 2)
                return node.edge[1], node.edge[3], edge3
            elseif(ind == 3)
                return node.edge[1], node.edge[2], edge3
            end
        else
            return node.edge[1], node.edge[2], node.edge[3]
        end
    end
end

# function to get the other two edges of a node
# besides the one specified
# it is called hybridEdges, but it not restricted to hybrid
function hybridEdges(node::Node, edge::Edge)
    size(node.edge,1) == 3 || error("node $(node.number) has $(size(node.edge,1)) edges instead of 3")
    edge1 = nothing
    edge2 = nothing
    for(e in node.edge)
        if(!isequal(e,edge))
            isa(edge1,Nothing) ? edge1 = e : edge2 = e
        end
    end
    return edge1,edge2
end


# function to remove an edge from a node
# warning: deletion is final, you can only
#          have edge back by pushing it again
# warning: if the edge removed is hybrid and node is tree,
#          node.hasHybEdge is set to false
#          assuming any tree node can only have one
#          one hybrid edge
function removeEdge!(node::Node,edge::Edge)
    try
        index = getIndexEdge(edge,node);
    catch e
        if isa(e, ErrorException)
            error("edge $(edge.number) not in node $(node.number)")
        end
    end
    index = getIndexEdge(edge,node);
    deleteat!(node.edge,index);
    all([!e.hybrid for e in node.edge]) ? node.hasHybEdge = false : node.hasHybEdge = true;
end

# function to remove a node from a edge
# warning: deletion is final, you can only
#          have node back by pushing it again
# warning: only removes node from edge, edge might still
#          be in node.edge
function removeNode!(node::Node,edge::Edge)
    try
        index = getIndexNode(edge,node);
    catch e
        if isa(e, ErrorException)
            error("node $(node.number) not in edge or strange edge with more than 2 nodes")
        end
    end
    index = getIndexNode(edge,node);
    deleteat!(edge.node,index);
end


# ----------------------------------------------------------------------------------------

# setLength
# warning: allows to change edge length for istIdentifiable=false
#          but issues a warning
# negative=true means it allows negative branch lengths (useful in qnet typeHyb=4)
function setLength!(edge::Edge, new_length::Float64, negative::Bool)
    (negative || new_length >= 0) || error("length has to be nonnegative: $(new_length), cannot set to edge $(edge.number)")
    new_length >= -log(1.5) || error("length can be negative, but not too negative (greater than $(-log(1.5))) or majorCF<0: new length is $(new_length)")
    #println("setting length $(new_length) to edge $(edge.number)")
    if(new_length > 10.0)
        new_length = 10.0;
    end
    edge.length = new_length;
    edge.y = exp(-new_length);
    edge.z = 1 - edge.y;
    #edge.istIdentifiable || warn("set edge length for edge $(edge.number) that is not identifiable")
end

setLength!(edge::Edge, new_length::Float64) = setLength!(edge, new_length, false)


# setGamma
# warning in the bad diamond/triangle cases because gamma is not identifiable
# updates isMajor according to gamma value
# changeOther = true, looks for the other hybrid edge and changes gamma too
# read = true, function called in readSubtree, needs to be the old one
function setGamma!(edge::Edge, new_gamma::Float64, changeOther::Bool, read::Bool)
    new_gamma >= 0 || error("gamma has to be positive: $(new_gamma)")
    new_gamma <= 1 || error("gamma has to be less than 1: $(new_gamma)")
    edge.hybrid || error("cannot change gamma in a tree edge");
    edge.isChild1 ? ind = 1 : ind = 2 ; # hybrid edge pointing at node 1 or 2
    node = edge.node[ind]
    node.hybrid || warn("hybrid edge $(edge.number) not pointing at hybrid node")
    !node.isBadDiamondI || warn("bad diamond situation: gamma not identifiable")
    if(!read)
        edges = hybridEdges(node,edge)
        length(edges) == 2 || error("strange here: node $(node.number) should have 3 edges and it has $(length(edges)+1).")
        if(edges[1].hybrid && !edges[2].hybrid)
            ind = 1
        elseif(edges[2].hybrid && !edges[1].hybrid)
            ind = 2
        else
            error("strange hybrid node $(node.number) with only one hybrid edge or with three hybrid edges")
        end
        if(changeOther)
            if(!approxEq(new_gamma,0.5))
                edge.gamma = new_gamma;
                edge.isMajor = (new_gamma>0.5) ? true : false
                edges[ind].gamma = 1 - new_gamma;
                edges[ind].isMajor = (new_gamma<0.5) ? true : false
            else #new gamma is 0.5
                edge.gamma = new_gamma
                edge.isMajor = true
                edges[ind].gamma = 1 - new_gamma
                edges[ind].isMajor = false
            end
        else
            if(!approxEq(new_gamma,0.5))
                edge.gamma = new_gamma;
                edge.isMajor = (new_gamma>0.5) ? true : false
            else #new gamma is 0.5
                edge.gamma = new_gamma
                edge.isMajor = !edges[ind].isMajor
            end
        end
    else # comes from readSubtree
        edge.gamma = new_gamma;
        edge.isMajor = (new_gamma>=0.5) ? true : false
    end
end

setGamma!(edge::Edge, new_gamma::Float64, changeOther::Bool) = setGamma!(edge, new_gamma, changeOther, false)

function numTreeEdges(net::HybridNetwork)
    2*net.numTaxa - 3 + net.numHybrids
end

function numIntTreeEdges(net::HybridNetwork)
    2*net.numTaxa - 3 + net.numHybrids - net.numTaxa
end
