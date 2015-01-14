# functions written in classes.jl and moved here after tested
# for pseudolikelihood implementation (Stage2)
# Claudia August 2014
#
# in julia: include("functions.jl")

# tests of functions in examples_classes.jl outside git_laptop

#------------- EDGE functions --------------------#

# warning: node needs to be defined as hybrid before adding to a
#          hybrid edge. First, an edge is defined as hybrid, and then
#          the nodes are added to it. If the node added is leaf, the
#          edge length is set unidentifiable (as it is external edge)
function setNode!(edge::Edge, node::Node)
  if(size(edge.node,1)  ==  2)
    error("vector of nodes already has 2 values");
  else
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
            edge.istIdentifiable=false;
        end
    else
        if(edge.hybrid)
	    if(node.hybrid)
                if(edge.node[1].hybrid)
                    error("hybrid edge has two hybrid nodes");
                else
                    edge.isChild1 = false;
	        end
	    else
	        if(!edge.node[1].hybrid)
	            Error("hybrid edge has no hybrid nodes");
	        else
	            edge.isChild1 = true;
	        end
	    end
        end
        if(node.leaf)
            if(edge.node[1].leaf)
                error("edge has two leaves")
            else
                edge.istIdentifiable=false;
            end
        end
    end
  end
end

# warning: node needs to be defined as hybrid before adding to a hybrid edge.
#          First, an edge is defined as hybrid, and then the nodes are added to it.
#          If there is a leaf in node, the edge.istIdentifiable=false
function setNode!(edge::Edge,node::Array{Node,1})
    size(node,1) !=  2 ?
    error("vector of nodes must have exactly 2 values") :
    edge.node = node;
    if(edge.hybrid)
      if(node[1].hybrid)
          edge.isChild1 = true;
      else
          if(node[2].hybrid)
              edge.isChild1 = false;
          else
              error("hybrid edge without hybrid node");
          end
      end
    end
    if(edge.node[1].leaf || edge.node[2].leaf)
        edge.istIdentifiable=false;
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
    while(i<= size(net.node,1) && !isequal(node,net.node[i]))
        i = i+1;
    end
    i>size(net.node,1)?error("node not in network"):return i;
end

function getIndex(edge::Edge, net::Network)
    i = 1;
    while(i<= size(net.edge,1) && !isequal(edge,net.edge[i]))
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
    if(size(edge.node,1) == 2)
        if(isequal(node,edge.node[1]))
            return 1
        elseif(isequal(node,edge.node[2]))
            return 2
        else
            error("node not in edge.node")
        end
    else
        error("this edge has more than 2 nodes")
    end
end

# function to find hybrid index in net.hybrid
function getIndexHybrid(node::Node, net::Network)
    if(node.hybrid)
        i = 1;
        while(i<= size(net.hybrid,1) && !isequal(node,net.hybrid[i]))
            i = i+1;
        end
        i>size(net.hybrid,1)?error("hybrid node not in network"):return i;
    else
        error("node $(node.number) is not hybrid so it cannot be in net.hybrid")
    end
end

# function to find leaf index in qnet.leaf
function getIndexLeaf(node::Node, qnet::QuartetNetwork)
    if(node.leaf)
        i = 1;
        while(i<= size(qnet.leaf,1) && !isequal(node,qnet.leaf[i]))
            i = i+1;
        end
        i>size(qnet.leaf,1)?error("leaf node not in q.network"):return i;
    else
        error("node $(node.number) is not leaf so it cannot be in qnet.leaf")
    end
end


# function that given a hybrid node, it gives you the minor hybrid edge
function getHybridEdge(node::Node)
    if(node.hybrid)
        a = nothing;
        for(e in node.edge)
            (e.hybrid && !e.isMajor) ? a = e : nothing;
        end
        isa(a,Nothing) ? error("hybrid node does not have minor hybrid edge") : return a
    else
        error("node is not hybrid node")
    end
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

# function to determine if two edges are equal
# used when comparing edges in QuartetNetwork
# and HybridNetwork (in particular in
# updateHasEdge)
# compares only the numbers of the edges
# warning: if number not uniquely determined,
#          it fails
function isequalEdge(ed1::Edge, ed2::Edge)
    if(ed1.number == ed2.number)
        return true
    else
        return false
    end
end

# function to determine if two nodes are equal
# compares only the numbers of the nodes
# warning: if number not uniquely determined,
#          it fails
function isequalNode(ed1::Node, ed2::Node)
    if(ed1.number == ed2.number)
        return true
    else
        return false
    end
end

# function to check in an edge is in an array by comparing
# the edges numbers (uses isequalEdge)
# needed for updateHasEdge
function isEdgeNumIn(edge::Edge,array::Array{Edge,1})
    return all([!isequalEdge(edge,e) for e in array]) ? false : true
end

# function to check in a leaf is in an array by comparing
# the numbers (uses isequalNode)
# needed for updateHasEdge
function isNodeNumIn(node::Node,array::Array{Node,1})
    return all([!isequalNode(node,e) for e in array]) ? false : true
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
    net.numTaxa -= n.leaf ? 1 : 0;
    if(net.root == index)
        warn("Root node deleted")
    end
    if(n.hybrid)
       removeHybrid!(net,n)
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
    if(n.hybrid)
        try
            index = getIndexHybrid(n,net);
        catch
            error("Hybrid Node $(n.number) not in network");
        end
        index = getIndexHybrid(n,net);
        deleteat!(net.hybrid,index);
        net.numHybrids -= 1;
    else
        error("cannot delete node $(n.number) from net.hybrid because it is not hybrid")
    end
end

# function to delete an internal node with only 2 edges
function deleteIntNode!(net::Network, n::Node)
    if(size(n.edge,1) == 2)
        index = n.edge[1].number < n.edge[2].number ? 1 : 2;
        edge1 = n.edge[index];
        edge2 = n.edge[index==1?2:1];
        node1 = getOtherNode(edge1,n);
        node2 = getOtherNode(edge2,n);
        removeEdge!(node2,edge2);
        removeNode!(n,edge1);
        setEdge!(node2,edge1);
        setNode!(edge1,node2);
        deleteNode!(net,n);
        deleteEdge!(net,edge2);
        edge1.hybrid = false;
    else
        error("node does not have only two edges")
    end
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
function printEdges(net::Network)
    println("Edge\tNode1\tNode2\tInCycle\tcontainRoot\tistIdentitiable\tLength\tisHybrid\tGamma")
    for e in net.edge
        println("$(e.number)\t$(e.node[1].number)\t$(e.node[2].number)\t$(e.inCycle)\t$(e.containRoot)\t\t$(e.istIdentifiable)\t\t$(e.length)\t$(e.hybrid)\t$(e.gamma)")
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
   if(node.hybrid)
       if(size(node.edge,1) == 3)
           hybmajor = nothing;
           hybminor = nothing;
           tree = nothing;
           for(e in node.edge)
               (e.hybrid && e.isMajor) ? hybmajor = e : nothing
               (e.hybrid && !e.isMajor) ? hybminor = e : nothing
               !e.hybrid ? tree = e : nothing
           end
           return hybmajor, hybminor, tree
       else
           error("node $(node.number) has $(size(node.edge,1)) edges instead of 3");
       end
   elseif(node.hasHybEdge)
       if(size(node.edge,1) == 3)
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
           error("node $(node.number) has more or less than 3 edges");
       end
   else
       warn("node $(node.number) is not hybrid $(node.hybrid) nor tree with hybrid edges (hasHybEdge) $(node.hasHybEdge), return the node.edge in order, unless a leaf is attached, then the edge attached to leaf is last");
       if(size(node.edge,1) == 3)
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
       else
           error("node $(node.number) has more or less than 3 edges");
       end
   end
end

# function to get the other two edges of a node
# besides the one specified
function hybridEdges(node::Node, edge::Edge)
    if(size(node.edge,1) == 3)
        edge1 = nothing
        edge2 = nothing
        for(e in node.edge)
            if(!isequal(e,edge))
                isa(edge1,Nothing) ? edge1 = e : edge2 = e
            end
        end
        return edge1,edge2
    else
        error("node $(node.number) has $(size(node.edge,1)) edges instead of 3");
    end
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
            error("edge not in node")
        end
    end
    index = getIndexEdge(edge,node);
    deleteat!(node.edge,index);
    all([!e.hybrid for e in node.edge]) ? node.hasHybEdge = false : node.hasHybEdge = true;
end

# function to remove an node from a edge
# warning: deletion is final, you can only
#          have node back by pushing it again
# warning: only removes node from edge, edge might still
#          be in node.edge
function removeNode!(node::Node,edge::Edge)
    try
        index = getIndexNode(edge,node);
    catch e
        if isa(e, ErrorException)
            error("node not in edge or strange edge with more than 2 nodes")
        end
    end
    index = getIndexNode(edge,node);
    deleteat!(edge.node,index);
end


# --------------------------------------- update incycle, root, gammaz -------------------------------------------

# function to update inCycle (with priority queue) after becoming part of a network
# based on program 3 CS367 with priority queue
# expected to be much faster than the other two udpateInCycle (queue and recursive)
# input: hybrid node around which we want to update inCycle
# needs module "Base.Collections"
# returns tuple: flag, nocycle, array of edges changed, array of nodes changed
#   flag: false if cycle intersects existing cycle or number of nodes in cycle < 3
#         (there is the possibility of returning edges in intersection: path)
#         true if cycle does not intersect existing cycle
#   nocycle: true if there is no cycle (instead of error). it is used in addHybridization
# calculates also the number of nodes in the cycle and put as hybrid node attribute "k"
# warning: it is not checking if hybrid node or minor hybrid edge
#          were already part of a cycle (inCycle!= -1)
#          But it is checking so for the other edges in cycle
# warning: it needs extra things: visited attribute, prev attribute
#          unlike updateInCycle recursive, but it is expected
#          to be much faster
function updateInCycle!(net::HybridNetwork,node::Node)
    if(node.hybrid)
        start = node;
        node.inCycle = node.number;
        node.k = 1;
        hybedge = getHybridEdge(node);
        hybedge.inCycle = node.number;
        last = getOtherNode(hybedge,node);
        dist = 0;
        queue = PriorityQueue();
        path = Node[];
        net.edges_changed = Edge[];
        net.nodes_changed = Node[];
        push!(net.edges_changed,hybedge);
        push!(net.nodes_changed,node);
        found = false;
        net.visited = [false for i = 1:size(net.node,1)];
        enqueue!(queue,node,dist);
        while(!found)
            if(isempty(queue))
                return false, true, net.edges_changed, net.nodes_changed
            else
                curr = dequeue!(queue);
                if(isequal(curr,last))
                    found = true;
                    push!(path,curr);
                else
                    if(!net.visited[getIndex(curr,net)])
                        net.visited[getIndex(curr,net)] = true;
                        if(isequal(curr,start))
                            for(e in curr.edge)
                                if(!e.hybrid || e.isMajor)
                                    other = getOtherNode(e,curr);
                                    other.prev = curr;
                                    dist = dist+1;
                                    enqueue!(queue,other,dist);
                                end
                            end
                        else
                            for(e in curr.edge)
                                other = getOtherNode(e,curr);
                                if(!other.leaf && !net.visited[getIndex(other,net)])
                                    other.prev = curr;
                                    dist = dist+1;
                                    enqueue!(queue,other,dist);
                                end
                            end
                        end
                    end
                end
            end
        end # end while
        curr = pop!(path);
        while(!isequal(curr, start))
            if(curr.inCycle!= -1)
                Data.Structures.enqueue!(path,curr);
                curr = curr.prev;
            else
                curr.inCycle = start.number;
                push!(net.nodes_changed, curr);
                node.k  =  node.k + 1;
                edge = getConnectingEdge(curr,curr.prev);
                edge.inCycle = start.number;
                push!(net.edges_changed, edge);
                curr = curr.prev;
            end
        end
        if(!isempty(path) || node.k<3)
            #error("new cycle intersects existing cycle")
            return false, false, net.edges_changed, net.nodes_changed
        else
            return true, false, net.edges_changed, net.nodes_changed
        end
    else
        error("node is not hybrid")
    end
end


# aux function to traverse the network for updateContainRoot
# it changes the containRoot argument to false
# of all the edges visited
# changed to recursive after Cecile's idea
# warning: it does not go accross hybrid node, minor hybrid edge
function traverseContainRoot(node::Node, edge::Edge, edges_changed::Array{Edge,1})
    if(!node.leaf && !node.hybrid)
        for(e in node.edge)
            if(!isequal(edge,e) && e.isMajor)
                other = getOtherNode(e,node);
                if(e.containRoot) # only considered changed those that were true and not hybrid
                    e.containRoot = false;
                    push!(edges_changed, e);
                end
                traverseContainRoot(other,e, edges_changed);
            end
        end
    end
end


# function to update containRoot (by default true)
# depending on the network
# input: hybrid node (can come from searchHybridNode)
# return flag, array of edges changed
#        flag: false if the set of edges to place the root is empty
function updateContainRoot!(net::HybridNetwork, node::Node)
    if(node.hybrid)
        net.edges_changed = Edge[];
        for (e in node.edge)
            if(!e.hybrid)
                other = getOtherNode(e,node);
                e.containRoot = false;
                push!(net.edges_changed,e);
                traverseContainRoot(other,e, net.edges_changed);
            end
        end
        if(all([!e.containRoot for e in net.edge]))
            return false,net.edges_changed
        else
            return true,net.edges_changed
        end
    else
        error("node is not hybrid")
    end
end

# function to identify if the network is one of the pathological cases
# see ipad notes: k = 0 (nonidentifiable), k = 1 (nonidentifiable, bad triangle I,II)
# k = 2 (bad diamond I,II)
# also checks if hybrid node has leaf child, in which case, major edge is non identifiable
# input: hybrid node around which to check (can come from searchHybridNode)
# updates gammaz with whatever edge lengths are originally in the network
# returns flag, array of edges changed (istIdentifiable)
#         false if the network is not identifiable for k=3
# warning: needs to have updateInCycle already done as it needs
#          inCycle, and k attributes
# check: assume any tree node that has hybrid Edge has only one tree edge in cycle (true?)
# fixit: we still need to be certain of generalization n> = 5 of bad diamond case,
#        here assumed: 1,1,1,> = 2 (see ipad figure)
# fixit: still unknown treatment for bad diamond II
function updateGammaz!(net::HybridNetwork, node::Node)
    if(node.hybrid)
        node.isBadTriangleI = false
        node.isBadTriangleII = false
        node.isBadDiamondI = false
        node.isBadDiamondII = false
        net.edges_changed = Edge[];
        edge_maj, edge_min, tree_edge2 = hybridEdges(node);
        other_maj = getOtherNode(edge_maj,node);
        other_min = getOtherNode(edge_min,node);
        if(node.k == 4) # could be bad diamond I,II
            edgebla,edge_min2,tree_edge3 = hybridEdges(other_min);
            edgebla,edge_maj2,tree_edge1 = hybridEdges(other_maj);
            other_min2 = getOtherNode(edge_min2,other_min);
            isLeaf1 = getOtherNode(tree_edge1,other_maj);
            isLeaf2 = getOtherNode(tree_edge2,node);
            isLeaf3 = getOtherNode(tree_edge3,other_min);
            tree_edge4 = nothing;
            for(e in other_min2.edge)
                if(isa(tree_edge4,Nothing) && e.inCycle == -1 && !e.hybrid)
                    tree_edge4 = e;
                end
            end
            if(isequal(other_min2,getOtherNode(edge_maj2,other_maj)) && isLeaf1.leaf && isLeaf2.leaf && isLeaf3.leaf) # bad diamond I
                warn("bad diamond I found")
                other_min.gammaz = edge_min.gamma*edge_min2.z;
                other_maj.gammaz = edge_maj.gamma*edge_maj2.z;
                edge_min2.istIdentifiable = false;
                edge_maj2.istIdentifiable = false;
                edge_maj.istIdentifiable = false;
                edge_min.istIdentifiable = false;
                push!(net.edges_changed,edge_min2);
                push!(net.edges_changed,edge_min);
                push!(net.edges_changed,edge_maj2);
                push!(net.edges_changed,edge_maj);
                node.isBadDiamondI = true;
            elseif(isequal(other_min2,getOtherNode(edge_maj2,other_maj)) && isLeaf1.leaf && !isLeaf2.leaf && isLeaf3.leaf && getOtherNode(tree_edge4,other_min2).leaf) # bad diamond II
                warn("bad diamond II found")
                node.isBadDiamondII = true;
                setLength!(tree_edge2,tree_edge2.length+edge_maj.length)
                setLength!(edge_maj,0.0)
                setLength!(edge_min,0.0)
                warn("bad diamond II: minor hybrid edge length set to zero and lost")
                push!(net.edges_changed,edge_maj)
                push!(net.edges_changed,edge_min)
                push!(net.edges_changed,tree_edge2)
                edge_maj.istIdentifiable = false
                edge_min.istIdentifiable = false
            end
        elseif(node.k == 3) # could be bad triangle I, II or non identifiable cases
            edgebla,tree_edge_incycle,tree_edge1 = hybridEdges(other_min);
            edgebla,edgebla,tree_edge3 = hybridEdges(other_maj);
            isLeaf1 = getOtherNode(tree_edge1,other_min);
            isLeaf2 = getOtherNode(tree_edge2,node);
            isLeaf3 = getOtherNode(tree_edge3,other_maj);
            if(isLeaf1.leaf && !isLeaf2.leaf && !isLeaf3.leaf) # bad triangle I
                warn("bad triangle I found")
                setLength!(edge_maj,edge_maj.length+tree_edge2.length)
                setLength!(tree_edge2, 0.0);
                setLength!(edge_min, 0.0);
                node.gammaz = edge_maj.gamma*edge_maj.gamma*edge_maj.z+edge_min.gamma*edge_min.gamma*tree_edge_incycle.z;
                other_min.gammaz = edge_min.gamma*tree_edge_incycle.z;
                tree_edge_incycle.istIdentifiable = false;
                edge_maj.istIdentifiable = false;
                edge_min.istIdentifiable = false;
                tree_edge2.istIdentifiable = false;
                node.isBadTriangleI = true;
                push!(net.edges_changed,tree_edge_incycle);
                push!(net.edges_changed,tree_edge2);
                push!(net.edges_changed,edge_maj);
                push!(net.edges_changed,edge_min);
            elseif(!isLeaf1.leaf && !isLeaf2.leaf && isLeaf3.leaf) # bad triangle I
                warn("bad triangle I found")
                setLength!(edge_min,edge_min.length + tree_edge2.length);
                setLength!(tree_edge2, 0.0);
                setLength!(edge_maj, 0.0);
                node.gammaz = edge_min.gamma*edge_min.gamma*edge_min.z+edge_maj.gamma*edge_maj.gamma*tree_edge_incycle.z;
                other_maj.gammaz = edge_maj.gamma*tree_edge_incycle.z;
                tree_edge_incycle.istIdentifiable = false;
                edge_maj.istIdentifiable = false;
                edge_min.istIdentifiable = false;
                tree_edge2.istIdentifiable = false;
                node.isBadTriangleI = true;
                push!(net.edges_changed,tree_edge_incycle);
                push!(net.edges_changed,tree_edge2);
                push!(net.edges_changed,edge_maj);
                push!(net.edges_changed,edge_min);
            elseif(!isLeaf1.leaf && isLeaf2.leaf && !isLeaf3.leaf) # bad triangle II
                warn("bad triangle II found")
                tree_edge3.istIdentifiable = false;
                tree_edge1.istIdentifiable = false;
                tree_edge_incycle.istIdentifiable = false;
                edge_maj.istIdentifiable = false;
                edge_min.istIdentifiable = false;
                push!(net.edges_changed,edge_maj);
                push!(net.edges_changed,edge_min);
                other_min.gammaz = tree_edge1.y*(1-(1-edge_min.gamma)*tree_edge_incycle.z)
                other_maj.gammaz = tree_edge3.y*(1-(1-edge_maj.gamma)*tree_edge_incycle.z)
                node.gammaz = tree_edge1.y*tree_edge3.y*(edge_min.gamma*tree_edge_incycle.z*tree_edge_incycle.z*(edge_min.gamma-1))
                node.isBadTriangleII = true;
                push!(net.edges_changed,tree_edge1);
                push!(net.edges_changed,tree_edge3);
                push!(net.edges_changed,tree_edge_incycle);
            elseif(sum([isLeaf1.leaf?1:0, isLeaf2.leaf?1:0, isLeaf3.leaf?1:0]) == 2) # non identifiable network
                return false, net.edges_changed
            end
        else
            if(getOtherNode(tree_edge2,node).leaf)
                edge_min.istIdentifiable = false;
                edge_maj.istIdentifiable = false;
                push!(net.edges_changed, edge_min);
                push!(net.edges_changed, edge_maj);
            end
        end
        if(node.k > 3 && !node.isBadDiamondI && !node.isBadDiamondII)
            #println("si entra el ultimo if de k>3 y no bad diamondI,II")
            edgebla,tree_edge_incycle,tree_edge1 = hybridEdges(other_min);
            if(!tree_edge_incycle.istIdentifiable)
                tree_edge_incycle.istIdentifiable = true;
                push!(net.edges_changed,tree_edge_incycle);
            end
            if(getOtherNode(tree_edge2,node).leaf)
                #println("si se da cuenta de q hybrid node tiene child leaf. here edge min is $(edge_min.number) y su istId is $(edge_min.istIdentifiable)")
                if(edge_maj.istIdentifiable)
                    edge_maj.istIdentifiable = false;
                    push!(net.edges_changed, edge_maj);
                end
                if(edge_min.istIdentifiable)
                    #println("si entra a cambiar edge min $(edge_min.number)")
                    edge_min.istIdentifiable = false;
                    push!(net.edges_changed, edge_min);
                end
            end
        end
        return true, net.edges_changed
    else
        error("node is not hybrid")
    end
end



# ---------------------------------------- undo update of new hybridization --------------------------------


# function to undo updateInCycle which returns an array
# of edges/nodes changed
function undoInCycle!(edges::Array{Edge,1},nodes::Array{Node,1})
    for(e in edges)
        e.inCycle = -1;
    end
    for(n in nodes)
        n.inCycle = -1;
    end
end

# function to undo updateContainRoot which returns an array
# of edges changed
function undoContainRoot!(edges::Array{Edge,1})
    for(e in edges)
        !e.containRoot ? e.containRoot = true : e.containRoot = false;
    end
end

# function to undo updateGammaz which returns an array
# of edges changed
# it only changes the status of istIdentifiable to true
function undoistIdentifiable!(edges::Array{Edge,1})
    for(e in edges)
        !e.istIdentifiable ? e.istIdentifiable = true : e.istIdentifiable = false;
    end
end


# function to undo updategammaz for the 4 cases:
# bad triangle I,II and bad diamond I,II
# input: hybrid node
# set length to edges that were not identifiable and
# changed the gammaz to -1
# recalculates the branch lengths in terms of gammaz
# warning: needs to know incycle attributes
# fixit: still unknown treatment for bad diamond II
function undoGammaz!(node::Node)
    if(node.hybrid)
        if(node.isBadTriangleI)
            edge_maj, edge_min, tree_edge2 = hybridEdges(node);
            other_maj = getOtherNode(edge_maj,node);
            other_min = getOtherNode(edge_min,node);
            edgebla,tree_edge_incycle,tree_edge1 = hybridEdges(other_min);
            if(node.gammaz != -1)
                setLength!(tree_edge2,-log(1-node.gammaz))
            else
                error("bad triangle I in node $(node.number) but no gammaz updated correctly")
            end
            if(other_maj.gammaz != -1)
                setLength!(tree_edge_incycle,-log(1-other_maj.gammaz))
            elseif(other_min.gammaz != -1)
                setLength!(tree_edge_incycle,-log(1-other_min.gammaz))
            else
                error("bad triangle I in node $(node.number) but no gammaz updated correctly")
            end
            setLength!(edge_maj,0.0) # t12 (tree_edge2) already has the right length
            setLength!(edge_min,0.0)
            node.gammaz = -1.0
            other_maj.gammaz = -1.0
            other_min.gammaz = -1.0
            edge_maj.istIdentifiable = true;
            edge_min.istIdentifiable = true;
            tree_edge_incycle.istIdentifiable = true;
            tree_edge2.istIdentifiable= true;
        elseif(node.isBadTriangleII)
            edge_maj, edge_min, tree_edge2 = hybridEdges(node);
            other_maj = getOtherNode(edge_maj,node);
            other_min = getOtherNode(edge_min,node);
            edgebla,tree_edge_incycle,tree_edge1 = hybridEdges(other_min);
            edgebla,edgebla,tree_edge3 = hybridEdges(other_maj)
            if(other_min.gammaz != -1)
                setLength!(tree_edge1,-log(other_min.gammaz))
            else
                error("bad triangle II in node $(node.number) but no gammaz updated correctly")
            end
            if(other_maj.gammaz != -1)
                setLength!(tree_edge3,-log(other_maj.gammaz))
            else
                error("bad triangle II in node $(node.number) but no gammaz updated correctly")
            end
            if(node.gammaz != -1)
                setLength!(tree_edge_incycle,-log(node.gammaz+other_min.gammaz*other_maj.gammaz)+log(other_min.gammaz)+log(other_maj.gammaz))
            else
                error("bad triangle II in node $(node.number) but no gammaz updated correctly")
            end
            node.gammaz = -1.0
            other_maj.gammaz = -1.0
            other_min.gammaz = -1.0
            tree_edge_incycle.istIdentifiable = true;
            tree_edge2.istIdentifiable = true;
            tree_edge1.istIdentifiable = true;
            tree_edge3.istIdentifiable = true;
         elseif(node.isBadDiamondI)
            edge_maj, edge_min, tree_edge2 = hybridEdges(node);
            other_maj = getOtherNode(edge_maj,node);
            other_min = getOtherNode(edge_min,node);
            edgebla,tree_edge_incycle1,tree_edge = hybridEdges(other_min);
            edgebla,tree_edge_incycle2,tree_edge = hybridEdges(other_maj);
            if(other_min.gammaz != -1)
                setLength!(tree_edge_incycle1,-log(1-other_min.gammaz))
            else
                error("bad diamond I in node $(node.number) but no gammaz updated correctly")
            end
            if(other_maj.gammaz != -1)
                setLength!(tree_edge_incycle2,-log(1-other_maj.gammaz))
            else
                error("bad diamond I in node $(node.number) but no gammaz updated correctly")
            end
            edge_maj.gamma = other_maj.gammaz / (other_maj.gammaz+other_min.gammaz)
            edge_min.gamma = other_min.gammaz / (other_maj.gammaz+other_min.gammaz)
            other_min.gammaz = -1.0
            other_maj.gammaz = -1.0
            tree_edge_incycle1.istIdentifiable = true;
            tree_edge_incycle2.istIdentifiable = true;
            edge_maj.istIdentifiable = true;
            edge_min.istIdentifiable = true;
        elseif(node.isBadDiamondII)
            edge_maj, edge_min, tree_edge2 = hybridEdges(node);
            edge_maj.istIdentifiable = true
        end
    else
        error("cannot undo gammaz if starting node is not hybrid")
    end
end



# ------------------------------------------------ add new hybridization------------------------------------

# function to add hybridization event
# input: edge1, edge2 are the edges to remove
#        edge3, edge4, edge5 are the new tree edges to add
#        net is the network
#        gamma is the gamma for the hybridization
# warning: assumes that edge1, edge2 are tree edges with inCycle=-1
#          assumes the hybrid edge goes from edge1 to edge2
#          sets minor hybrid edge length to zero
# this function create the hybrid node/edges and connect everything
# and deletes edge1,2 from the nodes, and removes the nodes from edge1,2
# returns the hybrid node to start future updates there
function createHybrid!(edge1::Edge, edge2::Edge, edge3::Edge, edge4::Edge, net::HybridNetwork, gamma::Float64)
    if(0 < gamma < 1)
        (edge1.hybrid || edge2.hybrid) ? error("edges to delete must be tree edges") : nothing
        (edge3.hybrid || edge4.hybrid) ? error("edges to add must be tree edges") : nothing
        pushEdge!(net,edge3);
        pushEdge!(net,edge4);
        # create hybridization
        max_node = maximum([e.number for e in net.node]);
        max_edge = maximum([e.number for e in net.edge]);
        hybrid_edge = Edge(max_edge+1,0.0,true,gamma,false);
        pushEdge!(net,hybrid_edge);
        hybrid_node = Node(max_node+1,false,true,[edge2,hybrid_edge,edge4]);
        tree_node = Node(max_node+2,false,false,[edge1,edge3,hybrid_edge]);
        setNode!(hybrid_edge,[tree_node,hybrid_node]);
        setNode!(edge3,[tree_node,edge1.node[2]]);
        setNode!(edge4,[hybrid_node,edge2.node[2]]);
        setEdge!(edge1.node[2],edge3);
        setEdge!(edge2.node[2],edge4);
        removeEdge!(edge2.node[2],edge2);
        removeEdge!(edge1.node[2],edge1);
        removeNode!(edge1.node[2],edge1);
        setNode!(edge1,tree_node);
        removeNode!(edge2.node[2],edge2);
        #[n.number for n in edge2.node]
        setNode!(edge2,hybrid_node)
        pushNode!(net,hybrid_node);
        pushNode!(net,tree_node);
        #pushHybrid!(net,hybrid_node);
        return hybrid_node
    else
        error("gamma must be between 0 and 1")
    end
end

# aux function for chooseEdgesGamma to identify
# if two edges are sisters and if they are cherry
# (have leaves)
# returns true/false for sisters, true/false for cherry
#         true/false for nonidentifiable (two leaves, k=1 node crossed by hybridization)
function sisterOrCherry(edge1::Edge,edge2::Edge)
    sisters = false
    cherry = false
    nonidentifiable = false
    node = nothing;
    if(isequal(edge1.node[1],edge2.node[1]) || isequal(edge1.node[1],edge2.node[2]))
        node = edge1.node[1];
    elseif(isequal(edge1.node[2],edge2.node[1]) || isequal(edge1.node[2],edge2.node[2]))
        node = edge1.node[2];
    end
    if(!isa(node,Nothing))
        if(size(node.edge,1) == 3)
            sisters = true
            if(getOtherNode(edge1,node).leaf && getOtherNode(edge2,node).leaf)
                cherry = true
            elseif(getOtherNode(edge1,node).leaf || getOtherNode(edge2,node).leaf)
                edge = nothing
                for(e in node.edge)
                    if(!isequal(e,edge1) && !isequal(e,edge2))
                        edge = e
                    end
                end
                if(getOtherNode(edge,node).leaf)
                    nonidentifiable = true
                end
            end
        else
            error("node found $(node.number) that does not have exactly 3 edges, it has $(size(node.edge,1)) edges instead.")
        end
    end
    return sisters, cherry, nonidentifiable
end

# aux function to addHybridization
# it chooses the edges in the network and the gamma value
# warning: chooses edge1, edge2, gamma randomly, but
#          we could do better later
# check: gamma is uniform(0,1/2) to avoid big gammas
function chooseEdgesGamma(net::HybridNetwork)
    index1 = 1;
    index2 = 1;
    while(index1 == index2 || index1 == 0 || index2 == 0 || index1 > size(net.edge,1) || index2 > size(net.edge,1) || net.edge[index1].inCycle != -1 || net.edge[index2].inCycle != -1 || cherry || nonidentifiable)
        index1 = iround(rand()*size(net.edge,1));
        index2 = iround(rand()*size(net.edge,1));
        sisters, cherry, nonidentifiable = sisterOrCherry(net.edge[index1],net.edge[index2]);
    end
    gamma = rand()*0.5;
    println("from $(edge1.number) to $(edge2.number), $(gamma)");
    return net.edge[index1],net.edge[index2],gamma
end


# aux function for addHybridization
# that takes the output edge1, edge2, gamma from
# chooseEdgesGamma and created necessary edges
# returns edge3, edge4, and adjusts edge1, edge2 to shorter length
function parameters4createHybrid!(edge1::Edge, edge2::Edge,net::HybridNetwork)
    max_edge = maximum([e.number for e in net.edge]);
    t1 = rand()*edge1.length;
    t3 = edge1.length - t1;
    edge3 = Edge(max_edge+1,t3);
    edge1.length = t1;
    t1 = rand()*edge2.length;
    t3 = edge2.length - t1;
    edge4 = Edge(max_edge+2,t3);
    edge2.length = t1;
    return edge3, edge4
end

# aux function to add the hybridization
# without checking all the updates
# returns the hybrid node of the new hybridization
# calls chooseEdgesGamma, parameter4createHybrid and createHybrid
function addHybridization!(net::HybridNetwork)
    edge1, edge2, gamma = chooseEdgesGamma(net);
    println("edge1, $(edge1.number), edge2 $(edge2.number)")
    edge3, edge4 = parameters4createHybrid!(edge1,edge2,net);
    hybrid = createHybrid!(edge1, edge2, edge3, edge4, net, gamma);
    return hybrid
end


# function to update who is the major hybrid
# after a new hybridization is added and
# inCycle is updated
# warning: needs updateInCycle! run first
function updateMajorHybrid!(net::HybridNetwork, node::Node)
    if(node.hybrid)
        for(e in node.edge)
            if(e.inCycle != -1 && !e.hybrid)
                e.hybrid = true;
                e.isMajor = true;
                getOtherNode(e,node).hasHybEdge = true;
                e.gamma = 1-getHybridEdge(node).gamma;
                isequal(e.node[1],node) ? e.isChild1 = true : e.isChild1 = false
            end
        end
    else
        error("node is not hybrid")
    end
end

# function to update everything of a new hybridization
# it follows the flow diagram in ipad
# input: new added hybrid, network,
#        updatemajor (bool) to decide if we need to update major edge
#        only need to update if new hybrid added, if read from file not needed
# returns: success (bool), hybrid, flag, nocycle, flag2, flag3
function updateAllNewHybrid!(hybrid::Node,net::HybridNetwork, updatemajor::Bool)
    flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,hybrid);
    if(nocycle)
        return false, hybrid, flag, nocycle, true, true
    else
        if(flag)
            if(updatemajor)
                updateMajorHybrid!(net,hybrid);
            end
            flag2, edgesGammaz = updateGammaz!(net,hybrid);
            if(flag2)
                flag3, edgesRoot = updateContainRoot!(net,hybrid);
                if(flag3)
                    return true, hybrid, flag, nocycle, flag2, flag3
                else
                    undoContainRoot!(edgesRoot);
                    undoistIdentifiable!(edgesGammaz);
                    undoGammaz!(hybrid);
                    undoInCycle!(edgesInCycle, nodesInCycle);
                    return false, hybrid, flag, nocycle, flag2, flag3
                end
            else
                undoistIdentifiable!(edgesGammaz);
                undoGammaz!(hybrid);
                undoInCycle!(edgesInCycle, nodesInCycle);
                return false, hybrid, flag, nocycle, flag2, true
            end
        else
            undoInCycle!(edgesInCycle, nodesInCycle);
            return false, hybrid, flag, nocycle, true, true
        end
    end
end

# function to add a new hybridization event
# it calls chooseEdgesGamma and createHybrid!
# input: network
# check: assumes that one of the two possibilities for
#        major hybrid edge gives you a cycle, true?
# warning: "while" removed, it does not attempt to add until
#          success, it attempts to add once
# returns: success (bool), hybrid, flag, nocycle, flag2, flag3
function addHybridizationUpdate!(net::HybridNetwork)
    hybrid = addHybridization!(net);
    updateAllNewHybrid!(hybrid,net,true)
end



# --------------------------------- delete hybridization -------------------------------

# function to identify inCycle (with priority queue)
# based on updateInCycle as it is the exact same code
# only without changing the incycle attribute
# only returning array of edges/nodes affected by the hybrid
# used when attempting to delete
# input: hybrid node around which we want to identify inCycle
# needs module "Base.Collections"
# returns tuple: nocycle, array of edges changed, array of nodes changed
# check: is this traversal much faster than a simple loop over
#        all edges/nodes and check if incycle==hybrid.number?
function identifyInCycle(net::Network,node::Node)
    if(node.hybrid)
        start = node;
        hybedge = getHybridEdge(node);
        last = getOtherNode(hybedge,node);
        dist = 0;
        queue = PriorityQueue();
        path = Node[];
        net.edges_changed = Edge[];
        net.nodes_changed = Node[];
        push!(net.edges_changed,hybedge);
        push!(net.nodes_changed,node);
        found = false;
        net.visited = [false for i = 1:size(net.node,1)];
        enqueue!(queue,node,dist);
        while(!found)
            if(isempty(queue))
                return true, net.edges_changed, net.nodes_changed
            else
                curr = dequeue!(queue);
                if(isequal(curr,last))
                    found = true;
                    push!(path,curr);
                else
                    if(!net.visited[getIndex(curr,net)])
                        net.visited[getIndex(curr,net)] = true;
                        if(isequal(curr,start))
                            for(e in curr.edge)
                                if(!e.hybrid || e.isMajor)
                                    other = getOtherNode(e,curr);
                                    other.prev = curr;
                                    dist = dist+1;
                                    enqueue!(queue,other,dist);
                                end
                            end
                        else
                            for(e in curr.edge)
                                other = getOtherNode(e,curr);
                                if(!other.leaf && !net.visited[getIndex(other,net)])
                                    other.prev = curr;
                                    dist = dist+1;
                                    enqueue!(queue,other,dist);
                                end
                            end
                        end
                    end
                end
            end
        end # end while
        curr = pop!(path);
        while(!isequal(curr, start))
            if(curr.inCycle == start.number)
                push!(net.nodes_changed, curr);
                edge = getConnectingEdge(curr,curr.prev);
                if(edge.inCycle == start.number)
                    push!(net.edges_changed, edge);
                end
                curr = curr.prev;
            end
        end
        return false, net.edges_changed, net.nodes_changed
    else
        error("node is not hybrid")
    end
end



# aux function to traverse the network
# similar to traverseContainRoot but only
# identifying the edges that would be changed by a given
# hybridization
# warning: it does not go accross hybrid node/edge,
#          nor tree node with minor hybrid edge
function traverseIdentifyRoot(node::Node, edge::Edge, edges_changed::Array{Edge,1})
    if(!node.leaf && !node.hybrid)
        for(e in node.edge)
            if(!isequal(edge,e) && e.isMajor)
                other = getOtherNode(e,node);
                push!(edges_changed, e);
                if(!other.hybrid)
                    if(!other.hasHybEdge)
                        traverseIdentifyRoot(other,e, edges_changed);
                    else
                        if(hybridEdges(other)[1].isMajor)
                            traverseIdentifyRoot(other,e, edges_changed);
                        end
                    end
                end
            end
        end
    end
end


# function to identify containRoot
# depending on a hybrid node on the network
# input: hybrid node (can come from searchHybridNode)
# return array of edges affected by the hybrid node
function identifyContainRoot(net::HybridNetwork, node::Node)
    if(node.hybrid)
        net.edges_changed = Edge[];
        for (e in node.edge)
            if(!e.hybrid)
                other = getOtherNode(e,node);
                push!(net.edges_changed,e);
                traverseIdentifyRoot(other,e, net.edges_changed);
            end
        end
        return net.edges_changed
    else
        error("node is not hybrid")
    end
end

# function to undo the effect of a hybridization
# and then delete it
# input: network, hybrid node, random flag
#        random = true, deletes one hybrid egde at random
#        (minor with prob 1-gamma, major with prob gamma)
#        random = false, deletes the minor edge always
# fixit: it uses the gamma of the hybrid edges even if
#          it is not identifiable
function deleteHybridizationUpdate!(net::HybridNetwork, hybrid::Node, random::Bool)
    if(hybrid.hybrid)
        nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,hybrid);
        if(!nocycle)
            edgesRoot = identifyContainRoot(net,hybrid);
            edges = hybridEdges(hybrid);
            undoGammaz!(hybrid);
            undoInCycle!(edgesInCycle, nodesInCycle);
            undoContainRoot!(edgesRoot);
            if(hybrid.isBadTriangleI || hybrid.isBadTriangleII)
                limit = 0.5
            else
                if(edges[1].gamma > 0.5 && edges[1].gamma != 1.0)
                    limit = edges[1].gamma
                else
                    error("strange major hybrid edge $(edges[1].number) with gamma either less than 0.5 or equal to 1.0")
                end
            end
            if(random)
                minor = rand() < limit ? false : true
            else
                minor = true;
            end
            deleteHybrid!(hybrid,net,minor)
        else
            error("the hybrid does not create a cycle")
        end
    else
        error("node $(hybrid.number) is not hybrid, so we cannot delete hybridization event around it")
    end
end

# function to delete a hybridization event
# input: hybrid node and network
#        minor: true (deletes minor edge), false (deletes major)
# warning: it is meant after undoing the effect of the
#          hybridization in deleteHybridizationUpdate!
#          by itself, it leaves things as is
function deleteHybrid!(node::Node,net::HybridNetwork,minor::Bool)
    if(node.hybrid)
        if(minor)
            hybedge1,hybedge2,treeedge1 = hybridEdges(node);
            other1 = getOtherNode(hybedge1,node);
            other2 = getOtherNode(hybedge2,node);
            other3 =  getOtherNode(treeedge1,node);
            if(hybedge1.number > treeedge1.number)
                setLength!(treeedge1, treeedge1.length + hybedge1.length);
                removeNode!(node,treeedge1);
                setNode!(treeedge1,other1);
                setEdge!(other1,treeedge1);
                removeEdge!(other1, hybedge1);
                deleteEdge!(net,hybedge1);
                treeedge1.containRoot = (!treeedge1.containRoot || !hybedge1.containRoot) ? false : true
            else
                hybedge1.hybrid = false;
                hybedge1.gamma = 1.0;
                hybedge1.isMajor = true;
                other1.hasHybEdge = false;
                setLength!(hybedge1, hybedge1.length + treeedge1.length);
                removeNode!(node,hybedge1);
                setNode!(hybedge1,other3);
                setEdge!(other3,hybedge1);
                removeEdge!(other3,treeedge1);
                deleteEdge!(net,treeedge1);
                hybedge1.containRoot = (!treeedge1.containRoot || !hybedge1.containRoot) ? false : true
            end
            hybindex = getIndex(true,[e.hybrid for e in other2.edge]);
            if(hybindex == 1)
                treeedge1 = other2.edge[2];
                treeedge2 = other2.edge[3];
            elseif(hybindex == 2)
                treeedge1 = other2.edge[1];
                treeedge2 = other2.edge[3];
            elseif(hybindex == 3)
                treeedge1 = other2.edge[1];
                treeedge2 = other2.edge[2];
            else
                error("strange node has more than three edges")
            end
            treenode1 = getOtherNode(treeedge1,other2);
            treenode2 = getOtherNode(treeedge2,other2);
            if(abs(treeedge1.number) > abs(treeedge2.number))
                setLength!(treeedge2, treeedge2.length + treeedge1.length);
                removeNode!(other2,treeedge2);
                setNode!(treeedge2,treenode1);
                setEdge!(treenode1,treeedge2);
                removeEdge!(treenode1,treeedge1);
                deleteEdge!(net,treeedge1);
                treeedge2.containRoot = (!treeedge1.containRoot || !treeedge2.containRoot) ? false : true
            else
                setLength!(treeedge1, treeedge2.length + treeedge1.length);
                removeNode!(other2,treeedge1);
                setNode!(treeedge1,treenode2);
                setEdge!(treenode2,treeedge1);
                removeEdge!(treenode2,treeedge2);
                deleteEdge!(net,treeedge2);
                treeedge1.containRoot = (!treeedge1.containRoot || !treeedge2.containRoot) ? false : true
            end
            #removeHybrid!(net,node);
            deleteNode!(net,node);
            deleteNode!(net,other2);
            deleteEdge!(net,hybedge2);
        else
            hybedge1,hybedge2,treeedge1 = hybridEdges(node);
            other1 = getOtherNode(hybedge1,node);
            other2 = getOtherNode(hybedge2,node);
            setLength!(treeedge1, treeedge1.length + hybedge2.length)
            removeEdge!(other2,hybedge2)
            removeNode!(node,treeedge1)
            setEdge!(other2,treeedge1)
            setNode!(treeedge1,other2)
            #removeHybrid!(net,node)
            deleteNode!(net,node)
            deleteEdge!(net,hybedge1)
            deleteEdge!(net,hybedge2)
            removeEdge!(other1,hybedge1)
            if(size(other1.edge,1) != 2)
                error("strange node $(other1.number) had 4 edges")
            end
            if(abs(other1.edge[1].number) < abs(other1.edge[2].number))
                edge = other1.edge[1]
                otheredge = other1.edge[2]
            else
                edge = other1.edge[2]
                otheredge = other1.edge[1]
            end
            setLength!(other1.edge[1], other1.edge[1].length + other1.edge[2].length)
            other3 =  getOtherNode(otheredge,other1);
            removeNode!(other1,edge)
            removeEdge!(other3,otheredge)
            setEdge!(other3,edge)
            setNode!(edge,other3)
            deleteNode!(net,other1)
            deleteEdge!(net,otheredge)
        end
    else
        error("node has to be hybrid")
    end
end


# -------------------------- change direction of minor hybrid edge ---------------------------------

# aux function to transform tree edge into hybrid edge
# input: new edge, hybrid node (needs to be attached to new edge)
#        new gamma
function makeEdgeHybrid!(edge::Edge,node::Node,gamma::Float64)
    if(!edge.hybrid)
        if(node.hybrid)
            println("estamos en make edge hybrid en edge $(edge.number) y node $(node.number)")
            println("vamos a hacer hashybedge true para $(getOtherNode(edge,node).number)")
            getOtherNode(edge,node).hasHybEdge = true
            println("$(getOtherNode(edge,node).hasHybEdge) debe ser true")
            if(size(edge.node,1) == 2)
                if(isequal(edge.node[1],node))
                    edge.isChild1 = true
                elseif(isequal(edge.node[2],node))
                    edge.isChild1 = false
                else
                    error("node $(node.number) is not attached to edge $(edge.number)")
                end
                edge.hybrid = true
                setGamma!(edge,gamma)
            else
                error("strange edge $(edge.number) has $(size(edge.node,1)) nodes instead of 2")
            end
        else
            error("to make edge $(edge.number) hybrid, you need to give the hybrid node it is going to point to and node $(node.number) is not hybrid")
        end
    else
        error("edge $(edge.number) already hybrid, cannot make it hybrid")
    end
end

# aux function to exchange who is the hybrid node
# input: current hybrid, new hybrid
# returns false if there is no need to updategammaz after
#         true if there is need to updategammaz after
function exchangeHybridNode!(current::Node,new::Node)
    if(current.hybrid && !new.hybrid)
        update = true
        new.hybrid = true
        new.k = current.k
        removeHybrid!(net,current)
        pushHybrid!(net,new)
        current.hybrid = false
        current.k = -1
        if(current.isBadDiamondI || current.isBadDiamondII)
            current.isBadDiamondI = false
            current.isBadDiamondII = false
            update = false
        elseif(current.isBadTriangleII || current.isBadTriangleI)
            current.isBadTriangleII = false
            current.isBadTriangleI = false
            update = true
        end
        return update
    else
        error("either current node $(current.number) is not hybrid: current.hybrid $(current.hybrid) or new node $(new.number) is already hybrid: new.hybrid $(new.hybrid)")
    end
end


# function to change the direction of a hybrid edge
# input: hybrid node, network
# warning: it assumes that undoGammaz has been run before
#          and that the new hybridization has been
#          approved by the containRoot criterion
# warning: it can only move direction of minor edge
#          because ow the minor hybrid edge becomes
#          tree edge (gamma=1.0)
# warning: it needs incycle attributes
# returns flag of whether it needs update gammaz
#         and new hybrid node
function changeDirection!(node::Node, net::HybridNetwork)
    if(node.hybrid)
        major,minor,tree = hybridEdges(node);
        othermin = getOtherNode(minor,node);
        othermaj = getOtherNode(major,node);
        edgebla,treecycle,edgebla = hybridEdges(othermin)
        update = exchangeHybridNode!(node,othermin)
        makeEdgeHybrid!(treecycle,othermin,major.gamma)
        if(othermin.k > 3)
            othermaj.hasHybEdge = false
        end
        major.hybrid = false
        major.gamma = 1.0
        major.isMajor = true
        return update,othermin
    else
        error("node $(node.number) is not hybrid, so we cannot change the direction of the hybrid edge")
    end
end

# function to change the direction of the minor hybrid edge
# and update necessary stepts before and after
# input: hybrid node and network
function changeDirectionUpdate!(node::Node,net::HybridNetwork)
    if(node.hybrid)
        undoGammaz!(node)
        edgesRoot = identifyContainRoot(net,node);
        undoContainRoot!(edgesRoot);
        update,hybrid = changeDirection!(node,net)
        if(hybrid.k == 3)
            updateGammaz!(net,hybrid)
        elseif(hybrid.k == 4 && update)
            updateGammaz!(net,hybrid)
        end
        updateContainRoot!(net,hybrid);
    else
        error("cannot change the direction of minor hybrid edge since node $(node.number) is not hybrid")
    end
end

# ------------------------- move origin of hybrid edge ---------------------------------


# aux function to choose the edge to move the origin
# or target of a hybrid edge
# input: network and hybrid node, and some tree edges that
#        cannot be picked (see ipad notes)
# it needs those edges as arguments, because to find them
# we need incycle attributes, but to choose we need to have
# undone the incycle
# returns new edge
# warning: excludes edges with incycle != -1, but
#          assumes the incycle for this particular
#          hybrid node have been undone already
#          also excludes tree1 and tree2 around
#          othermin, and all hybrid edges
function chooseEdge(net::HybridNetwork,node::Node,tree::Edge,tree1::Edge,tree2::Edge)
    index1 = 1;
    while(index1 == 0 || index1 > size(net.edge,1) || net.edge[index1].inCycle != -1 || isequal(net.edge[index1],tree) || isequal(net.edge[index1],tree1) || isequal(net.edge[index1],tree2) || net.edge[index1].hybrid)
        index1 = iround(rand()*size(net.edge,1));
    end
    return net.edge[index1]
end

# function to move the origin of a hybrid edge
# warning: it does not do any update
# input: necessary nodes/edges, see moveOriginUpdate and
#        ipad notes
# warning: it changes the branch lengths of newedge, tree1, tree2 to match the
#          optimum branch lengths in the corresponding other edge (see ipad notes)
function moveOrigin(node::Node, other1::Node, other2::Node, tree1::Edge, tree2::Edge,newedge::Edge)
    if(node.hybrid)
        if(size(newedge.node,1) == 2)
            node1 = newedge.node[1];
            node2 = newedge.node[2];
        else
            error("strange edge $(newedge.number) that has $(size(newedge.node,1)) nodes instead of 2")
        end
        removeEdge!(other1,tree1)
        removeNode!(other1,tree1)
        removeEdge!(node1,newedge)
        removeNode!(node1,newedge)
        setEdge!(other1,newedge)
        setNode!(newedge,other1)
        setEdge!(node1,tree1)
        setNode!(tree1,node1)
        removeEdge!(other2,tree2)
        removeNode!(other2,tree2)
        removeEdge!(node2,newedge)
        removeNode!(node2,newedge)
        setEdge!(other2,newedge)
        setNode!(newedge,other2)
        setEdge!(node2,tree2)
        setNode!(tree2,node2)
        t = newedge.length
        t1 = tree1.length
        t2 = tree2.length
        setLength!(newedge,t1+t2)
        setLength!(tree1,(t1/(t1+t2))*t)
        setLength!(tree2,(t2/(t1+t2))*t)
    else
        error("cannot move origin of hybridization because node $(node.number) is not hybrid")
    end
end

# function to move the origin of a hybrid edge
# and update everything that needs update: gammaz, incycle
# input: network, hybrid node, random flag
#        random = true, moves the origin of hybrid egde at random
#        (minor with prob 1-gamma, major with prob gamma)
#        random = false, moves the origin of the minor edge always
# returns: success (bool), flag, nocycle, flag2
function moveOriginUpdate!(net::HybridNetwork, node::Node, random::Bool)
    if(node.hybrid)
        nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,node);
        if(!nocycle)
            major,minor,tree = hybridEdges(node);
            if(random)
                if(major.gamma > 0.5 && major.gamma != 1.0)
                    othermin = rand() < major.gamma ? getOtherNode(major,node) : getOtherNode(minor,node)
                else
                    error("strange major hybrid edge $(major.number) with gamma $(major.gamma) either less than 0.5 or equal to 1.0")
                end
            else
                othermin = getOtherNode(minor,node);
            end
            edgebla, tree1, tree2 = hybridEdges(othermin);
            other1 = getOtherNode(tree1,othermin);
            other2 = getOtherNode(tree2,othermin);
            undoGammaz!(node);
            undoInCycle!(edgesInCycle, nodesInCycle);
            newedge = chooseEdge(net,node,tree,tree1,tree2)
            println("chosen edge $(newedge.number)")
            moveOrigin(node,other1,other2,tree1,tree2,newedge)
            flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,node);
            if(nocycle)
                return false,flag,nocycle,true
            else
                if(flag)
                    flag2, edgesGammaz = updateGammaz!(net,node)
                    if(flag2)
                        return true,flag,nocycle,flag2
                    else
                        undoistIdentifiable!(edgesGammaz);
                        undoGammaz!(node);
                        undoInCycle!(edgesInCycle, nodesInCycle);
                        return false, flag, nocycle, flag2
                    end
                else
                    undoInCycle!(edgesInCycle, nodesInCycle);
                    return false, flag, nocycle, true
                end
            end
        else
            error("the hybrid $(node.number) does not create a cycle")
        end
    else
        error("node $(hybrid.number) is not hybrid, so we cannot delete hybridization event around it")
    end
end


# ------------------------------- Read topology ------------------------------------------


# aux function to advance stream in readSubtree
# warning: s IOStream needs to be a file, not a stream converted
#          from string
function advance!(s::IOStream, c::Char, numLeft::Array{Int64,1})
    c = Base.peekchar(s)
    if(Base.eof(s))
        error("Tree ends prematurely while reading subtree after left parenthesis $(numLeft[1]).")
    end
    return read(s,Char)
end


# aux function to read all digits of taxon name
# it allows names with letters and numbers
# it also reads # as part of the name and returns pound=true
# it returns the node name as string as well to check if it exists already (as hybrid)
# warning: treats digit taxon numbers as strings to avoid repeated node numbers
function readNum(s::IOStream, c::Char, net::HybridNetwork, numLeft::Array{Int64,1})
    pound = 0;
    if(isdigit(c) || isalpha(c) || c == '#')
        pound += (c == '#') ? 1 : 0
        num = read(s,Char)
        c = Base.peekchar(s)
        while(isdigit(c) || isalpha(c) || c == '#')
            d = read(s,Char)
            num = string(num,d)
            if(d == '#')
                pound += 1;
               c = Base.peekchar(s);
               if(isdigit(c) || isalpha(c))
                   if(c != 'L' && c != 'H' && c != 'R')
                       warn("Expected H, R or LGT after # but received $(c) in left parenthesis $(numLeft[1]).")
                   end
               else
                   a = readall(s);
                   error("Expected name after # but received $(c) in left parenthesis $(numLeft[1]). Remaining is $(a).")
               end
            end
            c = Base.peekchar(s);
        end
        if(pound == 0)
            return size(net.names,1)+1, num, false
        elseif(pound == 1)
            return size(net.names,1)+1, num, true
        else
            a = readall(s);
            error("strange node name with $(pound) # signs. remaining is $(a).")
        end
    else
        a = readall(s);
        error("Expected int digit, alphanum or # but received $(c). remaining is $(a).");
    end
end


# aux function to read floats like length
function readFloat(s::IOStream, c::Char)
    if(isdigit(c))
        num = read(s,Char);
        c = Base.peekchar(s);
        while(isdigit(c) || c == '.')
            d = read(s,Char);
            num = string(num,d);
            c = Base.peekchar(s);
        end
        return float(num)
    else
        a = readall(s);
        error("Expected float digit after : but received $(c). remaining is $(a).");
    end
end



# aux function to read subtree
# warning: s IOStream needs to be a file, not a stream converted
#          from string
# warning: reads additional info :length:bootstrap:gamma
# warning: allows for name of internal nodes without # after: (1,2)A,...
# warning: warning if hybrid edge without gamma value, warning if gamma value (ignored) without hybrid edge

function readSubtree!(s::IOStream, parent::Node, numLeft::Array{Int64,1}, net::HybridNetwork, hybrids::Array{ASCIIString,1}, index::Array{Int64,1})
    c = Base.peekchar(s)
    e = nothing;
    hasname = false; # to know if the current node has name
    pound = false;
    if(c =='(')
       numLeft[1] += 1
       #println(numLeft)
       n = Node(-1*numLeft[1],false);
       c = read(s,Char)
       bl = readSubtree!(s,n,numLeft,net,hybrids,index)
       c = advance!(s,c,numLeft)
       br = false;
       if(c == ',')
           br = readSubtree!(s,n,numLeft,net,hybrids,index);
           c = advance!(s,c,numLeft)
       end
       if(c != ')')
           a = readall(s);
           error("Expected right parenthesis after left parenthesis $(numLeft[1]) but read $(c). The remainder of line is $(a).")
       end
        c = Base.peekchar(s);
        if(isdigit(c) || isalpha(c) || c == '#') # internal node has name
            hasname = true;
            num,name,pound = readNum(s,c,net,numLeft);
            n.number = num;
            c = Base.peekchar(s);
            if(!pound)
                warn("internal node with name without it being a hybrid node. node name might be meaningless after tree modifications.")
            end
        end
    elseif(isdigit(c) || isalpha(c) || c == '#')
        hasname = true;
        bl = true;
        num,name,pound = readNum(s,c,net,numLeft)
        n = Node(num,true);
    else
        a = readall(s);
        error("Expected beginning of subtree but read $(c), remaining is $(a).");
    end
    if(pound) # found pound sign in name
        n.hybrid = true;
        #println("encontro un hybrid $(name).")
        #println("hybrids list tiene size $(size(hybrids,1))")
        if(in(name,hybrids))
            #println("dice que $(name) esta en hybrids")
            ind = getIndex(name,hybrids);
            other = net.node[index[ind]];
            #println("other is leaf? $(other.leaf), n is leaf? $(n.leaf)")
            if(!n.leaf && !other.leaf)
                error("both hybrid nodes are internal nodes: successors of the hybrid node must only be included in the node list of a single occurrence of the hybrid node.")
            elseif(n.leaf)
                e = Edge(net.numEdges+1);
                e.hybrid = true
                e.isMajor = false;
                pushEdge!(net,e);
                setNode!(e,[other,parent]);
                setEdge!(other,e);
                setEdge!(parent,e);
            else # !n.leaf
                if(size(other.edge,1) == 1) #other should be a leaf
                    #println("other is $(other.number), n is $(n.number), edge of other is $(other.edge[1].number)")
                    otheredge = other.edge[1];
                    otherparent = getOtherNode(otheredge,other);
                    #println("parent of other is $(otherparent.number)")
                    removeNode!(other,otheredge);
                    deleteNode!(net,other);
                    setNode!(otheredge,n);
                    setEdge!(n,otheredge);
                    e = Edge(net.numEdges+1);
                    e.hybrid = true
                    setNode!(e,[n,parent]);
                    setEdge!(n,e);
                    setEdge!(parent,e);
                    pushNode!(net,n);
                    pushEdge!(net,e);
                    n.number = other.number;
                else
                    error("strange: node $(other.number) is a leaf hybrid node so it should have only one edge and it has $(size(other.edge,1))")
                end
            end
        else
            #println("dice que $(name) no esta en hybrids")
            if(bl || br)
                n.hybrid = true;
                push!(net.names,string(name));
                #println("aqui vamos a meter a $(name) en hybrids")
                push!(hybrids,string(name));
                pushNode!(net,n);
                push!(index,size(net.node,1));
                e = Edge(net.numEdges+1);
                e.hybrid = true
                n.leaf ? e.isMajor = false : e.isMajor = true
                pushEdge!(net,e);
                setNode!(e,[n,parent]);
                setEdge!(n,e);
                setEdge!(parent,e);
            end
        end
    else
        if(bl || br)
            if(hasname)
                push!(net.names,string(name));
            end
            pushNode!(net,n);
            e = Edge(net.numEdges+1);
            pushEdge!(net,e);
            setNode!(e,[n,parent]);
            setEdge!(n,e);
            setEdge!(parent,e);
        end
    end
    c = Base.peekchar(s);
    if(isa(e,Nothing))
        return false
    end
    #n.hybrid ? e.hybrid = true : e.hybrid =false
    #println("parent is $(parent.number) and hasHybEdge is $(parent.hasHybEdge) before reading :")
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
                            warn("third colon : without gamma value after in $(numLeft[1]) left parenthesis, ignored")
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
                        warn("third colon : without gamma value after in $(numLeft[1]) left parenthesis, ignored.")
                    end
                else
                    warn("second colon : read without any double in left parenthesis $(numLeft[1]), ignored.")
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
                        warn("third colon : without gamma value after in $(numLeft[1]) left parenthesis, ignored")
                    end
                else
                    e.hybrid ? warn("hybrid edge $(e.number) read but without gamma value in left parenthesis $(numLeft[1])") : nothing
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
                    warn("third colon : without gamma value after in left parenthesis number $(numLeft[1]), ignored")
                end
            else
                warn("second colon : read without any double in left parenthesis $(numLeft[1]), ignored.")
            end
        else
            warn("one colon read without double in left parenthesis $(numLeft[1]), ignored.")
        end
    end
    return true
end


# function to read topology from parenthetical format
# input: file name
# warning: allows numbers and/or letters (or #) in taxon names
#          but not other characters
function readTopology(file::String)
    net = HybridNetwork()
    try
        s = open(file)
    catch
        error("Could not find or open $(file) file");
    end
    s = open(file)
    line = readuntil(s,";");
    if(line[end] != ';')
        error("file does not end in ;")
    end
    seekstart(s)
    c = Base.peekchar(s)
    numLeft = [0]; # made Array to make it mutable
    hybrids = ASCIIString[];
    index = Int64[];
    if(c == '(')
       numLeft[1] += 1;
       #println(numLeft)
       n = Node(-1*numLeft[1],false);
       c = read(s,Char)
       b = false;
       while(c != ';')
           b |= readSubtree!(s,n,numLeft,net,hybrids,index)
           c = read(s,Char);
           if(eof(s))
               error("Tree ended while reading in subtree beginning with left parenthesis number $(numLeft[1]).")
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
           deleteEdge!(net,edge);
       else
           pushNode!(net,n);
           net.root = getIndex(n,net);
       end
    else
       error("Expected beginning of tree with ( but received $(c) instead")
    end
    cleanAfterRead!(net)
    storeHybrids!(net)
    return net
end

# aux function to solve a polytomy
# warning: chooses one resolution at random
# warning: new nodes have the same number as the node with polytomy
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
        n1 = Node(n.number,false,false,[edge3,edge4,ednew]);
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
    if(!n.hybrid)
        while(size(n.edge,1) > 3)
            solvePolytomyRecursive!(net,n);
        end
    else
        error("cannot solve polytomy in a hybrid node $(n.number).")
    end
end

# aux function to add a child to a leaf hybrid
function addChild!(net::HybridNetwork, n::Node)
    if(n.hybrid)
        ed1 = Edge(net.numEdges+1,0.0);
        n1 = Node(size(net.names,1)+1,true,false,[ed1]);
        setEdge!(n,ed1);
        setNode!(ed1,[n,n1]);
        pushNode!(net,n1);
        pushEdge!(net,ed1);
    else
        error("cannot add child to tree node.")
    end
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
        for(i in 1:size(n.edge,1))
            !n.edge[i].hybrid ? push!(hyb,n.edge[i]) : nothing
        end
        #println("hyb tiene $([e.number for e in hyb])")
        for(e in hyb)
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
#   default values of 0.9,0.1 if not present
# fixit: for only one hybrid, better error than warning
function cleanAfterRead!(net::HybridNetwork)
    mod(sum([!e.hybrid?e.gamma:0 for e in net.edge]),1) == 0 ? nothing : error("tree (not network) read and some tree edge has gamma different than 1")
    for(n in net.node)
        if(size(n.edge,1) == 2)
            if(!n.hybrid)
                deleteIntNode!(net,n);
            else
                hyb = sum([e.hybrid?1:0 for e in n.edge]);
                if(hyb == 1)
                    deleteIntNode!(net,n);
                end
            end
        end
        if(!n.hybrid)
            if(size(n.edge,1) > 3)
                warn("polytomy found in node $(n.number), random resolution chosen")
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
                    error("current hybrid node $(n.number) with name S(net.names[n.number]) has only one hybrid edge attached. there are other $(hybnodes-1) hybrids out there but this one remained unmatched")
                end
            elseif(hyb == 0)
                warn("hybrid node $(n.number) is not connected to any hybrid edges, it was transformed to tree edge")
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
                if(suma == 2)
                    warn("hybrid edges for hybrid node $(n.number) do not contain gamma value, set default: 0.9,0.1")
                    for(e in n.edge)
                        if(e.hybrid)
                            (!e.isMajor) ? setGamma!(e,0.1) : setGamma!(e,0.9)
                        end
                    end
                elseif(suma != 1)
                    ed1 = nothing
                    ed2 = nothing
                    for(e in n.edge)
                        if(e.hybrid)
                            isa(ed1,Nothing) ? ed1=e : ed2=e
                        end
                    end
                    if(ed1.gamma < 1 && ed2.gamma < 1) #both gammas were set, but contradictory
                        error("hybrid edges for hybrid node $(n.number) have gammas that do not sum up to one: $(ed1.gamma),$(ed2.gamma)")
                    elseif(ed1.gamma < 1)
                        warn("only one hybrid edge of hybrid node $(n.number) has gamma value $(ed1.gamma) set, the other edge will be assigned $(1-ed1.gamma).")
                        setGamma!(ed2,1-ed1.gamma);
                    else
                        warn("only one hybrid edge of hybrid node $(n.number) has gamma value $(ed2.gamma) set, the other edge will be assigned $(1-ed2.gamma).")
                        setGamma!(ed1,1-ed2.gamma);
                    end
                end
            end
        end
    end
end


# function to search for the hybrid nodes in a read network after cleaning it
# and store this information as a network's attribute
function storeHybrids!(net::HybridNetwork)
    flag = true;
    try
        hybrid = searchHybridNode(net)
    catch
        warn("topology read not a network, but a tree as it has no hybrid nodes")
        flag = false;
    end
    if(flag)
        hybrid = searchHybridNode(net);
        net.hybrid = hybrid;
        net.numHybrids = size(hybrid,1);
    end
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
    if(isempty(net.hybrid))
        warn("not a network read, but a tree as it does not have hybrid nodes")
        all([e.containRoot for e in net.edge]) ? nothing : error("some tree edge has contain root as false")
        all([!e.isHybrid for e in net.edge]) ? nothing : error("some edge is hybrid and should be all tree edges in a tree")
        all([!n.hasHybEdge for n in net.node]) ? nothing : error("some tree node has hybrid edge true, but it is a tree, there are no hybrid edges")
    else
        for(n in net.hybrid)
            success,hyb,flag,nocycle,flag2,flag3 = updateAllNewHybrid!(n,net,false)
            if(!success)
                error("current hybrid $(n.number) conflicts with previous hybrid by intersecting cycles: $(!flag), nonidentifiable topology: $(!flag2), empty space for contain root: $(!flag3), or does not create a cycle (probably problem with the root placement): $(nocycle).")
            end
        end
    end
end

# function to read a topology from file name and update it
# by calling updateAllReadTopology after
function readTopologyUpdate(file::String)
    net = readTopology(file)
    updateAllReadTopology!(net)
    parameters!(net)
    return net
end

# ----------------------------------------------------------------------------------------

# setLength
# warning: allows to change edge length for istIdentifiable=false
#          but issues a warning
function setLength!(edge::Edge, new_length::Float64)
  if(new_length<0)
      error("length has to be nonnegative: $(new_length)");
  else
      edge.length = new_length;
      edge.y = exp(-new_length);
      edge.z = 1 - edge.y;
      if(!edge.istIdentifiable)
          warn("set edge length for edge $(edge.number) that is not identifiable")
      end
  end
end


# setGamma
# warning: does not allow to change gamma in the bad diamond/triangle cases
# because gamma is not identifiable
# updates isMajor according to gamma value
function setGamma!(edge::Edge, new_gamma::Float64)
 if(edge.hybrid)
	if(0 < new_gamma < 1)
            edge.isChild1 ? ind = 1 : ind = 2 ; # hybrid edge pointing at node 1 or 2
            if(edge.node[ind].hybrid)
                if(edge.node[ind].isBadDiamondI || edge.node[ind].isBadDiamondII || edge.node[ind].isBadTriangleI || edge.node[ind].isBadTriangleII)
                    error("bad diamond or triangle situation: gamma not identifiable")
                else
                    edge.gamma = new_gamma;
                    edge.isMajor = (new_gamma>0.5) ? true : false
                end
            else
                warn("hybrid edge $(edge.number) not pointing at hybrid node")
                edge.gamma = new_gamma;
                edge.isMajor = (new_gamma>0.5) ? true : false
            end
	else
	     error("gamma has to be between 0 and 1: $(new_gamma)");
        end
  else
	error("cannot change gamma in a tree edge");
  end
end

# -------------------- delete Leaf-------------------------------------


# aux function to make a hybrid node a tree node
# used in deleteLeaf
# input: hybrid node
function makeNodeTree!(net::Network, hybrid::Node)
    if(hybrid.hybrid)
        warn("we make node $(hybrid.number) a tree node, but it can still have hybrid edges pointing at it")
        hybrid.gammaz = -1
        hybrid.isBadDiamondI = false
        hybrid.isBadDiamondII = false
        hybrid.isBadTriangleI = false
        hybrid.isBadTriangleII = false
        hybrid.k = -1
        removeHybrid!(net,hybrid)
        hybrid.hybrid = false
    else
        error("cannot make node $(hybrid.number) tree node because it already is")
    end
end

# aux function to make a hybrid edge tree edge
# used in deleteLeaf
# input: edge and hybrid node it pointed to
function makeEdgeTree!(edge::Edge, node::Node)
    if(edge.hybrid)
        if(node.hybrid)
            warn("we make edge $(edge.number) a tree edge, but it will still point to hybrid node $(node.number)")
            edge.hybrid = false
            edge.isMajor = true
            edge.gamma = 1.0
            getOtherNode(edge,node).hasHybEdge = false
        else
            error("need the hybrid node at which edge $(edge.number) is pointing to, node $(node.number) is tree node")
        end
    else
        error("cannot make edge $(edge.number) tree because it is tree already")
    end
end

# function to delete internal nodes until there is no
# more internal node with only two edges
# calls deleteIntLeaf! inside a while
# returns middle node (see ipad drawing)
function deleteIntLeafWhile!(net::Network, middle::Node, leaf::Node)
    while(size(middle.edge,1) == 2)
        middle = deleteIntLeaf!(net,middle,leaf)
    end
    return middle
end

# function to delete internal nodes until there is no
# more internal node with only two edges
# calls deleteIntLeaf! inside a while
function deleteIntLeafWhile!(net::Network, leafedge::Edge, leaf::Node)
    middle = getOtherNode(leafedge, leaf)
    while(size(middle.edge,1) == 2)
        middle = deleteIntLeaf!(net,middle,leaf)
    end
end


# function to delete an internal node in an external edge
# input: network, internal node and leaf
# returns the new middle
# note: leafedge is the edge that survives
function deleteIntLeaf!(net::Network, middle::Node, leaf::Node)
    println("calling deleteIntLeaf for middle $(middle.number) and leaf $(leaf.number)")
    if(size(middle.edge,1) == 2)
        if(isequal(getOtherNode(middle.edge[1],middle),leaf))
            leafedge = middle.edge[1]
            otheredge = middle.edge[2]
        elseif(isequal(getOtherNode(middle.edge[2],middle),leaf))
            leafedge = middle.edge[2]
            otheredge = middle.edge[1]
        else
            error("leaf $(leaf.number) is not attached to internal node $(middle.number) by any edge")
        end
        othernode = getOtherNode(otheredge,middle)
        setLength!(leafedge,leafedge.length + otheredge.length)
        removeNode!(middle,leafedge)
        removeEdge!(othernode,otheredge)
        setNode!(leafedge,othernode)
        setEdge!(othernode,leafedge)
        deleteNode!(net,middle)
        deleteEdge!(net,otheredge)
        return othernode
    else
        warn("internal node $(middle.number) does not have two edges only, it has $(size(middle.edge,1))")
    end
end

# function to delete an internal node in an external edge
# input: network, edge (from leaf.edge) to move in that direction and leaf
# returns the new middle
# note: leafedge is the edge that survives
function deleteIntLeaf!(net::Network, leafedge::Edge, leaf::Node)
    middle = getOtherNode(leafedge,leaf)
    if(size(middle.edge,1) == 2)
        if(isequal(getOtherNode(middle.edge[1],middle),leaf))
            otheredge = middle.edge[2]
        elseif(isequal(getOtherNode(middle.edge[2],middle),leaf))
            otheredge = middle.edge[1]
        else
            error("leaf $(leaf.number) is not attached to internal node $(middle.number) by any edge")
        end
        othernode = getOtherNode(otheredge,middle)
        setLength!(leafedge,leafedge.length + otheredge.length)
        removeNode!(middle,leafedge)
        removeEdge!(othernode,otheredge)
        setNode!(leafedge,othernode)
        setEdge!(othernode,leafedge)
        deleteNode!(net,middle)
        deleteEdge!(net,otheredge)
        return othernode
    else
        warn("internal node $(middle.number) does not have two edges only, it has $(size(middle.edge,1))")
    end
end

# function to delete a leaf from a network
# input: network, leaf node
# warning: it will delete from the actual network
#          need to create a copy before calling this
#          function
# fixit: still missing to test hybrid-> bad triangle II and
#        normal case (not bad triangle/diamond) of hasHybEdge
#        and normal case, delete not hybrid leaf bad triangle II
function deleteLeaf!(net::Network, leaf::Node)
    if(leaf.leaf && isNodeNumIn(leaf,net.leaf))
        if(size(leaf.edge,1) == 1)
            other = getOtherNode(leaf.edge[1],leaf);
            if(other.hybrid)
                edge1,edge2 = hybridEdges(other,leaf.edge[1]);
                if(!edge1.hybrid || !edge2.hybrid)
                    error("hybrid node $(other.node) does not have two hybrid edges, they are tree edges: $(edge1.number), $(edge2.number)")
                end
                other1 = getOtherNode(edge1,other);
                other2 = getOtherNode(edge2,other);
                removeEdge!(other1,edge1)
                removeEdge!(other2,edge2)
                deleteEdge!(net,edge1)
                deleteEdge!(net,edge2)
                deleteEdge!(net,leaf.edge[1])
                deleteNode!(net,other)
                deleteNode!(net,leaf)
                if(other.isBadTriangleII)
                    edge11 = isequal(getOtherNode(other1.edge[1],other1),other2) ? other1.edge[2] : other1.edge[1]
                    if(isequal(getOtherNode(other2.edge[1],other2),other1))
                        edge12 = other2.edge[2]
                        other3 = getOtherNode(other2.edge[2],other2)
                        deleteEdge!(net,other2.edge[1])
                    else
                        edge12 = other2.edge[1]
                        other3 = getOtherNode(other2.edge[1],other2)
                        deleteEdge!(net,other2.edge[2])
                    end
                    removeNode!(other1,edge11)
                    removeEdge!(other3,edge12)
                    setEdge!(other3,edge11)
                    setNode!(edge11,other3)
                    deleteEdge!(net,edge12)
                    deleteNode!(net,other2)
                    deleteNode!(net,other1)
                    setLength!(edge11,-log(other.gammaz+other1.gammaz*other2.gammaz))
                else
                    if(size(other1.edge,1) == 2)
                        node = getOtherNode(other1.edge[1],other1)
                        leaf1 = node.leaf ? node : getOtherNode(other1.edge[2],other1)
                        if(leaf1.leaf)
                            deleteIntLeafWhile!(net,other1,leaf1);
                        end
                    end
                    if(size(other2.edge,1) == 2)
                        node = getOtherNode(other2.edge[1],other1)
                        leaf1 = node.leaf ? node : getOtherNode(other2.edge[2],other2)
                        if(leaf1.leaf)
                            deleteIntLeafWhile!(net,other2,leaf1);
                        end
                    end
                end
            else
                if(other.hasHybEdge)
                    #println("entro a caso has hyb edge")
                    edge1,edge2 = hybridEdges(other,leaf.edge[1]);
                    if(!edge1.hybrid && !edge2.hybrid)
                        error("node $(other.number) has hybrid edge attribute true, but the edges $(edge1.number), $(edge2.number) are not hybrid (and the third edge has a leaf $(leaf.number)")
                    end
                    other1 = getOtherNode(edge1,other);
                    other2 = getOtherNode(edge2,other);
                    # warning: if something changed for other1, also need to do for other2
                    if(other1.hybrid)
                        if(other1.isBadDiamondI)
                            edgemaj,edgemin,edgebla = hybridEdges(other1)
                            edge4 = isequal(edge1,edgemaj) ? edgemin : edgemaj
                            other3 = getOtherNode(edge4,other1)
                            edgebla,edge3,edge5 = hybridEdges(other3)
                            leaf5 = getOtherNode(edge5,other3)
                            removeNode!(other,edge2)
                            removeEdge!(other1,edge1)
                            setNode!(edge2,other1)
                            setEdge!(other1,edge2)
                            if(other3.gammaz != -1)
                                #println("entro a cambiar length en edge $(edge2.number) con gammaz $(other3.gammaz)")
                                setLength!(edge2,-log(1-other3.gammaz))
                            else
                                error("hybrid node $(other1.number) is bad diamond, but for node $(other3.number), gammaz is not well updated, it is $(other3.gammaz)")
                            end
                            makeEdgeTree!(edge4,other1)
                            makeNodeTree!(net,other1)
                            removeEdge!(other2,edge3)
                            removeEdge!(other3,edge3)
                            deleteNode!(net,other)
                            deleteEdge!(net,edge1)
                            deleteEdge!(net,edge3)
                            removeNode!(other3,edge4)
                            removeEdge!(leaf5,edge5)
                            setNode!(edge4,leaf5)
                            setEdge!(leaf5,edge4)
                            setLength!(edge4,edge5.length)
                            deleteNode!(net,other3)
                            deleteEdge!(net,edge5)
                        elseif(other1.isBadTriangleI)
                            #println("identify node $(other1.number) as bad triangle")
                            edgemaj,edgemin,treeedge = hybridEdges(other1)
                            #println("treeedge is edge $(treeedge.number)")
                            edge3 = isequal(edge1,edgemaj) ? edgemin : edgemaj
                            removeNode!(other1,treeedge)
                            removeEdge!(other2,edge3)
                            removeEdge!(other1,edge1)
                            setEdge!(other2,treeedge)
                            setNode!(treeedge,other2)
                            if(other1.gammaz != -1)
                                setLength!(treeedge,-log(1-other1.gammaz))
                            else
                                error("hybrid node $(other1.number) is bad triangle I, but gammaz is not updated, it is $(other1.gammaz)")
                            end
                            deleteEdge!(net,edge2)
                            deleteEdge!(net,edge1)
                            deleteEdge!(net,edge3)
                            deleteNode!(net,other1)
                            deleteNode!(net,other)
                        else
                            removeEdge!(other,leaf.edge[1])
                            edgebla,edgebla,treeedge = hybridEdges(other1)
                            if(getOtherNode(treeedge,other1).leaf)
                                setLength!(edge1,edge2.length+edge1.length)
                                removeEdge!(other2,edge2)
                                removeNode!(other,edge1)
                                setNode!(edge1,other2)
                                setEdge!(other2,edge1)
                                deleteEdge!(net,edge2)
                                deleteNode!(net,other)
                            end
                        end
                        deleteNode!(net,leaf)
                        deleteEdge!(net,leaf.edge[1])
                    elseif(other2.hybrid)
                        if(other2.isBadDiamondI)
                            edgemaj,edgemin,edgebla = hybridEdges(other2)
                            edge4 = isequal(edge2,edgemaj) ? edgemin : edgemaj
                            other3 = getOtherNode(edge4,other2)
                            edgebla,edge3,edge5 = hybridEdges(other3)
                            leaf5 = getOtherNode(edge5,other3)
                            removeNode!(other,edge1)
                            removeEdge!(other2,edge2)
                            setNode!(edge1,other2)
                            setEdge!(other2,edge1)
                            if(other3.gammaz != -1)
                                setLength!(edge1,-log(1-other3.gammaz))
                            else
                                error("hybrid node $(other2.number) is bad diamond, but for node $(other3.number), gammaz is not well updated, it is $(other3.gammaz)")
                            end
                            makeEdgeTree!(edge4,other2)
                            makeNodeTree!(net,other2)
                            removeEdge!(other1,edge3)
                            removeEdge!(other3,edge3)
                            deleteNode!(net,other)
                            deleteEdge!(net,edge2)
                            deleteEdge!(net,edge3)
                            removeNode!(other3,edge4)
                            removeEdge!(leaf5,edge5)
                            setNode!(edge4,leaf5)
                            setEdge!(leaf5,edge4)
                            setLength!(edge4,edge5.length)
                            deleteNode!(net,other3)
                            deleteEdge!(net,edge5)
                        elseif(other2.isBadTriangleI)
                            edgemaj,edgemin,treeedge = hybridEdges(other2)
                            edge3 = isequal(edge2,edgemaj) ? edgemin : edgemaj
                            removeNode!(other2,treeedge)
                            removeEdge!(other1,edge3)
                            removeEdge!(other1,edge1)
                            setEdge!(other1,treeedge)
                            setNode!(treeedge,other1)
                            if(other2.gammaz != -1)
                                setLength!(treeedge,-log(1-other2.gammaz))
                            else
                                error("hybrid node $(other2.number) is bad triangle I, but gammaz is not updated, it is $(other2.gammaz)")
                            end
                            deleteEdge!(net,edge1)
                            deleteEdge!(net,edge2)
                            deleteEdge!(net,edge3)
                            deleteNode!(net,other2)
                            deleteNode!(net,other)
                        else
                            removeEdge!(other,leaf.edge[1])
                            edgebla,edgebla,treeedge = hybridEdges(other2)
                            if(getOtherNode(treeedge,other2).leaf)
                                setLength!(edge2,edge2.length+edge1.length)
                                removeEdge!(other1,edge1)
                                removeNode!(other,edge2)
                                setNode!(edge2,other1)
                                setEdge!(other1,edge2)
                                deleteEdge!(net,edge1)
                                deleteNode!(net,other)
                            end
                        end
                        deleteNode!(net,leaf)
                        deleteEdge!(net,leaf.edge[1])
                    else
                        error("node $(other.number) has hybrid edge, but neither of the other nodes $(other1.number), $(other2.number )are hybrid")
                    end
                else # other is tree node without hybrid edges
                    println("entra al caso (1)")
                    edge1,edge2 = hybridEdges(other,leaf.edge[1]);
                    other1 = getOtherNode(edge1,other);
                    other2 = getOtherNode(edge2,other);
                    removeEdge!(other,leaf.edge[1])
                    deleteNode!(net,leaf)
                    deleteEdge!(net,leaf.edge[1])
                    if(other1.leaf || other2.leaf)
                        if(other1.leaf && other2.leaf)
                            error("just deleted a leaf $(leaf.number) and its two attached nodes are leaves also $(other1.number), $(other2.number)")
                        end
                        newleaf = other1.leaf ? other1 : other2
                        middle = other
                        println("middle is $(middle.number), middle.hybrid $(middle.hybrid), middle.hasHybEdge $(middle.hasHybEdge)")
                        middle = deleteIntLeafWhile!(net,middle,newleaf)
                        println("middle is $(middle.number), middle.hybrid $(middle.hybrid), middle.hasHybEdge $(middle.hasHybEdge)")
                        if(middle.hybrid)
                            if(middle.isBadTriangleI)
                                edgemaj,edgemin,treeedge = hybridEdges(middle)
                                if(getOtherNode(edgemaj,middle).gammaz == -1)
                                    edge2 = edgemaj
                                    edge1 = edgemin
                                elseif(getOtherNode(edgemin,middle).gammaz == -1)
                                    edge2 = edgemin
                                    edge1 = edgemaj
                                else
                                    error("hybrid node $(middle.hybrid) is bad triangle I, so one of the two nodes in cycle must have gammaz = -1, and neither node has it")
                                end
                                other2 = getOtherNode(edge2,middle)
                                other1 = getOtherNode(edge1,middle)
                                edgebla,edge3,edgebla = hybridEdges(other2)
                                removeEdge!(other2,edge2)
                                removeEdge!(middle,edge2)
                                deleteEdge!(net,edge2)
                                deleteIntLeafWhile!(net,middle,getOtherNode(treeedge,middle));
                                if(other1.gammaz != -1)
                                    setLength!(edge3,-log(1-other1.gammaz));
                                else
                                    error("hybrid node $(middle.number) is bad triangle I, but node incycle $(other1.number) does not have gammaz updated: $(other1.gammaz)")
                                end
                            end
                        elseif(middle.hasHybEdge)
                            edge1,edge2,edgebla = hybridEdges(middle)
                            other2 = getOtherNode(edge1,middle)
                            if(other2.isBadTriangleII || other2.isBadTriangleI)
                                other3 = getOtherNode(edge2,middle)
                                edge3,edgebla,edge4 = hybridEdges(other3)
                                edgebla,edgebla,treeedge = hybridEdges(other2)
                                leaf2 = getOtherNode(treeedge,other2)
                                removeEdge!(middle,edge1)
                                println("delete edge1 which is $(edge1.number)")
                                deleteEdge!(net,edge1)
                                deleteIntLeafWhile!(net,middle,newleaf);
                                removeNode!(other2,treeedge)
                                removeEdge!(other3,edge3)
                                setNode!(treeedge,other3)
                                setEdge!(other3,treeedge)
                                deleteNode!(net,other2)
                                deleteEdge!(net,edge3)
                                if(other2.isBadTriangleII)
                                    if(other3.gammaz != -1)
                                        setLength!(edge4,-log(other3.gammaz))
                                    else
                                        error("hybrid node $(other2.number) is bad triangle II, but node $(other3.number) has no gammaz updated, it is $(other3.gammaz)")
                                    end
                                elseif(other2.isBadTriangleI)
                                    if(other2.gammaz != -1 && other3.gammaz != -1)
                                        setLength!(treeedge,-log(1-other2.gammaz+other3.gammaz))
                                    else
                                        error("hybrid node $(other2.number) is bad triangle I, but it does not have gammaz updated: $(other2.gammaz), nor node $(other3.number): $(other3.gammaz)")
                                    end
                                end
                            end
                        end
                    end
                end
            end
        else
            error("strange leaf with $(size(leaf.edge,1)) edges instead of 1")
        end
    else
        if(!leaf.leaf)
            error("node $(leaf.number) is not a leaf, cannot delete it")
        else
            error("node $(leaf.number) is not in net.leaf, cannot delete it")
        end
    end
end

# -------------------- extract quartet ---------------------------------------

# function to update hasEdge attribute in a
# QuartetNetwork after leaves deleted
# with deleteLeaf!
function updateHasEdge!(qnet::QuartetNetwork, net::HybridNetwork)
    warn("function to compare edges depends on edges number being unique")
    warn("assumes no bad scenario, that is, all gammas and internal t are identifiable")
    qnet.hasEdge = vcat([isNodeNumIn(n,qnet.hybrid) for n in net.hybrid],[e.istIdentifiable & isEdgeNumIn(e,qnet.edge) for e in net.edge])
end



# function to extract a quartet from a network
# input: QuartetNetwork (already created from HybridNetwork)
#        quartet: array with the 4 leaf nodes to keep
# return: QuartetNetwork with only 4 tips
# it updates qnet.hasEdge and qnet.indexht
function extractQuartet(net::HybridNetwork,quartet::Array{Node,1})
    if(size(quartet,1) != 4)
        error("quartet array should have 4 nodes, it has $(size(quartet,1))")
    else
        if(!quartet[1].leaf || !quartet[2].leaf || !quartet[3].leaf || !quartet[4].leaf)
            error("all four nodes to keep when extracting the quartet should be leaves: $([q.number for q in quartet])")
        else
            qnet = QuartetNetwork(net)
            leaves = copy(qnet.leaf)
            for(n in leaves)
                if(!isNodeNumIn(n,quartet))
                    println("delete leaf $(n.number)")
                    deleteLeaf!(qnet,n)
                end
            end
            updateHasEdge!(qnet,net)
            parameters!(qnet,net)
            return qnet
        end
    end
end

# function to extract a quartet from a Quartet object
# it calls the previous extractQuartet
# returns qnet (check: maybe not needed later) and assigns
# quartet.qnet = qnet
function extractQuartet!(net::HybridNetwork, quartet::Quartet)
    list = Node[]
    for(q in quartet.taxon)
        try
            getIndexNode(getIndex(q,net.names),net)
        catch
            error("taxon $(q) not in network")
        end
        push!(list, net.node[getIndexNode(getIndex(q,net.names),net)])
    end
    qnet = extractQuartet(net,list)
    qnet.quartetTaxon = quartet.taxon
    quartet.qnet = qnet
    #return qnet
end


# function to extract all quartets from net according
# to the array of quartets of a Data object
# it updates expCF, hasEdgem indexht
function extractQuartet!(net::HybridNetwork, d::DataCF)
    for(q in d.quartet)
        extractQuartet!(net,q)
        qnet = deepcopy(q.qnet);
        calculateExpCFAll!(qnet);
        q.qnet.expCF = qnet.expCF
    end
end


# ------------------------------- calculate expCF -------------------------------------

# ------- Identify Quartet

# function to identify the type of hybridization of a given hybrid node
# in a quartet network
# sets node.k to the updated count
# sets the node.typeHyb (1,2,3,4,5 see ipad notes) and pushes into qnet.typeHyb
# sets node.prev = other to keep track of the "other" node
# fixit: I think qnet.typeHyb is never used
function identifyQuartet!(qnet::QuartetNetwork, node::Node)
    if(node.hybrid)
        k = sum([(n.inCycle == node.number && size(n.edge,1) == 3) ? 1 : 0 for n in qnet.node])
        node.k = k
        if(k < 2)
            error("strange quartet network with a hybrid node $(node.number) but no cycle")
        elseif(k == 2)
            other = qnet.node[getIndex(true, [(n.inCycle == node.number && size(n.edge,1) == 3 && !isequal(n,node)) for n in qnet.node])]
            edgebla,edgebla,edge1 = hybridEdges(node)
            edgebla,edgebla,edge2 = hybridEdges(other)
            if(getOtherNode(edge1,node).leaf || getOtherNode(edge2,other).leaf)
                node.typeHyb = 1
                node.prev = other
                push!(qnet.typeHyb,1)
            else
                node.typeHyb = 3
                node.prev = other
                push!(qnet.typeHyb,3)
            end
        elseif(k == 3)
            edge1,edge2,edge3 = hybridEdges(node)
            if(getOtherNode(edge3,node).leaf)
                node.typeHyb = 2
                push!(qnet.typeHyb,2)
                for(n in qnet.node)
                    if(n.inCycle == node.number && size(n.edge,1) == 3 && !isequal(n,node))
                        edge1,edge2,edge3 = hybridEdges(n)
                        if(getOtherNode(edge3,n).leaf)
                            node.prev = n
                            break
                        end
                    end
                end
            else
                node.typeHyb = 4
                push!(qnet.typeHyb,4)
                for(n in qnet.node)
                    if(n.inCycle == node.number && size(n.edge,1) == 3 && !isequal(n,node))
                        node.prev = n
                        break
                    end
                end
            end
        elseif(k == 4)
            node.typeHyb = 5
            push!(qnet.typeHyb,5)
            node.prev = nothing
        else
            error("strange quartet network with $(k) nodes in cycle, maximum should be 4")
        end
    else
        error("cannot identify the hybridization around node $(node.number) because it is not hybrid node.")
    end
end

# function to identify the Quartet network as
# 1 (equivalent to tree), 2 (minor CF different)
function identifyQuartet!(qnet::QuartetNetwork)
    if(qnet.which == -1)
        if(qnet.numHybrids == 0)
            qnet.which = 1
        elseif(qnet.numHybrids == 1)
            identifyQuartet!(qnet, qnet.hybrid[1])
            if(qnet.typeHyb[1] == 5)
                qnet.which = 2
            else
                qnet.which = 1
            end
        else
            for(n in qnet.hybrid)
                identifyQuartet!(qnet, n)
            end
            if(all([(n.typeHyb != 5) for n in qnet.hybrid]))
                qnet.which = 1
            else
                if(all([(n.typeHyb == 5 || n.typeHyb == 1) for n in qnet.hybrid]))
                    warn("hybridization of type 5 found in quartet network along with other hybridizations of type 1. there is the possibility of overlapping cycles.")
                    qnet.which = 2
                else
                    error("found in the same quartet, two hybridizations with overlapping cycles: types of hybridizations are $([n.typeHyb for n in qnet.hybrid])")
                end
            end
        end
    end
end

# ---------------------- branch length optimization ---------------------------------


# function to get the branch lengths/gammas to optimize for a given
# network
# fixit: it is ignoring "bad" cases, listing all the parameters
# warning: order of parameters (h,t)
# updates net.numht also with the number of hybrid nodes and number of identifiable edges (n2,n)
function parameters(net::Network)
    warn("ignores bad cases, listing all the parameters")
    t = Float64[]
    h = Float64[]
    n = Int64[]
    n2 = Int64[]
    for(e in net.edge)
        if(e.istIdentifiable)
            push!(t,e.length)
            push!(n,e.number)
        end
        if(e.hybrid && !e.isMajor)
            push!(h,e.gamma)
            if(e.node[e.isChild1 ? 1 : 2].hybrid)
                push!(n2,e.node[e.isChild1 ? 1 : 2].number)
            else
                error("strange thing, hybrid edge $(e.number) pointing at tree node $(e.node[isChild1?1:2].number)")
            end
        end
    end
    size(t,1) == 0 ? error("net does not have identifiable branch lengths") : nothing
    return vcat(h,t),vcat(n2,n)
end

function parameters!(net::Network)
    warn("ignores bad cases, listing all the parameters")
    warn("deleting net.ht,net.numht and updating with current edge lengths (numbers)")
    net.ht,net.numht = parameters(net)
end

# function to update qnet.indexht based on net.numht
# warning: assumes net.numht is updated already with parameters!(net)
function parameters!(qnet::QuartetNetwork, net::HybridNetwork)
    if(size(net.numht,1) > 0)
        size(qnet.indexht,1) > 0 ? warn("deleting qnet.indexht to replace with info in net") : nothing
        n2 = net.numht[1:net.numHybrids]
        n = net.numht[net.numHybrids + 1 : length(net.numht)]
        qn2 = Int64[]
        qn = Int64[]
        for(e in qnet.edge)
            if(e.istIdentifiable)
                push!(qn, getIndex(e.number,n)+net.numHybrids)
            end
            if(e.hybrid && !e.isMajor)
                push!(qn2, getIndex(e.node[e.isChild1 ? 1 : 2].number,n2))
            end
        end
        qnet.indexht = vcat(qn2,qn)
    else
        error("net.numht not correctly updated, need to run parameters first")
    end
end




# ----------------------- Eliminate Hybridizations

# aux function to eliminate hybridizations
# in quartet
# input: quartet network,
#        node correspond to the hybrid node
#        internal: true if loop is in internal edge
function eliminateLoop!(qnet::QuartetNetwork, node::Node, internal::Bool)
    if(node.hybrid)
        edge1,edge2,edge3 = hybridEdges(node)
        deleteIntLeafWhile!(qnet, edge1, node)
        deleteIntLeafWhile!(qnet, edge2, node)
        other = getOtherNode(edge1,node)
        if(!isequal(getOtherNode(edge2,node),other))
            error("node $(node.number) and other $(other.number) are not the two nodes in a cycle with k=2")
        end
        removeEdge!(node,edge2)
        removeEdge!(other,edge2)
        deleteEdge!(qnet,edge2)
        if(internal)
            setLength!(edge1, -log(1-edge1.gamma*edge1.gamma*edge1.z-edge2.gamma*edge2.gamma*edge2.z))
        else
            leaf = getOtherNode(edge3,node)
            if(leaf.leaf)
                deleteIntLeafWhile!(qnet,node,leaf)
            else
                edge1,edge2,edge3 = hybridEdges(other)
                deleteIntLeafWhile!(qnet,other,getOtherNode(edge3,other))
            end
        end
    else
        error("cannot eliminate loop around node $(node.number) since it is not hybrid")
    end
end

# aux function to identify intermediate edge between two nodes
function intermediateEdge(node1::Node,node2::Node)
    edge = nothing
    for(e in node1.edge)
        if(isequal(getOtherNode(e,node1),node2))
            edge = e
        end
    end
    if(isa(edge, Nothing))
        error("nodes $(node1.number), $(node2.number) are not connected by an edge")
    end
    return edge
end

# function to eliminate a triangle hybridization
# input: quartet network,
#        node, other nodes in the hybridization
#        case: 1 (global case 2),2 (global case 4), 1 (global case 5)
function eliminateTriangle!(qnet::QuartetNetwork, node::Node, other::Node, case::Int64)
    println("start eliminateTriangle----")
    if(node.hybrid)
        edgemaj, edgemin, treeedge = hybridEdges(node)
        deleteIntLeafWhile!(qnet, edgemaj, node)
        deleteIntLeafWhile!(qnet, edgemin, node)
        if(isequal(getOtherNode(edgemaj,node),other))
            hybedge = edgemaj
            otheredge = edgemin
        elseif(isequal(getOtherNode(edgemin,node),other))
            hybedge = edgemin
            otheredge = edgemaj
        else
            error("node $(node.number) and other node $(other.number) are not connected by an edge")
        end
        println("hybedge is $(hybedge.number), otheredge is $(otheredge.number)")
        middle = qnet.node[getIndex(true, [(n.inCycle == node.number && size(n.edge,1) == 3 && !isequal(n,other) && !isequal(n,node)) for n in qnet.node])]
        println("middle node is $(middle.number) in eliminateTriangle")
        ind = getIndex(true,[(e.inCycle == node.number && !isequal(getOtherNode(e,middle),node)) for e in middle.edge])
        edge = middle.edge[ind]
        println("edge is $(edge.number) with length $(edge.length) in eliminateTriangle, will do deleteIntLeaf from middle through edge")
        deleteIntLeafWhile!(qnet,edge,middle)
        println("after deleteIntLeaf, edge $(edge.number) has length $(edge.length)")
        if(!isequal(getOtherNode(edge,middle),other))
            error("middle node $(middle.number) and other node $(other.number) are not connected by an edge")
        end
        if(case == 1)
            setLength!(edge,-log(1 - hybedge.gamma*edge.z))
            println("edge $(edge.number) length is $(edge.length) after updating")
            removeEdge!(middle,otheredge)
            removeEdge!(other,hybedge)
            removeEdge!(node,treeedge)
            setEdge!(other, treeedge)
            deleteEdge!(qnet,otheredge)
            deleteEdge!(qnet,hybedge)
            deleteNode!(qnet,node)
        elseif(case == 2)
            setLength!(edge, -log(otheredge.gamma*otheredge.gamma*otheredge.y + hybedge.gamma*otheredge.gamma*(3-edge.y) + hybedge.gamma*hybedge.gamma*hybedge.y))
            println("edge $(edge.number) length is $(edge.length) after updating")
            removeEdge!(middle,otheredge)
            removeEdge!(node,otheredge)
            deleteEdge!(qnet,otheredge)
            deleteIntLeafWhile!(qnet, node, getOtherNode(treeedge,node))
        else
            error("unknown case $(case), should be 1 or 2")
        end
    else
        error("cannot eliminate triangle around node $(node.number) since it is not hybrid")
    end
    println("end eliminateTriangle ---")
end

# function to polish quartet with hybridization type 5
# (2 minor CF different) and calculate the expCF
# CF calculated in the order 12|34, 13|24, 14|23 of the qnet.quartet.taxon
function quartetType5!(qnet::QuartetNetwork, node::Node)
    if(node.hybrid && node.typeHyb == 5)
        edge1,edge2,edge3 = hybridEdges(node)
        deleteIntLeafWhile!(qnet,edge1,node)
        deleteIntLeafWhile!(qnet,edge2,node)
        other1 = getOtherNode(edge1,node)
        other2 = getOtherNode(edge2,node)
        edgebla,edge5, edgetree1 = hybridEdges(other1)
        edgebla,edge6, edgetree2 = hybridEdges(other2)
        deleteIntLeafWhile!(qnet,edge5,other1)
        deleteIntLeafWhile!(qnet,edge6,other2)
        cf1 = edge1.gamma*(1-2/3*edge5.y) + edge2.gamma*1/3*edge6.y
        cf2 = edge1.gamma*1/3*edge5.y + edge2.gamma*(1-2/3*edge6.y)
        cf3 = edge1.gamma*1/3*edge5.y + edge2.gamma*1/3*edge6.y
        println("cf1,cf2,cf3: $(cf1),$(cf2),$(cf3)")
        leaf1 = getOtherNode(edge3,node)
        leaf2 = getOtherNode(edgetree1,other1)
        leaf3 = getOtherNode(edgetree2, other2)
        leaf4 = qnet.leaf[getIndex(true,[(!isequal(n,leaf1) && !isequal(n,leaf2) && !isequal(n,leaf3)) for n in qnet.leaf])]
        println("leaf1 is $(leaf1.number)")
        println("leaf2 is $(leaf2.number)")
        println("leaf3 is $(leaf3.number)")
        println("leaf4 is $(leaf4.number)")
        tx = whichLeaves(qnet,qnet.quartetTaxon[1],qnet.quartetTaxon[2], leaf1,leaf2,leaf3,leaf4)
        println("tx is $(tx)")
        if(tx == (1,2) || tx == (2,1) || tx == (3,4) || tx == (4,3))
            qnet.expCF[1] = cf1
            tx = whichLeaves(qnet,qnet.quartetTaxon[1],qnet.quartetTaxon[3], leaf1,leaf2,leaf3,leaf4)
            if(tx == (1,3) || tx == (3,1) || tx == (2,4) || tx == (4,2))
                qnet.expCF[2] = cf2
                qnet.expCF[3] = cf3
            elseif(tx == (1,4) || tx == (4,1) || tx == (3,2) || tx == (2,3))
                qnet.expCF[2] = cf3
                qnet.expCF[3] = cf2
            else
                error("strange quartet network, could not find which leaves correspond to taxon1, taxon3")
            end
        elseif(tx == (1,3) || tx == (3,1) || tx == (2,4) || tx == (4,2))
            qnet.expCF[1] = cf2
            tx = whichLeaves(qnet,qnet.quartetTaxon[1],qnet.quartetTaxon[3], leaf1,leaf2,leaf3,leaf4)
            if(tx == (1,2) || tx == (2,1) || tx == (3,4) || tx == (4,3))
                qnet.expCF[2] = cf1
                qnet.expCF[3] = cf3
            elseif(tx == (1,4) || tx == (4,1) || tx == (3,2) || tx == (2,3))
                qnet.expCF[2] = cf3
                qnet.expCF[3] = cf1
            else
                error("strange quartet network, could not find which leaves correspond to taxon1, taxon3")
            end
        elseif(tx == (1,4) || tx == (4,1) || tx == (2,4) || tx == (2,3))
            qnet.expCF[1] = cf3
            tx = whichLeaves(qnet,qnet.quartetTaxon[1],qnet.quartetTaxon[3], leaf1,leaf2,leaf3,leaf4)
            if(tx == (1,3) || tx == (3,1) || tx == (2,4) || tx == (4,2))
                qnet.expCF[2] = cf2
                qnet.expCF[3] = cf1
            elseif(tx == (1,2) || tx == (2,1) || tx == (3,4) || tx == (4,3))
                qnet.expCF[2] = cf1
                qnet.expCF[3] = cf2
            else
                error("strange quartet network, could not find which leaves correspond to taxon1, taxon3")
            end
        else
            error("strange quartet network, could not find which leaves correspond to taxon1, taxon2")
        end
        if(sum(qnet.expCF) < 0.99999 || sum(qnet.expCF) > 1.0000001)
            error("strange quartet network with hybridization in node $(node.number) of type 5: expCF do not add up to 1")
        end
    else
        error("cannot polish the quartet type 5 hybridization since either the node is not hybrid: $(!node.hybrid) or it has type $(node.typeHyb), different than 5")
    end
end

# function to eliminate a hybridization around a given
# hybrid node
function eliminateHybridization!(qnet::QuartetNetwork, node::Node)
    if(node.hybrid)
        if(node.typeHyb == 1)
            eliminateLoop!(qnet,node,false)
        elseif(node.typeHyb == 3)
            eliminateLoop!(qnet,node,true)
        elseif(node.typeHyb == 4)
            println("node is $(node.number), other node is $(node.prev.number)")
            eliminateTriangle!(qnet,node,node.prev,2)
        elseif(node.typeHyb == 2)
            println("node is $(node.number), other node is $(node.prev.number)")
            eliminateTriangle!(qnet,node,node.prev,1)
        elseif(node.typeHyb != 5)
            error("node type of hybridization should be 1,2,3,4 or 5, but for node $(node.number), it is $(node.typeHyb)")
        end
    else
        error("cannot eliminate hybridization around node $(node.number) since it is not hybrid node")
    end
end

# aux function to eliminate all internal nodes with only
# two edges in a quartet network with qnet.which=1
# fixit: no hay forma de saber quien es node sin tener que encontrarlo cada vez?
function internalLength!(qnet::QuartetNetwork)
    if(qnet.which == 1)
        node = qnet.node[getIndex(true,[size(n.edge,1) == 3 for n in qnet.node])]
        println("node is $(node.number)")
        edge = nothing
        for(e in node.edge)
            if(!getOtherNode(e,node).leaf)
                edge = e
            end
        end
        println("edge $(edge.number) has length $(edge.length) before lumping all internal edges into it")
        deleteIntLeafWhile!(qnet,edge,node)
        println("edge $(edge.number) has length $(edge.length) after lumping all internal edges into it, and set to qnet.t1")
        qnet.t1 = edge.length
    end
end

# function to eliminate hybridizations in a quartet network
# first step to later calculate expCF
# input: quartet network
function eliminateHybridization!(qnet::QuartetNetwork)
    if(qnet.which != -1)
        if(qnet.numHybrids > 0)
            for(n in qnet.hybrid)
                eliminateHybridization!(qnet,n)
            end
        end
        if(qnet.which == 1)
            internalLength!(qnet)
        end
    else
        error("qnet which has to be updated by now to 1 or 2, and it is $(qnet.which)")
    end
end

# ------------------------- update qnet.formula

# function to identify to which of the 4 leaves
# two taxa correspond. this is to identify later
# if the two taxa correspond to major/minor CF
# input: qnet, taxon1,taxon2
# returns leaf for taxon1, leaf for taxon2 (i.e. 12)
# warning: assumes that the numbers for the taxon in the output.csv table are the names
function whichLeaves(qnet::QuartetNetwork, taxon1::ASCIIString, taxon2::ASCIIString, leaf1::Node, leaf2::Node, leaf3::Node, leaf4::Node)
    warn("assume the numbers for the taxon read from the observed CF table match the numbers given to the taxon when creating the object network")
    if(taxon1 == qnet.names[leaf1.number])
        if(taxon2 == qnet.names[leaf2.number])
            return 1,2
        elseif(taxon2 == qnet.names[leaf3.number])
            return 1,3
        elseif(taxon2 == qnet.names[leaf4.number])
            return 1,4
        else
            error("taxon2 $(taxon2) is not one of the three remaining leaves: $(leaf2.number), $(leaf3.number), $(leaf4.number)")
        end
    elseif(taxon1 == qnet.names[leaf2.number])
        if(taxon2 == qnet.names[leaf1.number])
            return 2,1
        elseif(taxon2 == qnet.names[leaf3.number])
            return 2,3
        elseif(taxon2 == qnet.names[leaf4.number])
            return 2,4
        else
            error("taxon2 $(taxon2) is not one of the three remaining leaves: $(leaf1.number), $(leaf3.number), $(leaf4.number)")
        end
    elseif(taxon1 == qnet.names[leaf3.number])
        if(taxon2 == qnet.names[leaf1.number])
            return 3,1
        elseif(taxon2 == qnet.names[leaf2.number])
            return 3,2
        elseif(taxon2 == qnet.names[leaf4.number])
            return 3,4
        else
            error("taxon2 $(taxon2) is not one of the three remaining leaves: $(leaf2.number), $(leaf1.number), $(leaf4.number)")
        end
    elseif(taxon1 == qnet.names[leaf4.number])
        if(taxon2 == qnet.names[leaf2.number])
            return 4,2
        elseif(taxon2 == qnet.names[leaf3.number])
            return 4,3
        elseif(taxon2 == qnet.names[leaf1.number])
            return 4,1
        else
            error("taxon2 $(taxon2) is not one of the three remaining leaves: $(leaf2.number), $(leaf3.number), $(leaf4.number)")
        end
    else
        error("taxon1: $(taxon1) is not one of the 4 leaves in quartet: $(leaf1.number), $(leaf2.number), $(leaf3.number), $(leaf4.number)")
    end
end

# function to update the attribute split in qnet
# qnet.leaf is a vector [x,y,z,w] and qnet.split is a vector
# [1,1,2,2] that says in which side of the split is each leaf
# warning: it needs to be run after eliminating hybridization and uniting
# internal edge
function updateSplit!(qnet::QuartetNetwork)
    if(qnet.which == 1)
        if(qnet.split == [-1,-1,-1,-1])
            qnet.split = [2,2,2,2]
            middle = qnet.node[getIndex(true,[size(n.edge,1) == 3 for n in qnet.node])]
            leaf1 = middle.edge[getIndex(true,[getOtherNode(e,middle).leaf for e in middle.edge])]
            leaf2 = middle.edge[getIndex(true,[(getOtherNode(e,middle).leaf && !isequal(leaf1,e)) for e in middle.edge])]
            leaf1 = getOtherNode(leaf1,middle) #leaf1 was edge, now it is node
            leaf2 = getOtherNode(leaf2,middle)
            ind1 = getIndex(leaf1,qnet.leaf)
            ind2 = getIndex(leaf2,qnet.leaf)
            qnet.split[ind1] = 1
            qnet.split[ind2] = 1
        end
    elseif(qnet.which == -1)
        error("cannot update split in quartet network if it has not been identified, and eliminated hybridizations")
    end
end

# function to know which formula (minor/major) to compute
# for qnet.expCF[1,2,3] depending on the order of taxa in
# qnet.quartet
# warning: needs qnet.split updated already
function updateFormula!(qnet::QuartetNetwork)
    if(qnet.which == 1)
        if(qnet.formula == [-1,-1,-1])
            if(qnet.split != [-1,-1,-1,-1])
                qnet.formula = [2,2,2]
                for(i in 2:4)
                    if(size(qnet.leaf,1) != 4)
                        error("strange quartet with $(size(qnet.leaf,1)) leaves instead of 4")
                    end
                    tx1,tx2 = whichLeaves(qnet,qnet.quartetTaxon[1],qnet.quartetTaxon[i], qnet.leaf[1], qnet.leaf[2], qnet.leaf[3], qnet.leaf[4]) # index of leaf in qnet.leaf
                    if(qnet.split[tx1] == qnet.split[tx2])
                        qnet.formula[i-1] = 1
                        break
                    end
                end
            else
                error("cannot update qnet.formula if qnet.split is not updated: $(qnet.split)")
            end
        end
    elseif(qnet.which != 2)
        error("qnet.which should be updated already to 1 or 2, and it is $(qnet.which)")
    end
end


# --------------- calculate exp CF ----------------------


# function to calculate expCF for a quartet network
# warning: needs qnet.formula and qnet.t1 already updated
function calculateExpCF!(qnet::QuartetNetwork)
    if(qnet.which == 1)
        if(qnet.formula != [-1,-1,-1,-1] && qnet.t1 != -1)
            for(i in 1:3)
                qnet.expCF[i] = qnet.formula[i] == 1 ? 1-2/3*exp(-qnet.t1) : 1/3*exp(-qnet.t1)
            end
        else
            error("quartet network needs to have updated formula and t1 before computing the expCF")
        end
    elseif(qnet.which == 2)
        if(qnet.numHybrids == 1)
            if(qnet.hybrid[1].typeHyb == 5)
                quartetType5!(qnet,qnet.hybrid[1])
            else
                error("strange quartet network type $(qnet.which) with one hybrid node $(qnet.hybrid[1].number) but it is not type 5, it is type $(qnet.hybrid[1].typeHyb)")
            end
        else
            error("quartet network with type $(qnet.which) but with $(qnet.numHybrids) hybrid nodes. all hybridizations type 1 (not identifiable) have been eliminated already, so there should only be one hybrid.")
        end
    end
end

# function to compute all the process of calculating the expCF
# for a given qnet
function calculateExpCFAll!(qnet::QuartetNetwork)
    identifyQuartet!(qnet)
    eliminateHybridization!(qnet)
    updateSplit!(qnet)
    updateFormula!(qnet)
    calculateExpCF!(qnet)
end

# function to calculate expCF for all the quartets in data
# after extractQuartet(net,data) that updates quartet.qnet
# warning: only updates expCF for quartet.changed=true
function calculateExpCFAll!(data::DataCF)
    !all([q.qnet.numTaxa != 0 for q in data.quartet]) ? error("qnet in quartets on data are not correctly updated with extractQuartet") : nothing
    for(q in data.quartet)
        if(q.changed)
            qnet = deepcopy(q.qnet);
            calculateExpCFAll!(qnet);
            q.qnet.expCF = qnet.expCF
        end
    end
end


# ---------------------------- Pseudolik for a quartet -------------------------

# function to calculate the log pseudolikelihood function for a single
# quartet
# sum_i=1,2,3 (obsCF_i)*log(expCF_i)
# warning: assumes that quartet.qnet is already updated with extractQuartet and
#          calculateExpCF
function logPseudoLik(quartet::Quartet)
    if(sum(quartet.qnet.expCF) != 0.0)
        suma = 0
        for(i in 1:3)
            suma += quartet.obsCF[i]*log(quartet.qnet.expCF[i])
        end
        return suma
    else
        error("expCF not updated for quartet $(quartet.number)")
    end
end

# function to calculate the log pseudolikelihood function for array of
# quartets
# sum_q [sum_i=1,2,3 (obsCF_i)*log(expCF_i)]
# warning: assumes that quartet.qnet is already updated with extractQuartet and
#          calculateExpCF for all quartets
function logPseudoLik(quartet::Array{Quartet,1})
    suma = 0
    for(q in quartet)
        suma += logPseudoLik(q)
    end
    return suma
end

logPseudoLik(d::DataCF) = logPseudoLik(d.quartet)

