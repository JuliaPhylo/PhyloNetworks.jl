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

function getIndex(node::Node, net::HybridNetwork)
    i = 1;
    while(i<= size(net.node,1) && !isequal(node,net.node[i]))
        i = i+1;
    end
    i>size(net.node,1)?error("node not in network"):return i;
end

function getIndex(edge::Edge, net::HybridNetwork)
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


function getIndexNode(number::Int64,net::HybridNetwork)
    try
        getIndex(true,[number==n.number for n in net.node])
    catch
        error("node number not in net.node")
    end
    return getIndex(true,[number==n.number for n in net.node])
end

function getIndexEdge(number::Int64,net::HybridNetwork)
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
function getIndexHybrid(node::Node, net::HybridNetwork)
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


# function to push a Node in net.node and
# update numNodes and numTaxa
function pushNode!(net::HybridNetwork, n::Node)
    push!(net.node,n);
    net.numNodes += 1;
    net.numTaxa += n.leaf ? 1 : 0;
end

# function to push an Edge in net.edge and
# update numEdges
function pushEdge!(net::HybridNetwork, e::Edge)
    push!(net.edge,e);
    net.numEdges += 1;
end


# function to push a hybrid Node in net.hybrid and
# update numHybrids
function pushHybrid!(net::HybridNetwork, n::Node)
    if(n.hybrid)
        push!(net.hybrid,n);
        net.numHybrids += 1;
    else
        error("node $(n.number) is not hybrid, so cannot be pushed in net.hybrid")
    end
end


# function to delete a Node in net.node and
# update numNodes and numTaxa
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
end

# function to delete an Edge in net.edge and
# update numEdges
function deleteEdge!(net::HybridNetwork, e::Edge)
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
function removeHybrid!(net::HybridNetwork, n::Node)
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
function deleteIntNode!(net::HybridNetwork, n::Node)
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
function searchHybridNode(net::HybridNetwork)
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
function searchHybridEdge(net::HybridNetwork)
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
function printEdges(net::HybridNetwork)
    println("Edge\tNode1\tNode2\tInCycle\tcontainRoot\tistIdentitiable\tLength\tisHybrid\tGamma")
    for e in net.edge
        println("$(e.number)\t$(e.node[1].number)\t$(e.node[2].number)\t$(e.inCycle)\t$(e.containRoot)\t\t$(e.istIdentifiable)\t\t$(e.length)\t$(e.hybrid)\t$(e.gamma)")
    end
end

# print for every node, inCycle and edges
function printNodes(net::HybridNetwork)
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
       error("node $(node.number) is not hybrid $(node.hybrid) nor tree with hybrid edges (hasHybEdge) $(node.hasHybEdge)");
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
                warn("bad diamond II found, treatment not decided yet")
                node.isBadDiamondII = true;
            end
        elseif(node.k == 3) # could be bad triangle I, II or non identifiable cases
            edgebla,tree_edge_incycle,tree_edge1 = hybridEdges(other_min);
            edgebla,edgebla,tree_edge3 = hybridEdges(other_maj);
            isLeaf1 = getOtherNode(tree_edge1,other_min);
            isLeaf2 = getOtherNode(tree_edge2,node);
            isLeaf3 = getOtherNode(tree_edge3,other_maj);
            if(isLeaf1.leaf && !isLeaf2.leaf && !isLeaf3.leaf) # bad triangle I
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
                setLength!(edge_maj,edge_maj.length+tree_edge2.length)
                setLength!(tree_edge2, 0.0);
                setLength!(edge_min, 0.0);
            elseif(!isLeaf1.leaf && !isLeaf2.leaf && isLeaf3.leaf) # bad triangle I
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
                setLength!(edge_min,edge_min.length + tree_edge2.length);
                setLength!(tree_edge2, 0.0);
                setLength!(edge_maj, 0.0);
            elseif(!isLeaf1.leaf && isLeaf2.leaf && !isLeaf3.leaf) # bad triangle II
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
                setLength!(tree_edge1,-log(1-other_min.gammaz))
            else
                error("bad triangle II in node $(node.number) but no gammaz updated correctly")
            end
            if(other_maj.gammaz != -1)
                setLength!(tree_edge3,-log(1-other_maj.gammaz))
            else
                error("bad triangle II in node $(node.number) but no gammaz updated correctly")
            end
            setLength!(tree_edge_incycle,0.0) #t11 and t10 already have the right length
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
            other_min.gammaz = -1.0
            other_maj.gammaz = -1.0
            tree_edge_incycle1.istIdentifiable = true;
            tree_edge_incycle2.istIdentifiable = true;
            edge_maj.istIdentifiable = true;
            edge_min.istIdentifiable = true;
        elseif(node.isBadDiamondII)
            warn("bad diamond II detected, still no treatment for bad diamond II")
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
        pushHybrid!(net,hybrid_node);
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
# fixit: update net.hybrid and numHybrids (push hybrid function?),

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
function identifyInCycle(net::HybridNetwork,node::Node)
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
            if(random)
                if(edges[1].gamma > 0.5 && edges[1].gamma != 1.0)
                    minor = rand() < edges[1].gamma ? false : true
                else
                    error("strange major hybrid edge $(edges[1].number) with gamma either less than 0.5 or equal to 1.0")
                end
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
#          by itself, it leaves things as if
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
            removeHybrid!(net,node);
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
            removeHybrid!(net,node);
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
        removeHybrid!(net,current)
        pushHybrid!(net,new)
        current.hybrid = false
        new.k = current.k
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
    return net
end

# ----------------------------------------------------------------------------------------

# setLength
# warning: allows to change edge length for istIdentifiable=false
#          but issues a warning
function setLength!(edge::Edge, new_length::Float64)
  if(new_length<0)
      error("length has to be nonnegative");
  else
      edge.length = new_length;
      edge.y = exp(-new_length);
      edge.z = 1 - edge.y;
      if(!edge.istIdentifiable)
          warn("set edge length that is not identifiable")
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
	     error("gamma has to be between 0 and 1");
        end
  else
	error("cannot change gamma in a tree edge");
  end
end

