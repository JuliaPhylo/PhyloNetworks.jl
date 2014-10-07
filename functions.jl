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
    if(size(edge.node) == 1)
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
	            error("hybrid edge has no hybrid nodes");
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
   edge.hybrid ? node.hasHybEdge = true : node.hasHybEdge = false;
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

function getIndex(bool::Bool, array::Array{Any,1})
    i = 1;
    while(i<= size(array,1) && !isequal(bool,array[i]))
        i = i+1;
    end
    i>size(array,1)?error("$(bool) not in array"):return i;
end

# find the index of an edge in node.edge
function getIndexEdge(edge::Edge,node::Node)
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

# find the index in net.node for a node with given number
# warning: assumes number uniquely determined
function getIndexNumNode(number::Int64,net::HybridNetwork)
    if(sum([net.node[i].number == number?1:0 for i = 1:size(net.node,1)]) == 0)
        error("node number $(number) not in network")
    else
        getIndex(true,[net.node[i].number == number for i = 1:size(net.node,1)])
    end
end

# find the index in net.edge for a edge with given number
# warning: assumes number uniquely determined
function getIndexNumEdge(number::Int64,net::HybridNetwork)
    if(sum([net.edge[i].number == number?1:0 for i = 1:size(net.edge,1)]) == 0)
        error("edge number $(number) not in network")
    else
        getIndex(true,[net.edge[i].number == number for i = 1:size(net.edge,1)])
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

# function to delete a Node in net.node and
# update numNodes and numTaxa
function deleteNode!(net::HybridNetwork, n::Node)
    try
        index = getIndex(n,net);
    catch
        error("Node not in network");
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


# search the hybrid node(s) in network: returns the hybrid node(s)
# if more than one hybrid, return an array of nodes
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
        return net.node[k]
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
    println("Edge\tNode1\tNode2\tInCycle\tcontainRoot\tistIdentitiable\tLength\tisHybrid")
    for e in net.edge
        println("$(e.number)\t$(e.node[1].number)\t$(e.node[2].number)\t$(e.inCycle)\t$(e.containRoot)\t\t$(e.istIdentifiable)\t\t$(e.length)\t$(e.hybrid)")
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
           error("node has more or less than 3 edges");
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
           error("node has more or less than 3 edges");
       end
   else
       error("node is not hybrid nor tree with hybrid edges");
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
    if(edge.hybrid && !node.hybrid)
        node.hasHybEdge = false;
    end
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
# see ipad notes: k = 0 (nonidentifiable), k = 1 (nonidentifiable, set to zero, bad triangle)
# k = 2 (bad diamond)
# input: hybrid node around which to check (can come from searchHybridNode)
# updates gammaz with whatever edge lengths are originally in the network
# returns flag, array of edges changed (istIdentifiable)
#               without the major/minor hybrid edges
#               that would be the eliminated ones
#         false if the network is not identifiable for k=3
# warning: needs to have updateInCycle already done as it needs
#          inCycle, and k attributes
# check: assume any tree node that has hybrid Edge has only one tree edge in cycle (true?)
# fixit: we still need to be certain of generalization n> = 5 of bad diamond case,
#        here assumed: 1,1,1,> = 2 (see ipad figure)
function updateGammaz!(net::HybridNetwork, node::Node)
    if(node.hybrid)
        net.edges_changed = Edge[];
        edge_maj, edge_min, tree_edge2 = hybridEdges(node);
        other_maj = getOtherNode(edge_maj,node);
        other_min = getOtherNode(edge_min,node);
        edge_min.istIdentifiable = false;
        push!(net.edges_changed, edge_min);
        if(node.k == 4) # could be bad diamond
            edgebla,edge_min2,tree_edge3 = hybridEdges(other_min);
            edgebla,edge_maj2,tree_edge1 = hybridEdges(other_maj);
            if(isequal(getOtherNode(edge_min2,other_min),getOtherNode(edge_maj2,other_maj)) && getOtherNode(tree_edge1,other_maj).leaf && getOtherNode(tree_edge2,node).leaf && getOtherNode(tree_edge3,other_min).leaf) # if this, you have bad diamond
                other_min.gammaz = edge_min.gamma*edge_min2.z;
                other_maj.gammaz = edge_maj.gamma*edge_maj2.z;
                edge_min2.istIdentifiable = false;
                edge_maj2.istIdentifiable = false;
                edge_maj.istIdentifiable = false;
                push!(net.edges_changed,edge_min2);
                push!(net.edges_changed,edge_maj2);
                push!(net.edges_changed,edge_maj);
                node.isBadDiamond = true;
            else
                node.isBadDiamond = false;
            end
        elseif(node.k == 3) # could be bad triangle, non identifiable or set to zero cases
            edgebla,tree_edge_incycle,tree_edge1 = hybridEdges(other_min);
            edgebla,edgebla,tree_edge3 = hybridEdges(other_maj);
            isLeaf1 = getOtherNode(tree_edge1,other_min);
            isLeaf2 = getOtherNode(tree_edge2,node);
            isLeaf3 = getOtherNode(tree_edge3,other_maj);
            if(isLeaf1.leaf && !isLeaf2.leaf && !isLeaf3.leaf) # if this, you have bad triangle
                node.gammaz = edge_maj.gamma*edge_maj.gamma*edge_maj.z+edge_min.gamma*edge_min.gamma*tree_edge_incycle.z;
                other_min.gammaz = edge_min.gamma*tree_edge_incycle.z;
                tree_edge_incycle.istIdentifiable = false;
                edge_maj.istIdentifiable = false;
                tree_edge2.istIdentifiable = false;
                node.isBadTriangle = true;
                push!(net.edges_changed,tree_edge_incycle);
                push!(net.edges_changed,tree_edge2);
                push!(net.edges_changed,edge_maj);
                edge_maj.length = edge_maj.length + tree_edge2.length;
                tree_edge2.length = 0.0;
                edge_min.length = 0.0;
            elseif(!isLeaf1.leaf && !isLeaf2.leaf && isLeaf3.leaf) # set to zero case1
                edge_maj.length = edge_maj.length + tree_edge2.length;
                tree_edge2.length = 0.0;
                edge_min.length = 0.0;
                tree_edge_incycle.length = tree_edge_incycle.length + tree_edge1.length;
                tree_edge1.length = 0.0;
                tree_edge1.istIdentifiable = false;
                tree_edge2.istIdentifiable = false;
                push!(net.edges_changed,tree_edge2);
                push!(net.edges_changed,tree_edge1);
                node.isBadTriangle = false;
            elseif(!isLeaf1.leaf && isLeaf2.leaf && !isLeaf3.leaf) # set to zero case2
                tree_edge_incycle.length = tree_edge_incycle.length + tree_edge1.length;
                tree_edge1.length = 0.0;
                edge_min.length = 0.0;
                tree_edge1.istIdentifiable = false;
                push!(net.edges_changed,tree_edge1);
                node.isBadTriangle = false;
            elseif(sum([isLeaf1.leaf?1:0, isLeaf2.leaf?1:0, isLeaf3.leaf?1:0]) == 2) # non identifiable network
                return false, net.edges_changed
            end
        end
        if(node.k > 3 && !node.isBadDiamond)
            edgebla,tree_edge_incycle,tree_edge1 = hybridEdges(other_min);
            if(!tree_edge_incycle.istIdentifiable)
                tree_edge_incycle.istIdentifiable = true;
                push!(net.edges_changed,tree_edge_incycle);
            end
            if(getOtherNode(tree_edge2,node).leaf && edge_maj.istIdentifiable)
                edge_maj.istIdentifiable = false;
                push!(net.edges_changed,edge_maj);
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
        push!(net.edge,edge3);
        push!(net.edge,edge4);
        # create hybridization
        max_node = maximum([e.number for e in net.node]);
        max_edge = maximum([e.number for e in net.edge]);
        hybrid_edge = Edge(max_edge+1,0.0,true,gamma,false);
        push!(net.edge,hybrid_edge);
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
        push!(net.node,hybrid_node);
        push!(net.node,tree_node);
        return hybrid_node
    else
        error("gamma must be between 0 and 1")
    end
end

# aux function to addHybridization
# it chooses the edges in the network and the gamma value
# warning: chooses edge1, edge2, gamma randomly, but
#          we could do better later
# check: gamma is uniform(0,1/2) to avoid big gammas
function chooseEdgesGamma(net::HybridNetwork)
    index1 = 1;
    index2 = 1;
    while(index1 == index2 || index1 == 0 || index2 == 0 || index1 > size(net.edge,1) || index2 > size(net.edge,1) || net.edge[index1].inCycle != -1 || net.edge[index2].inCycle != -1)
        index1 = iround(rand()*size(net.edge,1));
        index2 = iround(rand()*size(net.edge,1));
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
# input: new added hybrid, network
# returns: success (bool), hybrid, flag, nocycle, flag2, flag3
function updateAllNewHybrid!(hybrid::Node,net::HybridNetwork)
    flag, nocycle, edgesInCycle, nodesInCycle = updateInCycle!(net,hybrid);
    if(nocycle)
        return false, hybrid, flag, nocycle, flag2, flag3
    else
        if(flag)
            updateMajorHybrid!(net,hybrid);
            flag2, edgesGammaz = updateGammaz!(net,hybrid);
            if(flag2)
                flag3, edgesRoot = updateContainRoot!(net,hybrid);
                if(flag3)
                    return true, hybrid, flag, nocycle, flag2, flag3
                else
                    undoContainRoot!(edgesRoot);
                    undoistIdentifiable!(edgesGammaz);
                    undoInCycle!(edgesInCycle, nodesInCycle);
                    return false, hybrid, flag, nocycle, flag2, flag3
                end
            else
                undoistIdentifiable!(edgesGammaz);
                undoInCycle!(edgesInCycle, nodesInCycle);
                return false, hybrid, flag, nocycle, flag2, flag3
            end
        else
            undoInCycle!(edgesInCycle, nodesInCycle);
            return false, hybrid, flag, nocycle, flag2, flag3
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
    updateAllNewHybrid!(hybrid,net)
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


# function to identify the edges that need to be changed to identifiable
# if the hybridization is deleted
# it also does the change the istIdentifiable attribute of such edges
# see ipad notes: only two cases, bad diamond/triangle
# input: hybrid node which will later be deleted
# returns true if changed something, false o.w.
# warning: the hybridization needs to be fully updated already:
#          inCycle, gammaz
# warning: this is the only of the "identify" functions that
#          does change the attributes
function identifyGammaz(net::HybridNetwork, node::Node)
    if(node.hybrid)
        edge_maj, edge_min, tree_edge2 = hybridEdges(node);
        other_maj = getOtherNode(edge_maj,node);
        if(node.isBadDiamond)
            edgebla,edge_maj2,tree_edge1 = hybridEdges(other_maj);
            edge_maj2.istIdentifiable = true;
            return true #, edge_maj2
        elseif(node.isBadTriangle)
            edge_maj.istIdentifiable = true;
            return true #, edge_maj
        end
        return false #, nothing
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
function deleteHybridizationUpdate!(net::HybridNetwork, hybrid::Node)
    nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,hybrid);
    if(!nocycle)
        edgesRoot = identifyContainRoot(net,hybrid);
        identifyGammaz(net,hybrid);
        undoInCycle!(edgesInCycle, nodesInCycle);
        undoContainRoot!(edgesRoot);
        deleteHybrid!(hybrid,net)
    else
        error("the hybrid does not create a cycle")
    end
end

# function to delete a hybridization event
# input: hybrid node and network
# warning: it is meant after undoing the effect of the
#          hybridization in deleteHybridizationUpdate!
#          by itself, it leaves things as if
function deleteHybrid!(node::Node,net::HybridNetwork)
    if(node.hybrid)
        hybedge1,hybedge2,treeedge1 = hybridEdges(node);
        other1 = getOtherNode(hybedge1,node);
        other2 = getOtherNode(hybedge2,node);
        other3 =  getOtherNode(treeedge1,node);
        if(hybedge1.number > treeedge1.number)
            treeedge1.length = treeedge1.length + hybedge1.length;
            removeNode!(node,treeedge1);
            setNode!(treeedge1,other1);
            setEdge!(other1,treeedge1);
            removeEdge!(other1, hybedge1);
            index = getIndex(hybedge1,net);
            deleteat!(net.edge,index);
            treeedge1.containRoot = (!treeedge1.containRoot || !hybedge1.containRoot) ? false : true
        else
            hybedge1.hybrid = false;
            hybedge1.gamma = 1.0;
            hybedge1.isMajor = true;
            other1.hasHybEdge = false;
            hybedge1.length = hybedge1.length + treeedge1.length;
            removeNode!(node,hybedge1);
            setNode!(hybedge1,other3);
            setEdge!(other3,hybedge1);
            removeEdge!(other3,treeedge1);
            index = getIndex(treeedge1,net);
            deleteat!(net.edge,index);
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
        if(treeedge1.number > treeedge2.number)
            treeedge2.length = treeedge2.length + treeedge1.length;
            removeNode!(other2,treeedge2);
            setNode!(treeedge2,treenode1);
            setEdge!(treenode1,treeedge2);
            removeEdge!(treenode1,treeedge1);
            index = getIndex(treeedge1,net);
            deleteat!(net.edge,index);
            treeedge2.containRoot = (!treeedge1.containRoot || !treeedge2.containRoot) ? false : true
        else
            treeedge1.length = treeedge2.length + treeedge1.length;
            removeNode!(other2,treeedge1);
            setNode!(treeedge1,treenode2);
            setEdge!(treenode2,treeedge1);
            removeEdge!(treenode2,treeedge2);
            index = getIndex(treeedge2,net);
            deleteat!(net.edge,index);
            treeedge1.containRoot = (!treeedge1.containRoot || !treeedge2.containRoot) ? false : true
        end
        index = getIndex(node,net);
        deleteat!(net.node,index);
        index = getIndex(other2,net);
        deleteat!(net.node,index);
        index = getIndex(hybedge2,net);
        deleteat!(net.edge,index);
    else
        error("node has to be hybrid")
    end
end

# ------------------------------- Read topology ------------------------------------------


# aux function to advance stream in readSubtree
# warning: s IOStream needs to be a file, not a stream converted
#          from string
function advance!(s::IOStream, c::Char, numLeft::Int64)
    c = Base.peekchar(s)
    if(Base.eof(s))
        error("Tree ends prematurely while reading subtree after left parenthesis $(numLeft).")
    end
    return read(s,Char)
end


# aux function to read all digits of taxon name
function readNum(s::IOStream, c::Char)
    if(isdigit(c))
        num = read(s,Char)
        c = Base.peekchar(s)
        while(isdigit(c))
            d = read(s,Char)
            num = string(num,d)
            c = Base.peekchar(s);
        end
        return int(num)
    else
        error("Expected int digit but received $(c)");
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
        error("Expected float digit after : but received $(c)");
    end
end


# aux function to read subtree
# warning: s IOStream needs to be a file, not a stream converted
#          from string
# warning: assumes taxon numbers, not taxon names
# fixit: allow taxon names also
function readSubtree!(s::IOStream, parent::Node, numLeft::Int64,net::HybridNetwork)
    c = Base.peekchar(s)
    e = nothing;
    if(c =='(')
       numLeft += 1
       println(numLeft)
       n = Node(numLeft,false);
       c = read(s,Char)
       bl = readSubtree!(s,n,numLeft,net)
       c = advance!(s,c,numLeft)
       br = false;
       if(c == ',')
           br = readSubtree!(s,n,numLeft,net);
       else
           a = readall(s);
           error("Expected comma after left parenthesis $(numLeft) but read $(c). The remainder of line is $(a).")
       end
       c = advance!(s,c,numLeft)
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
        num = readNum(s,c)
        println("creating node $(num)")
        n = Node(num,true);
        pushNode!(net,n);
        e = Edge(net.numEdges+1,1.0);
        pushEdge!(net,e);
        setNode!(e,[n,parent]);
        setEdge!(n,e);
        setEdge!(parent,e);
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
        length = readFloat(s,c);
        setLength!(e,length);
    end
    return true
end

# function to read topology from parenthetical format
# input: file name
# warning: at the moment, assumes a tree
# warning: assumes taxon numbers, not taxon names
function readTopology(file::String)
    net = HybridNetwork()
    try
        s = open(file)
    catch
        error("Could not find or open $(file) file");
    end
    s = open(file)
    c = Base.peekchar(s)
    numLeft = 0
    if(c == '(')
       numLeft += 1
       println(numLeft)
       n = Node(numLeft,false);
       c = read(s,Char)
       b = false;
       while(c != ';')
           numLeft += 1;
           b |= readSubtree!(s,n,numLeft,net)
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



# ----------------------------------------------------------------------------------------

# setLength
# warning: does not allow to change edge length for istIdentifiable=false
function setLength!(edge::Edge, new_length::Float64)
  if(new_length<0)
      error("length has to be nonnegative");
  else
      if(edge.istIdentifiable)
          edge.length = new_length;
          edge.y = exp(-new_length);
          edge.z = 1 - edge.y;
      else
          error("edge length is not identifiable")
      end
  end
end


# setGamma
# warning: does not allow to change gamma in the bad diamond/triangle cases
# because gamma is not identifiable
# warning: we are not updating the status of isMajor with the value of gamma
#          see ipad notes. isMajor is the original tree edge
function setGamma!(edge::Edge, new_gamma::Float64)
 if(edge.hybrid)
	if(0 < new_gamma < 1)
            edge.isChild1 ? ind = 1 : ind = 2 ; # hybrid edge pointing at node 1 or 2
            if(edge.node[ind].hybrid)
                if(edge.node[ind].isBadDiamond || edge.node[ind].isBadTriangle)
                    error("bad diamond or triangle situation: gamma not identifiable")
                else
                    edge.gamma = new_gamma;
                end
            else
                error("hybrid edge not pointing at hybrid node")
            end
	else
	     error("gamma has to be between 0 and 1");
        end
  else
	error("cannot change gamma in a tree edge");
  end
end

