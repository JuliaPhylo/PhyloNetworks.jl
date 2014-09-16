# functions written in classes.jl and moved here after tested
# for pseudolikelihood implementation (Stage2)
# Claudia August 2014
#
# in julia: include("functions.jl")

# tests of functions in examples_classes.jl outside git_laptop

#------------- EDGE functions --------------------#

# warning: node needs to be defined as hybrid before adding to a hybrid edge
#          First, an edge is defined as hybrid, and then the nodes are added to it
function setNode!(edge::Edge, node::Node)
  if(size(edge.node,1)  ==  2)
    error("vector of nodes already has 2 values");
  else
    push!(edge.node,node);
    n = size(edge.node,1);
    if(edge.hybrid)
	if(n == 1)
            if(node.hybrid)
               edge.isChild1 = true;
            else
               edge.isChild1 = false;
	    end
	else
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
    end
  end
end

# warning: node needs to be defined as hybrid before adding to a hybrid edge
#          First, an edge is defined as hybrid, and then the nodes are added to it
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
end


# -------------- NODE -------------------------#

# warning: not really used, read types.jl
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


# find the index in net.node for a node with given number
# warning: assumes number uniquely determined
function getIndexNumNode(number::Int64,net::HybridNetwork)
    if(sum([net.node[i].number == number?1:0 for i = 1:size(net.node,1)]) == 0)
        error("node number $(number) not in network")
    else
        getIndex(true,[net.node[i].number == number for i = 1:size(net.node,1)])
    end
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
    println("Edge\tNode1\tNode2\tInCycle\tcontainRoot\tistIdentitiable\tLength")
    for e in net.edge
        println("$(e.number)\t$(e.node[1].number)\t$(e.node[2].number)\t$(e.inCycle)\t$(e.containRoot)\t$(e.istIdentifiable)\t$(e.length)")
    end
end

# print for every node, inCycle and edges
function printNodes(net::HybridNetwork)
    println("Node\tIn Cycle\tEdges numbers")
    for n in net.node
        print("$(n.number)\t$(n.inCycle)")
        for e in n.edge
            print("\t$(e.number)")
        end
        print("\n")
    end
end

# find the hybrid edges for a given hybrid node
function hybridEdges(node::Node)
   if(node.hybrid)
      hybedges = Edge[];
      for(e in node.edge)
        if(e.hybrid)
	   push!(hybedges,e);
	end
      end
      if(size(hybedges,1) == 2)
        return hybedges
      else
        error("node has more or less than 2 hybrid edges");
      end
   else
      error("node is not hybrid");
   end
end


# function to update inCycle (with priority queue) after becoming part of a network
# based on program 3 CS367 with priority queue
# expected to be much faster than the other two udpateInCycle (queue and recursive)
# input: hybrid node around which we want to update inCycle
# needs module "Base.Collections"
# returns error if no cycle in the network
# returns tuple: flag, array of edges changed, array of nodes changed
#                      without including the minor hybrid edge and hybrid node
#   flag: false if cycle intersects existing cycle or number of nodes in cycle < 3
#         (there is the possibility of returning edges in intersection: path)
#         true if cycle does not intersect existing cycle
# calculates also the number of nodes in the cycle and put as hybrid node attribute "k"
# check: visited is an aux array, do we prefer as attribute in Node? other data structure more efficient?
# check: prev is attribute in Node, do we prefer as aux array?
# warning: it is not checking if hybrid node or minor hybrid edge
#          were already part of a cycle (inCycle!= -1)
#          But it is checking so for the other edges in cycle
# warning: it needs extra things: visited array, prev attribute
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
        found = false;
        visited = [false for i = 1:size(net.node,1)];
        enqueue!(queue,node,dist);
        while(!found)
            if(isempty(queue))
                error("no cycle in network")
            else
                curr = dequeue!(queue);
                if(isequal(curr,last))
                    found = true;
                    push!(path,curr);
                else
                    if(!visited[getIndex(curr,net)])
                        visited[getIndex(curr,net)] = true;
                        if(isequal(curr,start))
                            for(e in curr.edge)
                                if(e.isMajor && e.hybrid)
                                    other = getOtherNode(e,curr);
                                    other.prev = curr;
                                    dist = dist+1;
                                    enqueue!(queue,other,dist);
                                end
                            end
                        else
                            for(e in curr.edge)
                                other = getOtherNode(e,curr);
                                if(!other.leaf && !visited[getIndex(other,net)])
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
        edges_changed = Edge[];
        nodes_changed = Node[];
        while(!isequal(curr, start))
            if(curr.inCycle!= -1)
                Data.Structures.enqueue!(path,curr);
                curr = curr.prev;
            else
                curr.inCycle = start.number;
                push!(nodes_changed, curr);
                node.k  =  node.k + 1;
                edge = getConnectingEdge(curr,curr.prev);
                edge.inCycle = start.number;
                push!(edges_changed, edge);
                curr = curr.prev;
            end
        end
        if(!isempty(path) || node.k<3)
            #error("new cycle intersects existing cycle")
            return false, edges_changed, nodes_changed
        else
            return true, edges_changed, nodes_changed
        end
    else
        error("node is not hybrid")
    end
end


# aux function to traverse the network for updateContainRoot
# it changes the containRoot argument to false
# of all the edges visited
# changed to recursive after Cecile's idea
function traverseContainRoot(node::Node, edge::Edge, edges_changed::Array{Edge,1})
    if(!node.leaf)
        for(e in node.edge)
            if(!isequal(edge,e))
                other = getOtherNode(e,node);
                if(e.containRoot) # only considered changed those that were true
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
        edges_changed = Edge[];
        for (e in node.edge)
            if(!e.hybrid)
                other = getOtherNode(e,node);
                e.containRoot = false;
                push!(edges_changed,e);
                traverseContainRoot(other,e, edges_changed);
            end
        end
        if(all([!e.containRoot for e in net.edge]))
            return false,edges_changed
        else
            return true,edges_changed
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
# returns false if the network is not identifiable for k=3
#         in this case, we need to redo updateContainRoot only
# warning: needs to have updateInCycle already done as it needs
#          inCycle, and k attributes
# check: assume any tree node that has hybrid Edge has only one tree edge in cycle (true?)
# fixit: we still need to be certain of generalization n> = 5 of bad diamond case,
#        here assumed: 1,1,1,> = 2 (see ipad figure)
function updateGammaz!(net::HybridNetwork, node::Node)
    if(node.hybrid)
        if(node.k == 4) # could be bad diamond
            edge_maj = nothing;
            edge_min = nothing;
            edge_maj2 = nothing;
            edge_min2 = nothing;
            tree_edge1 = nothing;
            tree_edge2 = nothing;
            tree_edge3 = nothing;
            for(e in node.edge)
                if(isa(edge_maj,Nothing))
	            edge_maj = (e.hybrid && e.isMajor)? e : nothing;
                end
                if(isa(edge_min,Nothing))
                    edge_min = (e.hybrid && !e.isMajor)? e : nothing;
                end
                if(isa(tree_edge2,Nothing))
                    tree_edge2 = (!e.hybrid && e.inCycle  ==  -1)? e : nothing;
                end
            end
            other_maj = getOtherNode(edge_maj,node);
            other_min = getOtherNode(edge_min,node);
            for(e in other_min.edge)
                if(isa(edge_min2,Nothing))
                    edge_min2 = (!e.hybrid && e.inCycle !=  -1)? e : nothing;
                end
                if(isa(tree_edge3,Nothing))
                    tree_edge3 = (!e.hybrid && e.inCycle  ==  -1)? e : nothing;
                end
            end
            for(e in other_maj.edge)
                if(isa(edge_maj2,Nothing))
                    edge_maj2 = (!e.hybrid && e.inCycle !=  -1)? e : nothing;
                end
                if(isa(tree_edge1,Nothing))
                    tree_edge1 = (!e.hybrid && e.inCycle  ==  -1)? e : nothing;
                end
            end
            if(isequal(getOtherNode(edge_min2,other_min),getOtherNode(edge_maj2,other_maj)) && getOtherNode(tree_edge1,other_maj).leaf && getOtherNode(tree_edge2,node).leaf && getOtherNode(tree_edge3,other_min).leaf) # if this, you have bad diamond
                other_min.gammaz = edge_min.gamma*edge_min2.z;
                other_maj.gammaz = edge_maj.gamma*edge_maj2.z;
                edge_min2.istIdentifiable = false;
                edge_maj2.istIdentifiable = false;
                node.isBadDiamond = true;
            else
                node.isBadDiamond = false;
            end
        elseif(node.k == 3) # could be bad triangle, non identifiable or set to zero cases
            edge_maj = nothing;
            edge_min = nothing;
            edge_maj2 = nothing;
            edge_min2 = nothing;
            for(e in node.edge)
                if(isa(edge_maj,Nothing))
	            edge_maj = (e.hybrid && e.isMajor)? e:nothing;
                end
                if(isa(edge_min,Nothing))
                    edge_min = (e.hybrid && !e.isMajor)? e:nothing;
                end
            end
            other_maj = getOtherNode(edge_maj,node);
            other_min = getOtherNode(edge_min,node);
            tree_edge1 = nothing;
            tree_edge2 = nothing;
            tree_edge3 = nothing;
            tree_edge_incycle = nothing;
            for(e in other_min.edge)
                if(isa(tree_edge1,Nothing))
	            tree_edge1 = (!e.hybrid && e.inCycle == -1)? e:nothing;
                end
                if(isa(tree_edge_incycle,Nothing))
	            tree_edge_incycle = (!e.hybrid && e.inCycle!= -1)? e:nothing;
                end
            end
            for(e in other_maj.edge)
                if(isa(tree_edge2,Nothing))
	            tree_edge2 = (!e.hybrid && e.inCycle == -1)? e:nothing;
                end
            end
            for(e in node.edge)
                if(isa(tree_edge3,Nothing))
	            tree_edge3 = (!e.hybrid && e.inCycle == -1)? e:nothing;
                end
            end
            isLeaf1 = getOtherNode(tree_edge1,other_min);
            isLeaf2 = getOtherNode(tree_edge2,other_min);
            isLeaf3 = getOtherNode(tree_edge3,other_min);
            if(isLeaf1.leaf && !isLeaf2.leaf && !isLeaf3.leaf) # if this, you have bad triangle
                node.gammaz = edge_maj.gamma^2*edge_maj.z+edge_min.gamma^2*tree_edge_incycle.z;
                tree_edge3.istIdentifiable = false;
                tree_edge_incycle.istIdentifiable = false;
                edge_maj.istIdentifiable = false;
                tree_edge2.istIdentifiable = false;
                edge_min.istIdentifiable = false;
                node.isBadTriangle = true;
            elseif(!isLeaf1.leaf && !isLeaf2.leaf && isLeaf3.leaf) # set to zero case1
                edge_maj.length = edge_maj.length + tree_edge2.length;
                tree_edge2.length = 0;
                tree_edge_incycle.length = tree_edge_incycle.length + tree_edge1.length;
                tree_edge1.length = 0;
                tree_edge1.istIdentifiable = false;
                tree_edge2.istIdentifiable = false;
                edge_min.istIdentifiable = false;
                node.isBadTriangle = false;
            elseif(!isLeaf1.leaf && isLeaf2.leaf && !isLeaf3.leaf) # set to zero case2
                tree_edge_incycle.length = tree_edge_incycle.length + tree_edge1.length;
                tree_edge1.length = 0;
                tree_edge1.istIdentifiable = false;
                node.isBadTriangle = false;
            elseif(sum([isLeaf1?1:0, isLeaf2?1:0, isLeaf3?1:0]) == 2) # non identifiable network
                return false
            end
        end
        return true
    else
        error("node is not hybrid")
    end
end



# setLength using updateGammaz
# calls updateGammaz, updateGamma2z when needed
# warning: network needs to be updated to have isBadDiamond/Triangle set
function setLength!(edge::Edge, new_length::Float64,net::HybridNetwork)
  if(new_length<0)
      error("length has to be nonnegative");
  else
      edge.length = new_length;
      edge.y = exp(-new_length);
      edge.z = 1 - edge.y;
      if(edge.inCycle!= -1)
          node = net.node[getIndexNumNode(edge.inCycle,net)]; #hybrid node for this edge
          if(edge.hybrid && edge.isMajor && node.isBadTriangle)
              updateGamma2z!(net,node);
          elseif (!edge.hybrid && edge.inCycle!= -1 && !all([!edge.node[i].hasHybEdge for i = 1:size(edge.node,1)])) # tree edge with node.hasHybEdge true
              if(node.isBadDiamond)
                  updateGammaz!(net,node);
              elseif(node.isBadTriangle)
                  updateGamma2z!(net,node);
              end
          end
      end
  end
end


# setGamma using updateGammaz
# calls updateGammaz when needed, and updateGamma2z always
# warning: network needs to be updated to have isBadDiamond/Triangle set
# check: bad triangle status depends on value of gamma!
function setGamma!(edge::Edge, new_gamma::Float64, net::HybridNetwork)
 if(edge.hybrid)
	if(0 < new_gamma < 1)
	     edge.gamma = new_gamma;
	     new_gamma<0.5 ? edge.isMajor = false : edge.isMajor = true;
	     edge.isChild1 ? ind = 1 : ind = 2 ; # hybrid edge pointing at node 1 or 2
             hybedges = hybridEdges(edge.node[ind]); # hybrid edges for hybrid node
	     if(!isempty(hybedges) && !isa(hybedges,Nothing)) # change 1-gamma in other hybrid edge
	          isequal(hybedges[1],edge) ? other = 2 : other = 1;
		  hybedges[other].gamma = 1.0-new_gamma;
       	          new_gamma<0.5 ? hybedges[other].isMajor = true : hybedges[other].isMajor = false;
	     end
            if(edge.node[ind].isBadDiamond)
                updateGammaz!(net, edge.node[ind]);
            end
            #if(edge.node[ind].isBadTriangle), check: bad triangle status depends on gamma, it can become bad triangle
                updateGamma2z!(net, edge.node[ind]);
            #end
	else
	     error("gamma has to be between 0 and 1");
        end
  else
	error("cannot change gamma in a tree edge");
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
