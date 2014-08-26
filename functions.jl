# functions for pseudolikelihood implementation (Stage2)
# Claudia August 2014
#
# in julia: include("functions.jl")



#------------- EDGE functions --------------------#

#Node needs to be defined as hybrid before adding to a hybrid edge
#First, an edge is defined as hybrid, and then the nodes are added to it
function setNode!(edge::Edge, node::Node)
  if(size(edge.node,1) == 2)
    error("vector of nodes already has 2 values");
  else
    push!(edge.node,node);
    n=size(edge.node,1);
    if(edge.hybrid)
	if(n==1)
            if(node.hybrid)
               edge.isChild1=true;
            else
               edge.isChild1=false;
	    end
	else
	    if(node.hybrid)
               if(edge.node[1].hybrid)
                  error("hybrid edge has two hybrid nodes");
               else
                  edge.isChild1=false;
	       end
	    else
	       if(!edge.node[1].hybrid)
	          error("hybrid edge has no hybrid nodes");
	       else
	          edge.isChild1=true;
	       end
	    end
	end
    end
  end
end

function setNode!(edge::Edge,node::Array{Node,1})
  size(node,1) != 2 ?
  error("vector of nodes must have exactly 2 values") :
  edge.node=node;
  if(edge.hybrid)
      if(node[1].hybrid)
          edge.isChild1=true;
      else
          if(node[2].hybrid)
              edge.isChild1=false;
          else
              error("hybrid edge without hybrid node");
          end
      end
   end
end



# -------------- NODE -------------------------#

function setEdge!(node::Node,edge::Edge)
   push!(node.edge,edge);
   edge.hybrid ? node.hasHybEdge=true : node.hasHybEdge=false;
end

function getOtherNode(edge::Edge,node::Node)
  isequal(edge.node[1],node) ? edge.node[2] : edge.node[1]
end
# -------------- NETWORK ----------------------- #

function getIndex(node::Node, net::HybridNetwork)
    i=1;
    while(i<=size(net.node,1) && !isequal(node,net.node[i]))
        i=i+1;
    end
    i>size(net.node,1)?error("node not in network"):return i;
end

function getIndex(edge::Edge, net::HybridNetwork)
    i=1;
    while(i<=size(net.edge,1) && !isequal(edge,net.edge[i]))
        i=i+1;
    end
    i>size(net.edge,1)?error("edge not in network"):return i;
end

function getIndex(bool::Bool, array::Array{Bool,1})
    i=1;
    while(i<=size(array,1) && !isequal(bool,array[i]))
        i=i+1;
    end
    i>size(array,1)?error("$(bool) not in array"):return i;
end

function getIndex(bool::Bool, array::Array{Any,1})
    i=1;
    while(i<=size(array,1) && !isequal(bool,array[i]))
        i=i+1;
    end
    i>size(array,1)?error("$(bool) not in array"):return i;
end


# search the hybrid node in network. assumes only one hybrid node, or it returns the first one
# fixit: generalize to more hybrid nodes, need to change getIndex
function searchHybridNode(net::HybridNetwork)
     getIndex(true,[net.node[i].hybrid for i=1:size(net.node,1)])
end

# search the hybrid edge in network. assumes only one hybrid edge, or it returns the first one
# fixit: change to return the minor edge
function searchHybridEdge(net::HybridNetwork)
     getIndex(true,[net.edge[i].hybrid for i=1:size(net.edge,1)])
end

# function to update gammaz in a network. assumes: only one hybrid node
# needs to have inCycle attributes updated already
# check: assume any tree node that has hybrid Edge has only one tree edge in cycle (true?)
# fixit: add possibility of more hybrid nodes!
function updateGammaz!(net::HybridNetwork)
    index=searchHybridNode(net);
    node=net.node[index];
    if(node.hybrid)
        edge_maj=nothing;
        edge_min=nothing;
        edge_maj2=nothing;
        edge_min2=nothing;
        for(i=1:size(node.edge,1))
            if(isa(edge_maj,Nothing))
	       edge_maj=(node.edge[i].hybrid && node.edge[i].isMajor)? node.edge[i]:nothing;
            end
            if(isa(edge_min,Nothing))
               edge_min=(node.edge[i].hybrid && !node.edge[i].isMajor)? node.edge[i]:nothing;
            end
        end
        other_maj=getOtherNode(edge_maj,node);
        other_min=getOtherNode(edge_min,node);
        for(j=1:size(other_min.edge,1))
            if(isa(edge_min2,Nothing))
               edge_min2=(!other_min.edge[j].hybrid && other_min.edge[j].inCycle != -1)? other_min.edge[j]:nothing;
            end
        end
        for(j=1:size(other_maj.edge,1))
            if(isa(edge_maj2,Nothing))
              edge_maj2=(!other_maj.edge[j].hybrid && other_maj.edge[j].inCycle != -1)? other_maj.edge[j]:nothing;
            end
        end
        other_min.gammaz=edge_min.gamma*edge_min2.z;
        other_maj.gammaz=edge_maj.gamma*edge_maj2.z;
        node.gammaz=edge_maj.gamma^2*edge_maj.z+edge_min.gamma^2*edge_min2.z;
   end
end

function printEdges(net::HybridNetwork)
    println("Edge#\tNode1\tNode2")
    for i in (1:net.numEdges)
        println("$(net.edge[i].number)\t$(net.edge[i].node[1].number)\t$(net.edge[i].node[2].number)")
    end;
end;

function printNodes(net::HybridNetwork)
    println("Node#\tEdges numbers")
    for i in (1:net.numNodes)
        print(net.node[i].number)
        for j in (1:length(net.node[i].edge))
            print("\t$(net.node[i].edge[j].number)")
        end;
        print("\n")
    end;
end;

# find the hybrid edges for a given hybrid node
function hybridEdges(node::Node)
   if(node.hybrid)
      n=size(node.edge,1);
      hybedges=Edge[];
      for(i=1:n)
        if(node.edge[i].hybrid)
	   push!(hybedges,node.edge[i]);
	end
      end
      if(size(hybedges,1)==2)
        return hybedges
      else
        println("node has more or less than 2 hybrid edges");
      end
   else
      println("node is not hybrid");
   end
end


# setLength using updateGammaz
# check: do we want edge as parameter? or its number?
function setLength!(edge::Edge, new_length::Float64,net::HybridNetwork)
  if(new_length<0)
      error("length has to be nonnegative");
  else
      edge.length = new_length;
      edge.y = exp(-new_length);
      edge.z = 1 - edge.y;
      updateGammaz!(net);
  end
end


# setGamma using updateGammaz
# check: we need to put edge as parameter still, maybe we can have a function to search for the hybrid edge?
#        maybe we don't have the actual edge as parameter, only its number? the number of the hybrid node?
function setGamma!(edge::Edge, new_gamma::Float64, net::HybridNetwork)
 if(edge.hybrid)
	if(0 < new_gamma < 1)
	     edge.gamma = new_gamma;
	     new_gamma<0.5 ? edge.isMajor=false : edge.isMajor=true;
	     edge.isChild1 ? ind=1 : ind=2 ; # hybrid edge pointing at node 1 or 2
             hybedges=hybridEdges(edge.node[ind]); # hybrid edges for hybrid node
	     if(!isempty(hybedges) && !isa(hybedges,Nothing)) # change 1-gamma in other hybrid edge
	          isequal(hybedges[1],edge) ? other=2 : other=1;
		  hybedges[other].gamma=1.0-new_gamma;
       	          new_gamma<0.5 ? hybedges[other].isMajor=true : hybedges[other].isMajor=false;
	    end
        updateGammaz!(net);
	else
	     error("gamma has to be between 0 and 1");
        end
  else
	error("cannot change gamma in a tree edge");
  end
end
