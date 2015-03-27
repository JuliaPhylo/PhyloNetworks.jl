# example of 5 taxa tree in ipad notes
# used to test createHybrid!
# Claudia September 2014
#
# in julia: include("tree_example.jl")

ed11=Edge(11,1.5);
ed12=Edge(12,0.2);
ed3=Edge(3,0.9);
ed4=Edge(4,0.1);
ed5=Edge(5,0.2);
ed6=Edge(6,0.1);
ed7=Edge(7,0.1);
ed8=Edge(8,0.1);
ed9=Edge(9,0.1);
ed10=Edge(10,0.1);
ed13=Edge(13,0.1);


n1=Node(1,false,false,[ed11,ed5,ed6]);
n4=Node(4,true,false,[ed11]);
n5=Node(5,false,false,[ed12,ed5,ed9]);
n6=Node(6,true, false,[ed6]);
n7=Node(7,true, false,[ed12]);
n8=Node(8,true, false,[ed8]);
n9=Node(9,false,false,[ed8,ed9,ed10]);
n10=Node(10,true, false,[ed10]);


setNode!(ed5,[n5,n1]);
setNode!(ed6,[n1,n6]);
setNode!(ed8,[n8,n9]);
setNode!(ed9,[n5,n9]);
setNode!(ed10,[n9,n10]);
setNode!(ed11,[n1,n4]);
setNode!(ed12,[n5,n7]);

net=HybridNetwork([n1,n4,n5,n6,n7,n8,n9,n10],[ed5,ed6,ed8,ed9,ed10,ed11,ed12]);
