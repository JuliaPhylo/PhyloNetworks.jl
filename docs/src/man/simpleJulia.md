## Julia types

Objects in Julia are called *types*. You can manipulate objects in Julia quite easily.
For example, let's take an object *d* from DataCF type.
```julia
d=readTrees2CF("tableCF.txt");
```

Typing the following command will provide a list of objects saved in memory:
```julia
whos()
```

If we want to know the type of a particular object, we do:
```julia
typeof(d)
```

If we want to about the attributes the object has, we can type ? in Julia, followed by *DataCF* for a description.

For example, we can get the number of 4-taxon subsets in the data with
```julia
d.numQuartets
```

In particular, the attribute *quartet* contains a vector of Quartet objects inside *d*, so
```julia
d.quartet[1].taxon
```
will provide the list of taxon names for the first 4-taxon subset in the data. We can corroborate this is the firs 4-taxon subset by checking the file *tableCF.txt*.

To see the observed CF, we can type
```julia
d.quartet[1].obsCF
```

We can verify the type with
```julia
typeof(d.quartet[1])
```

You can also read a simple network in Julia and print the list of edges
```julia
str = "(A,((B,#H1),(C,(D)#H1)));"
net = readTopology(str)
printEdges(net)
```
You see that the edges do not have branch lengths, and the hybrid edges do not have gamma values, you can set them with
```julia
setLength!(net.edge[1],1.9)
setGamma!(net.edge[3],0.8)
printEdges(net)
```
where 1 and 3 correspond to the position of the given edge to modify in the list of edges.
You can only change the gamma value of hybrid edges.
Such an attempt below will cause an error with a message to explain that
the edge was a tree edge:
```julia
setGamma!(net.edge[4],0.7) # should return error message
```
