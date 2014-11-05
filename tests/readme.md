### Tests functions

#### add hybridization

tests_5taxon.jl runs all the tests for the eight 5taxon networks of
starting with tree_example.jl and adding one hybridization

It calls add_hybrid_caseC,D,E,F,G,H,I,J.jl

aux functions: print_add.jl and test_functions_5taxon.jl

#### delete hybridization

tests_5taxon_delete.jl runs all the tests for the eight 5taxon
networks of starting with tree_example.jl and adding one
hybridization, and then deleting it and comparing to the original tree
example

It calls delete_hybrid_caseC,D,E,F,G,H,I,J.jl

aux functions: test_functions_5taxon.jl

#### read topology

test_5taxon_readTopology.jl runs all the tests for the eight 5taxon
networks by reading from parenthetical format and then updating

aux functions: test_functions_5taxon_read.jl
