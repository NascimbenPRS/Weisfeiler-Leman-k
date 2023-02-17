## Weisfeiler-Leman-k
A C++ implementation of the k-dimensional Weisfeiler-Leman graph colouring algorithm to test isomorphism on a graph-pair `P`. The k-dimensional version of the algorithm iteratively colours k-tuples of both graphs, until either an identical stable colouring is reached (in which case the graphs **might** be isomorphic) or a mismatch happens (in which case the graphs are **definitely** not isomorphic). 

An alternative variant of WL-k, which colours k-subsets instead of k-tuples, is available by specifiying the "sets" selector string when launching the algorithm from command line. Such variant builds a new node- and edge-coloured graph-pair `P_k` from the original one `P`: we conjecture that WL-1 on `P_k` is equivalent to WL-k on `P`.   

# Compilation
Run `make` in the `source` folder.

# Execution
`WL_Test.exe [Graph1 file name] [Graph 2 file name] [k] [test selector string]`
- Graph files must be in DIMACS format
- Parameter `k` must be >= 1
- selector string must be "tuples" or "sets"

# Test
