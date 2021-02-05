# Getting started

This module accepts a list of DSGRN-computable networks in a text file and performs a user-specified DSGRN query on every parameter of every network.

__References:__ http://epubs.siam.org/doi/abs/10.1137/15M1052743, https://link.springer.com/chapter/10.1007/978-3-319-67471-1_19, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5975363/, https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006121, https://epubs.siam.org/doi/abs/10.1137/17M1134548, https://doi.org/10.1007/s00285-020-01471-4


__Dependencies:__ `Python 3.6/3.7`, `mpi4py 3.0.3`, `pandas`, `progressbar2`, `DSGRN` (https://github.com/shaunharker/DSGRN or https://github.com/marciogameiro/DSGRN), `min_interval_posets` (https://github.com/breecummins/min_interval_posets), and `dsgrn_utilities` (https://github.com/breecummins/dsgrn_utilities).

After installing all dependencies according to their instructions, do
```bash    
    cd dsgrn_net_query
    . install.sh
```

To run tests to confirm the code is working as expected, do

```bash    
    cd tests
    pytest
```

The recommended way to do a query is to use the `dsgrn_net_query/call_job.py` function.
```bash
    cd dsgrn_net_query
    python call_job.py <num_processes> <querymodule.py> <networks_file.txt> <params.json> <optional_results_directory>
```
This function saves commandline output to the log file `dsgrn_net_query.log`. 

NOTE: Calling `call_job.py` with the query module `CountStableFC_large_networks.py` means the number of processes will be doubly specified, since the number of processes is also a required argument in the parameter file for `CountStableFC_large_networks.py` only. The commandline argument will overwrite whatever number of processes is in the parameter file.

The argument `querymodule.py` is any module in `dsgrn_net_query/src/dsgrn_net_query/queries`.
Alternatively, direct calls on the command line use the full path to the query module. Here is the call for any module other than `CountStableFC_large_networks.py`.
```bash    
    cd dsgrn_net_query
    mpiexec -n <num_processes> python src/dsgrn_net_query/queries/<querymodule.py> <networks_file.txt> <params.json> <optional_results_directory>
```    
The call for `CountStableFC_large_networks.py` is 
```bash    
    cd dsgrn_net_query
    python src/dsgrn_net_query/queries/CountStableFC_large_networks.py <networks_file.txt> <params.json> <optional_results_directory>
```   
 Depending on the size and number of the networks, these computations can take a long time, and it is recommended to run via a scheduler or in the background.

All functions create a unique date-time stamped folder in which to store results, so that overwriting old results is not possible.


# Inputs 

`querymodule.py`           =   any module in dsgrn_net_query/queries; currently the following queries are available: `CountFPMatch.py`, `CountStableFC.py`, `CountStableFC_large_networks.py`, `CountPatternMatch.py`, and `CountPatternMatch_large_networks.py`.

`networks_file.txt`         =   path to a `.txt` file containing either a single DSGRN network specification
                            or a list of them (comma-separated and surrounded by square
                            brackets)
    
`params.json`    =     path to a `.json` file containing a dictionary with query module specific arguments. All queries require the key "count" which is 'true' or 'false' (no quotes) in the .json file. This determines whether to count all DSGRN parameters at which the desired query is true, or to check only for existence at at least one parameter. See individual query documentation strings for other arguments.
                            
`optional_results_directory`     =   optional path to a directory where results are to be stored; 
                            default is current directory
                            


# Outputs

Query output is saved as a dictionary to a `.json` file, for example
```
query_results.json
``` 
This can be imported as a Python dictionary using 
```python
import json
results = json.load(open("query_results.json"))
``` 
The keys are the DSGRN network specifications, and the values are usually `[# matches, param_graph_size]`, if `count = true`, and `[true_or_false, param_graph_size]`, if `count = false`. 

The module `CountPatternMatch.py` can save multiple files depending how many searches are indicated by the user. The example: `query_results_domain_wt_rnaseq_ts.json` indicates that queries were performed using the subroutine that checks for a match to the timeseries file `wt_rnaseq_ts.csv` (or `.tsv`) by searching through the whole domain graph. As before, the keys are DSGRN network specifications, and now the values are:
```
[(epsilon_1, # matches, param_graph_size), (epsilon_2, # matches, param_graph_size), ... ] 
``` 
If the output file has the form `query_results_stablefc_wt_rnaseq_ts.json`, then the query was performed only within stable full cycles of the Morse graph. The values for the network specifications are of the form:
```
[(epsilon_1, # matches, # stable full cycles, param_graph_size),   
 (epsilon_2, # matches, # stable full cycles, param_graph_size), ... ] 
``` 
Again, `# matches` will be `true` or `false` if `count` is set to `false`.

   # Troubleshooting

1. There are no DSGRN query matches.
    * This probably means you are in the wrong part of network space. Repeat the process with a different set of networks.
    * In the particular case of the query `CountPatternMatch.py`, there are additional parameters that affect search results. The choice of noise levels, `epsilons`, will strongly affect the number of matches. Generally speaking, `epsilons` in the range `[0,0.15]` are reasonable. An epsilon of 0.15 means 15% noise both above and below the curve; i.e., 30% of the difference between the global maximum and global minimum of each curve. Noise levels much higher than this are very permissive and may lead to spurious matches. Very low noise levels close or equal to zero may result in too few or no matches. Additionally, searches with `stablefc = true` are far more restrictive than searches for `domain = true`. For that reason it may be preferable to perform both searches when pattern matches within stable full cycles are desired, for the purpose of comparison.
    
2. You get an error when trying to query a specific network.
    * There are many causes of this. However, a common one is trying to query a network that is not DSGRN-computable.  To check if your network is computable, open ipython or a Jupyter notebook and do
        
        ``` 
        import DSGRN
        network = DSGRN.Network("networks_file.txt")
        pg = DSGRN.ParameterGraph(network) 
        ```
        if `"networks_file.txt"` has a single DSGRN network specification, or
        ```
        import ast, DSGRN
        networks = ast.literal_eval(open("networks_file.txt").read())
        network = DSGRN.Network(networks[i]) #For the i-th network.
        pg = DSGRN.ParameterGraph(network) 
        ```
        if `"networks_file.txt"` has list syntax. If you get a `Could not find logic resource` error, then your network is not DSGRN computable. The most likely cause of this is either that you have a version of DSGRN that will not compute any network that has a node without an out-edge, or that there are not too many in- or out-edges for any one node. At the time of this writing, 5 in-edges or 5 out-edges is likely too many (although not always). 
    * Another possibility when using `CountStableFC_large_networks.py` is that an out-of-memory error occurred. This can happen with very large networks when trying to compute the DSGRN database, which happens inside `CountStableFC_large_networks.py`, but in the other queries.
   
    
# Customizing 

The `queries` package is extensible. Users can add DSGRN query modules to the package `dsgrn_net_query.queries`. The required API is a function with the following signature:

      newmodule.query(networks_file.txt, params.json, optional_results_directory="")

  The query must be parallelized using `mpi4py` (see the function `query` in existing query modules). The results must be saved to a `.json` file with a dictionary keyed by network specifications within the `optional_results_directory` (see the function `record_results` in existing query modules). 

