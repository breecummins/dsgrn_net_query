# Getting started

This module accepts a list of DSGRN-computable networks in a text file and performs a user-specified DSGRN query on every parameter of every network.

__References:__ http://epubs.siam.org/doi/abs/10.1137/15M1052743, https://link.springer.com/chapter/10.1007/978-3-319-67471-1_19, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5975363/, https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006121

__Dependencies:__ Python 3.6/3.7, mpi4py 3.0.3, pandas, progressbar2, DSGRN (https://github.com/shaunharker/DSGRN or https://github.com/marciogameiro/DSGRN), and min_interval_posets (https://github.com/breecummins/min_interval_posets).

To install, do
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
python call_job.py <num_processes> <querymodule.py> <networks_file.txt> <params.json> <optional_results_directory>

```
The `querymodule.py` argument is any module in `dsgrn_net_query/queries` and the other inputs will be discussed below (see the files in the `tests` folder for examples of the input arguments to `call_jobs.py`.) The `call_job.py` function creates a unique date-time stamped folder in which to store results, so that overwriting old results is not possible. 

Alternatively, direct calls on the command line look like this (except for `CountStableFC_large_networks.py`).
```bash    
    mpiexec -n <num_processes> python <querymodule.py> <networks_file.txt> <params.json> <optional_results_directory>
```    
The call for `CountStableFC_large_networks.py` is 
```bash    
    python CountStableFC_large_networks.py <networks_file.txt> <params.json> <optional_results_directory>
```   
Neither of these last two calls creates a date-time stamped folder, so overwriting of results is possible.



# Inputs 

`querymodule.py`           =   any module in dsgrn_net_query/queries; Example: CountStableFC.py

`networks_file.txt`         =   path to a `.txt` file containing either a single DSGRN network specification
                            or a list of them (comma-separated and surrounded by square
                            brackets)
    
`params.json`    =     path to a `.json` file containing a dictionary containing query module specific arguments -- see 
                            individual query documentation. Can be empty for some queries. 
                            
*TODO: Insert arguments for each module here.*

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

However, the module `patternmatch.py` can save multiple files depending how many searches are indicated by the user. The example: `query_results_PathMatchInDomainGraph_wt_rnaseq_ts.json` indicates that queries were performed using the subroutine `PathMatchInDomainGraph` using the timeseries file `wt_rnaseq_ts.csv`. Both `.csv` and `.tsv` file types are acceptable. The file itself contains a list of results:
```
[(epsilon_1, # matches, param_graph_size), (epsilon_2, # matches, param_graph_size), ... ] 
``` 
If the subroutine `PathMatchInStableFullCycle` is specified, then the results are 
```
[(epsilon_1, # matches, # stable full cycles, param_graph_size),   
 (epsilon_2, # matches, # stable full cycles, param_graph_size), ... ] 
``` 
Again, `# matches` will be `true_or_false` if `count` is set to `false`.

# Customizing 

The `queries` package is extensible. Users can add DSGRN query modules to the package `dsgrn_net_query.queries` for inclusion in parameter files. The required API is a function with the following signature:

      newmodule.query(networks_file.txt, params.json, optional_results_directory="")

  The results must be saved to a `.json` file are saved to a file within the `optional_results_directory`. See the `queries` folder for already implemented queries with this API.

   # Troubleshooting

1. There are no DSGRN query matches.
    * This probably means you are in the wrong part of network space. Repeat the process with a different set of networks.
    
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

3. In the particular case of the query `patternmatch.py`, there are additional parameters that affect search results. The choice of noise levels, `epsilons`, will strongly affect the number of matches. Generally speaking, `epsilons` in the range `[0,0.15]` are reasonable. An epsilon of 0.15 means 15% noise both above and below the curve; i.e., 30% of the difference between the global maximum and global minimum of each curve. Noise levels much higher than this are very permissive and may lead to spurious matches. Very low noise levels close or equal to zero may result in too few matches.

    Additionally, the type of search may be more or less restrictive. `PathMatchinDomainGraph` is more permissive than `PathMatchinStableFullCycle`. Use of the `CycleMatch` functions are not recommended at this time, since they are extremely dependent on precisely the integer number of periods in the time series that all start at an extremum, rather than half-max. Because of this restriction, matches are extremely unlikely to occur strictly for spurious reasons.