import DSGRN
import json, os, sys, ast
from functools import partial
from inspect import getmembers, isfunction
from dsgrn_net_query.utilities.make_posets_from_time_series import calculate_posets_from_multiple_time_series
from dsgrn_net_query.utilities.file_utilities import read_networks
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor


def query(network_file,params_file,resultsdir=""):
    '''
    For each epsilon in a list of epsilons, a poset from time series data is created or accessed from the params dictionary.
    This poset is matched against the domain graph for each parameter in each network from a list of networks.
    The result is True if there is at least one match, False if not, or a count of the number of parameters with a match.

    :param network_file: a .txt file containing either a single DSGRN network specification or a list of network specification strings in DSGRN format
    :param params_file: A .json file containing a dictionary with the keys
        "matchingfunction" : a string or list of strings containing the name(s) of one of the matching functions in this module
        **NOTE** Cycle matches are not recommended. They will not work unless, for each variable, the first and last extrema have the same label.
        "count" : True or False, whether to count all params or shortcut at first success
        One can either specify posets directly, or extract posets from timeseries data.
        Include EITHER
        "posets" : a (quoted) dictionary of Python tuples of node names keying a list of tuples of epsilon with a DSGRN
        formatted poset:
                    '{ ("A","B","C") : [(eps11,poset11), (eps21,poset21),...], ("A","B","D") : [(eps12,poset12), (eps22,
                    poset22),...] }' (the quotes are to handle difficulties with the json format)
        OR the three keys
        "timeseriesfname" : path to a file containing the time series data from which to make posets
        "tsfile_is_row_format" : True if the time series file is in row format (times are in the first row); False if in
        column format (times are in the first column)
        "epsilons" : list of floats 0 <= x <= 0.5, one poset will be made for each x
                    Note that an epsilon of 0.10 means that the noise level is considered to be +/- 10% of the distance
                    between global maximum and global minimum for each time series. Thus all information on curve shape
                    is lost at epsilon = 0.5. It is recommended to stay far below that level
    :param resultsdir: optional path to directory where results will be written, default is current directory

    :return: Writes True (pattern match for the poset) or False (no pattern match) or
        parameter count (# successful matches) plus the number of parameters for each
         epsilon to a dictionary keyed by network spec, which is dumped to a json file:
         { networkspec : [(eps, result, num params)] }
         For PathMatchinStableFullCycle with count = True, the count of stable full cycles is added:
         { networkspec : [(eps, result, num stable full cycles, num params)] }
    '''

    networks = read_networks(network_file)
    params = json.load(open(params_file))

    if "posets" not in params:
        posets,networks = calculate_posets_from_multiple_time_series(params,networks)
    else:
        lit_posets = ast.literal_eval(params["posets"])
        posets = {}
        for names,pos in lit_posets.items():
            # make sure variables are in canonical order
            sort_names = tuple(sorted(list(names)))
            params["timeseriesfname"] = "no_time_series_file"
            posets[sort_names] = {"no_time_series_file" : pos}
    extract_queries(params)
    work_function = partial(search_over_networks, params, posets,len(networks))
    with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
        if executor is not None:
            print("Querying networks.")
            output=list(executor.map(work_function, enumerate(networks)))
            results = dict(output)
            record_results(results,resultsdir,params)


def record_results(results,resultsdir,params):
    if results:
        for name in results[next(iter(results.keys()))]:
            reparse = dict((k, []) for k in results.keys())
            for netspec in results.keys():
                reparse[netspec] = results[netspec][name]
            rname = os.path.join(resultsdir,"query_results_{}_{}.json".format(name[0],name[1].split(".")[0]))
            if os.path.exists(rname):
                os.rename(rname,rname+".old")
            json.dump(reparse,open(rname,'w'))
    else:
        tsf = params["timeseriesfname"]
        if isinstance(tsf,str):
            tsf = [tsf]
        for name in params["matchingfunction"]:
            for ts in tsf:
                rname = os.path.join(resultsdir, "query_results_{}_{}.json".format(name,ts.split(".")[0]))
                if os.path.exists(rname):
                    os.rename(rname, rname + ".old")
                json.dump(results, open(rname, 'w'))


def extract_queries(params):
    pmf = params['matchingfunction']
    matchingnames = [pmf] if isinstance(pmf,str) else pmf
    funcs = dict([o for o in getmembers(sys.modules[__name__]) if isfunction(o[1])])
    not_implemented = set(matchingnames).difference([o for o in funcs])
    if not_implemented:
        raise ValueError(
            "\nMatching function(s) {} not implemented in patternmatch.py.\n".format(not_implemented))
    params['matchingfunction'] = dict([o for o in funcs.items() if o[0] in matchingnames])


def search_over_networks(params,posets,N,enum_netspec):
    (k,netspec) = enum_netspec
    matchingfuncs = params["matchingfunction"]
    network = DSGRN.Network(netspec)
    names = tuple(sorted([network.name(k) for k in range(network.size())]))
    ER = dict([((k,ts),[]) for k in matchingfuncs.keys() for ts in posets[names]])
    for tsfile,pos_list in posets[names].items():
        for (eps, (events, event_ordering)) in pos_list:
            paramgraph, patterngraph = getGraphs(events, event_ordering, network)
            for name,mf in matchingfuncs.items():
                R = mf(paramgraph, patterngraph, params['count'])
                if name == "PathMatchInStableFullCycle" and params["count"]:
                    ER[(name,tsfile)].append((eps, R[0], R[1], paramgraph.size()))
                else:
                    ER[(name,tsfile)].append((eps, R, paramgraph.size()))
    print("Network {} of {} complete.".format(k+1,N))
    return (netspec, ER)


def getGraphs(events,event_ordering,network):
    '''
    Make pattern graph and parameter graph for the network and poset.
    '''
    poe = DSGRN.PosetOfExtrema(network,events,event_ordering)
    patterngraph = DSGRN.PatternGraph(poe)
    paramgraph = DSGRN.ParameterGraph(network)
    return paramgraph,patterngraph


def CycleMatchInStableMorseSet(paramgraph, patterngraph, count):
    '''
    Search for cycle matches in stable Morse sets only.
    :return: Integer count of parameters if count = True; if count = False return True if at least one match, False otherwise.
    '''
    # TODO: In order for cycle matches to work correctly, the last extremum on each time series with an even number of extrema must be removed
    numparams = 0
    for paramind in range(paramgraph.size()):
        domaingraph = DSGRN.DomainGraph(paramgraph.parameter(paramind))
        morsedecomposition = DSGRN.MorseDecomposition(domaingraph.digraph())
        morsegraph = DSGRN.MorseGraph(domaingraph,morsedecomposition)
        for i in range(0,morsedecomposition.poset().size()):
             if morsegraph.annotation(i)[0] in ["FC", "XC"] and len(morsedecomposition.poset().children(i)) == 0:
                searchgraph = DSGRN.SearchGraph(domaingraph,i)
                matchinggraph = DSGRN.MatchingGraph(searchgraph,patterngraph)
                if DSGRN.CycleMatch(matchinggraph):
                    if count:
                        numparams +=1
                        break
                    else:
                        return True
    return numparams if count else False


def CycleMatchInDomainGraph(paramgraph, patterngraph, count):
    '''
    Search for cycle matches anywhere in the domain graph.
    :return: Integer count of parameters if count = True; if count = False return True if at least one match, False otherwise.
    '''
    # TODO: In order for cycle matches to work correctly, the last extremum on each time series with an even number of extrema must be removed
    numparams = 0
    for paramind in range(paramgraph.size()):
        domaingraph = DSGRN.DomainGraph(paramgraph.parameter(paramind))
        searchgraph = DSGRN.SearchGraph(domaingraph)
        matchinggraph = DSGRN.MatchingGraph(searchgraph, patterngraph)
        if DSGRN.CycleMatch(matchinggraph):
            if count:
                numparams += 1
            else:
                return True
    return numparams if count else False


def PathMatchInDomainGraph(paramgraph, patterngraph, count):
    '''
    Search for path matches anywhere in the domain graph.
    :return: Integer count of parameters if count = True; if count = False return True if at least one match, False otherwise.
    '''
    numparams = 0
    for paramind in range(paramgraph.size()):
        domaingraph = DSGRN.DomainGraph(paramgraph.parameter(paramind))
        searchgraph = DSGRN.SearchGraph(domaingraph)
        matchinggraph = DSGRN.MatchingGraph(searchgraph, patterngraph)
        if DSGRN.PathMatch(matchinggraph):
            if count:
                numparams += 1
            else:
                return True
    return numparams if count else False


def PathMatchInStableFullCycle(paramgraph, patterngraph, count):
    '''
    Search for path matches in stable full cycles only.
    :return: Integer count of parameters if count = True; if count = False return True if at least one match, False otherwise.
    '''
    numparams = 0
    numFC = 0
    for paramind in range(paramgraph.size()):
        FC = False
        domaingraph = DSGRN.DomainGraph(paramgraph.parameter(paramind))
        morsedecomposition = DSGRN.MorseDecomposition(domaingraph.digraph())
        morsegraph = DSGRN.MorseGraph(domaingraph,morsedecomposition)
        for i in range(0,morsedecomposition.poset().size()):
             if morsegraph.annotation(i)[0]  == "FC" and len(morsedecomposition.poset().children(i)) == 0:
                if not FC:
                    numFC += 1
                    FC = True
                searchgraph = DSGRN.SearchGraph(domaingraph,i)
                matchinggraph = DSGRN.MatchingGraph(searchgraph,patterngraph)
                if DSGRN.PathMatch(matchinggraph):
                    if count:
                        numparams +=1
                        break
                    else:
                        return True
    return (numparams,numFC) if count else False


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
        "Calling signature has two required arguments \n " \
        "mpiexec -n <num_processes> python patternmatch.py <path_to_network_file> <path_to_parameter_file>"
        )
        exit(1)
    network_file = sys.argv[1]
    params_file = sys.argv[2]
    if len(sys.argv)>3:
        resultsdir = sys.argv[3]
        query(network_file, params_file, resultsdir)
    else:
        query(network_file,params_file)
