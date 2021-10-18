import DSGRN
import json, os, sys, ast
from functools import partial
from dsgrn_net_query.utilities.poset_utilities import calculate_posets_from_multiple_time_series,check_posets
from dsgrn_net_query.utilities.file_utilities import read_networks, create_results_folder
from mpi4py.futures import MPICommExecutor


def query(network_file,params_file,resultsdir=""):
    '''
    For each epsilon in a list of epsilons, a partially ordered set (poset) of maxima and minima of time series data
    is created or accessed from the params dictionary.
    This poset is matched against the domain graph for each DSGRN parameter for each network in a list of networks.
    The result is the count of the number of DSGRN parameters with a match OR is True if there is at least one match
    and False if not, depending on the choice of the parameter "count".

    :param network_file: a .txt file containing a single DSGRN network specification
    specification strings in DSGRN format
    :param params_file: A .json file containing a dictionary with the keys
        "domain" : True or False (true or false in .json format), whether or not to perform a path search anywhere in the domain graph
        "stablefc" : True or False (true or false in .json format), whether or not to perform a path search within stable full cycles
        Both "domain" and "stablefc" are allowed to be True.
        "count" : True or False (true or false in .json format), whether to count all DSGRN parameters or shortcut at first success
            NOTE: Only count = True is currently implemented
        "datetime" : optional datetime string to append to subdirectories in resultsdir, default = system time
        "parameter_list" : optional sublist of the parameter graph

        One can either specify posets directly, or extract posets from timeseries data.
        Include EITHER the three keys
        "timeseriesfname" : path to a file or list of files containing the time series data from which to make posets
        "tsfile_is_row_format" : True if the time series file is in row format (times are in the first row); False if in
        column format (times are in the first column)
        "epsilons" : list of floats 0 <= x <= 0.5, one poset will be made for each x
                    Note that an epsilon of 0.10 means that the noise level is considered to be +/- 10% of the distance
                    between global maximum and global minimum for each time series. Thus all information on curve shape
                    is lost at epsilon = 0.5. It is recommended to stay below that level.
        OR the single key
        "posets" : a (quoted) dictionary of Python tuples of node names keying a list of tuples of epsilon with a DSGRN
        formatted poset:
                    '{ ("A","B","C") : [(eps11,poset11), (eps21,poset21),...], ("A","B","D") : [(eps12,poset12), (eps22,
                    poset22),...] }' (the quotes are to handle difficulties with the json format)
        This key takes specialized information and is not likely to be specified in general usage.

    :param resultsdir: optional path to directory where results will be written, default is current directory

    :return: Writes a .json file containing a dictionary keyed by DSGRN network specification with a list of results.
        The results are DSGRN parameter count with successful matches, or True (existence of at least one match)
        or False (no pattern matches exist), depending on the value of the parameter "count". The size of the DSGRN
        parameter graph for the network and the value of 'epsilon' for the attempted match to the time series are
        also recorded.
        { networkspec : [(eps, result, DSGRN param graph size)] }.
        When "stablefc" and "count" are both True, the number of stable full cycles (with or without a pattern match)
        is also recorded:
        { networkspec : [(eps, result, num stable full cycles, DSGRN param graph size)] }

        Separate files will be written for "domain", "stablefc", and each timeseries file name.
        When multiple time series are specified and "count" is True, a file with the string "all" in the name will be
        written with the aggregated results across time series datasets (this eliminates the double-counting of DSGRN
        parameters that can occur across multiple time series). The file ending in "all" records the number of DSGRN
        parameters with a match to at least one time series dataset.
    '''

    spec = read_networks(network_file)
    param_dict = json.load(open(params_file))

    sanity_check(param_dict)

    posets,networks = get_posets(spec,param_dict)
    print("Querying networks.\n")

    if not networks:
        print("No networks available for analysis. Make sure network file is in the correct format\nand make sure that every network node name is the time series data or 'poset' value.")
        return None
    else:
        results = {}
        if param_dict["count"]:
            with MPICommExecutor() as executor:
                spec = spec[0]
                results[spec] = {}
                network = DSGRN.Network(spec)
                param_graph = DSGRN.ParameterGraph(network)
                if "parameter_list" in param_dict:
                    dsgrn_params = [(p,param_graph.parameter(p)) for p in param_dict["parameter_list"]]
                else:
                    dsgrn_params = [(p,param_graph.parameter(p)) for p in range(param_graph.size())]
                names = tuple(sorted([network.name(k) for k in range(network.size())]))
                work_function = partial(PathMatch, network, posets[names], param_dict["domain"], param_dict["stablefc"])
                output=dict(executor.map(work_function, dsgrn_params))
                results[spec] = reformat_output(output, list(posets[names].keys()), param_dict, param_graph.size())
                print("Network {} of {} complete.".format(1, 1))
                sys.stdout.flush()
                record_results(network_file, params_file, results, resultsdir, param_dict)
        else:
            raise ValueError("Existence of path match without counting is not yet implemented for large networks. Use CountPatternMatch.py.")


def sanity_check(params):
    '''
    Checks to be sure the correct keys are in the dictionary params.
    :param params: dictionary
    :return: None, errors are raised.
    '''
    if any(["timeseriesfname" in params, "tsfile_is_row_format" in params, "epsilons" in params]) and "posets" in params:
        raise ValueError("Only one of 'posets' or the three keys 'timeseriesfname', 'tsfile_is_row_format' and 'epsilons' may be specified in the parameter file.")
    if any(["timeseriesfname" in params, "tsfile_is_row_format" in params, "epsilons" in params]) and not all(["timeseriesfname" in params, "tsfile_is_row_format" in params, "epsilons" in params]):
        raise ValueError("All of the three keys 'timeseriesfname', 'tsfile_is_row_format' and 'epsilons' must be specified in the parameter file.")
    if not "posets" in params and not "timeseriesfname" in params:
        raise ValueError("Either 'posets' or the three keys 'timeseriesfname', 'tsfile_is_row_format' and 'epsilons' must be specified in the parameter file.")
    if any(["domain" not in params, "stablefc" not in params, "count" not in params]):
        raise ValueError("All of the three keys 'domain', 'stablefc' and 'count' must be specified in the parameter file.")


def get_posets(networks,params):
    if "posets" not in params:
        posets,networks = calculate_posets_from_multiple_time_series(params,networks)
    else:
        lit_posets = ast.literal_eval(params["posets"])
        posets = {}
        for names,pos in lit_posets.items():
            # make sure variables are in canonical order
            sort_names = tuple(sorted(list(names)))
            params["timeseriesfname"] = "no_time_series_file"
            posets[sort_names] = {"no_time_series_file" : sorted(pos)}
            networks = check_posets(networks,posets)
            if "epsilons" not in params:
                # side effect on parameter dictionary
                params["epsilons"] = sorted([eps for eps,_ in pos])
    return posets,networks


def record_results(network_file, params_file,results,resultsdir,params):
    '''
    Record results in a .json file.
    :param network_file: The input .txt file containing the list of DSGRN network specifications.
    :param params_file: The input .json parameter file.
    :param results: The dictionary of results.
    :param resultsdir: The location to save the dictionary of results.
    :param params: The dictionary of parameters generated from the .json parameter file.
    :return: None. File is written.
    '''
    if "datetime" in params:
        resultsdir = create_results_folder(network_file, params_file, resultsdir,params["datetime"])
    else:
        resultsdir = create_results_folder(network_file, params_file, resultsdir)

    def savefile(rname,rdict):
        if os.path.exists(rname):
            os.rename(rname, rname + ".old")
        json.dump(rdict, open(rname, 'w'))

    reparse = {}
    reparse_params = {}
    for netspec,ER in results.items():
        for search,tsdict in ER.items():
            if search != "param_list" and params[search]:
                for ts, rlist in tsdict.items():
                    key = (search,ts)
                    if key in reparse:
                        reparse[key].append((netspec,rlist))
                    else:
                        reparse[key] = [(netspec,rlist)]
            elif search == "param_list":
                for p_search,p_tsdict in tsdict.items():
                    for pts, plist in p_tsdict.items():
                        key = (p_search,pts)
                        if key in reparse_params:
                            reparse_params[key].append((netspec,plist))
                        else:
                            reparse_params[key] = [(netspec,plist)]
    for key,list_of_tup in reparse.items():
        ts = os.path.splitext(os.path.basename(key[1]))[0]
        # ts = key[1].split("/")[-1].split(".")[0]
        rname = os.path.join(resultsdir, "query_results_{}_{}.json".format(key[0], ts))
        savefile(rname,dict(list_of_tup))
    for key,list_of_tup in reparse_params.items():
        ts = os.path.splitext(os.path.basename(key[1]))[0]
        rname = os.path.join(resultsdir, "query_results_{}_{}_param_list.json".format(key[0], ts))
        savefile(rname,dict(list_of_tup))
    open(".query_results.log","w").write(resultsdir)


def reformat_output(output, tsfiles, param_dict, pgsize):
        ts_keys = tsfiles[:]
        if len(tsfiles) > 1:
            ts_keys += ["all"]
            all_match = {}
        domain_bool = param_dict["domain"]
        stablefc_bool = param_dict["stablefc"]
        eps = param_dict["epsilons"]
        res = {"param_list": {} }
        if domain_bool:
            res["domain"] = {t: set([]) for t in ts_keys}
            res["param_list"]["domain"] = {}
            if len(tsfiles) > 1:
                all_match["domain"] = {e: set([]) for e in eps}
        if stablefc_bool:
            res["stablefc"] = {t: set([]) for t in ts_keys}
            res["param_list"]["stablefc"] = {}
            if len(tsfiles) > 1:
                all_match["stablefc"] = {e: set([]) for e in eps}
        for tsfile in tsfiles:
            domain_match = {e: set([]) for e in eps}
            stablefc = 0
            stablefc_match = {e: set([]) for e in eps}
            for param_index, data in output.items():
                if domain_bool:
                    for epsilon, booln in data["domain"][tsfile].items():
                        if booln:
                            domain_match[epsilon].add(param_index)
                            if len(tsfiles) > 1:
                                all_match["domain"][epsilon].add(param_index)
                if stablefc_bool:
                    stablefc += data["stablefc"]
                    for epsilon, booln in data["match"][tsfile].items():
                        if booln:
                            stablefc_match[epsilon].add(param_index)
                            if len(tsfiles) > 1:
                                all_match["stablefc"][epsilon].add(param_index)
            if domain_bool:
                res["domain"][tsfile] = [d + [pgsize] for d in sorted([[m[0],len(m[1])] for m in domain_match.items()])]
                res["param_list"]["domain"][tsfile] = sorted([[m[0],sorted(list(m[1]))] for m in domain_match.items()])
            if stablefc_bool:
                res["stablefc"][tsfile] = [d + [stablefc, pgsize] for d in sorted([[m[0],len(m[1])] for m in stablefc_match.items()])]
                res["param_list"]["stablefc"][tsfile] = sorted([[m[0],sorted(list(m[1]))] for m in stablefc_match.items()])
        if len(tsfiles) > 1:
            if domain_bool:
                res["domain"]["all"] = sorted([[e, len(plist), pgsize] for e, plist in all_match["domain"].items()])
                res["param_list"]["domain"]["all"] = sorted([[m[0],sorted(list(m[1]))] for m in all_match["domain"].items()])
            if stablefc_bool:
                res["stablefc"]["all"] = sorted([[e, len(plist), stablefc, pgsize] for e, plist in all_match["stablefc"].items()])
                res["param_list"]["stablefc"]["all"] = sorted([[m[0],sorted(list(m[1]))] for m in all_match["stablefc"].items()])
        return res


def PathMatch(network, posets, domain, stablefc,dsgrn_param):
    '''
    Test for the existence of at least one pattern match in the domain graph and/or stable full cycles.
    :param network: DSGRN.Network object
    :param dsgrn_params: list of tuples of parameter indices and DSGRN.Parameter objects.
    :param posets: The partially ordered sets that are to be matched at each epsilon in DSGRN format.
    :param domain: True or False, search over whole domain graph.
    :param stablefc: True or False search over stable full cycles only.
    :return: {"domain" : {tsfile : {str(eps) : True or False}}, "match" : {tsfile : {str(eps) : True or False}, "stablefc" : True or False}
    '''
    param_index,param = dsgrn_param
    DomMatch = { tsfile : {} for tsfile,_ in posets.items()}
    FCMatch = { tsfile : {} for tsfile,_ in posets.items()}
    FC = False
    domaingraph = DSGRN.DomainGraph(param)
    for tsfile, poset_list in posets.items():
        for (eps, (events, event_ordering)) in poset_list:
            patterngraph = DSGRN.PatternGraph(DSGRN.PosetOfExtrema(network,events,event_ordering))
            if domain:
                DomMatch[tsfile][eps] = domain_check(domaingraph,patterngraph)
            if stablefc:
                FCMatch[tsfile][eps], FC = stableFC_check(domaingraph,patterngraph)
    results = {}
    if domain:
        results["domain"] = DomMatch
    if stablefc:
        results["stablefc"] = FC
        results["match"] = FCMatch
    return (param_index,results)


def domain_check(domaingraph,patterngraph):
    '''
    Check for match in domain graph for one parameter
    :param domaingraph: DSGRN domain graph object
    :param patterngraph: DSGRN pattern graph object
    :return: True or False
    '''
    ismatch = False
    searchgraph = DSGRN.SearchGraph(domaingraph)
    matchinggraph = DSGRN.MatchingGraph(searchgraph, patterngraph)
    if DSGRN.PathMatch(matchinggraph):
        ismatch = True
    return ismatch


def stableFC_check(domaingraph,patterngraph):
    '''
    Check for match in any stable full cycle for one parameter
    :param domaingraph: DSGRN domain graph object
    :param patterngraph: DSGRN pattern graph object
    :return: True or False for the existence of a match and True or False for the existence of a stable full cycle
    '''
    FC = False
    ismatch = False
    morsegraph = DSGRN.MorseGraph(domaingraph)
    for i in range(0, morsegraph.poset().size()):
        if morsegraph.annotation(i)[0] == "FC" and len(morsegraph.poset().children(i)) == 0:
            FC = True
            searchgraph = DSGRN.SearchGraph(domaingraph, i)
            matchinggraph = DSGRN.MatchingGraph(searchgraph, patterngraph)
            if DSGRN.PathMatch(matchinggraph):
                ismatch = True
                break
    return ismatch,FC



if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
        "Calling signature has two required arguments \n " \
        "mpiexec -n <num_processes> python CountPatternMatch.py <path_to_network_file> <path_to_parameter_file>"
        )
        exit(1)
    network_file = sys.argv[1]
    params_file = sys.argv[2]
    if len(sys.argv)>3:
        resultsdir = sys.argv[3]
        query(network_file, params_file, resultsdir)
    else:
        query(network_file,params_file)
