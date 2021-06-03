import DSGRN
import json, os, sys, ast
from functools import partial
from dsgrn_net_query.utilities.poset_utilities import calculate_posets_from_multiple_time_series,check_posets
from dsgrn_net_query.utilities.file_utilities import read_networks, create_results_folder
from dsgrn_net_query.utilities.signatures_no_mpi import make_db
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor


def query(network_file,params_file,resultsdir=""):
    '''
    For each epsilon in a list of epsilons, a partially ordered set (poset) of maxima and minima of time series data
    is created or accessed from the params dictionary.
    This poset is matched against the domain graph for each DSGRN parameter for each network in a list of networks.
    The result is the count of the number of DSGRN parameters with a match OR is True if there is at least one match
    and False if not, depending on the choice of the parameter "count".

    :param network_file: a .txt file containing either a single DSGRN network specification or a list of network
    specification strings in DSGRN format
    :param params_file: A .json file containing a dictionary with the keys
        "domain" : True or False (true or false in .json format), whether or not to perform a path search anywhere in the domain graph
        "stablefc" : True or False (true or false in .json format), whether or not to perform a path search within stable full cycles
        Both "domain" and "stablefc" are allowed to be True.
        "count" : True or False (true or false in .json format), whether to count all DSGRN parameters or shortcut at first success
        "datetime" : optional datetime string to append to subdirectories in resultsdir, default = system time

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

    networks = read_networks(network_file)
    params = json.load(open(params_file))

    sanity_check(params)

    posets,networks = get_posets(networks,params)

    if not networks:
        print("No networks available for analysis. Make sure network file is in the correct format\nand make sure that every network node name is the time series data or 'poset' value.")
        return None
    else:
        work_function = partial(search_over_networks, params, posets,len(networks))
        with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
            if executor is not None:
                print("Querying networks.")
                output=list(executor.map(work_function, enumerate(networks)))
                results = dict(output)
                record_results(network_file, params_file,results,resultsdir,params)


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
            posets[sort_names] = {"no_time_series_file" : pos}
            networks = check_posets(networks,posets)
    return posets,networks



def search_over_networks(params,posets,N,enum_netspec):
    '''
    Work function for parallelization.
    :param params: dictionary
    :param posets: dictionary of partially ordered sets of extrema for each time series in DSGRN format, keyed by the
            names of the genes in the time series file that are to be matched.
    :param N: size of the parameter graph for the specified network
    :param enum_netspec: an (integer, DSGRN network specification) pair
    :return: (DSGRN network specification, results) pair
    '''
    (k,netspec) = enum_netspec
    domain = params["domain"]
    stablefc = params["stablefc"]
    network = DSGRN.Network(netspec)
    names = tuple(sorted([network.name(k) for k in range(network.size())]))
    newposets = posets[names]
    ER = {}
    if params["count"] and not domain:
        dmatches, fcmatches = PathMatches_with_count_stablefc_only(network,newposets,k+1)
    elif params["count"]:
        dmatches, fcmatches = PathMatches_with_count(network,newposets,domain,stablefc)
    else:
        dmatches, fcmatches = PathMatches_without_count(network,newposets,domain,stablefc)
    if domain:
        ER["domain"]= dmatches
    if stablefc:
        ER["stablefc"]= fcmatches
    print("Network {} of {} complete.".format(k+1,N))
    sys.stdout.flush()
    return (netspec, ER)


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
    for netspec,ER in results.items():
        for search,tsdict in ER.items():
            if params[search]:
                for ts, rlist in tsdict.items():
                    key = (search,ts)
                    if key in reparse:
                        reparse[key].append((netspec,rlist))
                    else:
                        reparse[key] = [(netspec,rlist)]
    for key,list_of_tup in reparse.items():
        ts = key[1].split("/")[-1].split(".")[0]
        rname = os.path.join(resultsdir, "query_results_{}_{}.json".format(key[0], ts))
        savefile(rname,dict(list_of_tup))
    print(resultsdir)


def PathMatches_with_count_stablefc_only(network, posets, N):
    '''
    Count the number of pattern matches in the domain graph and/or stable full cycles.
    :param network: DSGRN network object.
    :param posets: The partially ordered sets that are to be matched at each epsilon in DSGRN format.
    :param N: network identifying integer
    :return: dictionary of results
    '''
    if len(posets) > 1:
        totalFC = {"all": {str(eps[0]) : set() for eps in posets[next(iter(posets))]} }
    numFCMatch = { tsfile : {str(eps[0]) : 0 for eps in poset_list} for tsfile,poset_list in posets.items()}
    paramgraph = DSGRN.ParameterGraph(network)
    # make DSGRN database to pull out stable FC parameters
    specfile = "temp{}.txt".format(N)
    dbfile = "temp{}.db".format(N)
    open(specfile,"w").write(network.specification())
    make_db(specfile,dbfile)
    params = DSGRN.StableFCQuery(DSGRN.Database(dbfile)).matches()
    os.remove(specfile)
    os.remove(dbfile)
    for paramind in params:
        domaingraph = DSGRN.DomainGraph(paramgraph.parameter(paramind))
        for tsfile, poset_list in posets.items():
            for (eps, (events, event_ordering)) in poset_list:
                patterngraph = DSGRN.PatternGraph(DSGRN.PosetOfExtrema(network,events,event_ordering))
                stabmatch, _ = stableFC_check(domaingraph,patterngraph)
                if stabmatch:
                    numFCMatch[tsfile][str(eps)]+=1
                    if len(posets) > 1:
                        totalFC["all"][str(eps)].add(paramind)
    fcmatches = {tsfile : [(float(eps),count,len(params),paramgraph.size()) for eps,count in edict.items()] for tsfile,edict in numFCMatch.items()}
    if len(posets)>1:
        fcmatches.update({"all" : [(float(eps),len(totalFC["all"][eps]),len(params),paramgraph.size()) for eps in totalFC["all"]]})
    return {},fcmatches


def PathMatches_with_count(network, posets, domain, stablefc):
    '''
    Count the number of pattern matches in the domain graph and/or stable full cycles.
    :param network: DSGRN network object.
    :param posets: The partially ordered sets that are to be matched at each epsilon in DSGRN format.
    :param domain: True or False, search over whole domain graph.
    :param stablefc: True or False search over stable full cycles only.
    :return: dictionary of results
    '''
    if len(posets) > 1:
        totalDom= {"all": {str(eps[0]) : set() for eps in posets[next(iter(posets))]} }
        totalFC = {"all": {str(eps[0]) : set() for eps in posets[next(iter(posets))]} }
    numDomMatch = { tsfile : {str(eps[0]) : 0 for eps in poset_list} for tsfile,poset_list in posets.items()}
    numFCMatch = { tsfile : {str(eps[0]) : 0 for eps in poset_list} for tsfile,poset_list in posets.items()}
    numFC = 0
    paramgraph = DSGRN.ParameterGraph(network)
    for paramind in range(paramgraph.size()):
        FC = False
        domaingraph = DSGRN.DomainGraph(paramgraph.parameter(paramind))
        for tsfile, poset_list in posets.items():
            for (eps, (events, event_ordering)) in poset_list:
                patterngraph = DSGRN.PatternGraph(DSGRN.PosetOfExtrema(network,events,event_ordering))
                if stablefc:
                    stabmatch, newFC = stableFC_check(domaingraph,patterngraph)
                    if newFC and not FC:
                        numFC +=1
                        FC = True
                    if stabmatch:
                        numFCMatch[tsfile][str(eps)]+=1
                        if len(posets) > 1:
                            totalFC["all"][str(eps)].add(paramind)
                    if stabmatch and domain:
                        numDomMatch[tsfile][str(eps)] += 1
                        if len(posets) > 1:
                            totalDom["all"][str(eps)].add(paramind)
                    elif not stabmatch and domain:
                        dommatch = domain_check(domaingraph, patterngraph)
                        if dommatch:
                            numDomMatch[tsfile][str(eps)] += 1
                            if len(posets) > 1:
                                totalDom["all"][str(eps)].add(paramind)
                if domain and not stablefc:
                    dommatch = domain_check(domaingraph,patterngraph)
                    if dommatch:
                        numDomMatch[tsfile][str(eps)]+=1
                        if len(posets) > 1:
                            totalDom["all"][str(eps)].add(paramind)
    dommatches = {tsfile : [(float(eps),count,paramgraph.size()) for eps,count in edict.items()] for tsfile,edict in numDomMatch.items()}
    fcmatches = {tsfile : [(float(eps),count,numFC,paramgraph.size()) for eps,count in edict.items()] for tsfile,edict in numFCMatch.items()}
    if len(posets)>1:
        dommatches.update({"all": [(float(eps),len(totalDom["all"][eps]),paramgraph.size()) for eps in totalDom["all"]]})
        fcmatches.update({"all" : [(float(eps),len(totalFC["all"][eps]),numFC,paramgraph.size()) for eps in totalFC["all"]]})
    return dommatches,fcmatches


def PathMatches_without_count(network, posets, domain, stablefc):
    '''
    Test for the existence of at least one pattern match in the domain graph and/or stable full cycles.
    :param network: DSGRN network object.
    :param posets: The partially ordered sets that are to be matched at each epsilon in DSGRN format.
    :param domain: True or False, search over whole domain graph.
    :param stablefc: True or False search over stable full cycles only.
    :return: dictionary of results
    '''

    def format(numDomMatch,numFCMatch):
        dommatches = {tsfile: [(float(eps), b, paramgraph.size()) for eps, b in edict.items()] for tsfile, edict in
             numDomMatch.items()}
        fcmatches = {tsfile: [(float(eps), b, paramgraph.size()) for eps, b in edict.items()] for tsfile, edict
             in numFCMatch.items()}
        return dommatches, fcmatches


    numDomMatch = { tsfile : {str(eps[0]) : False for eps in poset_list} for tsfile,poset_list in posets.items()}
    numFCMatch = { tsfile : {str(eps[0]) : False for eps in poset_list} for tsfile,poset_list in posets.items()}
    paramgraph = DSGRN.ParameterGraph(network)
    for paramind in range(paramgraph.size()):
        domaingraph = DSGRN.DomainGraph(paramgraph.parameter(paramind))
        for tsfile, poset_list in posets.items():
            for (eps, (events, event_ordering)) in poset_list:
                patterngraph = DSGRN.PatternGraph(DSGRN.PosetOfExtrema(network,events,event_ordering))
                if domain and not numDomMatch[tsfile][str(eps)]:
                    dommatch = domain_check(domaingraph,patterngraph)
                    if dommatch:
                        numDomMatch[tsfile][str(eps)] = True
                if stablefc and not numFCMatch[tsfile][str(eps)]:
                    stabmatch, newFC = stableFC_check(domaingraph,patterngraph)
                    if stabmatch:
                        numFCMatch[tsfile][str(eps)] = True
        b = True
        for tsfile,poset_list in posets.items():
            b = bool(b*all([numDomMatch[tsfile][str(eps[0])] for eps in poset_list]))
            b = bool(b*all([numFCMatch[tsfile][str(eps[0])] for eps in poset_list]))
        if b:
            return format(numDomMatch,numFCMatch)
    return format(numDomMatch,numFCMatch)


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


# def stableFC_check_buggy(domaingraph,patterngraph,paramind):
#     # Note: morsegraph.poset().stringify != morsedecomposition.poset().stringify in all cases.
#     '''
#     Check for match in any stable full cycle for one parameter
#     :param domaingraph: DSGRN domain graph object
#     :param patterngraph: DSGRN pattern graph object
#     :return: True or False for the existence of a match and True or False for the existence of a stable full cycle
#     '''
#     FC = False
#     ismatch = False
#     morsedecomposition = DSGRN.MorseDecomposition(domaingraph.digraph())
#     morsegraph = DSGRN.MorseGraph(domaingraph, morsedecomposition)
#     for i in range(0, morsegraph.poset().size()):
#         if morsegraph.annotation(i)[0] == "FC" and len(morsedecomposition.poset().children(i)) == 0:
#             FC = True
#             searchgraph = DSGRN.SearchGraph(domaingraph, i)
#             matchinggraph = DSGRN.MatchingGraph(searchgraph, patterngraph)
#             if DSGRN.PathMatch(matchinggraph):
#                 ismatch = True
#                 break
#     return ismatch,FC


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
