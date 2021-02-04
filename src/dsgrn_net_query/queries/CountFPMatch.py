import DSGRN
import os, json, sys
from functools import partial
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from dsgrn_net_query.utilities.file_utilities import read_networks, create_results_folder
from dsgrn_net_query.utilities.parameter_utilities import get_neighbors


def query(network_file,params_file,resultsdir=""):
    '''
    Take the intersection of an arbitrary number of DSGRN fixed points in a list.

    :param network_file: a .txt file containing either a single DSGRN network specification or
                a list of network specification strings in DSGRN format
    :param params_file: A json file with a dictionary containing the keys
                    "included_bounds", "excluded_bounds", and "count".
                    The two "bounds" variables are each a list of dictionaries of variable names common to all network
                    specifications with an associated integer range.
                    Example: [{"X1":[2,2],"X2":[1,1],"X3":[0,1]},{"X1":[0,1],"X2":[1,1],"X3":[2,3]}]
                    The integer ranges are the matching conditions for an FP.
                    For example, if there are four variables X1, X2, X3, X4 in the network spec,
                    the FP (2,1,0,*) would be a match to the first fixed point for any value of *.
                    The "included_bounds" are those fixed points that must be present and
                    the "excluded_bounds" are those that must be absent. Either one or both may be empty. When they are both empty,
                    the count is just the number of parameters with at least one fixed point.
                    "count" : True or False (true or false in .json format);
                    whether to count all parameters with a match or shortcut at first success
                    "datetime" : optional datetime string to append to subdirectories in resultsdir, default = system time
                    "neighbors" : optional True or False (true or false in .json format) stating whether to query DSGRN parameters that
                    neighbor essential DSGRN parameters, neighbor-checking is computationally expensive, default = False

    :param resultsdir: optional path to directory where results will be written, default is current directory

    :return: Writes a .json file containing a dictionary keyed by DSGRN network specification with a list of results.
            The results are DSGRN parameter count with successful matches to the fixed point bounds, or True
            (existence of at least one match) or False (no matches exist), depending on the value of the parameter "count".
            The size of the DSGRN parameter graph for the network is also recorded.
            { networkspec : [result, DSGRN param graph size] }.
            If "hex_constraints" are specified, then the number of DSGRN parameters that satsify those hex constraints
            (with or without a bounds match) are also recorded.
            { networkspec : [result, num params with hex constraints, DSGRN param graph size] }.
    '''

    networks = read_networks(network_file)
    params = json.load(open(params_file))
    datetime = None if "datetime" not in params else params["datetime"]

    sanity_check(params)

    if not networks:
        raise ValueError("No networks available for analysis. Make sure network file is in the correct format.")
    else:
        work_function = partial(search_over_networks, params, len(networks))
        with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
            if executor is not None:
                print("Querying networks.")
                output=list(executor.map(work_function, enumerate(networks)))
                results = dict(output)
                record_results(network_file, params_file,results,resultsdir,datetime)


def sanity_check(params):
    '''
    Checks to be sure the correct keys are in the dictionary params.
    :param params: dictionary
    :return: None, errors are raised.
    '''
    if not all(["included_bounds" in params,"excluded_bounds" in params, "count" in params]):
        raise ValueError("The parameter file must contain keys 'included_bounds', 'excluded_bounds', and 'count'.")


def record_results(network_file, params_file,results,resultsdir,datetime):
    '''
    Record results in a .json file.
    :param network_file: The input .txt file containing the list of DSGRN network specifications.
    :param params_file: The input .json parameter file.
    :param results: The dictionary of results.
    :param resultsdir: The location to save the dictionary of results.
    :param datetime: None or string with datetime
    :return: None. File is written.
    '''
    resultsdir = create_results_folder(network_file, params_file, resultsdir,datetime)
    rname = os.path.join(resultsdir,"query_results.json")
    if os.path.exists(rname):
        os.rename(rname,rname+".old")
    json.dump(results,open(rname,'w'))
    print(resultsdir)


def search_over_networks(params,N,enum_network):
    '''
    Work function for parallelization.
    :param params: dictionary containing the keys "included_bounds", "excluded_bounds", and "count"
    :param N: Size of the DSGRN parameter graph
    :param enum_network: An (integer, DSGRN network specification) pair
    :return: (DSGRN network specification, results) pair
    '''
    (k, netspec) = enum_network
    if "neighbors" in params and params["neighbors"] is True:
        noness_netspec, paramlist = get_neighbors(netspec)
        network,parametergraph = getpg(noness_netspec)
    else:
        network, parametergraph = getpg(netspec)
        paramlist = list(range(parametergraph.size()))
    numparams = 0
    for p in paramlist:
        if have_match(network, parametergraph.parameter(p), params["included_bounds"], params["excluded_bounds"]):
            if params["count"]:
                numparams +=1
            else:
                print("Network {} of {} complete.".format(k + 1, N))
                sys.stdout.flush()
                return (netspec,(True, parametergraph.size()))
    print("Network {} of {} complete.".format(k + 1, N))
    sys.stdout.flush()
    if params["count"]:
        return netspec,(numparams,len(paramlist))
    else:
        return netspec,(False,len(paramlist))


def getpg(netspec):
    '''
    Calculate DSGRN parameter graph
    :param netspec: DSGRN network specification
    :return: (DSGRN.Network object, DSGRN.ParameterGraph object)
    '''
    network = DSGRN.Network(netspec)
    parametergraph = DSGRN.ParameterGraph(network)
    return network,parametergraph


def have_match(network,parameter,included_bounds,excluded_bounds):
    '''
    Check if both included and excluded bounds are satisfied for any fixed point of the specified DSGRN parameter.
    :param network: DSGRN.Network object
    :param parameter: DSGRN.Parameter object
    :param included_bounds: List of dictionaries of DSGRN fixed point bounds to include
    :param excluded_bounds: List of dictionaries of DSGRN fixed point bounds to exclude
    :return: True or False
    '''
    stable_FP_annotations = DSGRN_Computation(parameter)
    included = all_included(network, included_bounds, stable_FP_annotations)
    excluded = all_excluded(network, excluded_bounds, stable_FP_annotations)
    if included and excluded:
        return True
    return False


def DSGRN_Computation(parameter):
    '''
    Get DSGRN annotations for all Morse sets that are fixed points.
    :param parameter: DSGRN.Parameter object
    :return: list of DSGRN annotations
    '''
    dg = DSGRN.DomainGraph(parameter)
    mg = DSGRN.MorseGraph(dg)
    return [mg.annotation(i)[0] for i in range(0, mg.poset().size()) if is_FP(mg.annotation(i)[0]) and len(mg.poset(
        ).children(i)) == 0]


def is_FP(annotation):
    '''
    Specifies whether a Morse set is a fixed point.
    :param annotation: DSGRN annotation string
    :return: True or False
    '''
    return annotation.startswith("FP")


def is_FP_match(bounds_ind, annotation):
    '''
    Specifies whether a DSGRN fixed point satisfies the user-specified bounds on fixed point location.
    :param bounds_ind: User-specified fixed point bounds using DSGRN indices for the variables.
    :param annotation: DSGRN annotation string
    :return: True or False
    '''
    digits = [int(s) for s in annotation.replace(",", "").split() if s.isdigit()]
    return all(digits[k] >= bounds_ind[k][0] and digits[k] <= bounds_ind[k][1]
           for k in bounds_ind)


def is_multistable_match(network,b,stable_FP_annotations):
    '''
    Checks across all fixed points for a match to bound b.
    :param network: DSGRN.Network object
    :param b: DSGRN fixed point bounds
    :param stable_FP_annotations: list of DSGRN annotations
    :return: True or False
    '''
    bounds_ind = {network.index(str(k)): b[k] for k in b}
    return any([is_FP_match(bounds_ind, a) for a in stable_FP_annotations])


def all_included(network,included_bounds,stable_FP_annotations):
    '''
    Checks if all included bounds are present in the collection of Morse sets.
    :param network: DSGRN.Network object
    :param included_bounds: list of DSGRN fixed point bounds
    :param stable_FP_annotations: list of DSGRN annotations
    :return: True or False
    '''
    for b in included_bounds:
        if not is_multistable_match(network, b, stable_FP_annotations):
            return False
    return True


def all_excluded(network,excluded_bounds,stable_FP_annotations):
    '''
    Checks if all excluded bounds are absent in the collection of Morse sets.
    :param network: DSGRN.Network object
    :param excluded_bounds: list of DSGRN fixed point bounds
    :param stable_FP_annotations: list of DSGRN annotations
    :return: True or False
    '''
    for b in excluded_bounds:
        if is_multistable_match(network, b, stable_FP_annotations):
            return False
    return True


if __name__ == "__main__":
     if len(sys.argv) < 3:
        print(
        "Calling signature has two required arguments \n " \
        "mpiexec -n <num_processes> python CountFPMatch.py <path_to_network_file> <path_to_parameter_file>"
        )
        exit(1)
     network_file = sys.argv[1]
     params_file = sys.argv[2]
     if len(sys.argv)>3:
        resultsdir = sys.argv[3]
        query(network_file, params_file, resultsdir)
     else:
        query(network_file,params_file)
