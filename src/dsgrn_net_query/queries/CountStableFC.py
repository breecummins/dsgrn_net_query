import DSGRN
import os, json,sys
from functools import partial
from dsgrn_net_query.utilities.file_utilities import read_networks, create_results_folder
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor


def query(network_file,params_file="",resultsdir=""):
    '''
    :param network_file: a .txt file containing either a single DSGRN network specification or a list of network specification strings in DSGRN format
    :param params_file: A json file with the key
            "count" = True or False (true or false in .json format)
                        whether or not to return the number of matches (True) or just whether or not there is at least one match (False)
            "datetime" : optional datetime string to append to subdirectories in resultsdir, default = system time
    :param resultsdir: optional path to directory where results will be written, default is current directory

    :return:  Writes a .json file containing a dictionary keyed by DSGRN network specification with a list of results.
            The results are DSGRN parameter count that have at least one Morse set that is a stable full cycle,
            or True (existence of at least one stable full cycle) or False (none exist), depending on the value of the parameter "count".
            The size of the DSGRN parameter graph for the network is also recorded.
            { networkspec : [result, DSGRN param graph size] }.
    '''

    networks = read_networks(network_file)
    params = json.load(open(params_file))
    datetime = None if "datetime" not in params else params["datetime"]

    if not networks:
        raise ValueError("No networks available for analysis. Make sure network file is in the correct format.")
    else:
        count = sanity_check(params)
        work_function = partial(search_over_networks, count, len(networks))
        with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
            if executor is not None:
                print("Querying networks.")
                output=list(executor.map(work_function, enumerate(networks)))
                results = dict(output)
                record_results(network_file,params_file,results,resultsdir,datetime)


def sanity_check(params):
    '''
    Checks to be sure the correct keys are in the dictionary params.
    :param params: dictionary
    :return: Either the value of the key "count" in the parameter dictionary, or an error is raised.
    '''
    if "count" not in params:
        raise ValueError("The key 'count' must be specified in the parameter file.")
    return params["count"]


def record_results(network_file,params_file,results,resultsdir,datetime):
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


def is_FC(annotation):
    '''
    Specifies whether a Morse set is a full cycle.
    :param annotation: DSGRN annotation string
    :return: True or False
    '''
    return annotation.startswith("FC")


def search_over_networks(count,N,enum_network):
    '''
    Work function for parallelization.
    :param count: True or False, count DSGRN parameters or shortcut to existence
    :param N: Size of the DSGRN parameter graph
    :param enum_network: An (integer, DSGRN network specification) pair
    :return: (DSGRN network specification, results) pair
    '''
    k,netspec = enum_network
    numparams = 0
    network = DSGRN.Network(netspec)
    parametergraph = DSGRN.ParameterGraph(network)
    for p in range(parametergraph.size()):
        parameter = parametergraph.parameter(p)
        dg = DSGRN.DomainGraph(parameter)
        mg = DSGRN.MorseGraph(dg)
        stable_FC_annotations = [mg.annotation(i)[0] for i in range(0, mg.poset().size())
                                 if is_FC(mg.annotation(i)[0]) and len(mg.poset().children(i)) == 0]
        if count and len(stable_FC_annotations) > 0:
            numparams+=1
        elif len(stable_FC_annotations) > 0:
            print("Network {} of {} complete".format(k + 1, N))
            sys.stdout.flush()
            return netspec,(True,parametergraph.size())
    print("Network {} of {} complete".format(k+1, N))
    sys.stdout.flush()
    if count:
        return netspec,(numparams,parametergraph.size())
    else:
        return netspec, (False, parametergraph.size())


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Calling signature has one required argument \n " \
            "mpiexec -n <num_processes> python CountStableFC.py <path_to_network_file>"
        )
        exit(1)
    network_file = sys.argv[1]
    if len(sys.argv) > 2:
        params_file = sys.argv[2]
    else:
        params_file = ""
    if len(sys.argv) > 3:
        resultsdir = sys.argv[3]
    else:
       resultsdir = ""
    query(network_file, params_file, resultsdir)
