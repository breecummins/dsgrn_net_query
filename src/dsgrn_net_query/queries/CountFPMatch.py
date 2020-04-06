import DSGRN
import os, json,sys
from functools import partial
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from dsgrn_net_query.utilities.file_utilities import read_networks, create_results_folder

def query(network_file,params_file,resultsdir=""):
    '''
    :param network_file: a .txt file containing either a single DSGRN network specification or a list of network specification strings in DSGRN format
    :param params_file: A json file with a dictionary containing the keys "count" and "bounds".
        "bounds" is a dictionary of variable names common to all network specifications with a range of values
            assigned to each. Example: {"X1":[2,2],"X2":[1,1],"X3":[0,1]}. The integer ranges
            are the matching conditions for an FP. For example, if there are four variables
            X1, X2, X3, X4 in the network spec, the FP (2,1,0,*) would be a match for any value of *.
        "count" : True or False (true or false in .json format);
                whether to count all parameters with a match or shortcut at first success
    :param resultsdir: optional path to directory where results will be written, default is current directory

    :return: Writes a .json file containing a dictionary keyed by DSGRN network specification with a list of results.
        The results are DSGRN parameter count with successful matches to the fixed point bounds, or True
        (existence of at least one match) or False (no matches exist), depending on the value of the parameter "count".
        The size of the DSGRN parameter graph for the network is also recorded.
        { networkspec : [result, num DSGRN params] }.
    '''

    networks = read_networks(network_file)
    params = json.load(open(params_file))

    bounds = params["bounds"]
    count = params["count"]

    work_function = partial(check_FP, bounds, count, len(networks))
    with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
        if executor is not None:
            print("Querying networks.")
            output=list(executor.map(work_function, enumerate(networks)))
            results = dict(output)
            record_results(network_file, params_file,results,resultsdir)


def record_results(network_file, params_file,results,resultsdir):
    '''
    Record results in a .json file.
    :param network_file: The input .txt file containing the list of DSGRN network specifications.
    :param params_file: The input .json parameter file.
    :param results: The dictionary of results.
    :param resultsdir: The location to save the dictionary of results.
    :return: None. File is written.
    '''
    resultsdir = create_results_folder(network_file, params_file, resultsdir)
    rname = os.path.join(resultsdir,"query_results.json")
    if os.path.exists(rname):
        os.rename(rname,rname+".old")
    json.dump(results,open(rname,'w'))
    print(resultsdir)


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


def check_FP(bounds,count,N,enum_network):
    '''
    Work function for parallelization.
    :param bounds: User-specified fixed point bounds
    :param count: True or False, whether to count parameters or shortcut for existence only.
    :param N: Size of the DSGRN parameter graph
    :param enum_network: An (integer, DSGRN network specification) pair
    :return: (DSGRN network specification, results) pair
    '''
    (k, netspec) = enum_network
    numparams = 0
    network = DSGRN.Network(netspec)
    bounds_ind = {network.index(str(k)): bounds[k] for k in bounds}
    parametergraph = DSGRN.ParameterGraph(network)
    for p in range(parametergraph.size()):
        parameter = parametergraph.parameter(p)
        dg = DSGRN.DomainGraph(parameter)
        md = DSGRN.MorseDecomposition(dg.digraph())
        mg = DSGRN.MorseGraph(dg, md)
        stable_FP_annotations = [mg.annotation(i)[0] for i in range(0, mg.poset().size())
                                 if is_FP(mg.annotation(i)[0]) and len(mg.poset().children(i)) == 0]
        if any([is_FP_match(bounds_ind,a) for a in stable_FP_annotations]):
            if count:
                numparams+=1
            else:
                print("Network {} of {} complete.".format(k + 1, N))
                sys.stdout.flush()
                return netspec, [True, parametergraph.size()]
    print("Network {} of {} complete.".format(k + 1, N))
    sys.stdout.flush()
    if count:
        return netspec,[numparams,parametergraph.size()]
    else:
        return netspec, [False, parametergraph.size()]


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Calling signature has two required arguments \n " \
            "mpiexec -n <num_processes> python CountFPMatch.py <path_to_network_file> <path_to_parameter_file>"
        )
        exit(1)
    network_file = sys.argv[1]
    params_file = sys.argv[2]
    if len(sys.argv) > 3:
        resultsdir = sys.argv[3]
        query(network_file, params_file, resultsdir)
    else:
        query(network_file, params_file)
