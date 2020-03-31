import DSGRN
import os, json,sys
from functools import partial
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from dsgrn_net_query.utilities.file_utilities import read_networks

def query(network_file,params_file,resultsdir=""):
    '''
    :param network_file: a .txt file containing either a single DSGRN network specification or a list of network specification strings in DSGRN format
    :param params_file: A json file with a dictionary containing the key "bounds". bounds is a dictionary
    of variable names common to all network specifications with a range of values
    assigned to each. Example: {"X1":[2,2],"X2":[1,1],"X3":[0,1]}. The integer ranges
    are the matching conditions for an FP. For example, if there are four variables
    X1, X2, X3, X4 in the network spec, the FP (2,1,0,*) would be a match for any
    value of *.
    :param resultsdir: optional path to directory where results will be written, default is current directory

    :return: Writes count of parameters with an FP match to a dictionary keyed by
    network spec, which is dumped to a json file.
    '''

    networks = read_networks(network_file)
    params = json.load(open(params_file))

    bounds = params["bounds"]

    work_function = partial(check_FP, bounds, len(networks))
    with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
        if executor is not None:
            print("Querying networks.")
            output=list(executor.map(work_function, enumerate(networks)))
            results = dict(output)
            record_results(results,resultsdir)


def is_FP(annotation):
    return annotation.startswith("FP")


def is_FP_match(bounds_ind, annotation):
    digits = [int(s) for s in annotation.replace(",", "").split() if s.isdigit()]
    return all(digits[k] >= bounds_ind[k][0] and digits[k] <= bounds_ind[k][1]
           for k in bounds_ind)


def check_FP(bounds,N,tup):
    (k, netspec) = tup
    count = 0
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
            count+=1
    print("Network {} of {} complete.".format(k + 1, N))
    sys.stdout.flush()
    return netspec,[count,parametergraph.size()]


def record_results(results,resultsdir):
    rname = os.path.join(resultsdir,"query_results.json")
    if os.path.exists(rname):
        os.rename(rname,rname+".old")
    json.dump(results,open(rname,'w'))


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
