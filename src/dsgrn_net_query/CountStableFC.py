import DSGRN
import os, json,sys
from functools import partial
from NetworkPerturbations.queries.query_utilities import read_networks
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor

def query(network_file,resultsdir,params_file=""):
    '''
    :param network_file: a .txt file containing either a single DSGRN network specification or a list of network specification strings in DSGRN format
    :param resultsdir: path to an existing directory where results file(s) will be stored
    :param params_file: An unnecessary .json file containing an empty dictionary that's here for API consistency only.

    :return: Writes count of parameters with a stable FC to a dictionary keyed by
    network spec, which is dumped to a json file.
    '''

    networks = read_networks(network_file)

    work_function = partial(check_FC, len(networks))
    with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
        if executor is not None:
            print("Querying networks.")
            output=list(executor.map(work_function, enumerate(networks)))
            results = dict(output)
            record_results(results,resultsdir)


def record_results(results,resultsdir):
    rname = os.path.join(resultsdir,"query_results.json")
    if os.path.exists(rname):
        os.rename(rname,rname+".old")
    json.dump(results,open(rname,'w'))


def is_FC(annotation):
    return annotation.startswith("FC")


def check_FC(N,tup):
    k,netspec = tup
    count = 0
    network = DSGRN.Network(netspec)
    parametergraph = DSGRN.ParameterGraph(network)
    for p in range(parametergraph.size()):
        parameter = parametergraph.parameter(p)
        dg = DSGRN.DomainGraph(parameter)
        md = DSGRN.MorseDecomposition(dg.digraph())
        mg = DSGRN.MorseGraph(dg, md)
        stable_FC_annotations = [mg.annotation(i)[0] for i in range(0, mg.poset().size())
                                 if is_FC(mg.annotation(i)[0]) and len(mg.poset().children(i)) == 0]
        if len(stable_FC_annotations) > 0:
            count+=1
    print("Network {} of {} complete".format(k+1, N))
    sys.stdout.flush()
    return netspec,(count,parametergraph.size())


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
        "Calling signature is \n " \
        "mpiexec -n <num_processes> python CountStableFC.py <path_to_network_file> <path_to_results_directory>"
        )
        exit(1)
    network_file = sys.argv[1]
    resultsdir = sys.argv[2]
    query(network_file,resultsdir)