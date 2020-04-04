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
    :param resultsdir: optional path to directory where results will be written, default is current directory

    :return: Writes count of parameters with a stable FC to a dictionary keyed by
    network spec, which is dumped to a json file.
    '''

    networks = read_networks(network_file)
    params = json.load(open(params_file))
    count = params["count"]

    work_function = partial(check_FC, count, len(networks))
    with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
        if executor is not None:
            print("Querying networks.")
            output=list(executor.map(work_function, enumerate(networks)))
            results = dict(output)
            record_results(network_file,params_file,results,resultsdir)


def record_results(network_file,params_file,results,resultsdir):
    resultsdir = create_results_folder(network_file, params_file, resultsdir)
    rname = os.path.join(resultsdir,"query_results.json")
    if os.path.exists(rname):
        os.rename(rname,rname+".old")
    json.dump(results,open(rname,'w'))
    print(resultsdir)


def is_FC(annotation):
    return annotation.startswith("FC")


def check_FC(count,N,tup):
    k,netspec = tup
    numparams = 0
    network = DSGRN.Network(netspec)
    parametergraph = DSGRN.ParameterGraph(network)
    for p in range(parametergraph.size()):
        parameter = parametergraph.parameter(p)
        dg = DSGRN.DomainGraph(parameter)
        md = DSGRN.MorseDecomposition(dg.digraph())
        mg = DSGRN.MorseGraph(dg, md)
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
