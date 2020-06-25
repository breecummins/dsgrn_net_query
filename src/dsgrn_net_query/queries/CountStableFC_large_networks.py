import DSGRN
import os, json, sys,subprocess,ast,shutil


def query(network_file,params_file,resultsdir=""):
    '''
    :param network_file: a .txt file containing either a single DSGRN network specification or a list of network specification strings in DSGRN format
    :param params_file: A json file with the keys
            "num_proc" = number of processes to use for each database creation
            "count" = True or False (true or false in .json format)
                        whether or not to return the number of matches (True) or just whether or not there is at least one match (False)
            "datetime" : optional datetime string to append to subdirectories in resultsdir, default = system time
    :param resultsdir: optional path to directory where results will be written, default is current directory

    :return: Writes a .json file containing a dictionary keyed by DSGRN network specification with a list of results.
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
        num_proc, count = sanity_check(params)
        results = {}
        for k,netspec in enumerate(networks):
            netfile = "temp{}.txt".format(k)
            dbfile = "temp{}.db".format(k)
            if os.path.exists(dbfile):
                os.remove(dbfile)
            with open(netfile,"w") as f:
                f.write(netspec)
            subprocess.check_call("mpiexec -n {} Signatures {} {}".format(num_proc,netfile,dbfile),shell=True)
            db = DSGRN.Database(dbfile)
            N = db.parametergraph.size()
            matches = len(DSGRN.StableFCQuery(db).matches())
            if count:
                results[netspec] = (matches,N)
            else:
                results[netspec] = (matches > 0, N)
            subprocess.call(["rm",netfile])
            subprocess.call(["rm",dbfile])
            print("Network {} of {} complete".format(k + 1, len(networks)))
            sys.stdout.flush()
        record_results(network_file,params_file,results,resultsdir,datetime)


def sanity_check(params):
    '''
    Checks to be sure the correct keys are in the dictionary params.
    :param params: dictionary
    :return: Either the values of the keys "num_proc" and "count" in the parameter dictionary, or an error is raised.
    '''
    if "num_proc" not in params or "count" not in params:
        raise ValueError("The keys 'num_proc' and 'count' must be specified in the parameter file.")
    return params["num_proc"],params["count"]


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


def read_networks(network_file):
    '''
    NOTE: Forced to copy from file_utilities due to collision between import of MPI and the mpiexec call inside this file.

    Read a .txt network file that has either a single DSGRN network specification or a list of them
    :param networks: A .txt file containing a single DSGRN network specification or a list of network specifications,
    :return: list of DSGRN network specifications
    '''

    network_str = open(network_file).read()
    if not network_str:
        networks = []
    elif network_str[0] == "[":
        networks = ast.literal_eval(network_str)
    else:
        while network_str[-1] == '\n':
            network_str = network_str[:-1]
        networks = [network_str]
    return networks


def create_results_folder(network_file, params_file, resultsdir,datetime):
    '''
    NOTE: Forced to copy from file_utilities due to collision between import of MPI and the mpiexec call inside this file.

    Create a date-time stamped folder to save results. Copy over input files.
    :param network_file: A .txt file
    :param params_file: A .json file
    :param resultsdir: optional path to directory where results will be written
    :return: string containing path to date-time stamped directory to save results file
    '''
    if datetime is None:
        datetime = subprocess.check_output(['date +%Y_%m_%d_%H_%M_%S'], shell=True).decode(sys.stdout.encoding).strip()
    dirname = os.path.join(os.path.expanduser(resultsdir), "dsgrn_net_query_results" + datetime)
    queriesdir = os.path.join(dirname, "queries" + datetime)
    os.makedirs(queriesdir)
    sys.stdout.flush()
    inputfilesdir = os.path.join(dirname, "inputs" + datetime)
    os.makedirs(inputfilesdir)
    # save input files to computations folder
    shutil.copy(network_file, inputfilesdir)
    if params_file:
        shutil.copy(params_file, inputfilesdir)
    return queriesdir


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Calling signature has two required arguments \n " \
            "python CountStableFC_large_networks.py <path_to_network_file> <path_to_parameter_file>"
        )
        exit(1)
    network_file = sys.argv[1]
    params_file = sys.argv[2]
    if len(sys.argv) > 3:
        resultsdir = sys.argv[3]
        query(network_file, params_file, resultsdir)
    else:
        query(network_file, params_file)
