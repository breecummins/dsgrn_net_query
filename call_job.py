import sys,subprocess,os,shutil


def create_results_folder(network_file, resultsdir, params_file):
    datetime = subprocess.check_output(['date +%Y_%m_%d_%H_%M_%S'], shell=True).decode(sys.stdout.encoding).strip()
    dirname = os.path.join(os.path.expanduser(resultsdir), "computations" + datetime)
    queriesdir = os.path.join(dirname, "queries" + datetime)
    os.makedirs(queriesdir)
    inputfilesdir = os.path.join(dirname, "inputs" + datetime)
    os.makedirs(inputfilesdir)
    # save input files to computations folder
    shutil.copy(network_file, inputfilesdir)
    if params_file:
        shutil.copy(params_file, inputfilesdir)
    return queriesdir


helpstring = "Calling signature has four required arguments \n " \
            "python call_job.py <num_processes> <path_to_query_module.py> <path_to_network_file> <path_to_parameter_file> <optional_path_to_resultsdir>"

if len(sys.argv) < 5:
    print(helpstring)
    exit(1)

#TODO: CountStableFC_large_networks.py does not need num_proc because it's in the parameter file
num_proc = sys.argv[1]
query = sys.argv[2]
network = sys.argv[3]
params = sys.argv[4]

if len(sys.argv) > 5:
    resultsdir = sys.argv[5]
else:
    resultsdir = ""

resultsdir = create_results_folder(network,resultsdir,params)

if "CountStableFC_large_networks.py" not in query:
    command = ["mpiexec", "-n", num_proc, "python", query, network, params, resultsdir]
else:
    command = ["python", query, network, params, resultsdir]

subprocess.check_call(command)




