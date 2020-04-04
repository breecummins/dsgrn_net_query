import sys,os


helpstring = "Calling signature has four required arguments \n " \
            "python call_job.py <num_processes> <path_to_query_module.py> <path_to_network_file> <path_to_parameter_file> <optional_path_to_resultsdir>"

if len(sys.argv) < 5:
    print(helpstring)
    exit(1)

#TODO: CountStableFC_large_networks.py does not need num_proc
num_proc = sys.argv[1]
query = sys.argv[2]
network = sys.argv[3]
params = sys.argv[4]

if len(sys.argv) > 5:
    resultsdir = sys.argv[5]
else:
    resultsdir = ""


if "CountStableFC_large_networks.py" not in query:
    command = " ".join(["mpiexec", "-n", num_proc, "python", query, network, params, resultsdir,">dsgrn_net_query.log","2>&1"])
else:
    command = " ".join(["python", query, network, params, resultsdir,">dsgrn_net_query.log","2>&1"])

os.system(command)



