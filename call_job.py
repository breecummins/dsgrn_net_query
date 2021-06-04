import sys,os,json,datetime


helpstring = "Calling signature has four required arguments \n " \
            "python call_job.py <num_processes> <query_module.py> <path_to_network_file> <path_to_parameter_file> <optional_path_to_resultsdir>"

if len(sys.argv) < 5:
    print(helpstring)
    exit(1)

num_proc = sys.argv[1]
query = sys.argv[2]
network_file = sys.argv[3]
param_file = sys.argv[4]

if len(sys.argv) > 5:
    resultsdir = sys.argv[5]
else:
    resultsdir = ""

param_dict = json.load(open(param_file))
if "datetime" not in param_dict:
    datetimestr = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
    param_dict["datetime"] = datetimestr
    json.dump(param_dict,open(param_file,"w"))
else:
    datetimestr = param_dict["datetime"]


if "CountStableFC_large_networks.py" not in query:
    command = " ".join(["mpiexec", "-n", num_proc, "python", "src/dsgrn_net_query/queries/{}".format(query), network_file, param_file, resultsdir,">dsgrn_net_query{}.log".format(datetimestr),"2>&1"])
else:
    # overwrite number of processes in parameter file with the commandline argument
    params = json.load(open(param_file))
    params["num_proc"] = num_proc
    json.dump(params,open(param_file,"w"))
    command = " ".join(["python", "src/dsgrn_net_query/queries/{}".format(query), network_file, param_file, resultsdir,">dsgrn_net_query{}.log".format(datetimestr),"2>&1"])

os.system(command)



