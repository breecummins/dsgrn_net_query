import subprocess,json,os,time
from pathlib import Path

Path("temp_results").mkdir(exist_ok=True)

def test_count_stableFCln():
    command = "python ../src/dsgrn_net_query/queries/CountStableFC_large_networks.py mpi_networks_FCln.txt mpi_params_FCln.json temp_results >dsgrn_net_query.log 2>&1"
    os.system(command)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file = os.path.join(qdir,"query_results.json")
    results = json.load(open(output_file))
    assert(results == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [2, 14], "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [0, 4]})
    time.sleep(1)


def test_count_stableFC():
    command = " ".join(["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountStableFC.py", "mpi_networks_FCln.txt", "mpi_params_FCln.json", "temp_results",">dsgrn_net_query.log","2>&1"])
    os.system(command)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file = os.path.join(qdir,"query_results.json")
    results = json.load(open(output_file))
    assert(results == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [2, 14], "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [0, 4]})
    subprocess.call(["rm","-r", "temp_results/"])


if __name__ == "__main__":
    test_count_stableFCln()
