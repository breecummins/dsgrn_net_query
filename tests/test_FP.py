import subprocess,json,os,time
from pathlib import Path

Path("temp_results").mkdir(exist_ok=True)

def test_count_stableFP():
    command = " ".join(["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountFPMatch.py", "mpi_networks_FP.txt", "mpi_params_FP.json", "temp_results",">dsgrn_net_query.log","2>&1"])
    os.system(command)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file = os.path.join(qdir,"query_results.json")
    results = json.load(open(output_file))
    assert(results == {'X1 : (X1)(~X3) : E\nX2 : (~X1) : E\nX3 : (X1 + X2) : E\n': [16, 168], 'X1 : (X1)(~X3) : E\nX2 : (X1) : E\nX3 : (X1 + X2) : E\n': [8, 168], 'X1 : (X1 + X2) : E\nX2 : (~X3) : E\nX3 : (X2) : E\n': [0, 4]})
    time.sleep(1)


def test_stableFP2():
    command = " ".join(["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountFPMatch.py", "mpi_networks_FP.txt", "mpi_params_FP2.json", "temp_results",">dsgrn_net_query.log","2>&1"])
    os.system(command)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file = os.path.join(qdir,"query_results.json")
    results = json.load(open(output_file))
    assert(results == {'X1 : (X1)(~X3) : E\nX2 : (~X1) : E\nX3 : (X1 + X2) : E\n': [True, 168], 'X1 : (X1)(~X3) : E\nX2 : (X1) : E\nX3 : (X1 + X2) : E\n': [True, 168], 'X1 : (X1 + X2) : E\nX2 : (~X3) : E\nX3 : (X2) : E\n': [False, 4]})
    subprocess.call(["rm","-r", "temp_results/"])

