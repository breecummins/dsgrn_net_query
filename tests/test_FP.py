import subprocess,json,os,time
from pathlib import Path


def test_count_stableFP():
    Path("temp_results").mkdir(exist_ok=True)
    command = " ".join(["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountFPMatch.py", "mpi_networks_FP.txt", "mpi_params_FP.json", "temp_results",">dsgrn_net_query.log","2>&1"])
    os.system(command)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file = os.path.join(qdir,"query_results.json")
    results = json.load(open(output_file))
    assert(results == {'X1 : (X1)(~X3) : E\nX2 : (~X1) : E\nX3 : (X1 + X2) : E\n': [16, 168], 'X1 : (X1)(~X3) : E\nX2 : (X1) : E\nX3 : (X1 + X2) : E\n': [8, 168], 'X1 : (X1 + X2) : E\nX2 : (~X3) : E\nX3 : (X2) : E\n': [0, 4]})
    subprocess.call(["rm","-r", "temp_results/"])
    time.sleep(1)


def test_stableFP2():
    Path("temp_results").mkdir(exist_ok=True)
    command = " ".join(["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountFPMatch.py", "mpi_networks_FP.txt", "mpi_params_FP2.json", "temp_results",">dsgrn_net_query.log","2>&1"])
    os.system(command)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file = os.path.join(qdir,"query_results.json")
    results = json.load(open(output_file))
    assert(results == {'X1 : (X1)(~X3) : E\nX2 : (~X1) : E\nX3 : (X1 + X2) : E\n': [True, 168], 'X1 : (X1)(~X3) : E\nX2 : (X1) : E\nX3 : (X1 + X2) : E\n': [True, 168], 'X1 : (X1 + X2) : E\nX2 : (~X3) : E\nX3 : (X2) : E\n': [False, 4]})
    subprocess.call(["rm","-r", "temp_results/"])
    time.sleep(1)


def test_multistability():
    Path("temp_results").mkdir(exist_ok=True)
    command = " ".join(["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountFPMatch.py", "mpi_networks_ME.txt", "mpi_params_ME.json", "temp_results",">dsgrn_net_query.log","2>&1"])
    subprocess.check_call(command,shell=True)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file = os.path.join(qdir,"query_results.json")
    results = json.load(open(output_file))
    assert(results == {"inducer1 : (inducer1) : E\ninducer2 : (inducer2) : E\nreporter : (~x3) : E\nx1 : (x2)(~inducer1) : E\nx2 : (inducer2) : E\nx3 : (x1 + x3) : E\nD : (D + reporter) : E": [48, 224]})
    subprocess.call(["rm","-r", "temp_results/"])
    time.sleep(1)



def test_multistability2():
    Path("temp_results").mkdir(exist_ok=True)
    command = " ".join(["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountFPMatch.py", "mpi_networks_ME.txt", "mpi_params_ME2.json", "temp_results",">dsgrn_net_query.log","2>&1"])
    os.system(command)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file = os.path.join(qdir,"query_results.json")
    results = json.load(open(output_file))
    assert(results == {"inducer1 : (inducer1) : E\ninducer2 : (inducer2) : E\nreporter : (~x3) : E\nx1 : (x2)(~inducer1) : E\nx2 : (inducer2) : E\nx3 : (x1 + x3) : E\nD : (D + reporter) : E": [True, 224]})
    subprocess.call(["rm","-r", "temp_results/"])


if __name__ == "__main__":
    test_multistability2()