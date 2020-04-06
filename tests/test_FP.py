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


def test_multistability():
    command = " ".join(["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountFPMatch.py", "mpi_networks_ME.txt", "mpi_params_ME.json", "temp_results",">dsgrn_net_query.log","2>&1"])
    subprocess.check_call(command,shell=True)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file = os.path.join(qdir,"query_results.json")
    results = json.load(open(output_file))
    assert(results == {'inducer1 :  : E\ninducer2 :  : E\nreceiver1 : (inducer1 + reporter) : E\nreceiver2 : (inducer1 + inducer2)(~receiver1) : E\nreporter : (inducer1 + inducer2)(~receiver2) : E\n': [0, 23328], 'inducer1 : (inducer1) : E\ninducer2 : (inducer2) : E\nreporter : (x3) : E\nx1 : (~inducer1)(~x2) : E\nx2 : (x3)(~inducer2) : E\nx3 : (x2) : E\n': [32, 224]})



def test_multistability2():
    command = " ".join(["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountFPMatch.py", "mpi_networks_ME.txt", "mpi_params_ME2.json", "temp_results",">dsgrn_net_query.log","2>&1"])
    os.system(command)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file = os.path.join(qdir,"query_results.json")
    results = json.load(open(output_file))
    assert(results == {'inducer1 :  : E\ninducer2 :  : E\nreceiver1 : (inducer1 + reporter) : E\nreceiver2 : (inducer1 + inducer2)(~receiver1) : E\nreporter : (inducer1 + inducer2)(~receiver2) : E\n': [False, 23328], 'inducer1 : (inducer1) : E\ninducer2 : (inducer2) : E\nreporter : (x3) : E\nx1 : (~inducer1)(~x2) : E\nx2 : (x3)(~inducer2) : E\nx3 : (x2) : E\n': [True, 224]})
    subprocess.call(["rm","-r", "temp_results/"])
