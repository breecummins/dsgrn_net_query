import subprocess,json,os
from pathlib import Path

Path("temp_results").mkdir(exist_ok=True)


def test_multistability():
    command = ["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/MultistabilityExists.py", "mpi_networks_ME.txt", "mpi_params_ME.json", "temp_results"]
    subprocess.check_call(command)
    output_file = "temp_results/query_results.json"
    results = json.load(open(output_file))
    assert(results == {'inducer1 :  : E\ninducer2 :  : E\nreceiver1 : (inducer1 + reporter) : E\nreceiver2 : (inducer1 + inducer2)(~receiver1) : E\nreporter : (inducer1 + inducer2)(~receiver2) : E\n': [0, 23328], 'inducer1 : (inducer1) : E\ninducer2 : (inducer2) : E\nreporter : (x3) : E\nx1 : (~inducer1)(~x2) : E\nx2 : (x3)(~inducer2) : E\nx3 : (x2) : E\n': [32, 224]})


def test_multistability2():
    command = ["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/MultistabilityExists.py", "mpi_networks_ME.txt", "mpi_params_ME2.json", "temp_results"]
    subprocess.check_call(command)
    output_file = "temp_results/query_results.json"
    results = json.load(open(output_file))
    assert(results == {'inducer1 :  : E\ninducer2 :  : E\nreceiver1 : (inducer1 + reporter) : E\nreceiver2 : (inducer1 + inducer2)(~receiver1) : E\nreporter : (inducer1 + inducer2)(~receiver2) : E\n': [False, 23328], 'inducer1 : (inducer1) : E\ninducer2 : (inducer2) : E\nreporter : (x3) : E\nx1 : (~inducer1)(~x2) : E\nx2 : (x3)(~inducer2) : E\nx3 : (x2) : E\n': [True, 224]})


def test_patternmatch():
    command = ["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/patternmatch.py", "mpi_networks_pm.txt", "mpi_params_pm.json", "temp_results"]
    subprocess.check_call(command)
    output_file1 = "temp_results/query_results_domain_wt1_microarray_coregenes_lifepoints_interpol_trim.json"
    output_file2 = "temp_results/query_results_stablefc_wt1_microarray_coregenes_lifepoints_interpol_trim.json"
    output_file3 = "temp_results/query_results_domain_wt_rnaseq_ts.json"
    output_file4 = "temp_results/query_results_stablefc_wt_rnaseq_ts.json"
    output_file5 = "temp_results/query_results_domain_all.json"
    output_file6 = "temp_results/query_results_stablefc_all.json"
    results1 = json.load(open(output_file1))
    results2 = json.load(open(output_file2))
    results3 = json.load(open(output_file3))
    results4 = json.load(open(output_file4))
    results5 = json.load(open(output_file5))
    results6 = json.load(open(output_file6))
    assert(results1 == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [[0.0, 0, 14], [0.01, 8, 14], [0.05, 5, 14]], "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [[0.0, 0, 4], [0.01, 0, 4], [0.05, 2, 4]]})
    assert(results2 == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [[0.0, 0, 2, 14], [0.01, 2, 2, 14], [0.05, 1, 2, 14]], "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [[0.0, 0, 0, 4], [0.01, 0, 0, 4], [0.05, 0, 0, 4]]})
    assert(results3 == {'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E': [[0.0, 0, 4], [0.01, 0, 4], [0.05, 0, 4]], 'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E': [[0.0, 0, 14], [0.01, 0, 14], [0.05, 0, 14]]})
    assert(results4 == {'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E': [[0.0, 0, 0, 4], [0.01, 0, 0, 4], [0.05, 0, 0, 4]], 'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E': [[0.0, 0, 2, 14], [0.01, 0, 2, 14], [0.05, 0, 2, 14]]})
    assert(results5==results1)
    assert(results6==results2)


def test_patternmatch2():
    command = ["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/patternmatch.py", "mpi_networks_pm.txt", "mpi_params_pm2.json", "temp_results"]
    subprocess.check_call(command)
    output_file1 = "temp_results/query_results_domain_wt1_microarray_coregenes_lifepoints_interpol_trim.json"
    output_file2 = "temp_results/query_results_stablefc_wt1_microarray_coregenes_lifepoints_interpol_trim.json"
    output_file3 = "temp_results/query_results_domain_wt_rnaseq_ts.json"
    output_file4 = "temp_results/query_results_stablefc_wt_rnaseq_ts.json"
    results1 = json.load(open(output_file1))
    results2 = json.load(open(output_file2))
    results3 = json.load(open(output_file3))
    results4 = json.load(open(output_file4))
    assert(results1 == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [[0.0, False, 14], [0.01, True, 14], [0.05, True, 14]], "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [[0.0, False, 4], [0.01, False, 4], [0.05, True, 4]]})
    assert(results2 == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [[0.0, False, 14], [0.01, True, 14], [0.05,True, 14]], "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [[0.0, False, 4], [0.01, False, 4], [0.05, False, 4]]})
    assert(results3 == {'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E': [[0.0, False, 4], [0.01, False, 4], [0.05, False, 4]], 'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E': [[0.0, False, 14], [0.01, False, 14], [0.05, False, 14]]})
    assert(results4 == {'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E': [[0.0, False, 4], [0.01, False, 4], [0.05, False, 4]], 'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E': [[0.0, False, 14], [0.01, False, 14], [0.05, False, 14]]})


def test_count_stableFCln():
    command = "python ../src/dsgrn_net_query/queries/CountStableFC_large_networks.py mpi_networks_FCln.txt mpi_params_FCln.json temp_results"
    subprocess.check_call(command,shell=True)
    output_file = "temp_results/query_results.json"
    results = json.load(open(output_file))
    assert(results == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [2, 14], "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [0, 4]})


def test_count_stableFC():
    command = ["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountStableFC.py", "mpi_networks_FCln.txt", "", "temp_results"]
    subprocess.check_call(command)
    output_file = "temp_results/query_results.json"
    results = json.load(open(output_file))
    assert(results == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [2, 14], "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [0, 4]})


def test_count_stableFP():
    command = ["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountFPMatch.py", "mpi_networks_FP.txt", "mpi_params_FP.json", "temp_results"]
    subprocess.check_call(command)
    output_file = "temp_results/query_results.json"
    results = json.load(open(output_file))
    assert(results == {'X1 : (X1)(~X3) : E\nX2 : (~X1) : E\nX3 : (X1 + X2) : E\n': [16, 168], 'X1 : (X1)(~X3) : E\nX2 : (X1) : E\nX3 : (X1 + X2) : E\n': [8, 168], 'X1 : (X1 + X2) : E\nX2 : (~X3) : E\nX3 : (X2) : E\n': [0, 4]})


def test_stableFP2():
    command = ["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountFPMatch.py", "mpi_networks_FP.txt", "mpi_params_FP2.json", "temp_results"]
    subprocess.check_call(command)
    output_file = "temp_results/query_results.json"
    results = json.load(open(output_file))
    assert(results == {'X1 : (X1)(~X3) : E\nX2 : (~X1) : E\nX3 : (X1 + X2) : E\n': [True, 168], 'X1 : (X1)(~X3) : E\nX2 : (X1) : E\nX3 : (X1 + X2) : E\n': [True, 168], 'X1 : (X1 + X2) : E\nX2 : (~X3) : E\nX3 : (X2) : E\n': [False, 4]})
    subprocess.call(["rm","-r", "temp_results/"])


if __name__ == "__main__":
    test_patternmatch()
