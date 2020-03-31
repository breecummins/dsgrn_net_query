import subprocess,json,os
from pathlib import Path

Path("temp_results").mkdir(exist_ok=True)


def test_multistability():
    command = ["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/MultistabilityExists.py", "mpi_networks_ME.txt", "mpi_params_ME.json", "temp_results"]
    subprocess.check_call(command)
    output_file = "temp_results/query_results.json"
    results = json.load(open(output_file))
    assert(set(results) == set(["inducer1 : (inducer1) : E\ninducer2 : (inducer2) : E\nreporter : (~x3) : E\nx1 : (~inducer1)(~x3) : E\nx2 : (inducer2 + x1)(~x3) : E\nx3 : (x1 + x2) : E\n", "inducer1 : (inducer1) : E\ninducer2 : (inducer2) : E\nreporter : (x3) : E\nx1 : (inducer1 + x3) : E\nx2 : (x1)(~inducer2)(~x3) : E\nx3 : (x2 + x3) : E\n", "inducer1 : (inducer1) : E\ninducer2 : (inducer2) : E\nreporter : (x3) : E\nx1 : (~inducer1)(~x2) : E\nx2 : (x3)(~inducer2) : E\nx3 : (x2) : E\n"]))


def test_patternmatch():
    command = ["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/patternmatch.py", "mpi_networks_pm.txt", "mpi_params_pm.json", "temp_results"]
    subprocess.check_call(command)
    output_file1 = "temp_results/query_results_PathMatchInDomainGraph_wt1_microarray_coregenes_lifepoints_interpol_trim.json"
    output_file2 = "temp_results/query_results_PathMatchInStableFullCycle_wt1_microarray_coregenes_lifepoints_interpol_trim.json"
    output_file3 = "temp_results/query_results_PathMatchInDomainGraph_wt_rnaseq_ts.json"
    output_file4 = "temp_results/query_results_PathMatchInStableFullCycle_wt_rnaseq_ts.json"
    results1 = json.load(open(output_file1))
    results2 = json.load(open(output_file2))
    results3 = json.load(open(output_file3))
    results4 = json.load(open(output_file4))
    assert(results1 == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [[0.0, 0, 14], [0.01, 8, 14], [0.05, 5, 14]], "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [[0.0, 0, 4], [0.01, 0, 4], [0.05, 2, 4]]})
    assert(results2 == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [[0.0, 0, 2, 14], [0.01, 2, 2, 14], [0.05, 1, 2, 14]], "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [[0.0, 0, 0, 4], [0.01, 0, 0, 4], [0.05, 0, 0, 4]]})
    assert(results3 == {'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E': [[0.0, 0, 4], [0.01, 0, 4], [0.05, 0, 4]], 'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E': [[0.0, 0, 14], [0.01, 0, 14], [0.05, 0, 14]]})
    assert(results4 == {'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E': [[0.0, 0, 0, 4], [0.01, 0, 0, 4], [0.05, 0, 0, 4]], 'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E': [[0.0, 0, 2, 14], [0.01, 0, 2, 14], [0.05, 0, 2, 14]]})



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
    subprocess.call(["rm","-r", "temp_results/"])


if __name__ == "__main__":
    test_patternmatch()
