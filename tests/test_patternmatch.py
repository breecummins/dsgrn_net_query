import subprocess,json,os,time
from pathlib import Path

Path("temp_results").mkdir(exist_ok=True)


def test_patternmatch():
    command = " ".join(["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountPatternMatch.py", "mpi_networks_pm.txt", "mpi_params_pm.json", "temp_results",">dsgrn_net_query.log","2>&1"])
    os.system(command)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file1 = os.path.join(qdir,"query_results_domain_wt1_microarray_coregenes_lifepoints_interpol_trim.json")
    output_file2 = os.path.join(qdir,"query_results_stablefc_wt1_microarray_coregenes_lifepoints_interpol_trim.json")
    output_file3 = os.path.join(qdir,"query_results_domain_wt_rnaseq_ts.json")
    output_file4 = os.path.join(qdir,"query_results_stablefc_wt_rnaseq_ts.json")
    output_file5 = os.path.join(qdir,"query_results_domain_all.json")
    output_file6 = os.path.join(qdir,"query_results_stablefc_all.json")
    results1 = json.load(open(output_file1))
    results2 = json.load(open(output_file2))
    results3 = json.load(open(output_file3))
    results4 = json.load(open(output_file4))
    results5 = json.load(open(output_file5))
    results6 = json.load(open(output_file6))
    assert(results1 == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [[0.0, 0, 14], [0.01, 8, 14], [0.05, 5, 14]],
                        "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [[0.0, 0, 4], [0.01, 0, 4], [0.05, 2, 4]]})
    assert(results2 == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [[0.0, 0, 2, 14], [0.01, 2, 2, 14], [0.05, 1, 2, 14]],
                        "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [[0.0, 0, 0, 4], [0.01, 0, 0, 4], [0.05, 0, 0, 4]]})
    assert(results3 == {'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E': [[0.0, 0, 4], [0.01, 0, 4], [0.05, 0, 4]],
                        'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E': [[0.0, 0, 14], [0.01, 0, 14], [0.05, 0, 14]]})
    assert(results4 == {'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E': [[0.0, 0, 0, 4], [0.01, 0, 0, 4], [0.05, 0, 0, 4]],
                        'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E': [[0.0, 0, 2, 14], [0.01, 0, 2, 14], [0.05, 0, 2, 14]]})
    assert(results5==results1)
    assert(results6==results2)



def test_patternmatch2():
    command = " ".join(["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountPatternMatch.py", "mpi_networks_pm.txt", "mpi_params_pm2.json", "temp_results",">dsgrn_net_query.log","2>&1"])
    os.system(command)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file1 = os.path.join(qdir,"query_results_domain_wt1_microarray_coregenes_lifepoints_interpol_trim.json")
    output_file2 = os.path.join(qdir,"query_results_stablefc_wt1_microarray_coregenes_lifepoints_interpol_trim.json")
    output_file3 = os.path.join(qdir,"query_results_domain_wt_rnaseq_ts.json")
    output_file4 = os.path.join(qdir,"query_results_stablefc_wt_rnaseq_ts.json")
    results1 = json.load(open(output_file1))
    results2 = json.load(open(output_file2))
    results3 = json.load(open(output_file3))
    results4 = json.load(open(output_file4))
    assert(results1 == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [[0.0, False, 14], [0.01, True, 14], [0.05, True, 14]],
                        "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [[0.0, False, 4], [0.01, False, 4], [0.05, True, 4]]})
    assert(results2 == {"SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E": [[0.0, False, 14], [0.01, True, 14], [0.05,True, 14]],
                        "SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E": [[0.0, False, 4], [0.01, False, 4], [0.05, False, 4]]})
    assert(results3 == {'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E': [[0.0, False, 4], [0.01, False, 4], [0.05, False, 4]],
                        'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E': [[0.0, False, 14], [0.01, False, 14], [0.05, False, 14]]})
    assert(results4 == {'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : NDD1 : E': [[0.0, False, 4], [0.01, False, 4], [0.05, False, 4]],
                        'SWI4 : (NDD1)(~YOX1) : E\nHCM1 : SWI4 : E\nNDD1 : HCM1 : E\nYOX1 : SWI4 : E': [[0.0, False, 14], [0.01, False, 14], [0.05, False, 14]]})


def test_patternmatch3():
    command = " ".join(["mpiexec", "-n", "2", "python", "../src/dsgrn_net_query/queries/CountPatternMatch.py", "networks_stable_X1X2X3.txt", "params_patternmatch_stable_X1X2X3.json", "temp_results",">dsgrn_net_query.log","2>&1"])
    os.system(command)
    qdir = subprocess.check_output("tail -n 1 dsgrn_net_query.log",shell=True).strip().decode("utf-8")
    output_file2 = os.path.join(qdir,"query_results_stablefc_no_time_series_file.json")
    results2 = json.load(open(output_file2))
    assert (results2 == {
        'X1 : (X1)(~X3) : E\nX2 : (X1 + X3) : E\nX3 : (X1 + X2) : E\n': [[0.0, 205, 532, 2352], [0.1, 317, 532, 2352]],
        'X1 : (X1)(~X3) : E\nX2 : (X1) : E\nX3 : (X1 + X2) : E\n': [[0.0, 40, 74, 168], [0.1, 54, 74, 168]],
        'X1 : (X1)(~X3) : E\nX2 : (X3)(~X1) : E\nX3 : (X1 + X2) : E\n': [[0.0, 0, 299, 2352], [0.1, 0, 299, 2352]]})
    subprocess.call(["rm","-r", "temp_results/"])


if __name__ == "__main__":
    test_patternmatch()