import pandas as pd
import ast, subprocess, os, shutil, sys

def extractdata(filename):
    '''
    Read a .csv or .tsv

    :param filename: A string with a .csv or .tsv file
    :return: A list of row names and a numpy array with file values
    '''
    file_type = filename.split(".")[-1]
    if file_type == "tsv":
        df = pd.read_csv(open(filename),comment="#",sep="\t")
    elif file_type == "csv":
        df = pd.read_csv(open(filename),comment="#")
    else:
        raise ValueError("File type not recognized. Require .tsv or .csv.")
    return list(df)[1:],df.values


def readrow(filename):
    '''
    Read time series data where time series occur in rows.

    :param filename: A string with a .csv or .tsv file with time points in the first row
    :return: A dictionary keying gene names to individual time series, and a 1D array of times
    '''
    times,data = extractdata(filename)
    times = [float(n) for n in times]
    names = data[:,0]
    if len(set(names)) < len(names):
        raise ValueError("Non-unique names in time series file.")
    return dict(zip(names,[data[k,1:] for k in range(data.shape[0])])), times


def readcol(filename):
    '''
    Read time series data where time series occur in columns.

    :param filename: A string with a .csv or .tsv file with time points in the first column
    :return: A dictionary keying gene names to individual time series, and a 1D array of times
    '''
    names,data = extractdata(filename)
    if len(set(names)) < len(names):
        raise ValueError("Non-unique names in time series file.")
    times = data[:,0]
    return dict(zip(names,[data[:,k] for k in range(1,data.shape[1])])), times


def read_networks(network_object):
    '''
    Identify a list of networks or generate a list of DSGRN network specifications from a .txt file, such as that produced in makejobs.Job.run().

    :param networks: Either a list of network specifications or a .txt file containing a single DSGRN network specification or a list of network specifications,
    :return: list of DSGRN network specifications
    '''

    if isinstance(network_object,list):
        networks = network_object
    elif isinstance(network_object,str):
        # read network file
        network_str = open(network_object).read()
        if not network_str:
            networks = []
        elif network_str[0] == "[":
            networks = ast.literal_eval(network_str)
        else:
            while network_str[-1] == '\n':
                network_str = network_str[:-1]
            networks = [network_str]
    return networks


def create_results_folder(network_file, params_file, resultsdir, datetime=None):
    '''
    Create a date-time stamped folder to save results. Copy over input files.
    :param network_file: a .txt file containing either a single DSGRN network specification or a list of network specification strings in DSGRN format
    :param params_file: A .json file containing a parameter dictionary
    :param resultsdir: path to directory where results will be written
    :param datetime: optional datetime string, default = system time
    :return: string containing path to date-time stamped directory to save results file
    '''
    if datetime is None:
        datetime = subprocess.check_output(['date +%Y_%m_%d_%H_%M_%S'], shell=True).decode(sys.stdout.encoding).strip()
    dirname = os.path.join(os.path.expanduser(resultsdir), "dsgrn_net_query_results" + datetime)
    queriesdir = os.path.join(dirname, "queries" + datetime)
    os.makedirs(queriesdir,exist_ok=True)
    sys.stdout.flush()
    inputfilesdir = os.path.join(dirname, "inputs" + datetime)
    os.makedirs(inputfilesdir,exist_ok=True)
    # save input files to computations folder
    shutil.copy(network_file, inputfilesdir)
    if params_file:
        shutil.copy(params_file, inputfilesdir)
    return queriesdir

