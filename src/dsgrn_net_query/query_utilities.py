import ast, DSGRN
import pandas as pd
from min_interval_posets import curve
from min_interval_posets import posets as make_posets


def satisfies_hex_constraints(param, hex_constraints):
    # param is a DSGRN.Parameter object
    #
    # hex_constraints is a dictionary of lists keying a tuple of two integers to a list
    # of hex numbers.
    # The tuple key describes the node type: (num inedges, num outedges) and the list contains the
    # allowable hex codes for this node type. In the algorithm below, only those hex codes in the list
    # are permitted for the node type.
    # Example: {(1,2) : ["C", "8"], (3,1) : ["0"]} means that any node with 1 in-edge and 2 out-edges must have
    # hex code 0x0C or 8 and any node with 3 in-edges and 1 outedge must have hex code 0.
    #
    # Node types not listed in hex_constraints are ignored.
    # Node types with different algebraic forms are not distinguished. For example, a node with input A+B+C and one
    # out-edge and a node with (A+B)(~C) and one out-edge are both called (3,1).
    # TODO: Extend API to distinguish between logic and type of regulation (activating or repressing).

    node_types = [(v, (v.numInputs(), v.numOutputs())) for v in param.logic()]
    for v,n in node_types:
        if n in hex_constraints and v.hex() not in hex_constraints[n]:
            return False
    return True


def extractdata(filename):
    '''
    Read a .csv or .tsv

    :param filename: A string with a .csv or .tsv file
    :return: A list of row names and a numpy array with file values
    '''
    file_type = filename.split(".")[-1]
    if file_type == "tsv":
        df = pd.read_csv(open(filename),comment="#",delim_whitespace=True)
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


def calculate_poset(params, networks):
    '''

    :param params: A dictionary with the keys
                    'timeseriesfname' : a .csv or .tsv file name
                    'tsfile_is_row_format' : True if time series are in rows, False if they are in columns.
                    "epsilons" : a list of noise values between 0.0 and 0.5

    :param networks:
    :return:
    '''
    data, times = readrow(params['timeseriesfname']) if params['tsfile_is_row_format'] else readcol(
        params['timeseriesfname'])
    posets = {}
    new_networks = []
    missing_names = set([])
    for networkspec in networks:
        network = DSGRN.Network(networkspec)
        names = tuple(sorted([network.name(k) for k in range(network.size())]))
        missing_names = missing_names.union(set(name for name in names if name not in data))
        if set(names).intersection(missing_names):
            continue
        if names not in posets.keys():
            curves = [curve.Curve(data[name], times, True) for name in names]
            pos = make_posets.eps_posets(dict(zip(names, curves)), params["epsilons"])
            if pos is None:
                raise ValueError("poset is None!")
            posets[names] = pos
        new_networks.append(networkspec)
    return posets, new_networks, missing_names


def calculate_posets_from_multiple_time_series(params,networks):
    '''


    :param params: A dictionary with the keys
                    'timeseriesfname' : a .csv or .tsv file name or a list of such file names
                    'tsfile_is_row_format' : True if time series are in rows, False if they are in columns.
                                    **NOTE** : All time series must in the same format, row or column
                    "epsilons" : a list of noise values between 0.0 and 0.5
    :param networks: Either a list of network specifications or a .txt file containing a list of network specifications
    :return: A dictionary of partial orders for multiple time series and a set of networks for which at least one pattern match will be performed.
    '''
    if isinstance(params['timeseriesfname'],str):
        timeseries_files = [params['timeseriesfname']]
    elif isinstance(params['timeseriesfname'],list):
        timeseries_files = params['timeseriesfname']
    else:
        raise ValueError("Input time series file names must be a string or a list of strings.")
    params_single = {"tsfile_is_row_format" : params["tsfile_is_row_format"], "epsilons" : params["epsilons"].copy()}
    new_networks = set([])
    posets = {}
    missing_names = set()
    for ts_file in timeseries_files:
        params_single["timeseriesfname"] = ts_file
        pos, nets, msn = calculate_poset(params_single,networks)
        new_networks.update(nets)
        missing_names.update(msn)
        for name,val in pos.items():
            if name not in posets:
                posets[name] = {ts_file : val}
            else:
                posets[name][ts_file] = val
    if missing_names:
        print(
            "No time series data for node(s) {} in at least one time series file. \nSkipping pattern matches whenever there is a missing time series.\nContinuing with {} networks.".format(
                missing_names, len(new_networks)))
    return posets, new_networks








