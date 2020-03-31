import ast, DSGRN
import pandas as pd


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









