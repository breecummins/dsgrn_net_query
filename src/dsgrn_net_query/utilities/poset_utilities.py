import DSGRN
from min_interval_posets.curve import Curve
from min_interval_posets.posets import eps_posets
from dsgrn_net_query.utilities.file_utilities import readcol,readrow


def calculate_poset(params, networks):
    '''
    Construct partial orders of extrema from time series files
    :param params: A dictionary with the keys
                    'timeseriesfname' : a .csv or .tsv file name
                    'tsfile_is_row_format' : True if time series are in rows, False if they are in columns.
                    "epsilons" : a list of noise values between 0.0 and 0.5

    :param networks: a list of DSGRN network specifications
    :return: (1) A dictionary of partially ordered sets keyed by sets of node names in network files. Any input network
            that has a node name that is not in the time series file will not be analyzed. (2) The list of pruned networks.
            (3) The set of names that were missing from the time series files.
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
            curves = [Curve(data[name], times, True) for name in names]
            pos = eps_posets(dict(zip(names, curves)), params["epsilons"])
            if pos is None:
                raise ValueError("poset is None!")
            posets[names] = pos
        new_networks.append(networkspec)
    return posets, new_networks, missing_names


def calculate_posets_from_multiple_time_series(params,networks):
    '''
    Calculate partially ordered sets of extrema for each epsilon in a list and for each time series in a list.

    :param params: A dictionary with the keys
                    'timeseriesfname' : a .csv or .tsv file name or a list of such file names
                    'tsfile_is_row_format' : True if time series are in rows, False if they are in columns.
                                    **NOTE** : All time series must in the same format, row or column
                    "epsilons" : a list of floating point noise values between 0.0 and 0.5
    :param networks: Either a list of network specifications or a .txt file containing a list of network specifications
    :return: A dictionary of partial orders for multiple time series and a set of networks that will be pattern matched
            across all time series.
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


def check_posets(networks,posets):
    new_networks = []
    missing_names = set([])
    for networkspec in networks:
        network = DSGRN.Network(networkspec)
        names = tuple(sorted([network.name(k) for k in range(network.size())]))
        poset_names = set([n for nodes in posets for n in nodes])
        missing_names = missing_names.union(set(name for name in names if name not in poset_names))
        if set(names).intersection(missing_names):
            continue
        new_networks.append(networkspec)
    if missing_names:
        print(
            "No data for node(s) {} in at least one partially ordered set. \nSkipping pattern matches whenever there is a missing name.\nContinuing with {} networks.".format(
                sorted(missing_names), len(new_networks)))
    return new_networks
