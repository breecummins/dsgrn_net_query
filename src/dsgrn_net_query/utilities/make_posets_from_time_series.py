from min_interval_posets import curve
from min_interval_posets import posets as make_posets


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

