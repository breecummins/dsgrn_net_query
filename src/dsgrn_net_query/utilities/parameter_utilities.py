import DSGRN
from dsgrn_utilities import get_parameter_neighbors as neighbors


def get_neighbors(ess_netspec):
    ess, noness_net_spec = neighbors.make_nonessential(ess_netspec)
    noness_pg = DSGRN.ParameterGraph(DSGRN.Network(noness_net_spec))
    ess_params, nbrs = neighbors.get_essential_parameter_neighbors(noness_pg)
    paramlist = ess_params + nbrs
    return noness_net_spec, paramlist