from distutils.core import setup

setup(
    name='dsgrn_net_query',
    version='0.0.1',
    package_dir={'':'src'},
    packages = ['dsgrn_net_query',"dsgrn_net_query.queries","dsgrn_net_query.utilities"],
    install_requires=["pandas","mpi4py","progressbar2","DSGRN","min_interval_posets"],
    author="Bree Cummins",
    url='https://github.com/breecummins/NetworkPerturbations'
    )