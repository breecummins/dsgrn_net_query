import DSGRN
import os, json, sys,subprocess,progressbar,ast,shutil
# from dsgrn_net_query.utilities.file_utilities import read_networks,create_results_folder
# from mpi4py import MPI
# from mpi4py.futures import MPICommExecutor
import sqlite3


def query(network_file,params_file,resultsdir=""):
    '''
    :param network_file: a .txt file containing either a single DSGRN network specification or a list of network specification strings in DSGRN format
    :param params_file: A json file with an empty dictionary (here for uniform API)
    :param resultsdir: optional path to directory where results will be written, default is current directory

    :return: Writes count of parameters with a stable FC to a dictionary keyed by
    network spec, which is dumped to a json file.
    '''

    networks = read_networks(network_file)
    params = json.load(open(params_file))
    num_proc = params["num_proc"]

    resultsdict = {}
    for k,netspec in enumerate(networks):
        netfile = "temp{}.txt".format(k)
        dbfile = "temp{}.db".format(k)
        with open(netfile,"w") as f:
            f.write(netspec)
        # fromSignatures(netfile,dbfile)
        subprocess.check_call("mpiexec -n {} Signatures {} {}".format(num_proc,netfile,dbfile),shell=True)
        db = DSGRN.Database(dbfile)
        N = db.parametergraph.size()
        matches = len(DSGRN.StableFCQuery(db).matches())
        resultsdict[netspec] = (matches,N)
        subprocess.call(["rm",netfile])
        subprocess.call(["rm",dbfile])
        print("Network {} of {} complete".format(k + 1, len(networks)))
        sys.stdout.flush()
    record_results(resultsdir,resultsdict)


def record_results(resultsdir,resultsdict):
    resultsdir = create_results_folder(network_file, params_file, resultsdir)
    rname = os.path.join(resultsdir,"query_results.json")
    if os.path.exists(rname):
        os.rename(rname,rname+".old")
    json.dump(resultsdict,open(rname,'w'))
    print(resultsdir)


def read_networks(network_object):
    '''
    Forced to copy from file_utilities due to collision between import of MPI and the mpiexec call inside this file.

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


def create_results_folder(network_file, params_file, resultsdir):
    # Forced to copy from file_utilities due to collision between import of MPI and the mpiexec call inside this file.
    datetime = subprocess.check_output(['date +%Y_%m_%d_%H_%M_%S'], shell=True).decode(sys.stdout.encoding).strip()
    dirname = os.path.join(os.path.expanduser(resultsdir), "dsgrn_net_query_results" + datetime)
    queriesdir = os.path.join(dirname, "queries" + datetime)
    os.makedirs(queriesdir)
    sys.stdout.flush()
    inputfilesdir = os.path.join(dirname, "inputs" + datetime)
    os.makedirs(inputfilesdir)
    # save input files to computations folder
    shutil.copy(network_file, inputfilesdir)
    if params_file:
        shutil.copy(params_file, inputfilesdir)
    return queriesdir

# def fromSignatures(specfile,outfile):
#     '''
#     Lifted from Signatures.py because it can't be imported.
#
#     :param specfile: network specification file
#     :param outfile: output database file
#     :return: None, writes to file
#     '''
#     gpg = DSGRN.ParameterGraph(DSGRN.Network(specfile))
#
#     def work(pi):
#         return (pi, DSGRN.MorseGraph(DSGRN.DomainGraph(gpg.parameter(pi))).stringify())
#
#     # global gpg
#     with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
#         if executor is not None:
#             bar = progressbar.ProgressBar(max_value=gpg.size())
#             print("Computing Morse Graphs")
#             results = list(bar(executor.map(work, range(0, gpg.size()), chunksize=65536)))
#             SaveDatabase(outfile, results, gpg)
#
#
# def SaveDatabase(filename, data, pg):
#     '''
#     Lifted from Signatures.py because it can't be imported.
#
#     :param filename: output database name
#     :param data: output from parallel process
#     :param pg: parameter graph
#     :return: None, writes to file
#     '''
#     print("Save Database")
#     N = pg.size()
#     conn = sqlite3.connect(filename)
#     conn.executescript("""
#       create table if not exists Signatures (ParameterIndex INTEGER PRIMARY KEY, MorseGraphIndex INTEGER);
#       create table if not exists MorseGraphViz (MorseGraphIndex INTEGER PRIMARY KEY, Graphviz TEXT);
#       create table if not exists MorseGraphVertices (MorseGraphIndex INTEGER, Vertex INTEGER);
#       create table if not exists MorseGraphEdges (MorseGraphIndex INTEGER, Source INTEGER, Target INTEGER);
#       create table if not exists MorseGraphAnnotations (MorseGraphIndex INTEGER, Vertex INTEGER, Label TEXT);
#       create table if not exists Network ( Name TEXT, Dimension INTEGER, Specification TEXT, Graphviz TEXT);
#       """)
#
#     # Postprocessing to give Morse Graphs indices
#     morsegraphs = []
#
#     def signatures_table(data):
#         bar = progressbar.ProgressBar(max_value=N)
#         morsegraphindices = {}
#         for (pi, mg) in data:
#             bar.update(pi)
#             if mg in morsegraphindices:  # ideally I'd have a graph isomorphism check
#                 mgi = morsegraphindices[mg]
#             else:
#                 mgi = len(morsegraphindices)
#                 morsegraphindices[mg] = mgi
#                 morsegraphs.append(mg)
#             yield (pi, mgi)
#         bar.finish()
#
#     def MG(mgi):
#         return DSGRN.MorseGraph().parse(morsegraphs[mgi])
#
#     name = filename
#     if filename[-3:] == '.db':
#         name = filename[:-3]
#
#     print("Inserting Network table into Database", flush=True)
#     conn.execute("insert into Network ( Name, Dimension, Specification, Graphviz) values (?, ?, ?, ?);",
#                  (name, pg.network().size(), pg.network().specification(), pg.network().graphviz()))
#
#     print("Inserting Signatures table into Database", flush=True)
#     conn.executemany("insert into Signatures (ParameterIndex, MorseGraphIndex) values (?, ?);", signatures_table(data))
#
#     print("Inserting MorseGraphViz table into Database", flush=True)
#     conn.executemany("insert into MorseGraphViz (MorseGraphIndex, Graphviz) values (?, ?);",
#                      ((mgi, MG(mgi).graphviz()) for mgi in progressbar.ProgressBar()(range(0, len(morsegraphs)))))
#
#     print("Inserting MorseGraphVertices table into Database", flush=True)
#     conn.executemany("insert into MorseGraphVertices (MorseGraphIndex, Vertex) values (?, ?);",
#                      ((mgi, v) for mgi in progressbar.ProgressBar()(range(0, len(morsegraphs))) for v in
#                       range(0, MG(mgi).poset().size())))
#
#     print("Inserting MorseGraphEdges table into Database", flush=True)
#     conn.executemany("insert into MorseGraphEdges (MorseGraphIndex, Source, Target) values (?, ?, ?);",
#                      ((mgi, s, t) for mgi in progressbar.ProgressBar()(range(0, len(morsegraphs))) for s in
#                       range(0, MG(mgi).poset().size()) for t in MG(mgi).poset().children(s)))
#
#     print("Inserting MorseGraphAnnotations table into Database", flush=True)
#     conn.executemany("insert into MorseGraphAnnotations (MorseGraphIndex, Vertex, Label) values (?, ?, ?);",
#                      ((mgi, v, label) for mgi in progressbar.ProgressBar()(range(0, len(morsegraphs))) for v in
#                       range(0, MG(mgi).poset().size()) for label in MG(mgi).annotation(v)))
#
#     print("Indexing Database.", flush=True)
#     conn.executescript("""
#       create index if not exists Signatures2 on Signatures (MorseGraphIndex, ParameterIndex);
#       create index if not exists MorseGraphAnnotations3 on MorseGraphAnnotations (Label, MorseGraphIndex);
#       create index if not exists MorseGraphViz2 on MorseGraphViz (Graphviz, MorseGraphIndex);
#       create index if not exists MorseGraphVertices1 on MorseGraphVertices (MorseGraphIndex, Vertex);
#       create index if not exists MorseGraphVertices2 on MorseGraphVertices (Vertex, MorseGraphIndex);
#       create index if not exists MorseGraphEdges1 on MorseGraphEdges (MorseGraphIndex);
#       create index if not exists MorseGraphAnnotations1 on MorseGraphAnnotations (MorseGraphIndex);
#       """)
#     conn.commit()
#     conn.close()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Calling signature has two required arguments \n " \
            "python CountStableFC_large_networks.py <path_to_network_file> <path_to_parameter_file>"
        )
        exit(1)
    network_file = sys.argv[1]
    params_file = sys.argv[2]
    if len(sys.argv) > 3:
        resultsdir = sys.argv[3]
        query(network_file, params_file, resultsdir)
    else:
        query(network_file, params_file)
