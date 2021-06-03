# Signatures
from DSGRN import *
import sqlite3

def SaveDatabase(filename, data, pg):
    # print("Save Database")
    conn = sqlite3.connect(filename)
    conn.executescript("""
      create table if not exists Signatures (ParameterIndex INTEGER PRIMARY KEY, MorseGraphIndex INTEGER);
      create table if not exists MorseGraphViz (MorseGraphIndex INTEGER PRIMARY KEY, Graphviz TEXT);
      create table if not exists MorseGraphVertices (MorseGraphIndex INTEGER, Vertex INTEGER);
      create table if not exists MorseGraphEdges (MorseGraphIndex INTEGER, Source INTEGER, Target INTEGER);
      create table if not exists MorseGraphAnnotations (MorseGraphIndex INTEGER, Vertex INTEGER, Label TEXT);
      create table if not exists Network ( Name TEXT, Dimension INTEGER, Specification TEXT, Graphviz TEXT);
      """)

    # Postprocessing to give Morse Graphs indices
    morsegraphs = []
    def signatures_table(data):
        morsegraphindices = {}
        for (pi, mg) in data:
            if mg in morsegraphindices: # ideally I'd have a graph isomorphism check
                mgi = morsegraphindices[mg]
            else:
                mgi = len(morsegraphindices)
                morsegraphindices[mg] = mgi
                morsegraphs.append(mg)
            yield (pi,mgi)

    def MG(mgi):
        return MorseGraph().parse(morsegraphs[mgi])

    name = filename
    if filename[-3:] == '.db':
        name = filename[:-3]

    # print("Inserting Network table into Database", flush=True)
    conn.execute("insert into Network ( Name, Dimension, Specification, Graphviz) values (?, ?, ?, ?);", (name, pg.network().size(), pg.network().specification(), pg.network().graphviz()))

    # print("Inserting Signatures table into Database", flush=True)
    conn.executemany("insert into Signatures (ParameterIndex, MorseGraphIndex) values (?, ?);", signatures_table(data))

    # print("Inserting MorseGraphViz table into Database", flush=True)
    conn.executemany("insert into MorseGraphViz (MorseGraphIndex, Graphviz) values (?, ?);",
      ( (mgi, MG(mgi).graphviz()) for mgi in range(0, len(morsegraphs))) )

    # print("Inserting MorseGraphVertices table into Database", flush=True)
    conn.executemany("insert into MorseGraphVertices (MorseGraphIndex, Vertex) values (?, ?);",
      ( (mgi, v) for mgi in range(0, len(morsegraphs)) for v in range(0,MG(mgi).poset().size()) ))

    # print("Inserting MorseGraphEdges table into Database", flush=True)
    conn.executemany("insert into MorseGraphEdges (MorseGraphIndex, Source, Target) values (?, ?, ?);",
      ( (mgi, s, t) for mgi in range(0, len(morsegraphs)) for s in range(0,MG(mgi).poset().size()) for t in MG(mgi).poset().children(s) ))

    # print("Inserting MorseGraphAnnotations table into Database", flush=True)
    conn.executemany("insert into MorseGraphAnnotations (MorseGraphIndex, Vertex, Label) values (?, ?, ?);",
      ( (mgi, v, label) for mgi in range(0, len(morsegraphs)) for v in range(0,MG(mgi).poset().size()) for label in MG(mgi).annotation(v) ))

    # print("Indexing Database.", flush=True)
    conn.executescript("""
      create index if not exists Signatures2 on Signatures (MorseGraphIndex, ParameterIndex);
      create index if not exists MorseGraphAnnotations3 on MorseGraphAnnotations (Label, MorseGraphIndex);
      create index if not exists MorseGraphViz2 on MorseGraphViz (Graphviz, MorseGraphIndex);
      create index if not exists MorseGraphVertices1 on MorseGraphVertices (MorseGraphIndex, Vertex);
      create index if not exists MorseGraphVertices2 on MorseGraphVertices (Vertex, MorseGraphIndex);
      create index if not exists MorseGraphEdges1 on MorseGraphEdges (MorseGraphIndex);
      create index if not exists MorseGraphAnnotations1 on MorseGraphAnnotations (MorseGraphIndex);
      """)
    conn.commit()
    conn.close()


def make_db(specfile, outfile):
    gpg = ParameterGraph(Network(specfile))
    results = []
    # print("Computing Morse Graphs")
    for pi in range(gpg.size()):
        results.append((pi, MorseGraph(DomainGraph(gpg.parameter(pi))).stringify()))
    SaveDatabase(outfile, results, gpg)
