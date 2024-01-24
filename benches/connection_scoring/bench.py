import argparse
import os
import glob
import time
import statistics
import json
import sys
import platform

import tqdm

sys.path.insert(0, os.path.realpath(os.path.join(__file__, "..", "..", "..")))

from pyrodigal import lib, Nodes, Sequence
from pyrodigal.lib import METAGENOMIC_BINS, ConnectionScorer
from pyrodigal.tests.fasta import parse


parser = argparse.ArgumentParser()
parser.add_argument("-r", "--runs", default=10, type=int)
parser.add_argument("-d", "--data", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()

BACKENDS = ["generic", None]
if lib._AVX512_RUNTIME_SUPPORT:
    BACKENDS.append("avx512")
if lib._AVX2_RUNTIME_SUPPORT:
    BACKENDS.append("avx")
if lib._MMX_RUNTIME_SUPPORT:
    BACKENDS.append("mmx")
if lib._SSE2_RUNTIME_SUPPORT:
    BACKENDS.append("sse")
if lib._NEON_RUNTIME_SUPPORT:
    BACKENDS.append("neon")


def score_connections(nodes, scorer, tinf):
    scorer.index(nodes)
    scorer.score_connections(nodes, tinf, final=True)


results = dict(results=[])
for filename in tqdm.tqdm(glob.glob(os.path.join(args.data, "*.fna"))):

    # load sequence
    with open(filename) as f:
        record = next(parse(f))
    seq = Sequence(record.seq)
    tinf = METAGENOMIC_BINS[0].training_info

    # create nodes
    nodes = Nodes()
    nodes.extract(seq, translation_table=tinf.translation_table)

    # run connection scoring
    for backend in BACKENDS:
        times = []
        for run in tqdm.tqdm(range(args.runs), desc=str(backend), leave=False):
            # initialize scorer
            scorer = ConnectionScorer(backend=backend)
            scorer_nodes = nodes.copy()
            # time how long it takes to score connections
            t1 = time.time()
            score_connections(scorer_nodes, scorer, tinf)
            t2 = time.time()
            # record runtime
            times.append(t2 - t1)
        # store benchmark result
        results["results"].append(
            {
                "sequence": os.path.basename(filename),
                "backend": backend,
                "node_count": len(nodes),
                "nucleotide_count": len(seq),
                "times": times,
                "mean": statistics.mean(times),
                "stddev": statistics.stdev(times),
                "median": statistics.median(times),
                "min": min(times),
                "max": max(times),
            }
        )


with open(args.output, "w") as f:
    json.dump(results, f, sort_keys=True, indent=4)
