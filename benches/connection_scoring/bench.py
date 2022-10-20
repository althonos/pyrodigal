import argparse
import os
import glob
import time
import statistics
import json
import sys
import platform

import tqdm

sys.path.append(os.path.realpath(os.path.join(__file__, "..", "..", "..")))

from pyrodigal import _pyrodigal, Nodes, Sequence
from pyrodigal._pyrodigal import METAGENOMIC_BINS, ConnectionScorer
from pyrodigal.tests.fasta import parse


parser = argparse.ArgumentParser()
parser.add_argument("-r", "--runs", default=10, type=int)
parser.add_argument("-d", "--data", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()

BACKENDS = ["generic", None]
if _pyrodigal._AVX2_RUNTIME_SUPPORT:
    BACKENDS.append("avx")
if _pyrodigal._MMX_RUNTIME_SUPPORT:
    BACKENDS.append("mmx")
if _pyrodigal._SSE2_RUNTIME_SUPPORT:
    BACKENDS.append("sse")
if _pyrodigal._NEON_RUNTIME_SUPPORT:
    BACKENDS.append("neon")


def score_connections(nodes, scorer, tinf):
    scorer.index(nodes)
    for i in range(500, len(nodes)):
        # compute boundary
        j = 0 if i < 500 else i - 500
        # score connections without fast-indexing skippable nodes
        scorer.compute_skippable(j, i)
        scorer.score_connections(nodes, j, i, tinf, final=True)


results = dict(results=[])
for filename in tqdm.tqdm(glob.glob(os.path.join(args.data, "*.fna"))):

    # load sequence
    with open(filename) as f:
        record = next(parse(f))
    seq = Sequence.from_string(record.seq)
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
