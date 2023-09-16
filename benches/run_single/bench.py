import argparse
import os
import glob
import time
import statistics
import json
import sys
import platform
import subprocess

import tqdm

sys.path.insert(0, os.path.realpath(os.path.join(__file__, "..", "..", "..")))

from pyrodigal import lib, GeneFinder
from pyrodigal.lib import METAGENOMIC_BINS, ConnectionScorer
from pyrodigal.tests.fasta import parse


parser = argparse.ArgumentParser()
parser.add_argument("-r", "--runs", default=5, type=int)
parser.add_argument("-d", "--data", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()

BACKENDS = ["generic", None]
if lib._AVX2_RUNTIME_SUPPORT:
    BACKENDS.append("avx")
if lib._MMX_RUNTIME_SUPPORT:
    BACKENDS.append("mmx")
if lib._SSE2_RUNTIME_SUPPORT:
    BACKENDS.append("sse")
if lib._NEON_RUNTIME_SUPPORT:
    BACKENDS.append("neon")


def run_pyrodigal(sequences, backend):
    orf_finder = GeneFinder(backend=backend)
    orf_finder.train(*sequences)
    all_genes = [orf_finder.find_genes(seq) for seq in sequences]
    return sum(len(genes.nodes) for genes in all_genes)


def run_prodigal(filename):
    proc = subprocess.run(["prodigal", "-i", filename], capture_output=True)
    proc.check_returncode()


results = dict(results=[])
for filename in tqdm.tqdm(glob.glob(os.path.join(args.data, "*.fna"))):

    # load sequences
    with open(filename) as f:
        sequences = [record.seq for record in parse(f)]

    # run training + testing for Pyrodigal
    for backend in BACKENDS:
        times = []
        for run in tqdm.tqdm(range(args.runs), desc=str(backend), leave=False):
            # time how long it takes to train & predict
            t1 = time.time()
            node_count = run_pyrodigal(sequences, backend)
            t2 = time.time()
            # record runtime
            times.append(t2 - t1)
        # store benchmark result
        results["results"].append(
            {
                "sequences": os.path.basename(filename),
                "backend": backend,
                "node_count": node_count,
                "nucleotide_count": sum(len(seq) for seq in sequences),
                "times": times,
                "mean": statistics.mean(times),
                "stddev": statistics.stdev(times),
                "median": statistics.median(times),
                "min": min(times),
                "max": max(times),
            }
        )
    
    # run training + testing with Prodigal
    times = []
    for run in tqdm.tqdm(range(args.runs), desc="prodigal", leave=False):
        t1 = time.time()
        run_prodigal(filename)
        t2 = time.time()
        # record runtime
        times.append(t2 - t1)
    # store benchmark result
    results["results"].append(
        {
            "sequences": os.path.basename(filename),
            "backend": "prodigal",
            "node_count": node_count,
            "nucleotide_count": sum(len(seq) for seq in sequences),
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
