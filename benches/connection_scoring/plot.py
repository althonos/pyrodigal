import argparse
import itertools
import json
import os
import re
import math

import numpy
import matplotlib.pyplot as plt
import scipy.stats
from palettable.cartocolors.qualitative import Bold_3

plt.rcParams["svg.fonttype"] = "none"


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
parser.add_argument("-o", "--output")
parser.add_argument("-s", "--show", action="store_true")
args = parser.parse_args()


with open(args.input) as f:
    data = json.load(f)
for result in data["results"]:
    if result["backend"] is None:
        result["backend"] = "None"
    else:
        result["backend"] = result["backend"].upper()

plt.figure(1, figsize=(12, 6))

plt.subplot(1, 2, 1)
data["results"].sort(key=lambda r: (r["backend"], r["node_count"]))
for color, (backend, group) in zip(
    Bold_3.hex_colors, itertools.groupby(data["results"], key=lambda r: r["backend"])
):
    group = list(group)
    X = numpy.array([r["node_count"] for r in group])
    Y = numpy.array([r["mean"] for r in group])
    reg = scipy.stats.linregress(X, Y)
    plt.plot([ 0, max(X) ], [ reg.intercept, reg.slope*max(X) + reg.intercept ], color=color, linestyle="--", marker="")
    ci = [1.96 * r["stddev"] / math.sqrt(len(r["times"])) for r in group]
    plt.errorbar(X, Y, ci, linestyle='', marker="+", color=color, elinewidth=0.3, ecolor='black', label=f"{backend} (R={reg.rvalue:.3f})")
plt.legend()
plt.xlabel("Node count")
plt.ylabel("Time (s)")


plt.subplot(1, 2, 2)
data["results"].sort(key=lambda r: (r["backend"], r["nucleotide_count"]))
for color, (backend, group) in zip(
    Bold_3.hex_colors, itertools.groupby(data["results"], key=lambda r: r["backend"])
):
    group = list(group)
    X = numpy.array([r["nucleotide_count"] / 1_000_000 for r in group])
    Y = numpy.array([r["mean"] for r in group])
    reg = scipy.stats.linregress(X, Y)
    plt.plot([ 0, max(X) ], [ reg.intercept, reg.slope*max(X) + reg.intercept ], color=color, linestyle="--", marker="")
    ci = [1.96 * r["stddev"] / math.sqrt(len(r["times"])) for r in group]
    plt.errorbar(X, Y, ci, linestyle='', marker="+", color=color, elinewidth=0.3, ecolor='black', label=f"{backend} (R={reg.rvalue:.3f})")
plt.legend()
plt.xlabel("Nucleotide count (Mbp)")
plt.ylabel("Time (s)")


output = args.output or args.input.replace(".json", ".svg")
plt.savefig(output, transparent=True)
if args.show:
    plt.show()
