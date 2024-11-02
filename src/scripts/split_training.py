import argparse
import collections
import itertools
import os
import sys

_TrainingInfo = collections.namedtuple(
    "_TrainingInfo",
    [
        "gc",
        "tt",
        "st_wt",
        "bias",
        "type_wt",
        "uses_sd",
        "rbs_wt",
        "ups_comp",
        "mot_wt",
        "no_mot",
        "gene_dc"
    ]
)

def _write_source_split(index, lines, dst):
    tinf = _TrainingInfo(*eval(
        "[" 
        + ''.join(lines[2:-4])
            .replace("{", "[")
            .replace("}", "]")
            .rstrip(";")
            .strip() 
        + "]"
    ))
    lines.clear()

    counts = collections.Counter()
    for i,j,k in itertools.product(range(4), range(4), range(4096)):
        counts[ tinf.mot_wt[i][j][k] ] += 1
    background = counts.most_common(1)[0][0]

    dst.write("#include <training.h>\n")

    dst.write(f"void initialize_metagenome_{index}(struct _training *tptr) {{\n")
    dst.write(f"    tptr->gc = {tinf.gc};\n")
    dst.write(f"    tptr->trans_table = {tinf.tt};\n")
    dst.write(f"    tptr->st_wt = {tinf.st_wt};")
    for i,b in enumerate(tinf.bias):
        dst.write(f"    tptr->bias[{i}] = {b};\n")
    for i,w in enumerate(tinf.type_wt):
        dst.write(f"    tptr->type_wt[{i}] = {w};\n")
    dst.write(f"    tptr->uses_sd = {tinf.uses_sd};\n")

    dst.write(f"    double rbs_wt[28] = {{ {', '.join(map(str, tinf.rbs_wt))} }};\n")
    dst.write(f"    memcpy(tptr->rbs_wt, rbs_wt, 28 * sizeof(double));\n")

    for i,x in enumerate(tinf.ups_comp):
        for j,y in enumerate(x):
            dst.write(f"    tptr->ups_comp[{i}][{j}] = {y};\n")

    dst.write(f"    tptr->no_mot = {tinf.no_mot};\n")
    dst.write(f"    double gene_dc[4096] = {{ {', '.join(map(str, tinf.gene_dc))} }};\n")
    dst.write(f"    memcpy(tptr->gene_dc, gene_dc, 4096 * sizeof(double));\n")

    dst.write(f"    for (size_t i = 0; i < 4; i++)\n")
    dst.write(f"        for (size_t j = 0; j < 4; j++)\n")
    dst.write(f"            for (size_t k = 0; k < 4096; k++)\n")
    dst.write(f"                tptr->mot_wt[i][j][k] = {background};\n")

    for i,j,k in itertools.product(range(4), range(4), range(4096)):
        wt = tinf.mot_wt[i][j][k]
        if wt != background:
            dst.write(f"    tptr->mot_wt[{i}][{j}][{k}] = {wt};\n")
    dst.write("}\n")

def _write_header(lines, dst):
    dst.writelines(lines)
    lines.clear()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--index", required=True, type=int)
    args = parser.parse_args()

    if os.path.exists(args.output):
        mtime_input = os.stat(args.input)
        mtime_output = os.stat(args.output)
        if mtime_output > mtime_input:
            sys.exit(0)

    with open(args.input, "r") as src, open(args.output, "w") as dst:
        lines = []
        index = 0
        for line in src:
            if line.startswith("void initialize_metagenome"):
                if line.startswith("void initialize_metagenome_0"):
                    lines.clear()
                elif index == args.index:
                    _write_source_split(index, lines, dst)
                    index += 1
                else:
                    lines.clear()
                    index += 1
            lines.append(line)
        if index == args.index:
            _write_source_split(index, lines, dst)
