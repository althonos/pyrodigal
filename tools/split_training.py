import argparse
import collections
import itertools
import os

TrainingInfo = collections.namedtuple(
    "TrainingInfo",
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

def write_header(outdir, index, lines):
    filename = os.path.join(outdir, "training_header.c".format(index))
    with open(filename, "wb") as dst:
        dst.writelines(lines)
    lines.clear()

def write_source_split(outdir, index, lines):
    filename = os.path.join(outdir, "training{:}.c".format(index))

    tinf = TrainingInfo(*eval(
        b"[" 
        + b''.join(lines[2:-4])
            .replace(b"{", b"[")
            .replace(b"}", b"]")
            .rstrip(b";")
            .strip() 
        + b"]"
    ))
    lines.clear()

    counts = collections.Counter()
    for i,j,k in itertools.product(range(4), range(4), range(4096)):
        counts[ tinf.mot_wt[i][j][k] ] += 1
    background = counts.most_common(1)[0][0]

    with open(filename, "w") as dst:
        dst.write('#include "training.h"\n')

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


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()

    with open(args.input, "rb") as src:
        lines = []
        index = 0
        for line in src:
            if line.startswith(b"void initialize_metagenome"):
                if line.startswith(b"void initialize_metagenome_0"):
                    write_header(args.output_dir, index, lines)
                else:
                    write_source_split(args.output_dir, index, lines)
                    index += 1
            # if line.lstrip().startswith(b"struct _training"):
            #     line = line.replace(b"struct _training", b"static const struct _training")
            lines.append(line)
        write_source_split(args.output_dir, index, lines)