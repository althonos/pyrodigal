# coding: utf-8

"""The Pyrodigal CLI, emulating the original Prodigal command line.
"""

import argparse
import contextlib
import sys

from . import __name__, __author__, __version__, Pyrodigal
from ._pyrodigal import _TRANSLATION_TABLES
from .tests.fasta import parse

def argument_parser():
    parser = argparse.ArgumentParser(prog=__name__, add_help=False)
    parser.add_argument("-a", required=False, metavar="trans_file", help="Write protein translations to the selected file.")
    parser.add_argument("-c", required=False, action="store_true", help="Closed ends. Do not allow genes to run off edges.", default=False)
    parser.add_argument("-d", required=False, metavar="nuc_file", help="Write nucleotide sequences of genes to the selected file.")
    parser.add_argument("-f", required=False, metavar="output_type", help="Select output format.", choices={"gff"}, default="gff")
    parser.add_argument("-g", required=False, metavar="tr_table", type=int, choices=_TRANSLATION_TABLES, help="Specify a translation table to use.", default=11)
    parser.add_argument("-i", metavar="input_file", required=True, help="Specify FASTA input file.")
    parser.add_argument("-p", required=False, metavar="mode", help="Select procedure.", choices={"single", "meta"}, default="single")
    parser.add_argument("-V", "--version", help="Show version number and exit.", action="version", version="{} v{}".format(__name__, __version__))
    parser.add_argument("-o", metavar="output_file", required=False, help="Specify output file.")
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    return parser

def main(argv=None, stdout=sys.stdout, stderr=sys.stderr):
    parser = argument_parser()
    args = parser.parse_args(argv)

    with contextlib.ExitStack() as ctx:
        # open output files if required
        nuc_file = None if args.d is None else ctx.enter_context(open(args.d, "w"))
        prot_file = None if args.a is None else ctx.enter_context(open(args.a, "w"))
        out_file = stdout if args.o is None else ctx.enter_context(open(args.o, "w"))
        # initialize the ORF finder
        pyrodigal = Pyrodigal(meta=args.p == "meta", closed=args.c)
        # find genes
        for i, seq in enumerate(parse(args.i)):
            # train if not in meta mode and encountering the first sequence
            if args.p == "single" and i == 0:
                pyrodigal.train(seq.seq, translation_table=args.g)
            # find genes with Pyrodigal
            preds = pyrodigal.find_genes(seq.seq)
            # write output in GFF format
            if args.f == "gff":
                preds.write_gff(out_file, prefix="{}_".format(seq.id), tool="{}_v{}".format(__name__, __version__))
            # if asked, write nucleotide sequences of genes
            if nuc_file is not None:
                preds.write_genes(nuc_file, prefix="{}_".format(seq.id))
            # if asked, write amino acide sequences of proteins
            if prot_file is not None:
                preds.write_translations(prot_file, prefix="{}_".format(seq.id))

    return 0
