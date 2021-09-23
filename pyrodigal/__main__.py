import argparse

from . import __name__, __author__, __version__, Pyrodigal
from .tests.fasta import parse

parser = argparse.ArgumentParser(prog=__name__)
# parser.add_argument("-a", required=False, help="Write protein translations to the selected file.")
parser.add_argument("-c", required=False, action="store_true", help="Closed ends. Do not allow genes to run off edges.", default=False)
# parser.add_argument("-d", required=False, help="Write nucleotide sequences of genes to the selected file.")
# parser.add_argument("-f", required=False, help="Select output format.", choices={"gbk", "gff", "sco"}, default="gbk")
# parser.add_argument("-g", required=False, type=int, choices=_TRANSLATION_TABLES, help="Specify a translation table to use.", default=11)
parser.add_argument("-i", required=True, help="Specify FASTA input file.")
parser.add_argument("-p", required=False, help="Select procedure.", choices={"single", "meta"}, default="single")

try:
    import argcomplete
    argcomplete.autocomplete(parser)
except ImportError:
    pass

args = parser.parse_args()

pyrodigal = Pyrodigal(meta=args.p == "meta", closed=args.c)


for i, seq in enumerate(parse(args.i)):
    if args.p == "single" and i == 0:
        pyrodigal.train(seq.seq)
    for pred in pyrodigal.find_genes(seq.seq):
        print(
            seq.id,
            "{}_v{}".format(__name__, __version__),
            "CDS",
            pred.begin,
            pred.end,
            "{:.1f}".format(pred.sscore + pred.cscore),
            "+" if pred.strand > 0 else "-",
            "0",
            ";".join([pred._gene_data, pred._score_data]),
            sep="\t",
        )
