# coding: utf-8

"""The Pyrodigal CLI, emulating the original Prodigal command line.
"""
import argparse
import contextlib
import io
import multiprocessing.pool
import os
import sys
import typing

from . import __name__, __author__, __version__
from .lib import TRANSLATION_TABLES, GeneFinder, TrainingInfo
from .tests.fasta import parse


_BZ2_MAGIC = b"BZh"
_GZIP_MAGIC = b"\x1f\x8b"
_XZ_MAGIC = b"\xfd7zXZ"
_LZ4_MAGIC = b"\x04\x22\x4d\x18"
_ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"

@contextlib.contextmanager
def zopen(path, mode='r', encoding=None, errors=None, newline=None) -> typing.Iterator[typing.BinaryIO]:
    with contextlib.ExitStack() as ctx:
        file = ctx.enter_context(open(path, "rb"))
        peek = file.peek()
        if peek.startswith(_GZIP_MAGIC):
            import gzip
            file = ctx.enter_context(gzip.open(file, mode="rb"))
        elif peek.startswith(_BZ2_MAGIC):
            import bz2
            file = ctx.enter_context(bz2.open(file, mode="rb"))
        elif peek.startswith(_XZ_MAGIC):
            import lzma
            file = ctx.enter_context(lzma.open(file, mode="rb"))
        elif peek.startswith(_LZ4_MAGIC):
            try:
                import lz4.frame
            except ImportError as err:
                raise RuntimeError("File compression is LZ4 but lz4 is not installed") from err
            file = ctx.enter_context(lz4.frame.open(file))
        elif peek.startswith(_ZSTD_MAGIC):
            try:
                import zstandard
            except ImportError as err:
                raise RuntimeError("File compression is ZSTD but zstandard is not installed") from err
            decompressor = zstandard.ZstdDecompressor()
            file = decompressor.stream_reader(file)
        if mode == "r":
            file = io.TextIOWrapper(file, encoding=encoding, errors=errors, newline=newline)
        yield file


def argument_parser(
    prog: str = __name__,
    version: str = __version__
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=prog, add_help=False)
    parser.add_argument(
        "-a",
        required=False,
        metavar="trans_file",
        help="Write protein translations to the selected file.",
    )
    parser.add_argument(
        "-c",
        required=False,
        action="store_true",
        help="Closed ends. Do not allow genes to run off edges.",
        default=False,
    )
    parser.add_argument(
        "-d",
        required=False,
        metavar="nuc_file",
        help="Write nucleotide sequences of genes to the selected file.",
    )
    parser.add_argument(
        "-f",
        required=False,
        metavar="output_type",
        help="Select output format.",
        choices={"gff", "gbk"},
        default="gff",
    )
    parser.add_argument(
        "-g",
        required=False,
        metavar="tr_table",
        type=int,
        choices=TRANSLATION_TABLES,
        help="Specify a translation table to use.",
        default=11,
    )
    parser.add_argument(
        "-i", metavar="input_file", required=True, help="Specify FASTA input file."
    )
    parser.add_argument(
        "-m",
        action="store_true",
        help="Treat runs of N as masked sequence; don't build genes across them.",
        default=False,
    )
    parser.add_argument(
        "-n",
        action="store_true",
        help="Bypass Shine-Dalgarno trainer and force a full motif scan.",
        default=False,
    )
    parser.add_argument(
        "-o", metavar="output_file", required=False, help="Specify output file."
    )
    parser.add_argument(
        "-p",
        required=False,
        metavar="mode",
        help="Select procedure.",
        choices={"single", "meta"},
        default="single",
    )
    parser.add_argument(
        "-s",
        required=False,
        metavar="start_file",
        help="Write all potential genes (with scores) to the selected file.",
    )
    parser.add_argument(
        "-t",
        required=False,
        metavar="training_file",
        help="Write a training file (if none exists); otherwise, read and use the specified training file.",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        required=False,
        default=1,
        metavar="jobs",
        help="The number of threads to use if input contains multiple sequences."
    )
    parser.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit."
    )
    parser.add_argument(
        "-V",
        "--version",
        help="Show version number and exit.",
        action="version",
        version="{} v{}".format(prog, version),
    )
    parser.add_argument(
        "--min-gene",
        help="The minimum gene length.",
        required=False,
        type=int,
        default=90,
    )
    parser.add_argument(
        "--min-edge-gene",
        help="The minimum edge gene length.",
        required=False,
        type=int,
        default=60,
    )
    parser.add_argument(
        "--max-overlap",
        help="The maximum number of nucleotides that can overlap between two genes on the same strand. This must be lower or equal to the minimum gene length.",
        required=False,
        type=int,
        default=60,
    )
    parser.add_argument(
        "--no-stop-codon",
        required=False,
        action="store_false",
        help="Disables translation of stop codons into star characters (*) for complete genes.",
        default=True,
    )
    return parser


def main(
    argv: typing.Optional[typing.List[str]] = None,
    stdout: typing.TextIO = sys.stdout,
    stderr: typing.TextIO = sys.stderr,
    *,
    gene_finder_factory: typing.Callable[..., GeneFinder] = GeneFinder,
    argument_parser: typing.Callable[[], argparse.ArgumentParser] = argument_parser,
) -> int:
    parser = argument_parser()
    args = parser.parse_args(argv)

    with contextlib.ExitStack() as ctx:
        try:
            # open output files if required
            nuc_file = None if args.d is None else ctx.enter_context(open(args.d, "w"))
            prot_file = None if args.a is None else ctx.enter_context(open(args.a, "w"))
            scores_file = (
                None if args.s is None else ctx.enter_context(open(args.s, "w"))
            )
            out_file = (
                stdout if args.o is None else ctx.enter_context(open(args.o, "w"))
            )

            # load training info
            if args.t is not None:
                if args.p == "meta":
                    print(
                        "Error: cannot specify metagenomic sequence with a training file.",
                        file=stderr,
                    )
                    return 1
                elif os.path.exists(args.t):
                    with open(args.t, "rb") as f:
                        training_info = TrainingInfo.load(f)
                else:
                    training_info = None
            else:
                training_info = None

            # open input (with support for compressed files)
            input_file = ctx.enter_context(zopen(args.i))

            # initialize the ORF finder
            gene_finder = gene_finder_factory(
                meta=args.p == "meta",
                closed=args.c,
                mask=args.m,
                training_info=training_info,
                min_gene=args.min_gene,
                min_edge_gene=args.min_edge_gene,
                max_overlap=args.max_overlap,
            )

            # pre-train if in training mode
            if args.p == "single":
                # use the same interleaving logic as Prodigal
                sequences = list(parse(input_file))
                training_info = gene_finder.train(
                    *(seq.seq for seq in sequences),
                    force_nonsd=args.n,
                    translation_table=args.g
                )
                # save the training info is desired
                if args.t is not None and not os.path.exists(args.t):
                    with open(args.t, "wb") as f:
                        training_info.dump(f)
            else:
                sequences = parse(input_file)

            # get the number of jobs
            if args.jobs == 0:
                args.jobs = os.cpu_count() or 1
            if args.jobs > 1:
                pool = ctx.enter_context(multiprocessing.pool.ThreadPool(args.jobs))
                parallel_map = pool.map
            else:
                parallel_map = map

            # find genes in parallel
            def process(sequence):
                return (sequence.id, gene_finder.find_genes(sequence.seq))

            for seq_id, preds in parallel_map(process, sequences):
                # write output in GFF or GBK format
                if args.f == "gff":
                    preds.write_gff(out_file, seq_id)
                elif args.f == "gbk":
                    preds.write_genbank(out_file, seq_id)
                # if asked, write nucleotide sequences of genes
                if nuc_file is not None:
                    preds.write_genes(nuc_file, seq_id)
                # if asked, write amino acid sequences of proteins
                if prot_file is not None:
                    preds.write_translations(prot_file, seq_id, include_stop=args.no_stop_codon)
                # if asked, write scores
                if scores_file is not None:
                    preds.write_scores(scores_file, seq_id)

        except Exception as err:
            print("Error: {}".format(err))
            return getattr(err, "errno", 1)
        else:
            return 0
