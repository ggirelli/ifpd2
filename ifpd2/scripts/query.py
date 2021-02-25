"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from ifpd2.asserts import enable_rich_assert
from ifpd2.walker import Walker
from ifpd2.logging import add_log_file_handler
from ifpd2.probe import OligoProbeBuilder
from ifpd2.probe_set import OligoProbeSetBuilder
import logging
import numpy as np  # type: ignore
import os


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description="""
    Lorem ipsum dolor sit amet, consectetur adipisicing elit. Soluta, aspernatur,
    natus. Possimus recusandae distinctio, voluptatem fuga delectus laudantium ut,
    inventore culpa sit amet ullam officiis, tenetur nobis eius vitae dolore.
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Check integrity of a database.",
    )

    parser.add_argument(
        "database", metavar="database", type=str, help="Path to database folder."
    )
    parser.add_argument(
        "chrom",
        type=str,
        help="Database feature to query for a probe.",
    )

    query = parser.add_argument_group("query arguments")
    query.add_argument(
        "--region",
        metavar=("chromStart", "chromEnd"),
        type=int,
        nargs=2,
        default=(0, np.inf),
        help="""Start and end locations (space-separated) of the region of interest.
        When a region is not provided (or start/end coincide),
        the whole feature is queried.""",
    )
    query.add_argument(
        "-X",
        metavar="probes",
        type=int,
        default=1,
        help="""Number of probes to query for. Default: 1.""",
    )
    query.add_argument(
        "-N",
        metavar="oligos",
        type=int,
        default=48,
        help="""Number of oligos per probe. Default: 48.""",
    )
    query.add_argument(
        "-D",
        metavar="dist",
        type=float,
        default=2,
        help="""Minimum distance between consecutive oligos. Default: 2""",
    )
    query.add_argument(
        "--single",
        action="store_const",
        dest="single",
        const=True,
        default=False,
        help="""Useful flag for
        designing single probes. Same as "-X 1 -w 1". It overrides -X and -w, and
        deactivate -W, if specified.""",
    )

    window = parser.add_argument_group("window arguments")
    window.add_argument(
        "-W", metavar="size", type=int, default=None, help="""Window size."""
    )
    window.add_argument(
        "-w",
        metavar="shift",
        type=float,
        default=0.1,
        help="""Window shift as a percentage of window size. Default: 0.1""",
    )

    focus = parser.add_argument_group("focus region arguments")
    focus.add_argument(
        "-R",
        metavar="size",
        type=float,
        default=8000,
        help="""Central focus region size, either in nt or as fraction of window
        size. Default: 8000""",
    )
    focus.add_argument(
        "-r",
        metavar="step",
        type=float,
        default=1000,
        help="""Central focus region step, either in nt or as fraction of central
        focus region size. Default: 1000""",
    )

    filters = parser.add_argument_group("filter arguments")
    filters.add_argument(
        "-F",
        metavar="nOT",
        type=int,
        nargs=2,
        default=(0, 100),
        help="""Acceptable range of off-targets. Default: (0,100)""",
    )
    filters.add_argument(
        "-G",
        metavar="dG",
        type=float,
        nargs=2,
        default=(0.0, 0.5),
        help="""Acceptable range of secondary structure delta free energy. Either
        as absolute kcal/mol os as a fraction of delta free energy of hybridization.
        Default: (.0,.5)""",
    )
    filters.add_argument(
        "-o",
        metavar="step",
        type=float,
        default=0.1,
        help="""Step for oligo score relaxation. Default: .1""",
    )
    filters.add_argument(
        "-t",
        metavar="temp",
        type=float,
        default=10.0,
        help="""Melting temperature range half-width. Default: 10.""",
    )
    filters.add_argument(
        "-P",
        metavar="size",
        type=float,
        default=10000,
        help="""Maximum probe size in nt. Default: 10000""",
    )
    filters.add_argument(
        "-H",
        metavar="size",
        type=float,
        default=0.1,
        help="""Maximum hole size in a probe as fraction of probe size.
        Default: .1""",
    )

    probe_set = parser.add_argument_group("probe set arguments")
    probe_set.add_argument(
        "-I",
        metavar="oligos",
        type=float,
        default=0.5,
        help="""Probe oligo intersection threshold as shared oligo fraction.
        Default: .5""",
    )

    advanced = parser.add_argument_group("advanced arguments")
    advanced.add_argument(
        "-k",
        metavar="length",
        type=int,
        help="""Oligo length in nt. Use this when all oligos in the
        database have the same length, for optimization purposes.""",
    )
    advanced.add_argument(
        "-O",
        metavar="outpath",
        type=str,
        default=".",
        help="Path to output directory. Default to current folder.",
    )
    advanced.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for parallelization. Default: 1",
    )
    advanced.add_argument(
        "--reuse",
        action="store_const",
        dest="reuse",
        const=True,
        default=False,
        help="Avoid overwriting previous database walk, load instead.",
    )

    parser.set_defaults(parse=parse_arguments, run=run)

    return parser


def assert_reusable(args):
    assert_msg = " ".join(
        [
            "output path should NOT direct towards an existing directory",
            "or file. Use '--reuse' to load previous results.",
            f"Provided path '{args.O}'",
        ]
    )
    assert not os.path.isfile(args.O), assert_msg + " leads to a file"
    if not args.reuse:
        assert not os.path.isdir(args.O), assert_msg + " leads to a directory."


@enable_rich_assert
def parse_arguments(args: argparse.Namespace) -> argparse.Namespace:
    assert not all([isinstance(args.X, type(None)), isinstance(args.W, type(None))])
    args.start, args.end = args.region
    args.start = int(args.start)
    if not np.isfinite(args.end):
        args.X = None
        logging.info(
            " ".join(
                [
                    "Cannot design a specific number of probes for queries on",
                    "a whole chromosome/feature.",
                ]
            )
        )
    args.end = int(args.end) if np.isfinite(args.end) else 0

    if args.R > 1:
        args.R = int(args.R)
    if args.r > 1:
        args.r = int(args.r)

    if args.single:
        args.X = 1
        args.w = 1.0
        args.W = None

    assert_reusable(args)

    return args


@enable_rich_assert
def run(args: argparse.Namespace) -> None:

    if not os.path.isdir(args.O):
        os.mkdir(args.O)
    add_log_file_handler("{0}/{1}.log".format(args.O, "ifpd2-main"))

    opb = OligoProbeBuilder()
    opb.N = args.N
    opb.D = args.D
    opb.Tr = args.t
    opb.Ps = args.P
    opb.Ph = args.H
    opb.Po = args.I
    opb.k = args.k
    opb.F = args.F
    opb.Gs = (float(args.G[0]), float(args.G[1]))
    opb.Ot = args.o

    opsb = OligoProbeSetBuilder(os.path.join(args.O, "probe_sets"))
    logging.info(opb.get_prologue())

    ow = Walker(args.database)
    ow.C = args.chrom
    ow.S = args.start
    ow.E = args.end
    ow.X = args.X
    ow.Ws = args.W
    ow.Wh = args.w
    ow.Rs = args.R
    ow.Rt = args.r
    ow.out_path = args.O
    ow.reuse = args.reuse
    ow.threads = args.threads

    ow.start(N=args.N, opb=opb, cfr_step=ow.Rt)

    opsb.build(ow.walk_results)
    opsb.export()

    logging.shutdown()
