"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from ifpd2.asserts import enable_rich_assert
from ifpd2 import const
from ifpd2.database2 import DataBase
from ifpd2.scripts import arguments as ap
from tqdm import tqdm  # type: ignore


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "dump",
        description="""
Dumps a database to a tsv.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Dump a database to tsv format.",
    )

    parser.add_argument("input", type=str, help="Path to database folder (input).")

    parser.add_argument(
        "--chrom",
        type=str,
        help="Database feature to dump.",
    )
    parser.add_argument(
        "--region-start",
        type=int,
        help="""Start location (space-separated) of the region of interest.
        When a region is not provided, the whole feature is dumped.""",
    )
    parser.add_argument(
        "--region-end",
        type=int,
        help="""End location (space-separated) of the region of interest.
        When a region is not provided, the whole feature is dumped.""",
    )

    parser = ap.add_version_option(parser)
    parser.set_defaults(parse=parse_arguments, run=run)

    return parser


@enable_rich_assert
def parse_arguments(args: argparse.Namespace) -> argparse.Namespace:
    if args.chrom is None and (
        args.region_start is not None or args.region_end is not None
    ):
        raise Exception(
            "cannot use --region-start or --region-end without --chrom. Stopped."
        )
    if args.region_start is not None and args.region_end is not None:
        assert (
            args.region_start < args.region_end
        ), "region end should be greater than region start"

    args.region = [0, -1]
    if args.region_end is not None:
        args.region[1] = args.region_end
    return args


@enable_rich_assert
def run(args: argparse.Namespace) -> None:
    DB = DataBase(args.input)

    chromosome_list = DB.chromosome_list
    if args.chrom is not None:
        assert args.chrom.encode() in chromosome_list, f"'{args.chrom}' not found"
        chromosome_list = [args.chrom.encode()]

    if args.region_start is not None:
        chrom_size_nt = DB.chromosome_sizes_nt[args.chrom.encode()]
        assert args.region[0] < chrom_size_nt
        args.region[0] = args.region_start

    print("\t".join(const.database_columns))
    for chromosome in chromosome_list:
        walker = tqdm(
            DB.buffer(chromosome, args.region[0], args.region[1]),
            desc=f"dumping '{chromosome.decode()}'",
            total=DB.chromosome_recordnos[chromosome],
        )
        for record in walker:
            print(record.to_csv("\t"))
