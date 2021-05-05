"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from ifpd2.asserts import enable_rich_assert
from ifpd2 import const
from ifpd2.database import DataBase
from ifpd2.scripts import arguments as ap
from ifpd2.walker2 import ChromosomeWalker
from tqdm import tqdm  # type: ignore
from typing import List


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


def get_chromosome_list(args: argparse.Namespace, DB: DataBase) -> List[bytes]:
    chromosome_list = DB.chromosome_list
    if args.chrom is not None:
        assert args.chrom.encode() in chromosome_list, f"'{args.chrom}' not found"
        chromosome_list = [args.chrom.encode()]
    return chromosome_list


def dump(args: argparse.Namespace, DB: DataBase) -> None:
    print("\t".join(const.database_columns))
    for chromosome in get_chromosome_list(args, DB):
        walker = ChromosomeWalker(DB, chromosome)
        for record in tqdm(
            walker.buffer(args.region[0], args.region[1]),
            desc=f"dumping '{chromosome.decode()}'",
            total=DB.chromosome_recordnos[chromosome],
        ):
            if args.region[0] > record["start"]:
                continue
            print(record.to_csv("\t"))


@enable_rich_assert
def run(args: argparse.Namespace) -> None:
    DB = DataBase(args.input)

    if args.region_start is not None:
        chrom_size_nt = DB.chromosome_sizes_nt[args.chrom.encode()]
        assert args.region_start < chrom_size_nt, " ".join(
            [
                f"requested window starts ({args.region_start})",
                f"after the chromosome end ({chrom_size_nt})",
            ]
        )
        args.region[0] = args.region_start

    dump(args, DB)
