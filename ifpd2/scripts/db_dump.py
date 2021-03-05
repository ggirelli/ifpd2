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
        "--region",
        metavar=("chromStart", "chromEnd"),
        type=int,
        nargs=2,
        help="""Start and end locations (space-separated) of the region of interest.
        When a region is not provided, the whole feature is dumped.""",
    )

    parser = ap.add_version_option(parser)
    parser.set_defaults(parse=parse_arguments, run=run)

    return parser


@enable_rich_assert
def parse_arguments(args: argparse.Namespace) -> argparse.Namespace:
    if args.chrom is None and args.region is not None:
        raise Exception("cannot use --region without --chrom. Stopped.")
    if args.region is not None:
        assert (
            args.region[0] < args.region[1]
        ), "region end should be greater than region start"
    return args


@enable_rich_assert
def run(args: argparse.Namespace) -> None:
    DB = DataBase(args.input)

    chromosome_list = DB.chromosome_list
    if args.chrom is not None:
        assert args.chrom in chromosome_list, f"'{args.chrom}' not found"
        chromosome_list = [args.chrom]
        if args.region is not None:
            chrom_size_nt = DB.chromosome_sizes_nt[args.chrom.encode()]
            assert args.region[0] < chrom_size_nt

    print("\t".join(const.database_columns))
    for chromosome in DB.chromosome_list:
        for record in tqdm(
            DB.walk_chromosome(chromosome),
            desc=f"dumping '{chromosome.decode()}'",
            total=DB.chromosome_recordnos[chromosome],
        ):
            print(record.to_csv("\t"))
