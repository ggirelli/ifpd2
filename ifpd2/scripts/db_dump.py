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

    parser = ap.add_version_option(parser)
    parser.set_defaults(parse=parse_arguments, run=run)

    return parser


@enable_rich_assert
def parse_arguments(args: argparse.Namespace) -> argparse.Namespace:
    return args


@enable_rich_assert
def run(args: argparse.Namespace) -> None:
    DB = DataBase(args.input)
    print("\t".join(const.database_columns))
    for chromosome in DB.chromosome_list:
        for record in tqdm(
            DB.walk_chromosome(chromosome), total=DB.chromosome_recordnos[chromosome]
        ):
            print(record.to_csv("\t"))
