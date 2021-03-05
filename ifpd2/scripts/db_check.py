"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from ifpd2.asserts import enable_rich_assert
from ifpd2.database2 import DataBase
from ifpd2.scripts import arguments as ap
import logging
import os
from tqdm import tqdm  # type: ignore


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "check",
        description="""
Check database integrity.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Check database integrity.",
    )

    parser.add_argument("input", type=str, help="Path to database folder (input).")

    parser = ap.add_version_option(parser)
    parser.set_defaults(parse=parse_arguments, run=run)

    return parser


@enable_rich_assert
def parse_arguments(args: argparse.Namespace) -> argparse.Namespace:
    assert os.path.isdir(args.input), f"cannot find database folder: '{args.input}'"
    assert os.path.isfile(
        os.path.join(args.input, "db.pickle")
    ), "db.pickle file is missing"
    return args


@enable_rich_assert
def run(args: argparse.Namespace) -> None:
    DB = DataBase(args.input)
    for chromosome in DB.chromosome_list:
        previous_position = -1
        for record in tqdm(
            DB.walk_chromosome(chromosome),
            desc=f"Checking sorting '{chromosome.decode()}'",
            total=DB.chromosome_recordnos[chromosome],
        ):
            assert record["start"] > previous_position
            previous_position = record["start"]
    logging.info("That's all! :smiley:")
