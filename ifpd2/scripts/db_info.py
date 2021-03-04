"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from ifpd2.asserts import enable_rich_assert
from ifpd2.scripts import arguments as ap
import os
import pickle
import logging


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "info",
        description="""
Shows database details.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Show database details.",
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
    with open(os.path.join(args.input, "db.pickle"), "rb") as IH:
        db_details = pickle.load(IH)

    logging.info(f"Database name: {db_details['args'].output}")
    logging.info(f"Sequence max length: {db_details['dtype']['sequence'][2:]}")
    logging.info("")
    logging.info("[bold]## Input files[/bold]")
    logging.info(f"hush files: {db_details['args'].hush}")
    logging.info(f"oligo-melting files: {db_details['args'].melting}")
    logging.info(f"OligoArrayAux files: {db_details['args'].secondary}")
    logging.info("")
    logging.info("[bold]## Chromosome details[/bold]")
    logging.info(f"Expecting {len(db_details['chromosomes'])} chromosomes.")
    logging.info("Chromosome sizes:")
    for chromosome, size in db_details['chromosomes'].items():
        logging.info(f"\t{chromosome.decode()} => {size}")

    for chromosome in db_details['chromosomes'].keys():
        chromosome_path = os.path.join(args.input, f"{chromosome.decode()}.bin")
        if not os.path.isfile(chromosome_path):
            logging.warning(f"missing expected chromosome file: '{chromosome_path}'")
    logging.info("")
    logging.info("That's all! :smiley:")
