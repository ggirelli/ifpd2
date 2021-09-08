"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from ifpd2.asserts import enable_rich_assert
from ifpd2.scripts import arguments as ap
import logging


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "reindex",
        description="""
Re-build a database's index.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Reindex a database.",
    )

    parser = ap.add_version_option(parser)
    parser.set_defaults(parse=parse_arguments, run=run)

    return parser


@enable_rich_assert
def parse_arguments(args: argparse.Namespace) -> argparse.Namespace:
    return args


@enable_rich_assert
def run(args: argparse.Namespace) -> None:
    raise NotImplementedError
