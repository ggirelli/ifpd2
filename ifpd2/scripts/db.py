"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from ifpd2.scripts import arguments as ap
from ifpd2 import scripts
import sys


def default_parser(*args) -> None:
    print("ifpd db -h for usage details.")
    sys.exit()


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description="""Possible ifpd db queries.""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Possible ifpd db queries.",
    )
    parser.set_defaults(parse=default_parser)
    parser = ap.add_version_option(parser)

    sub_subparsers = parser.add_subparsers(
        title="sub-commands",
        help="Access the help page for a sub-command with: sub-command -h",
    )

    scripts.db_check.init_parser(sub_subparsers)
    scripts.db_dump.init_parser(sub_subparsers)
    scripts.db_info.init_parser(sub_subparsers)
    scripts.db_make.init_parser(sub_subparsers)
    scripts.db_merge.init_parser(sub_subparsers)
    scripts.db_reindex.init_parser(sub_subparsers)

    return parser
