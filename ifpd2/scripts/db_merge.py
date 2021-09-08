"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from ifpd2.asserts import enable_rich_assert
from ifpd2.scripts import arguments as ap


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "merge",
        description="""
Merge multiple databases into one.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Merge databases.",
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
