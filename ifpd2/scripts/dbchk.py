"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from ifpd2.asserts import enable_rich_assert
from ifpd2.scripts import arguments as ap
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
    parser = ap.add_version_option(parser)
    parser.set_defaults(parse=parse_arguments, run=run)

    return parser


@enable_rich_assert
def parse_arguments(args: argparse.Namespace) -> argparse.Namespace:
    assert os.path.isdir(args.database), f"database folder not found: {args.database}"
    return args


@enable_rich_assert
def run(args: argparse.Namespace) -> None:
    raise NotImplementedError
