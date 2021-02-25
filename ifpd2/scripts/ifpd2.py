"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
from ifpd2.const import __version__
from ifpd2.scripts import arguments as ap
from ifpd2 import scripts
import sys


def default_parser(*args) -> None:
    print("ifpd2 -h for usage details.")
    sys.exit()


def main():
    parser = argparse.ArgumentParser(
        description=f"""
Version:    {__version__}
Author:     Gabriele Girelli
Docs:       http://ggirelli.github.io/ifpd2
Code:       http://github.com/ggirelli/ifpd2

Another iFISH probe design pipeline (II).
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.set_defaults(parse=default_parser)
    parser = ap.add_version_option(parser)

    subparsers = parser.add_subparsers(
        title="sub-commands",
        help="Access the help page for a sub-command with: sub-command -h",
    )

    scripts.dbchk.init_parser(subparsers)
    scripts.query.init_parser(subparsers)

    args = parser.parse_args()
    args = args.parse(args)
    args.run(args)
